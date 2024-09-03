from extract_on_target_primers import get_primers_for_ritz, collapse_primer, extract_regions_by_getfasta
from extract_on_target_primers import get_on_target_info_n2, get_primers_bruteforce, get_primers_for_ritz, get_on_target_info_crispritz
from calc_off_target_score import calc_off_target_score, calc_off_target_score_bowtie
from select_primers import run_selection
from util import bowtie_filter, calc_score_for_primer, read_forbidden_strings
import time
import argparse
import os


USE_CRISPRITZ = False
ONLY_CONSIDER_TOP100_GRNAS = False
COLLAPSE_PRIMER = False


def get_parser():
    # basic options
    parser = argparse.ArgumentParser()
    parser.add_argument('--chr', type=str, help='chr of target region in the fasta file')
    parser.add_argument('--start', type=int,
        help='Starting position of target region in the fasta file. \
        The starting position is included in the target region.')
    parser.add_argument('--end', type=int,
        help='Ending position of target region in the fasta file. \
        The ending position is not included in the target region.')

    parser.add_argument('--len', type=int, default=20, help='The length of gRNA (not include pam).')
    parser.add_argument('--formulation', type=int,
                        help='1: Constrain is on-target activity score, object is the minimized gRNA. \n \
                             2: Constrain is the gRNA number, object is the maximized on-target activity score.')
    parser.add_argument('--constraint', type=float,
                        help='When formulation=1, constrain is on-target activity score.\n \
                            When formulation=2, constrain is the gRNA number.')
    parser.add_argument('--pam_seq', nargs='?', default="NGG", type=str,
                        help='Default pam is NGG.')
    # TODO support more genomes
    # parser.add_argument('--reference-genome', nargs='?', default="hg38", type=str)

    # advanced options
    parser.add_argument('--mismatch', type=int,
                        help='The number of mismatches allowed in each gRNA in the target region')
    parser.add_argument('--off_target_window', nargs='?', default=0, type=int,
                        help='The window size that CRIBAR uses to check the off-target. \n \
                              Default value is the length of the target region')
    parser.add_argument('--off_target_ratio', nargs='?', default=0.2, type=float,
                        help='The threshold of off-target binding site density / on-target binding site density')
    parser.add_argument('--min_gap', nargs='?', default=30, type=int,
                        help='The minimum distance between two binding sites')
    parser.add_argument('--excluded_substring', type=str,
                        help='Excluded substrings separated by a comma. Example: AAAA,TTTT')

    # internal parameters
    parser.add_argument('--crispritz_genome_index_dir', type=str)
    parser.add_argument('--genome_prefix', type=str)
    parser.add_argument('--src_dir', type=str)
    parser.add_argument('--work_dir', nargs='?', default=os.getcwd() + "/", type=str)
    # target seq dir, CRISPRitz searches all fa files in a folder.
    # We re-use the target sequence path for bruteforce for simplicity
    parser.add_argument('--tar_seq_dir', nargs='?', default="tar_seq/", type=str)
    parser.add_argument('--crispritz_dir', nargs='?', default="../lib/linux-64_crispritz-2.6.6-py39h68928f9_1/bin/", type=str)
    parser.add_argument('--verbose', action='store_true')

    return parser


if __name__ == '__main__':
    start_time = time.time()
    args = get_parser().parse_args()

    # basic options
    tar_chr, tar_start, tar_end = args.chr, args.start, args.end
    pam_seq = args.pam_seq
    grna_len = args.len
    formulation = args.formulation
    constraint = args.constraint

    # advanced options
    mismatch_num = str(args.mismatch)
    off_target_ratio = args.off_target_ratio
    off_target_window = args.off_target_window
    excluded_substring = None if args.excluded_substring is None else set(args.excluded_substring.split(","))
    min_gap = args.min_gap

    # internal parameters
    genome_prefix = args.genome_prefix
    work_dir = args.work_dir
    src_dir = args.src_dir
    tar_dir = args.tar_seq_dir

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    if not os.path.exists(tar_dir):
        os.makedirs(tar_dir)

    fasta_file = genome_prefix + ".fa"
    bowtie_index = genome_prefix

    verbose = args.verbose

    crispritz_dir = args.crispritz_dir

    if not off_target_window:
        off_target_window = tar_end - tar_start + 1
    seq = extract_regions_by_getfasta(fasta_file, tar_chr, tar_start, tar_end, work_dir, tar_dir).upper()

    if not excluded_substring:
        forbidden_strings = read_forbidden_strings(src_dir)
    else:
        forbidden_strings = excluded_substring

    if os.path.isfile(crispritz_dir + "crispritz.py"):
        USE_CRISPRITZ = True
        if verbose:
            print ("Use CRISPRitz for on-target search")
    else:
        if verbose:
            print ("Did not find CRISPRitz, use build-in on-target searching")

    if USE_CRISPRITZ:
        pam_file = work_dir+"crispritz_pam.txt"
        with open(pam_file,'w') as pam_file_handler:
            pam_file_handler.write("N"*grna_len + pam_seq + " " + str(len(pam_seq)) + "\n")

        crispritz_index_dir = args.crispritz_genome_index_dir

        primers = get_primers_for_ritz(seq, forbidden_strings, grna_len, pam_seq)
        on_target_primer2info = get_on_target_info_crispritz(primers, work_dir, tar_dir, pam_file, src_dir, crispritz_dir, mismatch_num)
    else:
        primers = get_primers_bruteforce(seq, tar_chr, tar_start, forbidden_strings, pam_seq, grna_len)
        on_target_primer2info = get_on_target_info_n2(primers, work_dir, seq, pam_seq, src_dir, mismatch_num)

    # if ONLY_CONSIDER_TOP100_GRNAS:
    #     for p in sorted(on_target_primer2info.keys(), reverse=True, key=lambda x: len(on_target_primer2info[x]))[:100]:
    #         print (p, len(on_target_primer2info[p]))

    if on_target_primer2info is None or len(on_target_primer2info) == 0:
        print("Failed to enumerate sgRNA candidate in the target region. Please consider to change or extend target region.")
        exit(0)

    print("Candidate sgRNAs: ", len(on_target_primer2info))
    if verbose:
        print ("Enumerated candidate sgRNAs. Time:", round(time.time() - start_time, 3), "seconds")
        print ("Parsed candidate sgRNA info. Time:", round(time.time() - start_time, 3), "seconds")

    # if COLLAPSE_PRIMER:
    #     collapse_primer(on_target_primer2info)
    #     print "collapsed primers", len(on_target_primer2info)
    #     bowtie_filter(on_target_primer2info, 100, work_dir)
    #     # print primer after bowtie filter
    #     for p in sorted(on_target_primer2info.keys(), reverse=True, key=lambda x: len(on_target_primer2info[x]))[:100]:
    #         print (p, len(on_target_primer2info[p]))
    #     print "primers after bowtie filter", len(on_target_primer2info)

    primer_to_score = calc_score_for_primer(on_target_primer2info)

    # ritz
    # off_target_primer2info = calc_off_target_score(work_dir, fasta_dir, on_target_primer2info, pam_file, mismatch_num)

    # bowtie
    off_target_primer2info = calc_off_target_score_bowtie(work_dir,
                                                          bowtie_index,
                                                          on_target_primer2info,
                                                          src_dir,
                                                          pam_seq,
                                                          (tar_chr, tar_start, tar_end),
                                                          mismatch_num,
                                                          verbose)

    # remove with too many hits
    to_remove = set()
    for p in on_target_primer2info.keys():
        if p not in off_target_primer2info.keys():
            to_remove.add(p)

    for p in to_remove:
        on_target_primer2info.pop(p)
        primer_to_score.pop(p)

    to_remove = set()
    top_primers = sorted(primer_to_score.keys(), key=lambda x: primer_to_score[x], reverse=True)[:200]
    for p in on_target_primer2info.keys():
        if p not in top_primers:
            to_remove.add(p)

    for p in to_remove:
        on_target_primer2info.pop(p)

    if verbose:
        print("Parsed off target info. Time:", round(time.time() - start_time, 3), "seconds")
        print("off_target_primer2info size:", len(off_target_primer2info))
        print("on_target_primer2info size after removal:", len(on_target_primer2info))

    success, primers, score, hits = run_selection(on_target_primer2info,
                                                  off_target_primer2info,
                                                  (tar_chr, tar_start, tar_end),
                                                  formulation,
                                                  constraint,
                                                  start_time,
                                                  verbose,
                                                  off_target_ratio,
                                                  off_target_window,
                                                  min_gap)

    if verbose:
        print("Selected optimal sgRNAs. Time:", round(time.time() - start_time, 3), "seconds")

    if success:
        print("")
        print("Num of selected sgRNAs:", len(primers))

        print("Selected sgRNAs:")

        print("\n".join(map(lambda p : p+pam_seq, primers)))

        print("")
        print("On-target activity score:", score)

        print("On-target hits:")

        print("chr:start:strand" + " "*14 + "hit sequence" + " "*18 + "activity_score")

        for hit in hits:
            _hit = hit[1]
            hit_seq = _hit[0]
            chr = _hit[1]
            if ":" in chr:
                chr = chr.split(":")[0]
            start = str(_hit[2])
            strand = _hit[3]
            score = str(_hit[4])
            end = str(_hit[2] + grna_len + len(pam_seq))

            pos = chr + ":" + start + "-" + end + ":" + strand
            pos = pos + " "*(30 - len(pos))
            hit_seq = hit_seq + " "*(30 - len(hit_seq))
            print(pos + hit_seq + score)


    else:
        print("sgRNA search failed. Please relax the constraints. Consider to:\n"
              "1. For formulation 1, reduce the constraint value. For formulation 2, increase the constraint value.\n"
              "2. Change or extend the target region. This may increase the candidate sgRNA number.\n"
              "3. Increase mismatches allowed value, to search the binding sites with lower similarity.\n"
              "4. Reduce the off-target window or increase off-target ratio value, to relax the off-target threshold.\n"
              "5. Reduce the minimal gap value.\n")

    # TODO remove intermediate files

