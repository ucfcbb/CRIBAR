from collections import OrderedDict
import subprocess
from collections import defaultdict
from util import parse_crispritz_res
from calc_cfd import calc_cfd, calc_score_batch
from util import revcom


def extract_regions_by_getfasta(fasta_file, chr, start, end, work_dir, tar_dir):
    bed_file = work_dir + "tar.bed"
    output_fasta = tar_dir+"target.fa"
    with open(bed_file, 'w') as b:
        b.write("\t".join(map(str, [chr, start, end, 0, 0, "+"]))+"\n")
    cmd = "bedtools getfasta -s -fi {fi} -bed {bed} -fo {fo}".format(fi=fasta_file, bed=bed_file, fo=output_fasta)
    process = subprocess.Popen([cmd], shell=True)
    process.wait()

    with open(output_fasta) as f:
        f.readline()
        seq = f.readline()
    return seq


def collapse_primer(primer2info):
    hit2primer = defaultdict(set)
    primer2score = defaultdict(float)
    for primer in primer2info:
        for h in primer2info[primer]:
            hit2primer["$".join(map(str, h[1:-1]))].add(primer)
            primer2score[primer] += h[-1]

    to_rm = set()
    for hit in hit2primer:
        if len(hit2primer[hit]) <= 1:
            continue
        primers = list(hit2primer[hit])

        scores = [primer2score[p] for p in primers]
        to_keep = primers[scores.index(max(scores))]

        to_rm = to_rm.union(set(primers) - set([to_keep]))

    for p in to_rm:
        primer2info.pop(p)


def is_valid(primer_pam, GG, CC):
    if primer_pam[-2:] != GG and primer_pam[:2] != CC:
        return False
    return True


def have_forbidden_string(primer_pam, forbidden_strings, strand):
    if (strand=="-" and primer_pam[3] == 'A') or (strand=="+" and primer_pam[-4] == 'T'):
        return True
    for s in forbidden_strings:
        if strand=="-" and (revcom(s) in primer_pam):
            return True
        if strand=="+" and (s in primer_pam):
            return True
    return False


def get_primers_for_ritz(_seq, forbidden_strings, k=20, pam='NGG', gap = 20):
    # get the primer_count
    primers = set()
    len_pam = len(pam)
    len_primer = k + len_pam
    GG = pam[-2:]
    CC = revcom(GG)

    for pos in range(len(_seq) - len_primer):
        primer_pam = _seq[pos:pos + len_primer]

        if not is_valid(primer_pam, GG, CC):
            continue

        # primer always end with GG, if start with CC get reverse_complement
        if primer_pam[-2:] == GG and not have_forbidden_string(primer_pam, forbidden_strings, "+"):
            primers.add(primer_pam)

        if primer_pam[:2] == CC and not have_forbidden_string(primer_pam, forbidden_strings, "-"):
            primers.add(revcom(primer_pam))

    return primers


def get_primers_bruteforce(_seq, tar_chr, tar_start, forbidden_strings, pam='NGG', grna_len = 20):
    # get the primer_count
    primers = set()
    len_pam = len(pam)
    len_primer = grna_len + len_pam
    GG = pam[-2:]
    CC = revcom(GG)

    for pos in range(len(_seq) - len_primer):
        primer_pam = _seq[pos:pos + len_primer]

        if not is_valid(primer_pam, GG, CC):
            continue

        primer_seq, strand = "", "+"
        # primer always end with GG, if start with CC get reverse_complement
        if primer_pam[-2:] == GG and not have_forbidden_string(primer_pam, forbidden_strings, "+"):
            primer_seq = primer_pam
            if primer_seq != "":
                primers.add((primer_seq, tar_chr, tar_start+pos+1, "+"))
        elif primer_pam[:2] == CC and not have_forbidden_string(primer_pam, forbidden_strings, "-"):
            primer_seq = primer_pam
            if primer_seq != "":
                primer_seq = revcom(primer_seq)
                primers.add((primer_seq, tar_chr, tar_start+pos+1, "-"))
    return primers


def get_on_target_info_n2(primers, work_dir, target_seq, pam, mm_scores, pam_scores, mismatch_num):
    n, primers = len(primers), list(primers)
    primer_res = OrderedDict()
    for i in range(n):
        primer_seq_i = primers[i][0]
        key = primer_seq_i[:-3]
        primer_res[primer_seq_i] = set()
        for j in range(n):
            primer_seq_j = primers[j][0]
            if sum(a!=b for a, b in zip(primer_seq_i, primer_seq_j)) > int(mismatch_num):
                continue
            score = calc_cfd(key, primer_seq_j[:-3], primer_seq_j[-2:], mm_scores, pam_scores)
            primer_res[primer_seq_i].add(tuple(list(primers[j])+[score]))

    new_primer_res = OrderedDict()
    for p in primer_res:
        # sorted first by chr, then by pos
        new_primer_res[p[:-3]] = sorted(list(primer_res[p]), key=lambda x : (x[1], x[2]))

    return new_primer_res


def get_on_target_info_crispritz(primers, work_dir, target_seq_path, pam_f, mm_scores, pam_scores, crispritz_dir, mismatch_num):
    _primer_file = work_dir + "cribar_on_target_primer"
    with open(_primer_file, "w") as _primer_file_handler:
        for p in primers:
            if "N" in p[:-3]:    continue
            _primer_file_handler.write(p[:-3] + "NNN" +"\n")

    on_target_pos_file = work_dir + "cribar_on_target_pos"
    run_crispritz(_primer_file, on_target_pos_file, target_seq_path, pam_f, crispritz_dir, mismatch_num, index=False)

    primer_res = parse_crispritz_res(on_target_pos_file)

    calc_score_batch(primer_res, mm_scores, pam_scores)

    return primer_res


def run_crispritz(input_f, output_f, target_seq_f, pam_f, crispritz_dir, mismatch_num=1, index = False):
    # if not index, build index
    if index:
        cmd = crispritz_dir + "crispritz.py search "+ target_seq_f + " " + pam_f + " " + input_f+ " "+ output_f +" -index -mm "+mismatch_num+ " -bDNA 0 -bRNA 0 -r"
    else:
        cmd = crispritz_dir + "crispritz.py search " + target_seq_f + " " + pam_f + " " + input_f + " " + output_f + " -mm "+str(mismatch_num) + " -r"
    print (cmd)
    process = subprocess.Popen([cmd], shell=True)
    process.wait()
