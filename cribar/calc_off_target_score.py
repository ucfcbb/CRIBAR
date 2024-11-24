from extract_on_target_primers import *
from select_primers import remove_primer_thres
from util import *
from calc_cfd import calc_score_batch
import os


def calc_off_target_score(work_dir, index_dir, on_target_primer2info, pam_f, mismatch_num):
    _primer_file = work_dir + "cribar_top_primer"
    all_target_pos_file = work_dir + "cribar_all_target_pos"
    write_primer(on_target_primer2info.keys(), _primer_file)

    run_crispritz(_primer_file, all_target_pos_file, index_dir, pam_f, mismatch_num, index=True)

    off_target_primer2info = parse_crispritz_res(all_target_pos_file)

    calc_score_batch(off_target_primer2info)

    remove_primer_thres(on_target_primer2info, off_target_primer2info)

    return off_target_primer2info


def calc_off_target_score_bowtie(work_dir,
                                 bowtie_index,
                                 on_target_primer2info,
                                 pam_seq,
                                 target_coor,
                                 mm_scores,
                                 pam_scores,
                                 mismatch_num=3,
                                 verbose=False,
                                 max_hits_on_genome=300):
    bowtie_in_file = work_dir + "bowtie_input.fa"
    bowtie_out_file = work_dir + "bowtie_out.sam"
    if os.path.isfile(bowtie_in_file):
        os.remove(bowtie_in_file)
    if os.path.isfile(bowtie_out_file):
        os.remove(bowtie_out_file)

    write_primer_to_fa(on_target_primer2info.keys(), pam_seq, bowtie_in_file)
    cmd = "bowtie -f -S -a -y -m " + str(max_hits_on_genome) + " -n " + str(mismatch_num) + " "+\
          bowtie_index + " " + bowtie_in_file + " " + bowtie_out_file
    if verbose:
        print (cmd)
    else:
        cmd += " >/dev/null 2>&1"

    process = subprocess.Popen([cmd], shell=True)
    process.wait()

    off_target_primer2info = parse_bowtie_result_by_pysam(bowtie_out_file, target_coor, pam_seq)

    calc_score_batch(off_target_primer2info, mm_scores, pam_scores)
    return off_target_primer2info
