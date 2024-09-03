import subprocess
from collections import OrderedDict
import os
import pysam


def write_primer(primer, file):
    with open(file, "w") as _primer_file_handler:
        for p in primer:
            _primer_file_handler.write(p+"\n")


def write_primer_to_fa(primer, pam_seq, file):
    count = 1
    with open(file, "w") as _primer_file_handler:
        for p in primer:
            _primer_file_handler.write(">"+str(count)+"\n")
            _primer_file_handler.write(p+pam_seq+"\n")
            count += 1


def parse_bowtie_result(_f):
    ret = set()
    with open(_f) as fin:
        for line in fin:
            if line[0] == "@": continue
            line = line.strip().split()
            if len(line) > 3 and line[2] != "*":
                primer = line[9]
                if primer[-3:] == "NGG":
                    ret.add(primer[:-3] + "NNN")
                else:
                    ret.add(revcom(primer[3:]) + "NNN")
    return ret


# def parse_bowtie_result_1(_f):
#     primer_res = OrderedDict()
#     with open(_f + ".targets.txt") as fin:
#         for line in fin:
#             if line[0] == "@": continue
#             line = line.strip().split()
#
#             if not (len(line) > 3 and line[2] != "*"): # line[2] != "*" are unmapped reads
#                 continue
#
#             primer, md, chr, pos, dir = line[9], line[12], line[2], int(line[3]), "-" if int(line[1]) & 16 else "+"
#             if dir == "-":
#                 primer = Seq(primer).reverse_complement()
#
#             chr, pos, dir = line[2], int(line[3]), "-" if int(line[1]) & 16 else "+"
#             # primer -> dna, chr, pos, dir,
#             if primer not in primer_res:
#                 primer_res[primer] = set()
#             primer_res.get(primer).add(dna, chr, pos, dir)
#
#     for p in primer_res:
#         # sorted first by chr, then by pos
#         primer_res[p] = sorted(list(primer_res[p]), key=lambda x: (x[1], x[2]))
#     return primer_res


def parse_bowtie_result_by_pysam(_f, target_coor, pam_seq):
    primer_res = OrderedDict()
    with pysam.AlignmentFile(_f) as samfile:
        for read in samfile:
            if read.is_unmapped:
                continue
            dir = "-" if read.is_reverse else "+"
            dna = read.get_reference_sequence()
            primer = read.query_alignment_sequence
            chr = read.reference_name
            pos = read.reference_start

            if dir == "-":
                primer = revcom(primer)
                dna = revcom(dna)

            primer = primer[:-3]
            if primer not in primer_res:
                primer_res[primer] = set()
            primer_res[primer].add((dna, chr, pos, dir))

    for p in primer_res:
        # sorted first by chr, then by pos
        primer_res[p] = sorted(list(primer_res[p]), key=lambda x: (x[1], x[2]))
    return primer_res


def parse_bowtie_result_by_pysam_3mis(_f, target_coor):
    primer_res = OrderedDict()
    with pysam.AlignmentFile(_f) as samfile:
        for read in samfile:
            if read.is_unmapped:
                continue
            dir = "-" if read.is_reverse else "+"
            dna = read.get_reference_sequence()
            primer = read.query_alignment_sequence
            chr = read.reference_name
            pos = read.reference_start

            if read.is_reverse:
                primer = revcom(primer)[:-3] + "NNN"
                dna = revcom(dna)
            else:
                primer = primer[:-3] + "NNN"
            if primer not in primer_res:
                primer_res[primer] = set()
            primer_res[primer].add((dna, chr, pos, dir))

    for p in primer_res:
        # sorted first by chr, then by pos
        primer_res[p] = sorted(list(primer_res[p]), key=lambda x: (x[1], x[2]))
    return primer_res


def bowtie_filter(on_target_primer2info, cutoff, work_dir):
    bowtie_in_file = work_dir + "bowtie_input.fa"
    bowtie_out_file = work_dir + "bowtie_out.sam"
    bowtie_index = "GRCh38_noalt_as/GRCh38_noalt_as"
    # bowtie_index = "/home/xiaoli/cribar/GRCh38_noalt_as/GRCh38_noalt_as"  # coombs
    # write_primer_to_fa(on_target_primer2info.keys(), bowtie_in_file)
    write_primer_to_fa(on_target_primer2info.keys(), bowtie_in_file)
    cmd = "bowtie -f -S -m " + str(cutoff) + " -n 3 "+bowtie_index +" "+ bowtie_in_file +" "+ bowtie_out_file
    process = subprocess.Popen([cmd], shell=True)
    process.wait()

    less_than_cutoff_primers = parse_bowtie_result(bowtie_out_file)
    for p in on_target_primer2info.keys():
        if p not in less_than_cutoff_primers:
            on_target_primer2info.pop(p)


def parse_crispritz_res(f):
    primer_res = OrderedDict()
    with open(f + ".targets.txt") as crispritz_res:
        crispritz_res.readline()
        for line in crispritz_res:
            if line[0] == "#": continue
            line = line.strip().split()
            if len(line) < 7: break
            # primer -> dna, chr, pos, dir,
            p = line[1][:-3]
            if p not in primer_res:
                primer_res[p] = set()
            primer_res.get(p).add((line[2], line[3], int(line[4]), line[6]))

    for p in primer_res:
        # sorted first by chr, then by pos
        primer_res[p] = sorted(list(primer_res[p]), key=lambda x : (x[1], x[2]))
    return primer_res


def calc_score_for_primer(primer2info):
    primer2score = {}
    for p in primer2info:
        primer2score[p] = sum(i[-1] for i in primer2info[p])
    return primer2score


def read_forbidden_strings(src_dir):
    ret = set()
    if not os.path.isfile(src_dir + "forbidden_strings"):
        print ("Cannot find forbidden_strings")
        return ret
    with open(src_dir + "forbidden_strings") as f:
        ret.add(f.readline().strip())
    return [x for x in ret if x]

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A', 'N':'N',
                'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'u': 'a', 'n': 'n'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)