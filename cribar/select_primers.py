from ILP import *
import operator
from collections import defaultdict
import time


def remove_primer_thres(on_tar_d, off_tar_d, max_off_target_score = 5000):
    for p in off_tar_d.keys():
        if len(off_tar_d[p]) > max_off_target_score:
            off_tar_d.pop(p, None)
            on_tar_d.pop(p, None)


def remove_worst_primer(ilpinput, primer_to_offtarget):
    p = max(primer_to_offtarget, key=lambda k: primer_to_offtarget[k])

    # on_target_primer2info.pop(p, None)
    keys = list(ilpinput.on_target_primer2info.keys())

    ilpinput.removed_primer.add(keys.index(p))


def get_on_target_positions(on_target_chr,on_target_start,on_target_end, primer_targets_dic):
    OnTarget_dict = defaultdict(dict)
    for primer in primer_targets_dic:
        OnTarget_dict[primer][on_target_chr] = [pos for pos in primer_targets_dic[primer][on_target_chr] if on_target_start <= pos <=on_target_end]
    return OnTarget_dict


# to-do cluster_len from tar_pos
def check_off_targets(primers, on_target_score, off_target_primer2info, tar_pos, verbose=False, off_target_ratio = 0.2, off_tar_win=0):
    if not off_tar_win:
        off_tar_win = tar_pos[2] - tar_pos[1]
    cluster_score_cutoff = on_target_score * off_target_ratio
    # print cluster_score_cutoff
    chr2off_targets = defaultdict(list)
    for primer in primers:
        # primer = primer[:-3] + 'NNN' # tmp, will remove
        for hit in off_target_primer2info.get(primer, list()):
            # hit: ('TCTcgAAACtGCAGAGGCAtTGG', 'chr22', 18121736, '+', 1)
            p_chr, p_pos, p_score = hit[1], hit[2], hit[-1]
            chr2off_targets[p_chr].append((p_pos, p_score, primer))

    if verbose:
        print ("start check off targets")

    for _chr in chr2off_targets:
        cur_chr_hits = sorted(chr2off_targets[_chr])
        primer_score = defaultdict(float)
        # two pointers, i is the first hit in the cluster_len, j is the last hit in the cluster_len
        i, score_sum = 0, 0
        for j in range(len(cur_chr_hits)):
            # ignore the target region
            if _chr == tar_pos[0] and tar_pos[1] <= cur_chr_hits[j][0] <= tar_pos[2]:
                continue
            score_sum += cur_chr_hits[j][1]
            primer_score[cur_chr_hits[i][2]] += cur_chr_hits[i][1]
            while cur_chr_hits[j][0] - cur_chr_hits[i][0] > off_tar_win and i <= j:
                score_sum -= cur_chr_hits[i][1]
                primer_score[cur_chr_hits[i][2]] -= cur_chr_hits[i][1]
                i += 1
            if score_sum > cluster_score_cutoff:
                # print "off_target_hits:", cur_chr_hits[i:j+1]
                return False, primer_score

    return True, None


def run_selection(on_target_primer2info, off_target_primer2info, tar_pos, formulation, constraint, start_time, verbose, off_target_ratio,
                  off_tar_win=0, min_gap = 30):
    check = False
    ilp_input = construct_data(on_target_primer2info, min_gap)

    if verbose:
        print("ilp_input:")
        print(ilp_input.__dict__)
        print("ILP input constructed")
    count = 0
    while not check:
        selected_primer_index, selected_primer, on_target_score, selected_hit = ilp(ilp_input, formulation, constraint)

        # print "selected_primer len:", len(selected_primer)

        count += 1
        if verbose:
            print ("ilp iteration:", count, selected_primer, "Time:", round(time.time() - start_time, 3), "seconds")

        if len(selected_primer) == 0:
            return False, None, None, None

        check, primer_to_offtarget = check_off_targets(selected_primer,
                                                       on_target_score,
                                                       off_target_primer2info,
                                                       tar_pos,
                                                       verbose,
                                                       off_target_ratio,
                                                       off_tar_win)
        # if len(off_target_copy_dict) == 0:
        #     print "check_off_targets_fail"
        #     exit(1)
        if verbose:
            print ("ILP iteration ", count, check)

        if check:
            return True, selected_primer, on_target_score, selected_hit
        elif count >= 30 or len(ilp_input.removed_primer) + 1 == ilp_input.primer_num:
            return False, None, None, None
        # print len(off_target_copy_dict)
        # remove primer with most contribution to check_off_targets
        remove_worst_primer(ilp_input, primer_to_offtarget)

        if verbose:
            print ("Removed candidate sgRNAs:", len(ilp_input.removed_primer))


