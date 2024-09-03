import pickle
import re
from util import revcom

# Unpickle mismatch scores and PAM scores
def get_mm_pam_scores(src_dir):
    try:
        mm_scores = pickle.load(open(src_dir + '/mismatch_score.pkl','rb'))
        pam_scores = pickle.load(open(src_dir + '/pam_score.pkl','rb'))
        return mm_scores, pam_scores
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


# Calculates CFD score
def calc_cfd(wt, sg, pam, src_dir):
    # cfd score only works for 20+3 nt sgRNA
    # for sgRNA with other length, return 1 to only count the number of hits
    if len(wt) != 20 or len(sg) != 20:
        return 1

    if "N" in wt:
        return 0
    if "N" in sg:
        return 0
    mm_scores, pam_scores = get_mm_pam_scores(src_dir)
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i, sl in enumerate(s_list):
        if wt_list[i] == sl:
            continue
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score *= mm_scores[key]
    score *= pam_scores[pam]
    return score


def calc_score_batch(primer_res, src_dir):
    for p in primer_res.keys():
        tmp = []
        for hit in primer_res[p]:
            # primer -> dna, chr, pos, dir, score
            hit_seq = hit[0].upper()
            pam = hit_seq[-2:]
            tmp.append(tuple(list(hit) + [calc_cfd(p, hit_seq[:-3], pam, src_dir)]))
        primer_res[p] = tmp

