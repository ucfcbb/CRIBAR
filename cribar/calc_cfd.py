from util import revcom


# Unpickle mismatch scores and PAM scores
def get_score_dict(input_text):
    try:
        result_dict = {}

        for line in input_text.strip().splitlines():
            key, value = line.split('=')  # Split only on the first '='
            try:
                result_dict[key.strip()] = float(value.strip())  # Convert value to float
            except:
                print(f"Failed to parse score file line: {line}")
                exit(0)
        return result_dict
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


# Calculates CFD score
# https://github.com/maximilianh/crisporWebsite/tree/master/CFD_Scoring
def calc_cfd(wt, sg, pam, mm_scores, pam_scores):
    # cfd score only works for 20+3 nt sgRNA
    # for sgRNA with other length, return 1 to only count the number of hits
    if len(wt) != 20 or len(sg) != 20:
        return 1

    if "N" in wt:
        return 0
    if "N" in sg:
        return 0
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


def calc_score_batch(primer_res, mm_scores, pam_scores):
    for p in primer_res.keys():
        tmp = []
        for hit in primer_res[p]:
            # primer -> dna, chr, pos, dir, score
            hit_seq = hit[0].upper()
            pam = hit_seq[-2:]
            tmp.append(tuple(list(hit) + [calc_cfd(p, hit_seq[:-3], pam, mm_scores, pam_scores)]))
        primer_res[p] = tmp

