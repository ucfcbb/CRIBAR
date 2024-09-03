from pulp import *


class ilpInput:
    def __init__(self, hit_num, primer_num, hits_list, conflict_hits, on_target_primer2info, removed_primer):
        self.hit_num = hit_num
        self.primer_num = primer_num
        self.hits_list = hits_list
        self.conflict_hits = conflict_hits
        self.on_target_primer2info = on_target_primer2info
        self.removed_primer = removed_primer


def is_conflict(hit1, hit2, min_gap=30):
    if hit1[1] != hit2[1]:
        return False
    return hit1[2] + min_gap > hit2[2]


def construct_data(on_target_primer2info, min_gap=30):
    hit_num = sum(len(on_target_primer2info[i]) for i in on_target_primer2info)
    primer_num = len(on_target_primer2info)

    hits_list = []
    for p in on_target_primer2info.keys():
        for hit in on_target_primer2info[p]:
            hits_list.append(hit)

    # conflict hits
    conflict_hits = set()
    hits_to_index = {h:i for i, h in enumerate(hits_list)}
    # sort by chr and pos
    sorted_hits = sorted(hits_list, key=lambda x : (x[1], x[2]))

    for i in range(len(sorted_hits)-1):
        for j in range(i+1, len(sorted_hits)):
            if is_conflict(sorted_hits[i], sorted_hits[j], min_gap):
                conflict_hits.add((hits_to_index[sorted_hits[i]], hits_to_index[sorted_hits[j]]))
            else:
                break

    return ilpInput(hit_num, primer_num, hits_list, conflict_hits, on_target_primer2info, set())


def ilp(ilp_input, problem_type, constraint=5, min_gap=30):
    ilp_res = None
    if problem_type == 1:
        # constraint is on-target activity score threshold
        ilp_res = ilp1(ilp_input, constraint)
    elif problem_type == 2:
        # constraint is primer number
        ilp_res = ilp2(ilp_input, constraint)

    selected_primer_index = []
    selected_primer = []
    selected_hit = []

    on_target_primer2info_keys = list(ilp_input.on_target_primer2info.keys())
    for v in ilp_res:
        if v.varValue == 1:
            vtype, vnum = v.name.split("_")
            if vtype == "P":
                p = int(vnum)
                selected_primer_index.append(p)
                selected_primer.append(on_target_primer2info_keys[p])
            else:
                selected_hit.append((v.name, ilp_input.hits_list[int(vnum)]))

    # print "selected_primer:", selected_primer
    # print "selected_hit", selected_hit

    on_target_score = sum([x[1][-1] for x in selected_hit])

    # print "on_target_score", on_target_score

    return selected_primer_index, selected_primer, on_target_score, selected_hit


def ilp1(ilp_input, min_score=5, min_gap=30):
    prob = LpProblem("problem_1", LpMinimize)

    M = sys.maxsize
    hit_val = LpVariable.dicts(name='H', indexs=range(ilp_input.hit_num), lowBound=0, upBound=1, cat="Integer")
    primer_val = LpVariable.dicts(name='P', indexs=range(ilp_input.primer_num), lowBound=0, upBound=1, cat="Integer")

    # obj
    prob += pulp.lpSum([primer_val[i] for i in range(ilp_input.primer_num)])
    # score constraint, score = hits_list[i][-1]
    prob += pulp.lpSum([hit_val[i] * ilp_input.hits_list[i][-1] for i in range(ilp_input.hit_num)]) >= min_score

    # conflict hit constraint
    for c in ilp_input.conflict_hits:
        prob += hit_val[c[0]] + hit_val[c[1]] <= 1

    # removed primer
    for i in ilp_input.removed_primer:
        prob += primer_val[i] <= 0

    # primer and hit constraint
    cur_hit, cur_primer = 0, 0
    for p in ilp_input.on_target_primer2info.keys():
        num_of_hit_for_p = len(ilp_input.on_target_primer2info[p])
        prob += pulp.lpSum([hit_val[i]] for i in range(cur_hit, cur_hit+num_of_hit_for_p)) <= M * primer_val[cur_primer]
        cur_hit += num_of_hit_for_p
        cur_primer += 1

    prob.solve(GLPK_CMD(msg=0, timeLimit=60))

    return prob.variables()


def ilp2(ilp_input, max_primers=3, min_gap=30):
    prob = LpProblem("problem_2", LpMaximize)

    M = sys.maxsize
    hit_val = LpVariable.dicts(name='H', indexs=range(ilp_input.hit_num), lowBound=0, upBound=1, cat="Integer")
    primer_val = LpVariable.dicts(name='P', indexs=range(ilp_input.primer_num), lowBound=0, upBound=1, cat="Integer")

    # primer constraint
    prob += pulp.lpSum([primer_val[i] for i in range(ilp_input.primer_num)]) <= max_primers
    # obj, score = hits_list[i][-1]
    prob += pulp.lpSum([hit_val[i] * ilp_input.hits_list[i][-1] for i in range(ilp_input.hit_num)])

    # conflict hit constraint
    for c in ilp_input.conflict_hits:
        prob += hit_val[c[0]] + hit_val[c[1]] <= 1

    # removed primer
    for i in ilp_input.removed_primer:
        prob += primer_val[i] <= 0

    # primer and hit constraint
    cur_hit, cur_primer = 0, 0
    for p in ilp_input.on_target_primer2info.keys():
        num_of_hit_for_p = len(ilp_input.on_target_primer2info[p])
        prob += pulp.lpSum([hit_val[i]] for i in range(cur_hit, cur_hit+num_of_hit_for_p)) <= M * primer_val[cur_primer]
        cur_hit += num_of_hit_for_p
        cur_primer += 1

    prob.solve(GLPK_CMD(msg=0,timeLimit=60))

    return prob.variables()
