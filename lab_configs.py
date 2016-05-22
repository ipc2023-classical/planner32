#!/usr/bin/env python

import sys 

# Merge and shrink for unsolvability stuff:                                                                                          
merge_order = {
    'lev' : 'level',
    'rnd' : 'random',
    'rlev' : 'reverse_level'
}
                  
def get_merge_criteria(crit):
    def get_subopt_factor(s):
        opt_factor = '1'
        opt_diff = '0'
        if len(crit) > 5:
            if 'f' in s:
                if s[0] == 'f':
                    factor = s[1:]
                else:
                    opt_diff, factor = s.split('f')
                opt_factor = '0.' + factor
            else:
                opt_diff = s
        return opt_factor, opt_diff

    if crit == 'goal': #goal                                                                                                                           
        return 'goal'
    if crit.startswith("cg"): #cg, cgrev, cgcom, cgrevcom                                                                                              
        return 'scc(reverse=%s, complete_cg=%s)' % (str('rev' in crit), str('com' in crit))
    elif crit.startswith('empty'): #empty, empty10f9, ....                                                                                             
        opt_factor, opt_diff = get_subopt_factor(crit[5:])
        return 'tr(goal=false, empty=false, opt_factor=%s, opt_diff=%s), tr(goal=false, empty=true, opt_factor=%s, opt_diff=%s), tr(goal=true, empty=t\
rue, opt_factor=%s, opt_diff=%s)' % (opt_factor, opt_diff, opt_factor, opt_diff, opt_factor, opt_diff)
    elif crit.startswith('trnum'): #trnum, trnum10f9, ....                                                                                             
        opt_factor, opt_diff = get_subopt_factor(crit[5:])
        return 'tr(goal=false, empty=false, opt_factor=%s, opt_diff=%s)' % (opt_factor, opt_diff)
    elif crit.startswith('trempty'): #trempty, trempty10f9, ....                                                                                       
        opt_factor, opt_diff = get_subopt_factor(crit[7:])
        return 'tr(goal=false, empty=true, opt_factor=%s, opt_diff=%s)' % (opt_factor, opt_diff)
    elif crit.startswith('trgoal'): #trgoal, trgoal10f9, ....                                                                                          
        opt_factor, opt_diff = get_subopt_factor(crit[6:])
        return 'tr(goal=true, empty=true, opt_factor=%s, opt_diff=%s)' % (opt_factor, opt_diff)
    return None

def get_new_linear_merge_strategy(criteria, var_order):
    return "merge_linear_criteria(criteria=[%s], var_order=%s)" % (",".join(map(lambda x : get_merge_criteria(x),  criteria)), merge_order[var_order])


def get_merge_strategy(merge):
    if merge == 'dfp':
        return 'merge_dfp()'
    else:
        if '_' in merge:
            sp = merge.split('_')
            return get_new_linear_merge_strategy(sp[0:-1], sp[-1])
        else:
            return get_new_linear_merge_strategy([], merge)

shrink_bisimulation="shrink_bisimulation(max_states=infinity, greedy=false, threshold=1)"
shrink_bisimulation_agggoals="shrink_bisimulation(max_states=infinity, greedy=false, threshold=1, aggregate_goals=true)"
shrink_own_label="shrink_composite(strategies=[shrink_own_labels(goal_shrinking=true), %s])" % shrink_bisimulation_agggoals

shrink_st = {
    'bisimulation' :  shrink_bisimulation,
    'ownlabel' : shrink_own_label
}

label_reduction_old="label_reduction_method=old"

def get_mas(merge_strategy, shrink_strategy, parameters):
    return "merge_and_shrink(shrink_after_merge=true, merge_strategy=%s, shrink_strategy=%s, %s)" % (get_merge_strategy(merge_strategy), shrink_st[shrink_strategy], ",".join(parameters))

def get_mas_unsat(merge_strategy, shrink_strategy, parameters=[]):
    return get_mas(merge_strategy, shrink_strategy, ["cost_type=3"]+parameters)




#pruning_type=generation, pruning_dd=bdd_map, limit_merge=10000,            
def get_simulation(merge_strategy=None, parameters=[]):
    fixed_parameters = "pruning_dd=bdd, limit_seconds=300, use_bisimulation=true,pruning_type=expansion, remove_spurious=true, min_insertions=1000, min_desactivation_ratio=0.0, limit_merge=infinity, limit_transitions_merge=100000, apply_simulation_shrinking=false, apply_label_dominance_reduction=false"

    parameters = fixed_parameters
    
    if "op" in parameters:
        parameters += ", prune_dead_operators=true, store_original_operators=true"
    else:
        parameters += ", prune_dead_operators=false, store_original_operators=false"

    if "trpr":
        parameters += ", apply_subsumed_transitions_pruning=true"
    else:
        parameters += ", apply_subsumed_transitions_pruning=false"

    if "inc" in parameters:
        parameters += ", intermediate_simulations=true,incremental_simulations=true"
    else:
        parameters += ", intermediate_simulations=false,incremental_simulations=false"
        

    merge = "merge_strategy=%s" % get_merge_strategy(merge_strategy) if merge_strategy else ""
    return "prune=simulation(%s, %s)" % (merge, ",".join(parameters))


def get_config(configuration):
    config = configuration.split("-")

    heuristic = "blind"
    pruningMethod = None

    if config[0] == 'mas':
        if len(config) < 3:
            print "Error: python lab_configs.py mas-merge-shrink[-params]" 
            exit()
        heuristic = get_mas(config[1], config[2], config[3:])
    elif config[0] == 'masunsat':
        heuristic = get_mas_unsat(config[1], config[2], config[3:])
    elif config[0] == 'dom':
        
        pruningMethod = "prune = %s" % get_simulation(config[1], config[2:])

    return "astar(%s, %s)" % (heuristic, pruningMethod if pruningMethod else "")
 
if len(sys.argv) >= 2: 
    print get_config(sys.argv[1])
else:

    default_configs = ["masunsat-dfp-ownlabel", ]
    #"masunsat-cgcom_goal_lev-ownlabel"]

    print "[\n" + ", \n".join(["['fd_simulation/%s', REPO_MAS, \"('%s', ['--search', '%s'])\"" %  (c, c, get_config(c)) for  c in default_configs]) + "]\n]"


