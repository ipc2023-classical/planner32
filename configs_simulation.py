heuristics = {"lmcut" : "lmcut()", "blind" : "blind()"}
pruning_dds = {"lmcut" : "bdd_map", "blind" : "bdd"}

pruning_types = {"expbelt" : "expansion",
                 "genbelt" : "generation", 
                 "par" : "parent, min_desactivation_ratio=0, min_insertions=infinity", 
                 "exp" : "expansion, min_desactivation_ratio=0, min_insertions=infinity" }

optionals_sim = {"inc" :  "intermediate_simulations=true, incremental_simulations=true", 
                 "irr" : "apply_subsumed_transitions_pruning=false, prune_dead_operators=true, store_original_operators=true", } 

optionals_prune = {"nospu" : "remove_spurious=false"}

simulation_type = {"sim" : "simulation_type=SIMPLE, label_dominance_type=NONE",
                   "bisim" : "simulation_type=NONE, label_dominance_type=NONE", 
                   "ldsim" : "simulation_type=SIMPLE, label_dominance_type=NORMAL", 
                   "noopsim" :  "simulation_type=SIMPLE, label_dominance_type=NOOP" }

shrinking = {
    "simsh" : "shrink_after_merge=true, shrink_strategy=shrink_bisimulation_perfect(),forbid_lr=false, apply_label_dominance_reduction=true,apply_simulation_shrinking=true",
    "bissh" : "shrink_after_merge=true, shrink_strategy=shrink_bisimulation_perfect(),forbid_lr=false",
    "nosh" : "shrink_after_merge=false, shrink_strategy=none(),forbid_lr=false", 
    "noshlr" : "shrink_after_merge=false, shrink_strategy=none(),forbid_lr=true"}


def get_optionals_sim(opt):
    res = []
    for o in opt:
        if o in optionals_sim:
            res.append(optionals_sim[o])
    return res

def get_optionals_prune(opt):
    res = []
    for o in opt:
        if o in optionals_prune:
            res.append(optionals_prune[o])
    return res

merge_strategies = { "dfp" : "merge_dfp()"
}

def get_merge (merge_params):
    if merge_params == "atomic": 
        return "limit_transitions_merge=1, limit_merge=1"
        
    merge = merge_strategies[merge_params[0:3]]
    limit = int(merge_params[3:].lower().replace("k", "000").replace("m", "000000"))
    return "merge_strategy=%s, limit_transitions_merge=%d, limit_merge=infinity" % (merge, limit)

def get_simulation_config (s):
    parts = s.split("-")
    h, simtype, merge, shrink, ptype, opt = parts[0], parts[1], parts[2], parts[3], parts[4],  parts[5:]

    heuristic = heuristics [h]
    pruning_dd = pruning_dds [h]
    pruning_type = pruning_types[ptype]
    merge = get_merge(merge)
    shrink = shrinking[shrink]
    default = "compute_final_simulation=true, switch_off_label_dominance=infinity"
    optional_sim = get_optionals_sim(opt)
    optional_pr = get_optionals_prune(opt)
    
    builder_params = ", ".join([default, merge] + optional_sim  )
    simulation_params = ", ".join(["pruning_dd=%s" % pruning_dd, "pruning_type=%s" %  pruning_type] + optional_pr  )
    config_pruning = "prune=simulation(%s, abs=builder_massim(%s))" % (simulation_params, builder_params)

    config = "astar(%s, %s)" % (heuristic, config_pruning)
    return config


def print_config(config): 
    print "%s %s" % (config, get_simulation_config(config))

# Experiment #1: simulation type and pruning types
merge_strategies_exp1 = ["atomic", "dfp10k", "dfp50k", "dfp100k", "dfp200k"]
for h in heuristics: 
    for sim in simulation_type:
        for mer in merge_strategies_exp1: 
            sh = "bissh"
            pr = "exp"
            config = "%s-%s-%s-%s-%s" %  (h, sim, mer, sh, pr)
            print_config (config)

for h in heuristics: 
    for sim in simulation_type:
        for pr in pruning_types: 
            mer = "dfp100k"
            sh = "bissh"
            config = "%s-%s-%s-%s-%s" %  (h, sim, mer, sh, pr)

            print_config(config)

