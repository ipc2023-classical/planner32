
heuristics = {"lmcut" : "lmcut()", "blind" : "blind()"}

def pruning_dds (h, ptype):
    if h == "blind" and ptype != "gen":
        return "bdd"
    return "bdd_map"

pruning_types = {"par" : "parent", 
                 "exp" : "expansion",
                 "gen" : "generation",
}

optionals_sim = {"inc" :  "intermediate_simulations=true, incremental_simulations=true", 
                 "irr" : "apply_subsumed_transitions_pruning=false, prune_dead_operators=true, store_original_operators=true",
                 "usedominatesin" : "num_labels_to_use_dominates_in=5000", } 

optionals_prune = {"nospu" : "remove_spurious=false"}

simulation_type = {"sim" : ["simulation_type=SIMPLE, label_dominance_type=NONE"],
                   "bisim" : ["simulation_type=NONE, label_dominance_type=NONE"], 
                   "ldsim" : ["simulation_type=SIMPLE, label_dominance_type=NORMAL"],
                   "ldsimalt" : ["simulation_type=SIMPLE, label_dominance_type=ALTERNATIVE"], 
                   "noopsim" :  ["simulation_type=SIMPLE, label_dominance_type=NOOP"]
}

numeric_pruning_types = {"parsucc" : "prune_successors=true, prune_dominated_by_parent=true, prune_dominated_by_closed=false, prune_dominated_by_open=false",
                         "par" : "prune_successors=false, prune_dominated_by_parent=true, prune_dominated_by_closed=false, prune_dominated_by_open=false", 
                         "exp" : "prune_successors=false, prune_dominated_by_parent=false, prune_dominated_by_closed=true, prune_dominated_by_open=false",
                         "gen" : "prune_successors=false, prune_dominated_by_parent=false, prune_dominated_by_closed=false, prune_dominated_by_open=true",
                         "gensucc" : "prune_successors=true, prune_dominated_by_parent=true, prune_dominated_by_closed=false, prune_dominated_by_open=true",
                         "succ" : "prune_successors=true, prune_dominated_by_parent=false, prune_dominated_by_closed=false, prune_dominated_by_open=false"     
}

numeric_simulation_type = {
    "qual"  : ["use_quantified_dominance=false"],
    "qtrade" : ["use_quantified_dominance=true, trade_off_dominance=true"],
    "qpos" : ["use_quantified_dominance=true, only_positive_dominance=true"],
    "qrel" : ["use_quantified_dominance=true"]
}




optionals_sim = {#"inc" :  "intermediate_simulations=true, incremental_simulations=true", 
                 #"irr" : "apply_subsumed_transitions_pruning=false, prune_dead_operators=true, store_original_operators=true", }
}

shrinking = {
    "simsh" : "shrink_after_merge=true, shrink_strategy=shrink_bisimulation_perfect(),forbid_lr=false, apply_label_dominance_reduction=true,apply_simulation_shrinking=true",
    "bissh" : "shrink_after_merge=true, shrink_strategy=shrink_bisimulation_perfect(),forbid_lr=false",
    "nosh" : "shrink_after_merge=false, shrink_strategy=none(),forbid_lr=false", 
    "noshlr" : "shrink_after_merge=false, shrink_strategy=none(),forbid_lr=true"
}

merge_strategies = {
    "dfp" : "merge_dfp()"
}


def get_merge (merge_params):
    if merge_params == "atomic": 
        return "limit_transitions_merge=1, limit_merge=1"
        
    merge = merge_strategies[merge_params[0:3]]
    if merge_params.endswith('states'):
        merge_params = merge_params.replace('states', '')
        limit_on = "limit_transitions_merge=infinity, limit_merge={}"
    else:
        limit_on = "limit_transitions_merge={}, limit_merge=infinity"
        
    limit = int(merge_params[3:].lower().replace("k", "000").replace("m", "000000"))
    return "merge_strategy={}, ".format(merge) + limit_on.format (limit)


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
    if "belt" in opt:
        res.append("min_desactivation_ratio=0, min_insertions=1000")
    else:
        res.append("min_desactivation_ratio=0, min_insertions=infinity")
    return res



def get_simulation_config (s):
    parts = s.split("-")
    h, simtype, merge, shrink, ptype, opt = parts[0], parts[1], parts[2], parts[3], parts[4],  parts[5:]

    heuristic = heuristics [h]
    pruning_dd = pruning_dds (h, ptype)
    pruning_type = pruning_types[ptype]
    merge = get_merge(merge)
    shrink = shrinking[shrink]
    default = "compute_final_simulation=true, switch_off_label_dominance=infinity"
    optional_sim = get_optionals_sim(opt)
    optional_pr = get_optionals_prune(opt)
    
    builder_params = ", ".join(simulation_type[simtype] + [default, merge, shrink] + optional_sim  )
    simulation_params = ", ".join( ["pruning_dd=%s" % pruning_dd, "pruning_type=%s" %  pruning_type] + optional_pr  )
    config_pruning = "prune=simulation(%s, abs=builder_massim(%s))" % (simulation_params, builder_params)

    config = "astar(%s, %s)" % (heuristic, config_pruning)
    return config



def get_numeric_simulation_config (s):
    parts = s.split("-")
    h, simtype, trval, merge, shrink, ptype, opt = parts[0], parts[1], parts[2], parts[3], parts[4],  parts[5], parts[6:]

    heuristic = heuristics [h]
    pruning_dd = "use_single_bdd=true" if h == "blind" and ptype != "gen" and  simtype == "qual"  else "use_single_bdd=false"
    
    pruning_type = numeric_pruning_types[ptype]
    merge = get_merge(merge)
    shrink = shrinking[shrink]
    default = "compute_final_simulation=false, switch_off_label_dominance=infinity"
    optional_sim = get_optionals_sim(opt)
    optional_pr = get_optionals_prune(opt)

    if "fulltau" in opt:
        optional_pr += ["tau_labels_self_loops=true", "tau_labels_recursive=true", "tau_labels_noop=true"]
    elif "recurtau" in opt:
        optional_pr += ["tau_labels_self_loops=true", "tau_labels_recursive=true", "tau_labels_noop=false"]   
    else:
        optional_pr += ["tau_labels_self_loops=true", "tau_labels_recursive=false", "tau_labels_noop=false"]   
        
    builder_params = ", ".join([default, merge] + optional_sim  )
    simulation_params = ", ".join(numeric_simulation_type[simtype] + [pruning_dd, "pruning_type=%s" %  pruning_type] + optional_pr  )
    
    config_pruning = "prune=num_simulation({simulation_params}, truncate_value={trval}, abs=builder_massim({builder_params}))".format(**locals())

    config = "astar(%s, %s)" % (heuristic, config_pruning)
    return config

