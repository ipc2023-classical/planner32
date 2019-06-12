import configs
from configs_simulation import *

import sys


REVISION = "b8314abb922f"
SERVERS = "new_servers" 

# Experiment #1: simulation type and pruning types
merge_strategies = ["atomic", "dfp50k"] #["atomic", "dfp10k", "dfp50k", "dfp100k", "dfp200k", "dfp100states", "dfp1kstates", "dfp10kstates"]

heuristic = "blind"
pruning_type = "gen"
sh = "bissh"
trval = 10

CONFIGS = []
for mer in merge_strategies: 
    for sim in ["sim", "bisim", "ldsim", "noopsim"]:
        config = "-".join(map(str, [heuristic, sim, mer, sh, pruning_type]))
        CONFIGS.append(configs.Config(config, config, get_simulation_config(config), 'optimal', REVISION, SERVERS))

    for opt in [["nooptau"], []]:
        for sim in ["qpos", "qtrade", "qrel", "qual"]:
            config = "-".join(map(str, [heuristic, sim, trval, mer, sh, pruning_type] + opt))
            CONFIGS.append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', REVISION, SERVERS))



        


#prune=num_simulation(prune_successors=False, prune_dominated_by_parent=False, prune_dominated_by_closed=True, prune_dominated_by_open=False, insert_dominated=true, use_quantified_dominance=False, exit_after_preprocessing=True, use_single_bdd=True, trade_off_dominance=False, only_positive_dominance=False, use_ADDs=False, compute_tau_labels_with_noop_dominance=False, truncate_value=100000, min_desactivation_ratio=0, min_insertions=infinity,  abs=builder_massim(compute_final_simulation=false, switch_off_label_dominance=infinity, merge_strategy=merge_dfp(), limit_transitions_merge=100000, limit_merge=infinity))'

# remove_spurious=true,



# 'prune=num_simulation(insert_dominated=true, use_quantified_dominance=False,
# exit_after_preprocessing=False, use_single_bdd=False, trade_off_dominance=False,
# only_positive_dominance=False, use_ADDs=False,
# compute_tau_labels_with_noop_dominance=False,
# compute_tau_labels_as_self_loops_everywhere=True, truncate_value=10,








# abs=builder_massim(compute_final_simulation=false, switch_off_label_dominance=infinity, merge_strategy=merge_dfp(), limit_transitions_merge=10000, min_limit_merge=0, limit_merge=infinity, , forbid_lr=false, shrink_strategy=shrink_bisimulation(max_states=infinity, greedy=false, threshold=1, aggregate_goals=true))))'])
