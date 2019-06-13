import configs
from configs_simulation import *

import sys

from collections import defaultdict

REVISION = "5a3b399a5ef0"
SERVERS = "new_servers" 

# Experiment #1: simulation type and pruning types
merge_strategies = ["atomic", "dfp50k"] #["atomic", "dfp10k", "dfp50k", "dfp100k", "dfp200k", "dfp100states", "dfp1kstates", "dfp10kstates"]

heuristic = "blind"
pruning_type = "exp"
sh = "bissh"
trval = 10

CONFIGS = defaultdict(list)
for mer in merge_strategies:
    CONFIGS["journal1-{}".format(mer)].append(configs.Config('blind', 'blind', "astar(blind())", 'optimal', 'd559b7e49375', SERVERS))
    
    for sim in ["sim", "bisim", "ldsim", "noopsim"]:
        config = "-".join(map(str, [heuristic, sim, mer, sh, pruning_type]))
        CONFIGS["journal1-{}".format(mer)].append(configs.Config(config, config, get_simulation_config(config), 'optimal', REVISION, SERVERS))

    for sim in ["qpos", "qtrade", "qrel", "qual"]:
        config = "-".join(map(str, [heuristic, sim, trval, mer, sh, pruning_type]))
        CONFIGS["journal1-{}".format(mer)].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', REVISION, SERVERS))



        

