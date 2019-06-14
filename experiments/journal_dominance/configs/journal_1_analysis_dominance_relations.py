import configs
from configs_simulation import *

import sys

from collections import defaultdict

REVISION = "6320039e08bb"
SERVERS = "new_servers" 

# Experiment #1: simulation type and pruning types
merge_strategies = ["atomic", "dfp50k"] #["atomic", "dfp10k", "dfp50k", "dfp100k", "dfp200k", "dfp100states", "dfp1kstates", "dfp10kstates"]

heuristic = "blind"
sh = "bissh"
trval = 10

CONFIGS = defaultdict(list)

blind_config = configs.Config('blind', 'blind', "astar(blind())", 'optimal', '6320039e08bb', SERVERS)


for pruning_type in ["exp", "gen"]:
    for mer in merge_strategies:
        CONFIG_NAME = "journal1-{}-{}".format(mer, pruning_type)
        CONFIGS[CONFIG_NAME].append(blind_config)
    
        for sim in ["sim", "bisim", "ldsim", "noopsim"]:
            config = "-".join(map(str, [heuristic, sim, mer, sh, pruning_type]))
            CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_simulation_config(config), 'optimal', REVISION, SERVERS))

        for sim in ["qpos", "qtrade", "qrel", "qual"]:
            config = "-".join(map(str, [heuristic, sim, trval, mer, sh, pruning_type]))
            CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', REVISION, SERVERS))
