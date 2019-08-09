import configs
from configs_simulation import *

import sys

from collections import defaultdict

REVISION = "ded5386d5ce3"
SERVERS = "new_servers" 



# Experiment #1: simulation type and pruning types
# merge_strategies = ["atomic", "dfp50k"]

heuristic = "blind"
sh = "nosh"
trval = 10

CONFIGS = defaultdict(list)

blind_config = configs.Config('blind', 'blind', "astar(blind())", 'optimal', '6320039e08bb', SERVERS)


pruning_type = "gen"
merge_strategy = "atomic"

sh = "nosh"
CONFIG_NAME = "journal1-{}".format(merge_strategy)
CONFIGS[CONFIG_NAME].append(blind_config)
        
for sim in ["sim", "bisim", "ldsim", "ldsimalt", "noopsim"]:
    config = "-".join(map(str, [heuristic, sim, merge_strategy, sh, pruning_type]))
    CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_simulation_config(config), 'optimal', REVISION, SERVERS))

for sim in ["qpos", "qtrade", "qrel", "qual"]:
    config = "-".join(map(str, [heuristic, sim, trval, merge_strategy, sh, pruning_type]))
    CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', REVISION, SERVERS))

# Extra numeric configurations 
config = "-".join(map(str, ["blind", "qrel", "10", "atomic", "nosh", "gensucc"]))
CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', REVISION, SERVERS))

config = "-".join(map(str, ["blind", "qrel", "10", "atomic", "nosh", "succ"]))
CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', REVISION, SERVERS))

for trval in [0, 1, 2, 5, 100, 1000]:
    config = "-".join(map(str, [heuristic, sim, trval, merge_strategy, sh, pruning_type]))
    CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', REVISION, SERVERS))

