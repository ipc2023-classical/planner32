import configs
from configs_simulation import *

import sys

from collections import defaultdict

REVISION = "05a2f1687ac7"
SERVERS = "old_servers" 

def add_config(CONFIG_NAME, config_list, revision = REVISION, servers = SERVERS):
    config = "-".join(map(str, config_list))
    if len(config.split("-")) == 5: 
        CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_simulation_config(config), 'optimal', revision, servers))
    else:
        CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', revision, servers))

# Experiment #1: simulation type and pruning types
# merge_strategies = ["atomic", "dfp50k"]

CONFIGS = defaultdict(list)

CONFIGS["journal1-atomic"].append(configs.Config('blind', 'blind', "astar(blind())", 'optimal', '6320039e08bb', SERVERS))


for sim in ["sim", "bisim", "ldsimalt", "noopsim", "qpos-10", "qtrade-10", "qrel-10", "qual-10"]:
    add_config("journal1-atomic", ["blind", sim, "atomic", "nosh", "gen"])

add_config("journal1-atomic", ["blind", "qrel", "10", "atomic", "nosh", "gensucc"])
add_config("journal1-atomic", ["blind", "qrel", "10", "atomic", "nosh", "succ"])

add_config("journal1-atomic", ["blind", "qrel", "10", "atomic", "nosh", "gensucc", "simpletau"])
add_config("journal1-atomic", ["blind", "qrel", "10", "atomic", "nosh", "gensucc", "recurtau"])

add_config("journal1-atomic", ["blind", "qrel", "10", "dfp50k", "bissh", "gensucc", "simpletau"])
add_config("journal1-atomic", ["blind", "qrel", "10", "dfp50k", "bissh", "gensucc", "recurtau"])

for trval in [0, 1, 100, 1000]:
    add_config("journal1-atomic", ["blind", "qrel", trval, "atomic", "nosh", pruning_type])

for sh in ["nosh", "bissh"]:
    for sim in ["sim", "bisim", "ldsimalt", "noopsim", "qpos-10", "qtrade-10", "qrel-10", "qual-10"]:
        add_config("journal1-nonatomic", ["blind", sim, "dfp50k", sh, "gen"])

        add_config("journal1-nonatomic", ["blind", "qrel", "10", "dfp50k", sh, "gensucc"])
        add_config("journal1-nonatomic", ["blind", "qrel", "10", "dfp50k", sh, "succ"])




for sim in ["bisim", "ldsimalt", "qrel-10"]:
    add_config("journal2-heuristic", ["lmcut", sim, "atomic", "bissh", "gen"])
    add_config("journal2-heuristic", ["lmcut", sim, "dfp50k", "bissh", "gen"])
    add_config("journal2-heuristic", ["lmcut", sim, "atomic", "bissh", "gensucc"])
    add_config("journal2-heuristic", ["lmcut", sim, "dfp50k", "bissh", "gensucc"])
 
