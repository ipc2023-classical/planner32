import configs
from configs_simulation import *

import sys

from collections import defaultdict

REVISION_OLD = "05a2f1687ac7"
REVISION = "ec11c019d0fc0f4a1ce8b3dd78e57afb2ac356f6"
REVISION_TIE_BREAKING = "284e65fbd1240c26db33c75773df4d85c3e60885"

SERVERS = "old_servers" 

def add_config(CONFIG_NAME, config_list, revision = REVISION, servers = SERVERS):
    config = "-".join(map(str, config_list))
    if len(config.split("-")) == 5: 
        CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_simulation_config(config), 'optimal', revision, servers))
    else:
        CONFIGS[CONFIG_NAME].append(configs.Config(config, config, get_numeric_simulation_config(config), 'optimal', revision, servers))

def add_config_tie_breaking(CONFIG_NAME, config_list, revision = REVISION, servers = SERVERS):
    config = "-".join(map(str, config_list))
    if len(config.split("-")) == 5: 
        CONFIGS[CONFIG_NAME].append(configs.Config(config + "-tie" , config+ "-tie", get_simulation_config(config, tie_breaking_by_g=True), 'optimal', revision, servers))
    else:
        CONFIGS[CONFIG_NAME].append(configs.Config(config + "-tie", config + "-tie", get_numeric_simulation_config(config, tie_breaking_by_g=True), 'optimal', revision, servers))

# Experiment #1: simulation type and pruning types
# merge_strategies = ["atomic", "dfp50k"]

CONFIGS = defaultdict(list)

CONFIGS["journal1-atomic"].append(configs.Config('blind', 'blind', "astar(blind())", 'optimal', '6320039e08bb', SERVERS))


for sim in ["sim", "bisim", "ldsimalt", "noopsim"]:
    add_config("journal1-atomic", ["blind", sim, "atomic", "nosh", "gen"], REVISION_OLD)

for sim in ["qpos-10", "qtrade-10", "qrel-10", "qual-10"]:
    add_config("journal1-atomic", ["blind", sim, "atomic", "nosh", "gen"])

add_config("journal1-atomic", ["blind", "qrel-10", "atomic", "nosh", "gensucc"])
add_config("journal1-atomic", ["blind", "qrel-10", "atomic", "nosh", "succ"])



# add_config("journal1-atomic", ["blind", "qrel-10", "atomic", "nosh", "gen", "fulltau"])
# add_config("journal1-atomic", ["blind", "qrel-10", "atomic", "nosh", "gen", "recurtau"])

# add_config("journal1-nonatomic", ["blind", "qrel-10", "dfp50k", "bissh", "gen", "fulltau"])
# add_config("journal1-nonatomic", ["blind", "qrel-10", "dfp50k", "bissh", "gen", "recurtau"])

# add_config("journal1-atomic", ["blind", "qrel-10", "atomic", "nosh", "gensucc", "fulltau"])
# add_config("journal1-atomic", ["blind", "qrel-10", "atomic", "nosh", "gensucc", "recurtau"])

# add_config("journal1-nonatomic", ["blind", "qrel-10", "dfp50k", "bissh", "gensucc", "fulltau"])
# add_config("journal1-nonatomic", ["blind", "qrel-10", "dfp50k", "bissh", "gensucc", "recurtau"])

for trval in [0, 1, 100, 1000]:
    add_config("journal1-atomic", ["blind", "qrel", trval, "atomic", "nosh", "gen"])

for sh in ["nosh", "bissh"]:
    for sim in ["sim", "bisim", "ldsimalt", "noopsim"]:
        add_config("journal1-nonatomic", ["blind", sim, "dfp50k", sh, "gen"], REVISION_OLD)

    for sim in ["qpos-10", "qtrade-10", "qrel-10", "qual-10"]:
        add_config("journal1-nonatomic", ["blind", sim, "dfp50k", sh, "gen"])


    add_config("journal1-nonatomic", ["blind", "qrel-10", "dfp50k", sh, "gensucc"])
    add_config("journal1-nonatomic", ["blind", "qrel-10", "dfp50k", sh, "succ"])

for sim in ["bisim", "ldsimalt", "qrel-10"]:
    add_config("journal2-heuristic", ["lmcut", sim, "atomic", "bissh", "gen"])
    add_config("journal2-heuristic", ["lmcut", sim, "dfp50k", "bissh", "gen"])

add_config("journal2-heuristic", ["lmcut", "qrel-10", "atomic", "bissh", "gensucc"])
add_config("journal2-heuristic", ["lmcut", "qrel-10", "dfp50k", "bissh", "gensucc"])


add_config("tie-breaking-analysis", ["hmax", "ldsimalt", "atomic", "nosh", "gen"], revision = REVISION_TIE_BREAKING)
add_config("tie-breaking-analysis", ["hmax", "ldsimalt", "dfp50k", "bissh", "gen"], revision = REVISION_TIE_BREAKING)
add_config_tie_breaking("tie-breaking-analysis", ["lmcut", "ldsimalt", "atomic", "nosh", "gen"], revision = REVISION_TIE_BREAKING)
add_config_tie_breaking("tie-breaking-analysis", ["lmcut", "ldsimalt", "dfp50k", "bissh", "gen"], revision = REVISION_TIE_BREAKING)
