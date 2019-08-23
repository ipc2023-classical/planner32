#! /usr/bin/env python

from lab.parser import Parser
from collections import defaultdict
import re

eval = Parser()


regexps = [re.compile("Compute LDSim on (?P<lts_num>(\d+)) LTSs. Total size: (?P<lts_total_size>(\d+)) Total trsize: (?P<lts_total_trsize>(\d+)) Max size: (?P<lts_max_size>(\d+)) Max trsize: (?P<lts_max_trsize>(\d+))"),
           re.compile(".*Init LDSim in (?P<time_init_ldsim>(.*))s:.*"),
           re.compile("LDSim computed (?P<time_ldsim>(.*))"),
           re.compile("Dead operators due to dead labels: (?P<dead_ops_by_labels>(\d+)) / (?P<orig_ops>(\d+)) \((?P<perc_dead_ops_by_labels>(\d*[.\d+]*))\%\)"), 
           re.compile("Dead operators detected by storing original operators: (?P<dead_ops_by_stored>(\d+)) / (?P<orig_ops>(\d+)) \((?P<perc_dead_ops_by_stored>(\d*[.\d+]*))\%\)"), 
           re.compile("Simulation pruning (?P<pruning_desactivated>(.*)): (?P<pruned_desactivated>(\d+)) pruned (?P<checked_desactivated>(\d+)) checked (?P<inserted_desactivated>(\d+)) inserted (?P<deadends_desactivated>(\d+)) deadends"),
           re.compile("Numeric LDSim computed(?P<time_ldsim> (.*))"),
           re.compile("Numeric LDSim outer iterations: (?P<outer_iterations_numeric_ldsimulation>(.*))"),
           re.compile("Numeric LDSim inner iterations: (?P<inner_iterations_numeric_ldsimulation>(.*))"),
           re.compile("First node pruned after checking (?P<dom_checked_before_first_pruned>(\d+)) and inserting (?P<dom_inserted_before_first_pruned>(\d+))"),
           re.compile("Done initializing simulation heuristic \[(?P<total_simulation_time>(.*))s\]"),
           re.compile('Done initializing merge-and-shrink heuristic \[(?P<total_abstraction_time>(.*))s\]'),
           re.compile('Final abstractions: (?P<final_abstractions>(\d+))'),
           re.compile('Useless vars: (?P<useless_vars>(\d+))'),
           re.compile('Total Simulations: (?P<total_simulations>(\d+))'),
           re.compile('Only Simulations: (?P<only_simulations>(\d+))'),
           re.compile('Similarity equivalences: (?P<similarity_equivalences>(\d+))'),
           re.compile('Completed preprocessing: (?P<time_completed_preprocessing>(.*))'),
]

type_atr = {'dead_ops_by_labels' : int, 'perc_dead_ops_by_labels' : float, 'orig_ops' : int, 
            'dead_ops_by_stored' : int, 'perc_dead_ops_by_stored' : float, 
            'lts_num' : int, 'lts_total_size' : int,  'lts_max_size' : int, 'lts_total_trsize' : int, 'lts_max_trsize' : int, 
            'pruning_desactivated' : (lambda x : 1 if "desactivated" == x else 0 ), 
            'pruned_desactivated' : int, 'checked_desactivated' : int, 'inserted_desactivated' : int, 'deadends_desactivated' : int,
            'time_init_ldsim' : lambda x : max(0.01, float(x)),
            'outer_iterations_numeric_ldsimulation' : int,
            'inner_iterations_numeric_ldsimulation' : int,
            "total_simulation_time" : lambda x : max(0.01, float(x)),
            "total_abstraction_time" : lambda x : max(0.01, float(x)), "final_abstractions" : int,
            "useless_vars" : int, "total_simulations" : int, "only_simulations" : int, "similarity_equivalences" : int,
            'dom_inserted_before_first_pruned' : int,  'dom_checked_before_first_pruned' : int,
            'time_ldsim' : lambda x : max(0.01, float(x)),
            'time_completed_preprocessing' : lambda x : max(0.01, float(x)),
        }

def parse_regexps (content, props):
    for l in content.split("\n"):
        for reg in regexps:
            mx = reg.match(l)
            if mx:
                data = mx.groupdict()
                for item in data:
                    props[item] = type_atr[item](data[item])
                break

    props["did_prune"] = 1 if "dom_checked_before_first_pruned" in props else 0



def parse_numeric_dominance (content, props):
    check = False
    min_val = 100000000
    max_val = -100000000

    for l in content.split("\n"):
        if check: 
            if l == "Init partitions" or l.startswith("Completed preprocessing"):
                if min_val == 100000000 and max_val == -100000000:
                    min_val = -100000000
                    max_val = 100000000
                props['min_negative_dominance'] = min_val 
                props['max_positive_dominance'] = max_val
                props["has_positive_dominance"] = 1 if (max_val > 0) else 0
                props["has_negative_dominance"] = 1 if (min_val < 0) and (min_val > -100000000) else 0
                return
            if ":" in l and not "infinity" in l:
                val = l.split(":")[0]
                if "(" in val:
                    val = val.split("(")[0]
                min_val = min(min_val, int(val))
                max_val = max(max_val, int(val))
        elif l == "------": 
            check = True


        
            
def fix_error (content, props):
    if not props["error"].startswith("unexplained"):
        return
    for l in content.split("\n"):
        if l.startswith("Peak memory: Failed to allocate memory. Released memory buffer.") or l.startswith("CUDD: out of memory allocating"):
            props["error"] = "out-of-memory"
            return


def set_algorithm_prop (content, props):
    if "algorithm" not in props:
        props["algorithm"] = props["experiment_name"]

# def real_search_time(content, props ):
#     if props['domain'][-4:] == "-por":
#         props['domain'] = props['domain'][:-4]
#         props['id'] [1] = props['domain']

#     if props['domain'][-1] == "-":
#         props['domain'] = props['domain'][:-1]
#         props['id'] [1] = props['domain']


#     finished_sim = False
#     regexp_mas = re.compile("Abstraction \(.*\): .*\[t=(?P<time>(\d*[.\d+]*))s\]")
#     regexp_mas2 = re.compile("Atomic abstraction .*: .*\[t=(?P<time>(\d*[.\d+]*))s\]")
#     regexp_sim = re.compile("Done initializing simulation heuristic \[(?P<time>(\d*[.\d+]*))s\]")
#     regexp_prep2 = re.compile("Best heuristic value:.*t=(?P<time>(\d*[.\d+]*))s.*")
#     regexp_prep = re.compile("Completed preprocessing: (?P<time>(\d*[.\d+]*))")
#     regexp_desac = re.compile("Desactivation of .*: (?P<pruned>(\d*)) pruned (?P<inserted>(\d*)) inserted")
#     regexp_prep3 = re.compile("f = .*t=(?P<time>(\d*[.\d+]*))s.*")

#     time_mas = None
#     time_sim = None
#     time_prep = None

#     for l in content.split("\n"):
#         if not finished_sim and regexp_mas2.match(l):
#             time_mas = float(regexp_mas2.match(l).groupdict()["time"])
#         elif not finished_sim and regexp_mas.match(l):
#             time_mas = float(regexp_mas.match(l).groupdict()["time"])
#         elif not finished_sim and regexp_sim.match(l):
#             time_sim = float(regexp_sim.match(l).groupdict()["time"])
#             finished_sim = True
#         elif finished_sim and regexp_prep.match(l):
#             time_prep = float(regexp_prep.match(l).groupdict()["time"])
#             break
#         elif finished_sim and regexp_prep2.match(l):
#             time_prep = float(regexp_prep2.match(l).groupdict()["time"])
#             break
#         elif finished_sim and regexp_prep3.match(l):
#             time_prep = float(regexp_prep3.match(l).groupdict()["time"])
#             break


#     for l in content.split("\n"):
#         if regexp_desac.match(l):
#             props['desactivated'] = 1
#             props['desactivated_inserted'] = regexp_desac.match(l).groupdict()["inserted"]
#             props['desactivated_pruned'] = regexp_desac.match(l).groupdict()["pruned"]
            
            
#     #print time_mas, time_sim, time_prep
#     if time_mas:
#         props['time_mas'] = time_mas

#     try:
#         if time_sim:
#             props['time_massim'] = time_sim
#             props['time_sim'] = time_sim - time_mas
#     except:
#         pass


#     try:
#         if time_prep:
#             props['time_prep'] = time_prep
#             props['time_bdd'] = time_prep - time_sim
#     except:
#         pass
            
#     try: 
#         if 'search_time' in props:
#             if time_prep:        
#                 props['real_search_time'] = max(0.1, float(props['search_time']) - time_prep)
#                 props['total_time_minus_mas'] = max(0.1, float(props['search_time']) - time_mas)
#                 props['total_time_minus_massim'] = max(0.1, float(props['search_time']) - time_sim)
#             else:
#                 props['real_search_time'] = float(props['search_time'])

#         if 'real_search_time' in props and 'expansions' in props and float(props['expansions']) > 0.0:
#             props['time_per_node'] = float(props['real_search_time'])/float(props['expansions'])
#             props['nodes_per_second'] = float(props['expansions'])/float(props['real_search_time'])
#     except: 
#         pass

# eval.add_function(real_search_time)

eval.add_function(parse_regexps)

#eval.add_function(desactivation)

eval.add_function(parse_numeric_dominance)

eval.add_function(fix_error)
eval.add_function(set_algorithm_prop)


eval.parse()
