#! /usr/bin/env python

from lab.parser import Parser
from collections import defaultdict
import re

eval = Parser()

eval.add_pattern('total_simulation_time', r'Done initializing simulation heuristic [(.+)s]', type=float, required=False)
eval.add_pattern('total_search_time', r'Search time: (.+)s', type=float, required=False)
eval.add_pattern('total_abstraction_time', r'Done initializing merge-and-shrink heuristic [(.+)s]', type=float, required=False)
eval.add_pattern('final_abstractions', r'Final abstractions: (d+)]', type=int, required=False)
eval.add_pattern('useless_vars', r'Useless vars: (d+)]', type=int, required=False)

def desactivation(content, props):
    regexp_desactivation = re.compile("Simulation pruning (.*)activated: (?P<pruned>(\d+)) pruned (?P<checked>(\d+)) checked (?P<inserted>(\d+)) inserted (?P<deadends>(\d+)) deadends")

    for l in content.split("\n"):
        mx = regexp_desactivation.match(l)
        if mx:
            props['pruning_desactivated'] =  1 if "desactivated" in l else 0 
            props['checked_desactivated'] = int(mx.groupdict()["checked"])
            props['pruned_desactivated'] = int(mx.groupdict()["pruned"])
            props['inserted_desactivated'] = int(mx.groupdict()["inserted"])
            props['deadends_desactivated'] = int(mx.groupdict()["deadends"])
            break
    

def dead_ops(content, props):
    regexp_dead_labels = re.compile("Dead operators due to dead labels: (?P<dead>(\d+)) / (?P<ops>(\d+)) \((?P<perc>(\d*[.\d+]*))\%\)")
    regexp_dead_orig = re.compile("Dead operators detected by storing original operators: (?P<dead>(\d+)) / (?P<ops>(\d+)) \((?P<perc>(\d*[.\d+]*))\%\)")
    
    for l in content.split("\n"):
        mx = regexp_dead_labels.match(l)
        if mx:
            props['dead_ops_by_labels'] = int(mx.groupdict()["dead"])
            props['perc_dead_ops_by_labels'] = float(mx.groupdict()["perc"])
            props['orig_ops'] = int(mx.groupdict()["ops"])
            continue
        mx = regexp_dead_orig.match(l)
        if mx:
            props['dead_ops_by_stored'] = int(mx.groupdict()["dead"])
            props['perc_dead_ops_by_stored'] = float(mx.groupdict()["perc"])
            props['orig_ops'] = int(mx.groupdict()["ops"])
            break


eval.add_function(dead_ops)

eval.add_function(desactivation)

    


eval.parse()

eval2 = Parser()

def real_search_time(content, props ):
    if props['domain'][-4:] == "-por":
        props['domain'] = props['domain'][:-4]
        props['id'] [1] = props['domain']

    if props['domain'][-1] == "-":
        props['domain'] = props['domain'][:-1]
        props['id'] [1] = props['domain']


    finished_sim = False
    regexp_mas = re.compile("Abstraction \(.*\): .*\[t=(?P<time>(\d*[.\d+]*))s\]")
    regexp_mas2 = re.compile("Atomic abstraction .*: .*\[t=(?P<time>(\d*[.\d+]*))s\]")
    regexp_sim = re.compile("Done initializing simulation heuristic \[(?P<time>(\d*[.\d+]*))s\]")
    regexp_prep2 = re.compile("Best heuristic value:.*t=(?P<time>(\d*[.\d+]*))s.*")
    regexp_prep = re.compile("Completed preprocessing: (?P<time>(\d*[.\d+]*))")
    regexp_desac = re.compile("Desactivation of .*: (?P<pruned>(\d*)) pruned (?P<inserted>(\d*)) inserted")
    regexp_prep3 = re.compile("f = .*t=(?P<time>(\d*[.\d+]*))s.*")

    time_mas = None
    time_sim = None
    time_prep = None

    for l in content.split("\n"):
        if not finished_sim and regexp_mas2.match(l):
            time_mas = float(regexp_mas2.match(l).groupdict()["time"])
        elif not finished_sim and regexp_mas.match(l):
            time_mas = float(regexp_mas.match(l).groupdict()["time"])
        elif not finished_sim and regexp_sim.match(l):
            time_sim = float(regexp_sim.match(l).groupdict()["time"])
            finished_sim = True
        elif finished_sim and regexp_prep.match(l):
            time_prep = float(regexp_prep.match(l).groupdict()["time"])
            break
        elif finished_sim and regexp_prep2.match(l):
            time_prep = float(regexp_prep2.match(l).groupdict()["time"])
            break
        elif finished_sim and regexp_prep3.match(l):
            time_prep = float(regexp_prep3.match(l).groupdict()["time"])
            break


    for l in content.split("\n"):
        if regexp_desac.match(l):
            props['desactivated'] = 1
            props['desactivated_inserted'] = regexp_desac.match(l).groupdict()["inserted"]
            props['desactivated_pruned'] = regexp_desac.match(l).groupdict()["pruned"]
            
            
    #print time_mas, time_sim, time_prep
    if time_mas:
        props['time_mas'] = time_mas

    try:
        if time_sim:
            props['time_massim'] = time_sim
            props['time_sim'] = time_sim - time_mas
    except:
        pass


    try:
        if time_prep:
            props['time_prep'] = time_prep
            props['time_bdd'] = time_prep - time_sim
    except:
        pass
            
    try: 
        if 'search_time' in props:
            if time_prep:        
                props['real_search_time'] = max(0.1, float(props['search_time']) - time_prep)
                props['total_time_minus_mas'] = max(0.1, float(props['search_time']) - time_mas)
                props['total_time_minus_massim'] = max(0.1, float(props['search_time']) - time_sim)
            else:
                props['real_search_time'] = float(props['search_time'])

        if 'real_search_time' in props and 'expansions' in props and float(props['expansions']) > 0.0:
            props['time_per_node'] = float(props['real_search_time'])/float(props['expansions'])
            props['nodes_per_second'] = float(props['expansions'])/float(props['real_search_time'])
    except: 
        pass

eval2.add_function(real_search_time)

eval2.parse()
