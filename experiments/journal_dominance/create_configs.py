#!/usr/bin/env python
# encoding: utf-8

import sys
sys.path.append('configs')

import errno
import os

from subprocess import call

import configs

def get_suite(suite_type):
    if suite_type == 'optimal':
        return 'suites.suite_optimal_strips()'
    else:
        return 'suites.suite_satisficing_strips()'

def get_queue(machines):
    if (machines == 'old_servers'):
        return "OracleGridEngineEnvironment(queue='all.q@fai01.cs.uni-saarland.de,all.q@fai02.cs.uni-saarland.de,all.q@fai03.cs.uni-saarland.de,all.q@fai04.cs.uni-saarland.de,all.q@fai05.cs.uni-saarland.de,all.q@fai06.cs.uni-saarland.de,all.q@fai07.cs.uni-saarland.de,all.q@fai08.cs.uni-saarland.de')"
    elif (machines == 'new_servers'):
        return "OracleGridEngineEnvironment(queue='all.q@fai11.cs.uni-saarland.de,all.q@fai12.cs.uni-saarland.de,all.q@fai13.cs.uni-saarland.de,all.q@fai14.cs.uni-saarland.de,all.q@fai14.cs.uni-saarland.de')"
    elif (machines == 'all_servers'):
        return "OracleGridEngineEnvironment(queue='all.q@fai01.cs.uni-saarland.de,all.q@fai02.cs.uni-saarland.de,all.q@fai03.cs.uni-saarland.de,all.q@fai04.cs.uni-saarland.de,all.q@fai05.cs.uni-saarland.de,all.q@fai06.cs.uni-saarland.de,all.q@fai07.cs.uni-saarland.de,all.q@fai08.cs.uni-saarland.de,all.q@fai11.cs.uni-saarland.de,all.q@fai12.cs.uni-saarland.de,all.q@fai13.cs.uni-saarland.de,all.q@fai14.cs.uni-saarland.de')"
    else:
        print("ERROR: please specify on which machines the experiments should be running!", machines)
        exit(1)


def get_script(config):
    QUEUES = get_queue(config.machines)
    SUITE = get_suite(config.type)

    return "#! /usr/bin/env python\n\
\n\
import os\n\
import subprocess\n\
\n\
import suites\n\
\n\
from lab.steps import Step\n\
from downward.checkouts import Translator, Preprocessor, Planner\n\
from downward.experiment import DownwardExperiment\n\
from downward.reports.absolute import AbsoluteReport\n\
\n\
from lab.environments import OracleGridEngineEnvironment\n\
\n\
REVISION = '" + config.revision + "'\n\
\n\
EXPPATH = '/mnt/data_server/torralba/dominance-journal/results/" + config.machines + "/' + REVISION  + '/" + config.folder + "/'\n\
\n\
REPO = '/mnt/data_server/torralba/dominance-journal/fd_simulation/'\n\
\n\
ENV = " + QUEUES + "\n\
\n\
ATTRIBUTES = ['search_time', 'reopened', 'memory', 'evaluations', 'total_time', 'expansions', 'error', 'generated', 'initial_h_values', 'coverage', 'plan_length', 'max_number_dups', 'cost']\n\
\n\
LIMITS={'search_time': 1800,  'search_memory' : 4096}\n\
\n\
COMBINATIONS = [(Translator(repo=REPO), Preprocessor(repo=REPO), Planner(repo=REPO, rev=REVISION))]\n\
\n\
exp = DownwardExperiment(path=EXPPATH, repo=REPO, environment=ENV, combinations=COMBINATIONS, limits=LIMITS, cache_dir='/mnt/data_server/torralba/dominance-journal/lab-data/')\n\
\n\
exp.add_search_parser(REPO + '/lab_parser.py')\n\
\n\
exp.add_config('" + config.nick + "' + REVISION, " + str(config.config) + ")\n\
\n\
exp.add_suite(" + SUITE + ", benchmark_dir='/mnt/data_server/torralba/downward-benchmarks/')\n\
\n\
def remove_work_tag(run):\n\
    \"\"\"Remove \"WORK-\" from the configs.\"\"\"\n\
    config = run['config']\n\
    config = config[5:] if config.startswith('WORK-') else config\n\
    # Shorten long config name.\n\
    config = config.replace('downward-', '')\n\
    run['config'] = config\n\
    return run\n\
\n\
# Make a report containing absolute numbers (this is the most common report).\n\
exp.add_report(AbsoluteReport(attributes=ATTRIBUTES), '/mnt/data_server/torralba/dominance-journal/reports/" + config.folder + "-rev=' + REVISION)\n\
\n\
# Parse the commandline and show or run experiment steps.\n\
exp()"


if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print ("please specify experiment and output folder")
        print ("\n  ".join(configs.CONFIGS.keys()))
        
        exit()
    
    experiment = sys.argv[1]
    folder = sys.argv[2]
  
    if (not folder.endswith('/')):
        folder = folder + '/'
    print "writing lab scripts to this subfolder: " + folder
    try:
        os.makedirs(folder)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(folder):
            pass
        else: raise
        
    for config in configs.CONFIGS[experiment]:
        data = get_script(config)

        EXPPATH = '/mnt/data_server/torralba/dominance-journal/results/{}/{}/{}'.format(config.machines, config.revision, config.folder)
        if os.path.isdir(EXPPATH):
            continue
        name = folder + config.folder + ".py"
        with open(name, "w") as file:
            file.write(data)
            call(["chmod", "+x", name])
            print "successfully written " + name 
