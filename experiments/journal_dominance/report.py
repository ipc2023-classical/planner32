#! /usr/bin/env python

"""Solve some tasks with A* and the LM-Cut heuristic."""

import os
import platform
import subprocess
import sys

from lab.environments import OracleGridEngineEnvironment
from lab.steps import Step
from downward.suites import *
from downward.configs import *
from downward.experiment import DownwardExperiment
from downward.reports.absolute import AbsoluteReport


NAME = "all"

sys.path.append('configs')
import configs

ENV = OracleGridEngineEnvironment(queue='all.q@fai01.cs.uni-saarland.de,all.q@fai02.cs.uni-saarland.de,all.q@fai03.cs.uni-saarland.de,all.q@fai04.cs.uni-saarland.de,all.q@fai05.cs.uni-saarland.de,all.q@fai06.cs.uni-saarland.de,all.q@fai07.cs.uni-saarland.de,all.q@fai08.cs.uni-saarland.de')

REPO = '/mnt/data_server/torralba/dominance-journal/fd_simulation/'

PARSER = REPO + '/lab_parser.py'

LIMITS={'search_time': 1800,  'search_memory' : 4096}
ATTRIBUTES = ['cost', 'unsolvable', 'coverage', 'expansions', 'expansions_until_last_jump', 'evaluations', 'memory', 'last_logged_time', 'total_simulation_time', 'total_time', 'search_time', 'pruned', 'time_numeric_ldsimulation', 'action_selection_rules', 'search_restarts', 'pruned_states', 'action_selection_rules']
SUITE = suite_satisficing_strips_with_ipc14()

REP_NAME = 'report-{}'.format(NAME)

EXPPATH = 'reports/' + REP_NAME

exp = DownwardExperiment(path=EXPPATH, repo=REPO, environment=ENV, limits=LIMITS)
exp.add_suite(SUITE)

for config in configs.get_configs(NAME):
#         config.revision
        EXPPATH = '/mnt/data_server/torralba/dominance-journal/results/{}/{}/{}'.format(config.machines, config.revision, config.folder)
        exp.add_fetcher(EXPPATH, parsers=[PARSER])
        #exp.add_fetcher(EXPPATH)

# Make a report containing absolute numbers (this is the most common report).
report = os.path.join(exp.eval_dir, 'report.html')
exp.add_report(AbsoluteReport(attributes=ATTRIBUTES), outfile=report)



# Parse the commandline and show or run experiment steps.
exp()
