#! /usr/bin/env python

import os
import subprocess

import suites

from lab.steps import Step
from downward.checkouts import Translator, Preprocessor, Planner
from downward.experiment import DownwardExperiment
from downward.reports.absolute import AbsoluteReport

from lab.environments import OracleGridEngineEnvironment

REVISION = 'd559b7e49375'

EXPPATH = '/mnt/data_server/torralba/dominance-journal/results/new_servers/' + REVISION  + '/blind/'

REPO = '/mnt/data_server/torralba/dominance-journal/fd_simulation/'

ENV = OracleGridEngineEnvironment(queue='all.q@fai11.cs.uni-saarland.de,all.q@fai12.cs.uni-saarland.de,all.q@fai13.cs.uni-saarland.de,all.q@fai14.cs.uni-saarland.de,all.q@fai14.cs.uni-saarland.de')

ATTRIBUTES = ['search_time', 'reopened', 'memory', 'evaluations', 'total_time', 'expansions', 'error', 'generated', 'initial_h_values', 'coverage', 'plan_length', 'max_number_dups', 'cost']

LIMITS={'search_time': 1800,  'search_memory' : 4096}

COMBINATIONS = [(Translator(repo=REPO), Preprocessor(repo=REPO), Planner(repo=REPO, rev=REVISION))]

exp = DownwardExperiment(path=EXPPATH, repo=REPO, environment=ENV, combinations=COMBINATIONS, limits=LIMITS, cache_dir='/mnt/data_server/torralba/dominance-journal/lab-data/')

exp.add_search_parser(REPO + '/lab_parser.py')

exp.add_config('blind' + REVISION, ['--search', 'astar(blind())'])

exp.add_suite(suites.suite_optimal_strips(), benchmark_dir='/mnt/data_server/torralba/downward-benchmarks/')

def remove_work_tag(run):
    """Remove "WORK-" from the configs."""
    config = run['config']
    config = config[5:] if config.startswith('WORK-') else config
    # Shorten long config name.
    config = config.replace('downward-', '')
    run['config'] = config
    return run

# Make a report containing absolute numbers (this is the most common report).
exp.add_report(AbsoluteReport(attributes=ATTRIBUTES), '/mnt/data_server/torralba/dominance-journal/reports/blind-rev=' + REVISION)

# Parse the commandline and show or run experiment steps.
exp()
