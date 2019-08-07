#! /usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import logging
import numpy
import os

from collections import defaultdict

from downward.experiment import FastDownwardExperiment
from downward.reports import PlanningReport
from downward.reports.absolute import AbsoluteReport
from downward.reports.compare import ComparativeReport
from downward.reports.scatter import ScatterPlotReport

from lab.reports import Attribute, geometric_mean, arithmetic_mean

from lab import tools

import sys
sys.path.insert(1, './scripts')

from common_setup import ReportExperiment, fix_algorithm, joint_domains
from domain_comparison import DomainComparisonReport
from personalized_table import PersonalizedTableReport, ColumnCompare
from cumulative import CumulativePlotReport


exp = ReportExperiment("report")

exp.add_fetcher('../properties/journal1-all/', postprocess_functions=[fix_algorithm, joint_domains])


all_configs=[
    'blind-ldsim-atomic-bissh-exp',
    'blind',
    'blind-ldsim-atomic-bissh-gen',
]


attributes = list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES)

exp.add_report(AbsoluteReport(attributes=attributes,filter_algorithm=all_configs))


exp.add_scatter_plot_step([
    # ("report-eval/scatter/expansions-base-vs-ldsim-atomic", ScatterPlotReport(filter_algorithm=["blind", "blind-ldsim-atomic-bissh-gen"], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"],)),
    # ("report-eval/scatter/simulationtime-ldsim-vs-qrel-atomic", ScatterPlotReport(filter_algorithm=["blind-qrel-10-atomic-bissh-gen", "blind-ldsim-atomic-bissh-gen"], attributes=["time_ldsim"], get_category=lambda run1, run2:  run1["domain_category"],)),
    ("report-eval/scatter/totalsimulationtime-ldsim-vs-qrel-atomic", ScatterPlotReport(filter_algorithm=["blind-qrel-10-atomic-bissh-gen", "blind-ldsim-atomic-bissh-gen"], attributes=["total_simulation_time"], get_category=lambda run1, run2:  run1["domain_category"],))

    
]
)
    
# , ,
#     _configpairs = [
#         (),
#         ("blind-ldsim-atomic-bissh-exp", "blind-ldsim-atomic-bissh-gen"),
#         ("blind", "blind-qrel-10-atomic-bissh-gen"),
#         ("blind-ldsim-atomic-bissh-gen", "blind-qrel-10-atomic-bissh-gen"),                              
#         ("blind", "blind-ldsim-dfp50k-bissh-gen"),
#         ("blind-ldsim-dfp50k-bissh-exp", "blind-ldsim-dfp50k-bissh-gen"),
#         ("blind", "blind-qrel-10-dfp50k-bissh-gen"),
#         ("blind-ldsim-dfp50k-bissh-gen", "blind-qrel-10-dfp50k-bissh-gen"),
#         ("blind-ldsim-dfp50k-bissh-exp", "blind-ldsim-dfp50k-bissh-gen"),
#         ("blind-ldsim-atomic-bissh-gen",  "blind-ldsim-dfp50k-bissh-gen"),
#     ],
# )


alg_list_atomic = [ 'blind',
             'blind-bisim-atomic-bissh-gen',
             'blind-sim-atomic-bissh-gen',
             'blind-noopsim-atomic-bissh-gen',
             'blind-ldsim-atomic-bissh-gen',
             "blind-qual-10-atomic-bissh-gen",
             "blind-qpos-10-atomic-bissh-gen",
             "blind-qtrade-10-atomic-bissh-gen",
             "blind-qrel-10-atomic-bissh-gen"
]
exp.add_report(
    PersonalizedTableReport(
        filter_algorithm=alg_list_atomic,
        #filter_run=(lambda x :  x["algorithm"] == "blind" and ("expansions_until_last_jump" not in x or x["expansions_until_last_jump"] < 1000)), 
        columns = [ ColumnCompare("$>$ blind", 'expansions_until_last_jump', lambda x : 'blind', lambda x, y : x < y),
                    ColumnCompare("$>$ blind x 2", 'expansions_until_last_jump', lambda x : 'blind', lambda x, y : x*2 < y),
                    ColumnCompare("$>$ blind x 10", 'expansions_until_last_jump', lambda x : 'blind', lambda x, y : x*10 < y),
                    ColumnCompare("$>$ -1", 'expansions_until_last_jump', lambda x : alg_list_atomic[alg_list_atomic.index(x) - 1 ], lambda x, y : x < y),
                    ColumnCompare("$>$ -1 x2", 'expansions_until_last_jump', lambda x : alg_list_atomic[alg_list_atomic.index(x) - 1 ], lambda x, y : x*2 < y),
                    ColumnCompare("$>$ -1 x10", 'expansions_until_last_jump', lambda x : alg_list_atomic[alg_list_atomic.index(x) - 1 ], lambda x, y : x*10 < y),

        ],
        algo_to_print= {
            'blind-bisim-atomic-bissh-gen':    'bisim',
            'blind-sim-atomic-bissh-gen':      'sim',
            'blind-noopsim-atomic-bissh-gen':  'noopsim', 
            'blind-ldsim-atomic-bissh-gen':    'ldsim',
            "blind-qual-10-atomic-bissh-gen":  "qual10",
            "blind-qpos-10-atomic-bissh-gen":  "qpos10",
            "blind-qrel-10-atomic-bissh-gen":  "qrel10",
            "blind-qtrade-10-atomic-bissh-gen": "qtrade10"
        },
        format='tex'
    ),
    outfile=os.path.join(exp.eval_dir, 'comparison-atomic-expansions_until_last_jump.tex'),
)




alg_list_dfp = [ 'blind',
             'blind-bisim-dfp50k-bissh-gen',
             'blind-sim-dfp50k-bissh-gen',
             'blind-noopsim-dfp50k-bissh-gen', 
             'blind-ldsim-dfp50k-bissh-gen',
             "blind-qual-10-dfp50k-bissh-gen",
             "blind-qpos-10-dfp50k-bissh-gen",
             "blind-qtrade-10-dfp50k-bissh-gen",
             "blind-qrel-10-dfp50k-bissh-gen"
]
exp.add_report(
    PersonalizedTableReport(
        filter_algorithm=alg_list_dfp,
#        filter_run=(lambda x :  x["algorithm"] == "blind" and ("expansions_until_last_jump" not in x or x["expansions_until_last_jump"] < 1000)), 
        columns = [ ColumnCompare("$>$ bisim", 'expansions_until_last_jump', lambda x : "blind-bisim-dfp50k-bissh-gen", lambda x, y : x < y),
                    ColumnCompare("$>$ bisim x 2", 'expansions_until_last_jump', lambda x : "blind-bisim-dfp50k-bissh-gen", lambda x, y : x*2 < y),
                    ColumnCompare("$>$ bisim x 10", 'expansions_until_last_jump', lambda x : "blind-bisim-dfp50k-bissh-gen", lambda x, y : x*10 < y, _printList=True),
                    ColumnCompare("$>$ -1", 'expansions_until_last_jump', lambda x : alg_list_dfp[alg_list_dfp.index(x) - 1 ], lambda x, y : x < y),
                    ColumnCompare("$>$ -1 x2", 'expansions_until_last_jump', lambda x : alg_list_dfp[alg_list_dfp.index(x) - 1 ], lambda x, y : x*2 < y),
                    ColumnCompare("$>$ -1 x10", 'expansions_until_last_jump', lambda x : alg_list_dfp[alg_list_dfp.index(x) - 1 ], lambda x, y : x*10 < y),

        ],
        algo_to_print= {
            'blind-bisim-dfp50k-bissh-gen':    'bisim',
            'blind-sim-dfp50k-bissh-gen':      'sim',
            'blind-noopsim-dfp50k-bissh-gen':  'noopsim', 
            'blind-ldsim-dfp50k-bissh-gen':    'ldsim',
            "blind-qual-10-dfp50k-bissh-gen":  "qual10",
            "blind-qpos-10-dfp50k-bissh-gen":  "qpos10",
            "blind-qrel-10-dfp50k-bissh-gen":  "qrel10",
            "blind-qtrade-10-dfp50k-bissh-gen": "qtrade10"
        },
        format='tex'
    ),
    outfile=os.path.join(exp.eval_dir, 'comparison-dfp50k-expansions_until_last_jump.tex'),
)


exp.add_cumulative_plot_step('search_time', ['blind', 'blind-bisim-dfp50k-bissh-gen'])

exp.run_steps()

