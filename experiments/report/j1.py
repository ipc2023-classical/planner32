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

def invert_min_negative_dominance(props):
    for p in props:
        if "min_negative_dominance" in props[p]:
            props[p] ["min_negative_dominance_inverted"] = -props[p] ["min_negative_dominance"]

exp.add_fetcher('../properties/journal1-all/', postprocess_functions=[fix_algorithm, joint_domains, invert_min_negative_dominance])


exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations"],filter_algorithm=[
    'blind',
    #'blind-ldsim-atomic-bissh-exp',
    'blind-ldsim-atomic-nosh-gen',
    'blind-ldsimalt-atomic-nosh-gen',
#    'blind-qrel-10-atomic-bissh-gen',
]))


exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations"],filter_algorithm=[
    'blind-ldsimalt-atomic-nosh-gen',
    'blind-qual-10-atomic-nosh-gen',
]), outfile='report-ldsimalt-vs-qual.html')


exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations", "min_negative_dominance", "max_positive_dominance", "has_positive_dominance", "has_negative_dominance"],filter_algorithm=[
    'blind-qrel-0-atomic-nosh-gen',
    'blind-qrel-1-atomic-nosh-gen',
    'blind-qrel-2-atomic-nosh-gen',
    'blind-qrel-5-atomic-nosh-gen',
    'blind-qrel-10-atomic-nosh-gen',
    'blind-qrel-100-atomic-nosh-gen',
    'blind-qrel-1000-atomic-nosh-gen',
]), outfile='report-k.html')


# exp.add_scatter_plot_step([
#     ("report-eval/scatter/expansions-base-vs-ldsimalt-atomic", ScatterPlotReport(filter_algorithm=["blind", "blind-ldsimalt-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     ("report-eval/scatter/expansions-ldsimalt-vs-qrel10-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-qrel-10-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"],)), 
#     # ("report-eval/scatter/expansions-base-vs-ldsim-atomic", ScatterPlotReport(filter_algorithm=["blind", "blind-ldsim-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     # ("report-eval/scatter/simulationtime-ldsim-vs-qrel-atomic", ScatterPlotReport(filter_algorithm=["blind-qrel-10-atomic-nosh-gen", "blind-ldsim-atomic-nosh-gen"], attributes=["time_ldsim"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     # ("report-eval/scatter/preprocessingtime-ldsim-vs-qrel-atomic", ScatterPlotReport(filter_algorithm=["blind-ldsim-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], attributes=["time_completed_preprocessing"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     # ("report-eval/scatter/preprocessingtime-noopsim-vs-qrel-atomic", ScatterPlotReport(filter_algorithm=['blind-noopsim-atomic-nosh-gen', "blind-qrel-10-atomic-nosh-gen"], attributes=["time_completed_preprocessing"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     # ("report-eval/scatter/preprocessingtime-noopsim-vs-ldsim-atomic", ScatterPlotReport(filter_algorithm=['blind-noopsim-atomic-nosh-gen', "blind-ldsim-atomic-nosh-gen"], attributes=["time_completed_preprocessing"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     # ("report-eval/scatter/preprocessingtime-ldsimalt-vs-ldsim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-ldsim-atomic-nosh-gen"], attributes=["time_completed_preprocessing"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     # ("report-eval/scatter/ldsimtime-ldsimalt-vs-ldsim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-ldsim-atomic-nosh-gen"], attributes=["time_ldsim"], get_category=lambda run1, run2:  run1["domain_category"],)),
#     # ("report-eval/scatter/expansions-ldsimalt-vs-ldsim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-ldsim-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"],)), 
#     # ("report-eval/scatter/expansions-ldsimalt-vs-noopsim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-noopsim-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"],)), 
#     # ("report-eval/scatter/expansions-ldsimalt-vs-sim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-sim-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"],)), 
#     # ("report-eval/scatter/ldsimtime-ldsimalt-vs-noopsim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-noopsim-atomic-nosh-gen"], attributes=["time_ldsim"], get_category=lambda run1, run2:  run1["domain_category"],)), 
#     # ("report-eval/scatter/ldsimtime-ldsimalt-vs-sim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-sim-atomic-nosh-gen"], attributes=["time_ldsim"], get_category=lambda run1, run2:  run1["domain_category"],)), 
#     # # ("report-eval/scatter/preprocessingtimebylabel-ldsim-atomic", ScatterPlotReport(filter_algorithm=["blind-ldsim-atomic-nosh-gen"], attributes=["labels", "time_completed_preprocessing"], get_category=lambda run1, run2:  run1["domain_category"],))
#     # ("report-eval/scatter/ldsimtime-ldsimalt-vs-sim-atomic", ScatterPlotReport(filter_algorithm=['blind-ldsimalt-atomic-nosh-gen', "blind-sim-atomic-nosh-gen"], attributes=["time_ldsim"], get_category=lambda run1, run2:  run1["domain_category"],)), 
# ]
# )



scatter_plots_k_expansions = [("report-eval/scatter/expansions-qrel0-vs-qrel{}-atomic".format(i), ScatterPlotReport(filter_algorithm=[ 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i)], attributes=["expansions_until_last_jump"], get_category=lambda run1, run2:  run1["domain_category"])) for i in [1, 2, 5, 10, 100, 1000]]

scatter_plots_k_time_preprocessing = [("report-eval/scatter/time-preprocessing-qrel0-vs-qrel{}-atomic".format(i), ScatterPlotReport(filter_algorithm=[ 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i)], attributes=["time_completed_preprocessing"], get_category=lambda run1, run2:  run1["domain_category"])) for i in [1, 2, 5, 10, 100, 1000]] 


scatter_plots_k_time_ldsim = [("report-eval/scatter/time-ldsim-qrel0-vs-qrel{}-atomic".format(i), ScatterPlotReport(filter_algorithm=[ 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i)], attributes=["time_ldsim"], get_category=lambda run1, run2:  run1["domain_category"])) for i in [1, 2, 5, 10, 100, 1000]] 


scatter_plots_k_min_dom = [("report-eval/scatter/min-dom-qrel0-vs-qrel{}-atomic".format(i), ScatterPlotReport(filter_algorithm=[ 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i)], attributes=["min_negative_dominance_inverted"], get_category=lambda run1, run2:  run1["domain_category"])) for i in [1, 2, 5, 10, 100, 1000]]

scatter_plots_k_dom_vars = [("report-eval/scatter/dom-vars-qrel0-vs-qrel{}-atomic".format(i), ScatterPlotReport(filter_algorithm=[ 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i)], attributes=["num_variables_with_dominance"], get_category=lambda run1, run2:  run1["domain_category"])) for i in [1, 2, 5, 10, 100, 1000]]

scatter_plots_k_dom_vars0 = [("report-eval/scatter/dom-vars0-qrel0-vs-qrel{}-atomic".format(i), ScatterPlotReport(filter_algorithm=[ 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i)], attributes=["num_variables_with_dominance_geq0"], get_category=lambda run1, run2:  run1["domain_category"])) for i in [1, 2, 5, 10, 100, 1000]]

scatter_plots_k_dom_vars1 = [("report-eval/scatter/dom-vars1-qrel0-vs-qrel{}-atomic".format(i), ScatterPlotReport(filter_algorithm=[ 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i)], attributes=["num_variables_with_dominance_geq1"], get_category=lambda run1, run2:  run1["domain_category"])) for i in [1, 2, 5, 10, 100, 1000]]

exp.add_scatter_plot_step(# scatter_plots_k + scatter_plots_k_time_preprocessing
    #scatter_plots_k_time_ldsim)    scatter_plots_k_min_dom)  
    scatter_plots_k_dom_vars + scatter_plots_k_dom_vars1 + scatter_plots_k_dom_vars0)
    
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
             'blind-bisim-atomic-nosh-gen',
             'blind-sim-atomic-nosh-gen',
             'blind-noopsim-atomic-nosh-gen',
             'blind-ldsim-atomic-nosh-gen',
             "blind-qual-10-atomic-nosh-gen",
             "blind-qpos-10-atomic-nosh-gen",
             "blind-qtrade-10-atomic-nosh-gen",
             "blind-qrel-10-atomic-nosh-gen"
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
            'blind-bisim-atomic-nosh-gen':    'bisim',
            'blind-sim-atomic-nosh-gen':      'sim',
            'blind-noopsim-atomic-nosh-gen':  'noopsim', 
            'blind-ldsim-atomic-nosh-gen':    'ldsim',
            "blind-qual-10-atomic-nosh-gen":  "qual10",
            "blind-qpos-10-atomic-nosh-gen":  "qpos10",
            "blind-qrel-10-atomic-nosh-gen":  "qrel10",
            "blind-qtrade-10-atomic-nosh-gen": "qtrade10"
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

