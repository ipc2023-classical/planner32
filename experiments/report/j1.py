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
from downward.reports.plot_attributes import AttributeScatterPlotReport

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

exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations",  "time_ldsim"],filter_algorithm=[
    "blind-qrel-100-atomic-nosh-gen",
    "blind-qrel-100-atomic-nosh-gen-test2",
    'blind-qrel-100-atomic-nosh-gen-usedominatesin-test2'
]), outfile='report-test-qrel.html')

exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations", "min_negative_dominance", "max_positive_dominance", "has_positive_dominance", "has_negative_dominance", "percentage_variables_with_dominance", "num_variables_with_dominance", "total_num_variables", "time_ldsim"],filter_algorithm=[
    'blind-qrel-100-atomic-nosh-gen',
    'blind-qrel-100-atomic-nosh-gen-test2',
]), outfile='report-test.html')


exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations",  "time_ldsim"],filter_algorithm=[
    'blind-ldsimalt-atomic-nosh-gen-test2',
    'blind-ldsimalt-atomic-nosh-gen',
]), outfile='report-test-ldsim.html')


exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations"],filter_algorithm=[
    'blind',
    'blind-ldsim-atomic-nosh-gen',
    'blind-ldsimalt-atomic-nosh-gen',
]))


exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations"],filter_algorithm=[
    'blind-ldsimalt-atomic-nosh-gen',
    'blind-qual-10-atomic-nosh-gen',
]), outfile='report-ldsimalt-vs-qual.html')


exp.add_report(AbsoluteReport(attributes=list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES) + ["time_completed_preprocessing", "total_simulations", "only_simulations", "min_negative_dominance", "max_positive_dominance", "has_positive_dominance", "has_negative_dominance", "percentage_variables_with_dominance", "num_variables_with_dominance", "total_num_variables"],filter_algorithm=[
    'blind-qrel-0-atomic-nosh-gen',
    'blind-qrel-1-atomic-nosh-gen',
    'blind-qrel-2-atomic-nosh-gen',
    'blind-qrel-5-atomic-nosh-gen',
    'blind-qrel-10-atomic-nosh-gen',
    'blind-qrel-100-atomic-nosh-gen',
    'blind-qrel-1000-atomic-nosh-gen',
]), outfile='report-k.html')



def scatter_alg(name, alg1, alg2, atr, domain_category):
    return ("report-eval/scatter/{}".format(name), ScatterPlotReport(filter_algorithm=[alg1, alg2], attributes=[atr], get_category=lambda run1, run2:  run1[domain_category]))


def scatter_atr(name, alg_list, atr1, atr2, domain_category):
    return ("report-eval/scatter/{}".format(name), AttributeScatterPlotReport(format='png', xscale='log', filter_algorithm=alg_list, attributes=[atr2, atr1], get_category=lambda run:  run[domain_category]))



scatter_plots_tex = [
    # ("/home/alvaro/projects/journal_dominance/AIJ15-DominancePruning/figures/scatter/expansions-base-vs-ldsim-atomic",
    #  ScatterPlotReport(format='tex', get_category=lambda run1, run2:  run1["domain_category_2"], filter_algorithm=["blind", "blind-ldsimalt-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], xlabel="Expansions base", ylabel="Expansions LD simulation", xmin=1, ymin=1, title="Expansions until last jump")),
    # ("/home/alvaro/projects/journal_dominance/AIJ15-DominancePruning/figures/scatter/expansions-ldsim-vs-qrel10-atomic",
    #  ScatterPlotReport(format='tex', get_category=lambda run1, run2:  run1["domain_category_2"], filter_algorithm=["blind-ldsimalt-atomic-nosh-gen","blind-qrel-10-atomic-nosh-gen"], attributes=["expansions_until_last_jump"], ylabel="Expansions Qrel 10", xlabel="Expansions LD simulation", xmin=1, ymin=1, title="Expansions until last jump")),
    # ("/home/alvaro/projects/journal_dominance/AIJ15-DominancePruning/figures/scatter/timesim-ldsim-vs-qrel10-atomic",
    #  ScatterPlotReport(format='tex', get_category=lambda run1, run2:  run1["domain_category_timesim"], filter_algorithm=["blind-ldsimalt-atomic-nosh-gen","blind-qrel-10-atomic-nosh-gen"], attributes=["time_ldsim"], ylabel="Time Qrel 10 (s)", xlabel="Time LD simulation (s)", title="Time Simulation"))
    # ("/home/alvaro/projects/journal_dominance/AIJ15-DominancePruning/figures/scatter/ldsim-labels-vs-timesim-atomic",
    #  AttributeScatterPlotReport(format='png', get_category=lambda run :  run["algorithm"], filter_algorithm=["blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], attributes=["translator_operators", "time_ldsim"], yscale='log', xscale='log', ymin=100, ) ) 
]

scatter_plots_atomic = [scatter_alg("expansions-base-vs-ldsim-atomic", "blind", "blind-ldsimalt-atomic-nosh-gen", "expansions_until_last_jump", "domain_category_2"),
                        scatter_alg("expansions-ldsim-vs-qrel-atomic", "blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen", "expansions_until_last_jump", "domain_category_1"),
                        scatter_alg("timesim-ldsim-vs-qrel-atomic", "blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen", "time_ldsim", "domain_category"),
                        scatter_alg("timepreprocess-ldsim-vs-qrel-atomic", "blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen", "time_preprocessing", "domain_category"),
]


scatter_plots_k = [scatter_alg("report-eval/scatter/scatter_k_{}-qrel0-vs-qrel{}-atomic".format(atr, i), 'blind-qrel-0-atomic-nosh-gen',  'blind-qrel-{}-atomic-nosh-gen'.format(i), atr, "domain_category") for i in [1, 10, 100, 1000] for atr in ["expansions_until_last_jump", "time_completed_preprocessing", "time_ldsim", "min_negative_dominance_inverted", "percentage_variables_with_dominance","percentage_variables_with_dominance_geq0","percentage_variables_with_dominance_geq1" ]]


scatter_plots_timesim = [
    scatter_atr("operators-vs-timesim", ["blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], "translator_operators", "time_ldsim", "algorithm"),
    scatter_atr("variables-vs-timesim", ["blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], "translator_variables", "time_ldsim", "algorithm"),
    scatter_atr("facts-vs-timesim", ["blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], "translator_facts", "time_ldsim", "algorithm"),
    scatter_atr("inner_iterations-vs-timesim", ["blind-qrel-10-atomic-nosh-gen"], "inner_iterations_numeric_ldsimulation", "time_ldsim", "algorithm"),
    scatter_atr("operators-vs-timeprepro", ["blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], "translator_operators", "time_completed_preprocessing", "algorithm"),
    scatter_atr("variables-vs-timeprepro", ["blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], "translator_variables", "time_completed_preprocessing", "algorithm"),
    scatter_atr("facts-vs-timeprepro", ["blind-ldsimalt-atomic-nosh-gen", "blind-qrel-10-atomic-nosh-gen"], "translator_facts", "time_completed_preprocessing", "algorithm"),
    scatter_atr("inner_iterations-vs-timeprepro", ["blind-qrel-10-atomic-nosh-gen"], "inner_iterations_numeric_ldsimulation", "time_completed_preprocessing", "algorithm")
]


# for p in scatter_plots_k_dom_vars + scatter_plots_k_dom_vars0 + scatter_plots_k_dom_vars1:
#     p[1]._set_scales("linear", "linear")

    
exp.add_scatter_plot_step(scatter_plots_timesim)


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

