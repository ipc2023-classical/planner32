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

exp = ReportExperiment("report")

exp.add_fetcher('../properties/journal1-all/', postprocess_functions=[fix_algorithm, joint_domains])


all_configs=[
    'blind-ldsim-atomic-bissh-exp',
    'blind',
    'blind-ldsim-atomic-bissh-gen',
]


attributes = list(ReportExperiment.DEFAULT_TABLE_ATTRIBUTES)

exp.add_report(AbsoluteReport(attributes=attributes,filter_algorithm=all_configs))

exp.add_scatter_plot_step()

exp.run_steps()

