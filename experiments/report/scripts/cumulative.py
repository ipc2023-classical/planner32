# -*- coding: utf-8 -*-

from collections import defaultdict
import logging
import math
import os

from lab import tools

from downward.reports.plot import MatplotlibPlot, Matplotlib, PgfPlots, \
    PlotReport, MIN_AXIS


# class CumulativeMatplotlib(Matplotlib):
#     @classmethod
#     def _plot(cls, report, axes, categories, styles):
#         # Display grid
#         axes.grid(b=True, linestyle='-', color='0.75')

#         has_points = False
#         # Generate the scatter plots
#         for category, coords in sorted(categories.items()):
#             X, Y = zip(*coords)
#             axes.scatter(X, Y, s=42, label=category, **styles[category])
#             if X and Y:
#                 has_points = True

#         if report.xscale == 'linear' or report.yscale == 'linear':
#             plot_size = report.missing_val * 1.01
#         else:
#             plot_size = report.missing_val * 1.5

#         # Plot a diagonal black line. Starting at (0,0) often raises errors.
#         axes.plot([0.001, plot_size], [0.001, plot_size], 'k')

#         axes.set_xlim(report.xlim_left or -1, report.xlim_right or plot_size)
#         axes.set_ylim(report.ylim_bottom or -1, report.ylim_top or plot_size)

#         for axis in [axes.xaxis, axes.yaxis]:
#             MatplotlibPlot.change_axis_formatter(
#                 axis, report.missing_val if report.show_missing else None)
#         return has_points


class CumulativePgfPlots(PgfPlots):
    @classmethod
    def _format_coord(cls, coord):
        def format_value(v):
            return str(v) if isinstance(v, int) else '%f' % v
        return '(%s, %s)' % (format_value(coord[0]), format_value(coord[1]))

    @classmethod
    def _get_plot(cls, report):
        options = cls._get_axis_options(report)
        plots_data = []
        for algorithm, coords in sorted(report.categories.items()):
            plot_tex_line = """
            \\addplot+[nomarks] coordinates {{
           {}
            }};
            \\addlegendentry{{{}}}
            """
            
            plots_data.append(plot_tex_line.format(coords, algorithm))

        plot_tex = """
        \begin{{axis}}[{axis_options}
        {data}
        \end{{axis}}
        """.format(axis_options=cls._format_options(options), data="\n".join(plots_data))

        return [plot_tex]


class CumulativePlotReport(PlotReport):
    """
    Generate a scatter plot for a specific attribute.
    """
    def __init__(self, show_missing=True, get_category=None, **kwargs):
        """
        See :class:`.PlotReport` for inherited arguments.

        The keyword argument *attributes* must contain exactly one
        attribute.
        """
        # If the size has not been set explicitly, make it a square.
        matplotlib_options = kwargs.get('matplotlib_options', {})
        matplotlib_options.setdefault('figure.figsize', [8, 8])
        kwargs['matplotlib_options'] = matplotlib_options
        PlotReport.__init__(self, **kwargs)
        if not self.attribute:
            logging.critical('CumulativePlotReport needs exactly one attribute')
        # By default all values are in the same category.
        self.get_category = get_category or (lambda run1, run2: None)
        self.show_missing = show_missing
        self.xlim_left = self.xlim_left or MIN_AXIS
        self.ylim_bottom = self.ylim_bottom or MIN_AXIS
        if self.output_format == 'tex':
            self.writer = CumulativePgfPlots
        else:
            raise "not supported"
            # self.writer = CumulativeMatplotlib

    def _set_scales(self, xscale, yscale):
        PlotReport._set_scales(self, xscale or self.attribute.scale or 'log', yscale)

    def _fill_categories(self, runs):
        # We discard the *runs* parameter.
        # Map category names to value tuples
        categories = defaultdict(list)
        for (domain, problem), runs in self.problem_runs.items():
            if len(runs) != 2:
                continue
            run1, run2 = runs
            assert (run1['algorithm'] == self.algorithms[0] and
                    run2['algorithm'] == self.algorithms[1])
            val1 = run1.get(self.attribute)
            val2 = run2.get(self.attribute)
            if val1 is None and val2 is None:
                continue
            category = self.get_category(run1, run2)
            categories[category].append((val1, val2))
        return categories

    def write(self):
        if not len(self.algorithms) == 2:
            logging.critical(
                'Cumulative plots need exactly 2 algorithms: %s' % self.algorithms)
        self.xlabel = self.xlabel or self.algorithms[0]
        self.ylabel = self.ylabel or self.algorithms[1]

        suffix = '.' + self.output_format
        if not self.outfile.endswith(suffix):
            self.outfile += suffix
        tools.makedirs(os.path.dirname(self.outfile))
        self._write_plot(self.runs.values(), self.outfile)
