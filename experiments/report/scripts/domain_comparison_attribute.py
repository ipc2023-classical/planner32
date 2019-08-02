from downward.reports import PlanningReport

from collections import defaultdict

def domain_mapping(domain):
    if domain == 'openstacks':
        return 'openstacks-06'
    elif 'openstacks' in domain:
        return 'openstacks-08-11-14'
    elif 'logistics' in domain:
        return 'logistics'
    else:
        result = domain.replace('-adl', '')
        result = domain.replace('-strips', '')
        result = result.replace('-08', '')
        result = result.replace('-opt08', '')
        result = result.replace('-opt11', '')
        result = result.replace('-opt14', '')
        result = result.replace('-opt18', '')
        result = result.replace('-sat08', '')
        result = result.replace('-sat11', '')
        result = result.replace('-sat14', '')
        result = result.replace('-sat18', '')
        return result

# TODO: this currently uses coverage and sum as the hard-coded attribute and
# aggregation function.
# NOTE: hacked to consider same domains of different years as one domain.
class DomainComparisonAttributeReport(PlanningReport):
    def __init__(self, algo_to_print = {}, attribute = "expansions_until_last_jump", comparison = (lambda x, y : x<y), **kwargs):
        PlanningReport.__init__(self, **kwargs)
        self.algo_to_print = algo_to_print
        self.attribute = attribute
        self.comparison = comparison

    def get_text(self):
        def turn_list_into_table_row(line):
            result = ''
            for index, value in enumerate(line):
                result += '{}'.format(value)
                if index == len(line) - 1:
                    result += ' \\\\'
                else:
                    result += ' & '
            return result

        def format_algo(algo):
            if algo in self.algo_to_print:
                return self.algo_to_print[algo]
            return algo

        domains = set()
        domain_and_algo_to_data = defaultdict(int)
        instances_with_data = defaultdict(int)
        for (domain, problem), runs in self.problem_runs.items():
            for run in runs:
                if not run['algorithm'] in self.algorithms or not self.attribute in run or not run[self.attribute]:
                    continue
                assert domain == run['domain']
                domain_and_algo_to_data[(domain, problem, run['algorithm'])] = run[self.attribute]
                instances_with_data [(domain, problem)] += 1

        selected_instances = set([(domain, problem) for (domain, problem) in instances_with_data if instances_with_data [(domain, problem)] == len(self.algorithms)])
               
        # end for
        lines = [turn_list_into_table_row([''] + [format_algo(algo) for algo in self.algorithms])]

        for algo1 in self.algorithms:
            line = []
            line.append(format_algo(algo1))
            for algo2 in self.algorithms:
                num_algo1_better = 0
                num_algo2_better = 0
                domains_algo1_better = set()
                domains_algo2_better = set()
                for (domain, problem) in selected_instances:
                    mapped_domain = domain_mapping(domain)
                    value1 = domain_and_algo_to_data[(domain, problem, algo1)]
                    value2 = domain_and_algo_to_data[(domain, problem, algo2)]
                    if self.comparison(value1, value2):
                        domains_algo1_better.add(mapped_domain)
                        num_algo1_better += 1
                    elif self.comparison(value2, value1):
                        domains_algo2_better.add(mapped_domain)
                        num_algo2_better += 1

                num_domains_algo1_better = len(domains_algo1_better)
                num_domains_algo2_better = len(domains_algo2_better)

                if algo1 == algo2:
                    entry = "--"
                elif num_algo1_better >= num_algo2_better:
                    entry = "\\textbf{{{} ({})}}".format(num_algo1_better, num_domains_algo1_better)
                else:
                    entry = "{} ({})".format(num_algo1_better, num_domains_algo1_better)
                line.append(entry)
            # aggregated value column
            lines.append(turn_list_into_table_row(line))
        return '\n'.join(lines)
