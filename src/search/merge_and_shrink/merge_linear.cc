#include "merge_linear.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../utilities.h"
#include "abstraction.h"

#include <cassert>
#include <cstdlib>

using namespace std;

MergeLinear::MergeLinear(const Options &opts)
    : MergeStrategy(),
      order(VariableOrderType(opts.get_enum("variable_order"))),
      need_first_index(true) {
}

// Alvaro: Merge strategies have now a limit on the size of the merge.
// If specified (> 0), the pair returned should fit the constraint:
// a1.size()*a2.size()<=limit
pair<int, int> MergeLinear::get_next(const std::vector<Abstraction *> &all_abstractions, 
				     int limit_abstract_states_merge, int /*min_limit_abstract_states_merge*/, 
				     int /*limit_transitions_merge*/) {
    assert(!done() && !order.done());

    int first;
    if (need_first_index) {
        need_first_index = false;
        first = order.next();
        //This may happen if some variables are eliminated due to being irrelevant
        while (!all_abstractions[first] && !done() && !order.done()) {
            first = order.next();
            remaining_merges--;
        }
        if (!all_abstractions[first]) {
            return make_pair(-1, -1);
        }
        
        cout << "First variable: " << first << endl;
    } else {
        // The most recent composite abstraction is appended at the end of
        // all_abstractions in merge_and_shrink.cc
        first = all_abstractions.size() - 1;
    }
    int second = order.next();

    while (!all_abstractions[second] && !done() && !order.done()) {
        cout << "Skipping var " << second << endl;
        second = order.next();
        remaining_merges--;
    }
    if (!all_abstractions[second]) {
        return make_pair(-1, -1);
    }

    cout << "Next variable: " << second << endl;
    assert(all_abstractions[first]);
    assert(all_abstractions[second]);
    --remaining_merges;
    if (done() && !order.done()) {
        cerr << "Variable order finder not done, but no merges remaining" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    if(limit_abstract_states_merge && 
       all_abstractions[first]->size()*all_abstractions[second]->size() > limit_abstract_states_merge){
      if(!done()){
	return get_next(all_abstractions, limit_abstract_states_merge);
      }else{
	return make_pair(-1, -1);
      }
    }
    return make_pair(first, second);
}

void MergeLinear::dump_strategy_specific_options() const {
    cout << "Linear merge strategy: ";
    order.dump();
}

string MergeLinear::name() const {
    return "linear";
}

bool MergeLinear::is_linear() const {
    return true;
}

static MergeStrategy *_parse(OptionParser &parser) {
    vector<string> merge_strategies;
    //TODO: it's a bit dangerous that the merge strategies here
    // have to be specified exactly in the same order
    // as in the enum definition. Try to find a way around this,
    // or at least raise an error when the order is wrong.
    merge_strategies.push_back("CG_GOAL_LEVEL");
    merge_strategies.push_back("CG_GOAL_RANDOM");
    merge_strategies.push_back("GOAL_CG_LEVEL");
    merge_strategies.push_back("RANDOM");
    merge_strategies.push_back("LEVEL");
    merge_strategies.push_back("REVERSE_LEVEL");
    parser.add_enum_option("variable_order", merge_strategies,
                           "the order in which atomic abstractions are merged",
                           "CG_GOAL_LEVEL");

    Options opts = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return new MergeLinear(opts);
}

static Plugin<MergeStrategy> _plugin("merge_linear", _parse);

