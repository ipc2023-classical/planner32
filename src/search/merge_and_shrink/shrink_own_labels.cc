#include "shrink_own_labels.h"

#include "abstraction.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../globals.h"
#include "../operator.h"
#include "../scc.h"
#include "../equivalence_relation.h"

#include <cassert>
#include <iostream>
#include <limits>

ShrinkOwnLabels::ShrinkOwnLabels(const Options &opts) : 
    ShrinkStrategy(opts), 
    perform_sg_shrinking (opts.get<bool>("goal_shrinking"))  {
}

ShrinkOwnLabels::~ShrinkOwnLabels() {
}

string ShrinkOwnLabels::name() const {
    return "own labels (to identify unsol. tasks)";
}

void ShrinkOwnLabels::dump_strategy_specific_options() const {
    cout << "Aggregate with goal states: " << 
	(perform_sg_shrinking? "yes" : "no") << endl;
}

bool ShrinkOwnLabels::reduce_labels_before_shrinking() const {
    return true;
}

void ShrinkOwnLabels::shrink(Abstraction &abs, int target, bool /*force*/) {
    /*
     * apply two rules:
     * (1) aggregate all states s1, ..., sn if they lie on an own-label cycle
     *     idea: use Tarjan's algorithm to identify strongly connected components (SCCs)
     * (2) aggregate state s with goal state g if:
     *     (a) perform_sg_shrinking is activated
     *     (b) all goal variables are in abstraction and
     *     (c) there is an own-label path from s to g
     *     idea: set goal-status for all states already during SCC detection
     *     PETER: why only aggregate s with g and not all states that
     *     are marked as goal states?
     *     PETER: when we are at it, why not remove all outgoing
     *            transitions of goal states (we cannot leave them
     *            anyway any more)?
     *     PETER: actually, if we use the first optimization, we must
     *            use the second one as well, as otherwise unreachable
     *            states might suddenly become reachable resulting in
     *            an immense overhead in abstract states
     *     PETER: we might actually do this as a general optimization
     *     in Abstraction::normalize: whenever all goal variables are
     *     merged in, we can safely remove any transitions starting at
     *     goal states
     *     PETER: removing goal transitions is not safe.
     * In case of own-label shrinking we inherit the information of
     * being a goal state from the own-label reachable
     * successors. If this way the initial state is marked as a goal,
     * but some states are not marked as goal because their outgoing
     * actions are not yet own-labeled, then removing the outgoing
     * actions of the goal states will mean that those states will
     * become unreachable and thus they will be pruned
     * CONCLUSION PETER: aggregating all goal-states is added as an
     * option in bisimulation, cause the abstraction size may increase
     * due to less reachability pruning. Pruning goal transitions is
     * disabled because is not safe
     */

    int num_states = abs.size();
    std::vector<bool> is_goal (abs.get_goal_states());
    
    /* this is a rather memory-inefficient way of implementing
       Tarjan's algorithm, but it's the best I got for now */
    vector<vector<int> > adjacency_matrix(num_states);
    int num_labels = abs.get_num_labels();
    for (int label_no = 0; label_no < num_labels; ++label_no) {
	if (!abs.is_own_label(label_no) || abs.get_label_cost_by_index(label_no) > 0)
	    continue;
	const vector<AbstractTransition> &transitions =
	    abs.get_transitions_for_label(label_no);
	for (size_t i = 0; i < transitions.size(); ++i) {
	    const AbstractTransition &trans = transitions[i];
	    adjacency_matrix[trans.src].push_back(trans.target);
	}
    }

    /* PETER: Can we do better than this, i.e., prevent the
       sorting? */
    /* remove duplicates in adjacency matrix */
    for (size_t i = 0; i < num_states; i++) {
        ::sort(adjacency_matrix[i].begin(), adjacency_matrix[i].end());
        vector<int>::iterator it = unique(adjacency_matrix[i].begin(), adjacency_matrix[i].end());
        adjacency_matrix[i].erase(it, adjacency_matrix[i].end());
    }

    /* perform Tarjan's algorithm for finding SCCs */
    EquivalenceRelation final_sccs;
    SCC::compute_scc_equivalence (adjacency_matrix, is_goal, final_sccs);

    /*cout << "===========================================" << endl;
      for (int i = 0; i < num_states; i++) {
      cout << "edges from " << i << " to";
      for (int j = 0; j < adjacency_matrix[i].size(); j++)
      cout << " " << adjacency_matrix[i][j];
      cout << endl;
      }
      cout << "found SCCs:" << endl;
      for (int i = 0; i < final_sccs.size(); i++) {
      for (int j = 0; j < final_sccs[i].size(); j++)
      cout << final_sccs[i][j] << " ";
      cout << endl;
      }
      cout << "===========================================" << endl;*/

    /* free some memory */
    vector<vector<int> > ().swap(adjacency_matrix);

    int new_size = final_sccs.size();
    if (perform_sg_shrinking && abs.get_all_goal_vars_in()) {
        /* now bring those groups together that follow the second rule */
        cout << "also using second rule of own-label shrinking" << endl;
        int goal_scc = -1;
        for (int i = 0; i < final_sccs.size(); i++) {
	    if (is_goal[final_sccs[i].front()]) {
		if (goal_scc == -1)
		    goal_scc = i;
		else {
		    final_sccs[goal_scc].splice(final_sccs[goal_scc].end(), final_sccs[i]);
		    new_size--;
		}
	    }
        }
    }

    if (new_size < num_states) {
        // only need to apply abstraction if this actually changes anything
	EquivalenceRelation equivalence_relation(new_size);
	int counter = 0;
	for (size_t group = 0; group < final_sccs.size(); ++group) {
	    if (final_sccs[group].empty())
		continue;
	    equivalence_relation[counter].swap(final_sccs[group]);
	    counter++;
	}
	apply(abs, equivalence_relation, target);
	//vector<int> ().swap(state_to_group);
	EquivalenceRelation ().swap(equivalence_relation);
    } else {
	cout << "Own-label shrinking does not reduce states" << endl;
    }
    return;
}

void ShrinkOwnLabels::shrink_atomic(Abstraction &abs) {
    shrink(abs, abs.size(), true);
}

void ShrinkOwnLabels::shrink_before_merge(Abstraction &abs1, 
					    Abstraction &abs2) {
    shrink(abs1, abs1.size(), true);
    shrink(abs2, abs2.size(), true);
    return;
}

ShrinkOwnLabels *ShrinkOwnLabels::create_default() {
    const int infinity = numeric_limits<int>::max();
    Options opts;
    opts.set("max_states", infinity);
    opts.set("max_states_before_merge", infinity);
    return new ShrinkOwnLabels(opts);
}

static ShrinkStrategy *_parse(OptionParser &parser) {
    ShrinkStrategy::add_options_to_parser(parser);

    parser.add_option<bool>("goal_shrinking",
                            "performs goal shrinking. Aggregate state s with goal state g if:"
			    "   (a) this parameter is activated"
			    "   (b) all goal variables are in abstraction and"
			    "   (c) there is an own-label path from s to g",
                            "true");

    Options opts = parser.parse();
    ShrinkStrategy::handle_option_defaults(opts);



    if (!parser.dry_run())
	return new ShrinkOwnLabels(opts);
    else
	return 0;
}

static Plugin<ShrinkStrategy> _plugin("shrink_own_labels", _parse);
