#include "numeric_simulation_relation.h"

#include "abstraction.h" 
#include "numeric_label_relation.h" 
#include "labelled_transition_system.h" 
#include "../priority_queue.h"

using namespace std;

NumericSimulationRelation::NumericSimulationRelation(Abstraction * _abs) : abs(_abs) { 
}

void NumericSimulationRelation::init_goal_respecting() {
    assert(abs->are_distances_computed());
    int num_states = abs->size();
    //const std::vector <bool> & goal_states = abs->get_goal_states();
    const std::vector <int> & goal_distances = abs->get_goal_distances();
    relation.resize(num_states);
    for(int s = 0; s < num_states; s++){
        relation[s].resize(num_states);
	for (int t = 0; t < num_states; t++){
	    // qrel (t, s) = h*(t) - h*(s)
	    relation[s][t] = goal_distances[t] - goal_distances[s];
	}
    }
}


void  NumericSimulationRelation::update (int lts_id, const LabelledTransitionSystem * lts,
					 const NumericLabelRelation & label_dominance) {
    bool changes = true;
    while (changes) {
        changes = false;
        for (int s = 0; s < lts->size(); s++) {
            for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
		if (s != t && may_simulate(t, s)) {
                    int previous_value = q_simulates(t, s); 
		    int minimum_value = minus_shortest_path_with_tau(t,s);

		    assert(minimum_value <= previous_value);
		    if(minimum_value == previous_value) {
			continue;
		    }

                    //Check if really t simulates s
                    //for each transition s--l->s':
                    // a) with noop t >= s' and l dominated by noop?
                    // b) exist t--l'-->t', t' >= s' and l dominated by l'?
                    lts->applyPostSrc(s, [&](const LTSTransition & trs) {
			    int new_value = minimum_value;
			    //cout << "Checking transition " << trs.label << " " << g_operators[trs.label].get_name() << " to " << trs.target << endl;
                            if(may_simulate(t, trs.target) && label_dominance.may_dominated_by_noop(trs.label, lts_id)) {
				new_value = std::max(new_value, 
						     q_simulates(t, trs.target) + label_dominance.q_dominated_by_noop(trs.label, lts_id));
				if (new_value >= previous_value) {
				    //cout << "Dominated by noop!" << endl;
				    return false; // Go to next transition
				}
                            }
                            
                            lts->applyPostSrc(t,[&](const LTSTransition & trt) {
				    if(label_dominance.may_dominate (trt.label, trs.label, lts_id) &&
				       may_simulate(trt.target, trs.target)) {
					new_value = std::max(new_value, 
							     label_dominance.q_dominates(trt.label, trs.label, lts_id)
							     + q_simulates(trt.target, trs.target));

					return new_value >= previous_value; 
					//break if we have found a transition that simulates with the best result possible
				    }
				    return false;
				});

                            if(new_value < previous_value) {
				previous_value = new_value;
				update_value(s, t, new_value);
                                changes = true;
                            }
                            return minimum_value == new_value;
                        });
		}
            }
        }
    }
}



bool NumericSimulationRelation::pruned(const State & state) const {
    return abs->get_abstract_state(state) == -1;
}

bool NumericSimulationRelation::simulates (const State & t, const State & s) const{
    int tid = abs->get_abstract_state(t);
    int sid = abs->get_abstract_state(s);
    return simulates(tid, sid);
}


//Copied from abstraction.cc (move somewhere else?)
static void dijkstra_search(
    const vector<vector<pair<int, int> > > &graph,
    AdaptiveQueue<int> &queue,
    vector<int> &distances) {
    while (!queue.empty()) {
        pair<int, int> top_pair = queue.pop();
        int distance = top_pair.first;
        int state = top_pair.second;
        int state_distance = distances[state];
        assert(state_distance <= distance);
        if (state_distance < distance)
            continue;
        for (int i = 0; i < graph[state].size(); i++) {
            const pair<int, int> &transition = graph[state][i];
            int successor = transition.first;
            int cost = transition.second;
            int successor_cost = state_distance + cost;
            if (distances[successor] > successor_cost) {
                distances[successor] = successor_cost;
                queue.push(successor_cost, successor);
            }
        }
    }
}


void NumericSimulationRelation::precompute_shortest_path_with_tau(const LabelledTransitionSystem * lts,
								  const vector<int> & tau_labels_){
    if(tau_labels_ == tau_labels) return;
    tau_labels = tau_labels_;

    int num_states = lts->size();
    //Create copy of the graph only with tau transitions
    vector<vector<pair<int, int> > > tau_graph(num_states);
    for (int label_no : tau_labels) {
        int label_cost = abs->get_label_cost_by_index(label_no);
        for (const auto & trans : lts->get_transitions_label(label_no)) {
	    if(trans.src != trans.target) {
		tau_graph[trans.src].push_back(make_pair(trans.target, label_cost));
	    }
        }
    }

    AdaptiveQueue<int> queue;
    distances_with_tau.resize(num_states);
    for(int s = 0; s < num_states; ++s) {
	//Perform Dijkstra search from s
	auto & distances = distances_with_tau [s];
	distances.resize(num_states);
	std::fill(distances.begin(), distances.end(), std::numeric_limits<int>::max());
	distances[s] = 0;
        queue.clear();
	queue.push(0, s);
	dijkstra_search(tau_graph, queue, distances);
    }
}
