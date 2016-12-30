#include "numeric_simulation_relation.h"

#include "numeric_label_relation.h" 
#include "../merge_and_shrink/abstraction.h" 
#include "../merge_and_shrink/labelled_transition_system.h" 
#include "../priority_queue.h"
#include "../debug.h"

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
    precompute_shortest_path_with_tau(lts, lts_id, label_dominance);
    cout << "NumericSimulationRelation::update" << endl;
    // for (int s = 0; s < lts->size(); s++) {	
    // 	cout << g_fact_names[lts_id][s] << endl;
    // }

    while (changes) {
	cout << "updating numSim" << endl;	
	// dump();
        changes = false;
        for (int s = 0; s < lts->size(); s++) {
            for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
		// cout << "s: " << s << " t: " << t << endl;
		if (s != t && may_simulate(t, s)) {
                    int previous_value = q_simulates(t, s);
		    // cout << "prev: " << previous_value << endl;
		    int lower_bound = minus_shortest_path_with_tau(t,s);
		    // cout << "minimum value: " << previous_value << endl;		    

		    assert(lower_bound <= previous_value);
		    if(lower_bound == previous_value) {
			continue;
		    }

		    int min_value = previous_value;

		    if (lts->is_goal(t) || !lts->is_goal(s) ){
			//Check if really t simulates s
			//for each transition s--l->s':
			// a) with noop t >= s' and l dominated by noop?
			// b) exist t--l'-->t', t' >= s' and l dominated by l'?
			lts->applyPostSrc(s, [&](const LTSTransition & trs) {
				int max_value = std::numeric_limits<int>::lowest();
				// cout << "Checking transition " << trs.label << " "
				//      << g_operators[trs.label].get_name() << " to " << trs.target << endl;
				if(may_simulate(t, trs.target) && label_dominance.may_dominated_by_noop(trs.label, lts_id)) {
				    max_value = q_simulates(t, trs.target) +
					abs->get_label_cost_by_index(trs.label) +
					label_dominance.q_dominated_by_noop(trs.label, lts_id);

				    // cout << q_simulates(t,
				    // trs.target) << " " <<
				    // label_dominance.q_dominated_by_noop(trs.label,
				    // lts_id) << endl;
				
				    if (max_value >= min_value) {
					// cout << "Dominated by noop!" << endl;
					return false; // Go to next transition
				    }
				}
				
				lts->applyPostSrc(t,[&](const LTSTransition & trt) {
					if(label_dominance.may_dominate (trt.label, trs.label, lts_id) &&
					   may_simulate(trt.target, trs.target)) {
					    max_value = std::max(max_value, 
								 label_dominance.q_dominates(trt.label, trs.label, lts_id)
								 + abs->get_label_cost_by_index(trs.label)
								 - abs->get_label_cost_by_index(trt.label) +
								 + q_simulates(trt.target, trs.target));

					    return max_value >= min_value; 
					    //break if we have found a transition that simulates with the best result possible
					}
					return false;
				    });

				// if(min_value > 0 && max_value <= 0) {
				//     std::cout << lts->name(t) << " does not positevely simulate "
				// 	      << lts->name(s) << " because of "
				// 	      << lts->name(trs.src) << " => "
				// 	      << lts->name(trs.target) << " ("
				// 	      << trs.label << ")"; // << std::endl;
				//     std::cout << "  Simulates? "
				// 	      << q_simulates(t, trs.target);
				//     std::cout << "  domnoop? "
				// 	      << label_dominance.q_dominated_by_noop(trs.label, lts_id) << "   ";
				//     // label_dominance.dump(trs.label);
				//     for (auto trt : lts->get_transitions(t)) {
				// 	std::cout << "Tried with: "
				// 		  << lts->name(trt.src) << " => "
				// 		  << lts->name(trt.target) << " ("
				// 		  << trt.label << ")" << " label dom: "
				// 		  << label_dominance.q_dominates(trt.label,
				// 						 trs.label, lts_id)
				// 		  << " target sim "
				// 		  << q_simulates(trt.target, trs.target)
				// 		  << std::endl;
                                //     }
				// }
				min_value  = std::min(min_value, max_value);

				return min_value <= lower_bound;
			    });
			assert(min_value < std::numeric_limits<int>::max());
			
			min_value = std::max(min_value, lower_bound);
		    } else {
			min_value = lower_bound;
		    }

		    assert(min_value <= previous_value);
		    if(min_value < previous_value) {
			// cout << "Updating " << lts->get_names()[s] << " <= " << lts->get_names()[t]
			//      << " with " << new_value << " before " << previous_value
			//      << " minimum: " << minimum_value << endl;
			update_value(t, s,  min_value);
			changes = true;
		    }
		}
            }
        }
    }

    
    // for (int s = 0; s < lts->size(); s++) {	
    // 	cout << g_fact_names[lts_id][s] << endl;
    // }
    dump();

}



bool NumericSimulationRelation::pruned(const State & state) const {
    return abs->get_abstract_state(state) == -1;
}

int NumericSimulationRelation::q_simulates (const State & t, const State & s) const{
    int tid = abs->get_abstract_state(t);
    int sid = abs->get_abstract_state(s);
    return q_simulates(tid, sid);
}

int NumericSimulationRelation::q_simulates (const vector<int> & t, const vector<int> & s) const{
    int tid = abs->get_abstract_state(t);
    int sid = abs->get_abstract_state(s);
    return q_simulates(tid, sid);
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


void NumericSimulationRelation::precompute_shortest_path_with_tau(const LabelledTransitionSystem * lts, int lts_id,
								  const NumericLabelRelation & label_dominance) {
    const vector<int> & tau_labels = label_dominance.get_tau_labels(lts_id);
	
    // cout << "Number of tau labels: " << tau_labels_.size() << endl;
    // if(!distances_with_tau.empty() && tau_labels_ == tau_labels) return;
    // tau_labels = tau_labels_;

    int num_states = lts->size();
    //Create copy of the graph only with tau transitions
    vector<vector<pair<int, int> > > tau_graph(num_states);
    for (int label_no : tau_labels) {
        int label_cost = abs->get_label_cost_by_index(label_no) -
	    std::max(0, label_dominance.q_dominates_noop(label_no, lts_id));
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

void NumericSimulationRelation::dump(const vector<string> & names) const{ 
    cout << "SIMREL:" << endl;
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
            if(may_simulate(j, i) && i != j){
		cout << names[i] << " <= " << names[j] << " (" << q_simulates(j, i) << ")" << endl;
            }
        }
    }
}

void NumericSimulationRelation::dump() const{ 
    cout << "SIMREL:" << endl;
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
	    cout << q_simulates(j, i)  << " ";
	}
	cout << endl;
    }
}

