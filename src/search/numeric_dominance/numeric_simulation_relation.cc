#include "numeric_simulation_relation.h"

#include "numeric_label_relation.h" 
#include "../merge_and_shrink/abstraction.h" 
#include "../merge_and_shrink/labelled_transition_system.h" 
#include "../priority_queue.h"
#include "../debug.h"

using namespace std;

NumericSimulationRelation::NumericSimulationRelation(Abstraction * _abs, int truncate_value_) : abs(_abs), truncate_value(truncate_value_) { 
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

int NumericSimulationRelation::compare_transitions(int lts_id, const LTSTransition & trs, 
						   const LTSTransition & trt, 
						   int tau_distance, const NumericLabelRelation & label_dominance) const {
    if(label_dominance.may_dominate (trt.label, trs.label, lts_id) &&
       may_simulate(trt.target, trs.target)) {
	return tau_distance + label_dominance.q_dominates(trt.label, trs.label, lts_id)
	    + abs->get_label_cost_by_index(trs.label)
	    - abs->get_label_cost_by_index(trt.label)
	    + q_simulates(trt.target, trs.target);
    } else {
	return std::numeric_limits<int>::lowest();
    }     
}

int NumericSimulationRelation::compare_noop(int lts_id, const LTSTransition & trs, 
					    int t,  int tau_distance, 
					    const NumericLabelRelation & label_dominance) const {

    // Checking noop 
    if(may_simulate(t, trs.target) && label_dominance.may_dominated_by_noop(trs.label, lts_id)) {
	return tau_distance + 
	    q_simulates(t, trs.target) +
	    abs->get_label_cost_by_index(trs.label) +
	    label_dominance.q_dominated_by_noop(trs.label, lts_id);
    } else {
	return std::numeric_limits<int>::lowest();
    }    
}				


int NumericSimulationRelation::update (int lts_id, const LabelledTransitionSystem * lts,
				       const NumericLabelRelation & label_dominance) {
    int num_iterations = 0;
    bool changes = true;
    precompute_shortest_path_with_tau(lts, lts_id, label_dominance);
    // for (int s = 0; s < lts->size(); s++) {	
    // 	cout << g_fact_names[lts_id][s] << endl;
    // }

    while (changes) {
	num_iterations ++;
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

				for (int t2 : reachable_with_tau[t]) {
				    int tau_distance = minus_shortest_path_with_tau(t, t2);

				    max_value = max(max_value, 
						    compare_noop(lts_id, trs, t2, tau_distance, label_dominance));
				    if (max_value >= min_value) {
					return false; // Go to next transition
				    }
			       
				    lts->applyPostSrc(t2,[&](const LTSTransition & trt) {
					    max_value = max(max_value, 
							    compare_transitions(lts_id, trs, trt, tau_distance, label_dominance));
					    return max_value >= min_value; 
					    //break if we have found a transition that simulates with the best result possible
					});

				    if (max_value >= min_value) {
					return false; 
				    }
				}

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

    return num_iterations;
    
    // for (int s = 0; s < lts->size(); s++) {	
    // 	cout << g_fact_names[lts_id][s] << endl;
    // }
    //dump(g_fact_names[lts_id]);
}



bool NumericSimulationRelation::pruned(const State & state) const {
    return abs->get_abstract_state(state) == -1;
}

int NumericSimulationRelation::q_simulates (const State & t, const State & s) const{
    int tid = abs->get_abstract_state(t);
    int sid = abs->get_abstract_state(s);
    if(sid == -1 || tid == -1) {
	return std::numeric_limits<int>::lowest(); 
    }
    return q_simulates(tid, sid);
}

int NumericSimulationRelation::q_simulates (const vector<int> & t, const vector<int> & s) const{
    int tid = abs->get_abstract_state(t);
    int sid = abs->get_abstract_state(s);
    if(sid == -1 || tid == -1) {
	return std::numeric_limits<int>::lowest(); 
    }
    return q_simulates(tid, sid);
}

//Copied from abstraction.cc (move somewhere else?)
static void dijkstra_search(
    const vector<vector<pair<int, int> > > &graph,
    AdaptiveQueue<int> &queue,
    vector<int> &distances,
    vector<int> & states_reached) {
    while (!queue.empty()) {
        pair<int, int> top_pair = queue.pop();
        int distance = top_pair.first;
        int state = top_pair.second;
        int state_distance = distances[state];
	states_reached.push_back(state);
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

    if(!label_dominance.have_tau_labels_changed(lts_id) && 
       !distances_with_tau.empty()) return;
    const vector<int> & tau_labels = label_dominance.get_tau_labels(lts_id);
	
    //cout << "Number of tau labels: " << tau_labels.size() << endl;
    // if(!distances_with_tau.empty() && tau_labels_ == tau_labels) return;
    // tau_labels = tau_labels_;

    int num_states = lts->size();
    //Create copy of the graph only with tau transitions
    vector<vector<pair<int, int> > > tau_graph(num_states);
    for (int label_no : tau_labels) {
	if(label_dominance.may_dominate_noop_in(label_no, lts_id)) {
	    int label_cost = abs->get_label_cost_by_index(label_no) +
		std::min(0,  - label_dominance.q_dominates_noop(label_no, lts_id));
	    for (const auto & trans : lts->get_transitions_label(label_no)) {
		if(trans.src != trans.target) {
		    tau_graph[trans.src].push_back(make_pair(trans.target, label_cost));
		}
	    }
	}
    }

    AdaptiveQueue<int> queue;
    distances_with_tau.resize(num_states);
    reachable_with_tau.resize(num_states);
    for(int s = 0; s < num_states; ++s) {
	//Perform Dijkstra search from s
	auto & distances = distances_with_tau [s];
	distances.resize(num_states);
	std::fill(distances.begin(), distances.end(), std::numeric_limits<int>::max());
	distances[s] = 0;
        queue.clear();
	queue.push(0, s);
	dijkstra_search(tau_graph, queue, distances, reachable_with_tau[s]);
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



void NumericSimulationRelation::statistics() const {
    map<int, int> values;
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
	    if(i == j) continue;
	    // if (relation[i][j] > std::numeric_limits<int>::lowest()) {
	    if(values.count(relation[i][j])) {
		values[relation[i][j]] ++;
	    } else {
		values[relation[i][j]] = 1;
	    }
	    // }
	}
    }

    for (auto & it : values) {
	if(it.first == std::numeric_limits<int>::lowest()) {
	    cout << "-infinity"; 
	} else {
	    cout << it.first;
	}

	cout << ": " << it.second << endl;
    }
}


void NumericSimulationRelation::precompute_absstate_bdds(SymVariables * vars_){
    vars = vars_;
    abs->getAbsStateBDDs(vars, abs_bdds);
}


void NumericSimulationRelation::precompute_bdds(bool dominating, bool quantified, bool use_ADD) {
    if(!quantified) {
	if(dominating) {
	    cout << "Precomputing dominating_bdds ... " << endl;
	    dominating_bdds.resize(abs->size(), vars->zeroBDD());
	}else {
	    cout << "Precomputing dominated_bdds ... " << endl;
	    
	    // for (int i = 0; i < abs->size(); i++){
	    // 	dominated_bdds.push_back(vars->zeroBDD());
	    // }

	    dominated_bdds.resize(abs->size(), vars->zeroBDD());
	}
	
	for(int i = 0; i < relation.size(); ++i){
	    for(int j = 0; j < relation.size(); ++j){
		if(may_simulate(i, j) && q_simulates (i, j) >= 0){
		    if(dominating) {
			dominating_bdds[j] += abs_bdds[i];
		    } else { 
			dominated_bdds[i] += abs_bdds[j];
		    }
		}
	    }
	}
    } else if (use_ADD) {
	if(dominating) {
	    dominating_adds.resize(abs->size(), vars->getADD(0));
	    may_dominating_bdds.resize(abs->size(), vars->zeroBDD());
	}else {
	    dominated_adds.resize(abs->size(), vars->getADD(0));
	    may_dominated_bdds.resize(abs->size(), vars->zeroBDD());
	}
	
	for(int i = 0; i < relation.size(); ++i){
	    for(int j = 0; j < relation.size(); ++j){
		int value = may_simulate(i, j) ? q_simulates (i, j) : 
		    std::numeric_limits<int>::lowest(); 

		ADD val = vars->getADD(value);
		if(dominating) {
		    may_dominating_bdds[j] += abs_bdds[i];
		    dominating_adds[j] += abs_bdds[i].Add()*val;
		} else { 
		    may_dominated_bdds[i] +=  abs_bdds[j];
		    dominated_adds[i] += abs_bdds[j].Add()*val;
		}
	    }
	}	
    } else {
	if(dominating) {
	    cout << "Precomputing dominating_bdd_maps ... " << endl;

	    dominating_bdd_maps.resize(abs->size());
	    may_dominating_bdds.resize(abs->size(), vars->zeroBDD()); 
	}else {
	    cout << "Precomputing dominated_bdd_maps ... " << endl;

	    dominated_bdd_maps.resize(abs->size());
	    may_dominated_bdds.resize(abs->size(), vars->zeroBDD());
	}
	
	for(int i = 0; i < relation.size(); ++i){
	    for(int j = 0; j < relation.size(); ++j){
		if(may_simulate(i, j)) {
		    int val = q_simulates (i, j);
		    if(dominating) {
			may_dominating_bdds[j] += abs_bdds[i];
			
			if(dominating_bdd_maps[j].count(val)) {
			    dominating_bdd_maps[j][val] += abs_bdds[i];
			} else { 
			    dominating_bdd_maps[j][val] = abs_bdds[i];
			}
		    }else {
			may_dominated_bdds[i] +=  abs_bdds[j];
			if(dominated_bdd_maps[i].count(val)) {
			    dominated_bdd_maps[i][val] += abs_bdds[j];
			} else { 
			    dominated_bdd_maps[i][val] = abs_bdds[j];
			}
		    }  
		}
	    }
	}
    }
}

BDD NumericSimulationRelation::getSimulatedBDD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return dominated_bdds[absstate];
}

BDD NumericSimulationRelation::getSimulatingBDD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return dominating_bdds[absstate];
}


const map<int, BDD> & NumericSimulationRelation::getSimulatedBDDMap(const State & state) const{
    assert(!dominated_bdd_maps.empty());
    int absstate = abs->get_abstract_state(state);
    assert(absstate != -1); 
    return dominated_bdd_maps[absstate];
}

const map<int, BDD> & NumericSimulationRelation::getSimulatingBDDMap(const State & state) const{
    assert(!dominated_bdd_maps.empty());
    int absstate = abs->get_abstract_state(state);
    assert(absstate == -1);
    return dominating_bdd_maps[absstate];
}


BDD NumericSimulationRelation::getMaySimulatedBDD(const State & state) const{
    assert(!may_dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return may_dominated_bdds[absstate];
}

BDD NumericSimulationRelation::getMaySimulatingBDD(const State & state) const{
    assert(!may_dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return may_dominating_bdds[absstate];
}


ADD NumericSimulationRelation::getSimulatedADD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->getADD(std::numeric_limits<int>::lowest());
    else return dominated_adds[absstate];
}

ADD NumericSimulationRelation::getSimulatingADD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->getADD(std::numeric_limits<int>::lowest());
    else return dominating_adds[absstate];
}
