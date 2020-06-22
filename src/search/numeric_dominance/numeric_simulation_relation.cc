#include "numeric_simulation_relation.h"
#include "numeric_label_relation.h" 
#include "../merge_and_shrink/abstraction.h" 
#include "../merge_and_shrink/labelled_transition_system.h" 
#include "../priority_queue.h"
#include "../debug.h"
#include "dijkstra_search_epsilon.h"
#include "breadth_first_search.h"

using namespace std;

template <typename T>
NumericSimulationRelation<T>::NumericSimulationRelation(Abstraction * _abs, 
							int truncate_value_,
							std::shared_ptr<TauLabelManager<T>> tau_labels_mgr) : abs(_abs), 
													   truncate_value(truncate_value_), tau_labels(tau_labels_mgr),
													   max_relation_value(0), cancelled(false) { 
    assert(abs);
}

template <typename T>
void NumericSimulationRelation<T>::init_goal_respecting() {
   
    assert(abs->are_distances_computed());
    //assert(abs->no_dead_labels()); 
    int num_states = abs->size();
    //const std::vector <bool> & goal_states = abs->get_goal_states();
    const std::vector <int> & goal_distances = abs->get_goal_distances();

    relation.resize(num_states);
    // is_relation_stable.resize(num_states);
    for(int s = 0; s < num_states; s++){
        relation[s].resize(num_states);
	// is_relation_stable[s].resize(num_states, false);
	for (int t = 0; t < num_states; t++){
	    // qrel (t, s) = h*(t) - h*(s)
	    // if (!abs->is_goal_state(t) && abs->is_goal_state(s)) {
	    // 	relation[s][t] = goal_distances_with_tau[t] ;
	    // } else {
	    // 	relation[s][t] = goal_distances[t] - goal_distances[s];
	    // }

	    //Here we have not computed tau distances yet (because
	    //with dominated by noop version may change)
	    relation[s][t] = goal_distances[t] - goal_distances[s];
	}
    }
    tau_distances_id = 0;
}


template <>
void NumericSimulationRelation<IntEpsilon>::init_goal_respecting() {
    assert(abs->are_distances_computed());
    //assert(abs->no_dead_labels()); 
    int num_states = abs->size();
    cout << "Recompute distances with epsilon" << endl;
    std::vector <IntEpsilonSum> goal_distances = abs->recompute_goal_distances_with_epsilon();
    relation.resize(num_states);
    // is_relation_stable.resize(num_states);

    for(int s = 0; s < num_states; s++){
        relation[s].resize(num_states);
	// is_relation_stable[s].resize(num_states, false);

	for (int t = 0; t < num_states; t++) {
	    // qrel (t, s) = h*(t) - h*(s)
	    // if (!abs->is_goal_state(t) && abs->is_goal_state(s)) {
	    // 	relation[s][t] = goal_distances_with_tau[t;]
	    // } else {
		IntEpsilonSum rel = (goal_distances[t] - goal_distances[s]);
		relation[s][t] = rel.get_epsilon_negative();
	    // }
	}
    }
    tau_distances_id = 0;
}

template <typename T>
T NumericSimulationRelation<T>::compare_transitions(int lts_id, int tr_s_target, int tr_s_label,
						    int tr_t_target, int tr_t_label,
						    T tau_distance, 
						    const NumericLabelRelation<T> & label_dominance) const {
    if(label_dominance.may_dominate (tr_t_label, tr_s_label, lts_id) &&
       may_simulate(tr_t_target, tr_s_target)) {

	return tau_distance +
	    label_dominance.q_dominates(tr_t_label, tr_s_label, lts_id)
	    + label_dominance.get_label_cost(tr_s_label)
	    - label_dominance.get_label_cost(tr_t_label)
	    + q_simulates(tr_t_target, tr_s_target);
    } else {
	return std::numeric_limits<int>::lowest();
    }         
}

template <typename T>
T NumericSimulationRelation<T>::compare_noop(int lts_id, int tr_s_target, int tr_s_label,
					    int t,  T tau_distance, 
					    const NumericLabelRelation<T> & label_dominance) const {

    // Checking noop 
    if(may_simulate(t, tr_s_target) &&
       label_dominance.may_dominated_by_noop(tr_s_label, lts_id)) {
	return tau_distance + 
	    q_simulates(t, tr_s_target) +
	    label_dominance.get_label_cost(tr_s_label) +
	    label_dominance.q_dominated_by_noop(tr_s_label, lts_id);
    } else {
	return std::numeric_limits<int>::lowest();
    }    
}				



template <typename T> void NumericSimulationRelation<T>::
cancel_simulation_computation(int lts_id, const LabelledTransitionSystem * lts) {


    const auto & tau_distances = tau_labels->get_tau_distances(lts_id);
    int new_tau_distances_id = tau_distances.get_id();
    if(new_tau_distances_id != tau_distances_id || !cancelled) { 
    	const auto & tau_distances = tau_labels->get_tau_distances(lts_id);
	cancelled = true;
	for (int s = 0; s < lts->size(); s++) {
	    for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
		update_value(t, s,  tau_distances.minus_shortest_path (t, s));
	    }
	}
    }
}


template <typename T>
vector<int> NumericSimulationRelation<T>::get_dangerous_labels(const LabelledTransitionSystem * lts) const {
    vector<int> dangerous_labels;
    
    int num_states = lts->size();
    vector<bool> is_state_to_check(num_states);
    vector<bool> is_ok(num_states);
    for(LabelGroup lg(0); lg.group < lts->get_num_label_groups(); ++lg) {
	std::fill(is_ok.begin(), is_ok.end(), false);
	std::fill(is_state_to_check.begin(), is_state_to_check.end(), false);
	const std::vector<TSTransition> & trs = lts->get_transitions_label_group(lg);
	std::vector<int> states_to_check;

	for (const auto & tr : trs) {
	    int s = tr.src;
	    int t = tr.target;
	    if (is_ok[s]) {
		continue;
	    } else if(may_simulate(t, s) && q_simulates (t, s) >= 0){
		is_ok[s] = true;
	    } else if (!is_state_to_check[s]) {
		states_to_check.push_back (s);
		is_state_to_check[s] = true;
	    }
	}

	for (int s : states_to_check) {
	    if (!is_ok[s]) {
		// Add dangerous labels
		const std::vector<int> & labels = lts->get_labels(lg);
		dangerous_labels.insert(dangerous_labels.end(), labels.begin(), labels.end());
		break;
	    }
	}
    }
    //for (int i  : dangerous_labels) cout << g_operators[i].get_name() << endl;
    return dangerous_labels;   
}



template <typename T>
int NumericSimulationRelation<T>::update_pair_stable (int lts_id, const LabelledTransitionSystem * lts,
					       const NumericLabelRelation<T> & label_dominance,
					       const TauDistances<T> & tau_distances,
					       int s, int t) {
    assert (s != t // && !is_relation_stable[s][t]
            && may_simulate(t, s)) ;
    
    T lower_bound = tau_distances.minus_shortest_path(t,s);
    T previous_value = q_simulates(t, s);
    // cout << "prev: " << previous_value << endl;
		    
    assert(lower_bound <= previous_value);
    if(lower_bound == previous_value) {
	return false;
    } 

    T min_value = previous_value;
    // if (lts->is_goal(t) || !lts->is_goal(s)) {

    //Check if really t simulates s
    //for each transition s--l->s':
    // a) with noop t >= s' and l dominated by noop?
    // b) exist t--l'-->t', t' >= s' and l dominated by l'?
    lts->applyPostSrc(s, [&](const LTSTransition & trs) {
	    for(int tr_s_label : lts->get_labels(trs.label_group)) {
		T max_value = std::numeric_limits<int>::lowest();
		for (int t2 : tau_distances.states_reachable_from(t)) {
		    T tau_distance = tau_distances.minus_shortest_path(t, t2);			    

		    max_value = max(max_value, 
				    compare_noop(lts_id, trs.target, tr_s_label, t2, tau_distance, label_dominance));

		    if (max_value >= min_value) {
			continue; // Go to next transition
		    }

		    // is_relation_stable[s][t] = lts->applyPostSrc(t2,[&](const LTSTransition & trt) {
		    // 	    if(is_relation_stable[trs.target][trt.target] || (trs.target == s && trt.target == t)) {
		    // 		for(int tr_t_label : lts->get_labels(trt.label_group)) {
		    // 		    max_value_stable = max(max_value_stable, compare_transitions(lts_id, trs.target, tr_s_label, trt.target, tr_t_label, tau_distance, label_dominance));
		    // 		    if (max_value_stable == max_value) {
		    // 			//break if we have found a transition that simulates with the best result possible
		    // 			return true;
		    // 		    }
		    // 		}
		    // 	    }
		    // 	    return false;
		    // 	});

		    lts->applyPostSrc(t2,[&](const LTSTransition & trt) {
			    for(int tr_t_label : lts->get_labels(trt.label_group)) {
				max_value = max(max_value, compare_transitions(lts_id, trs.target, tr_s_label, trt.target, tr_t_label, tau_distance, label_dominance));
				if (max_value >= min_value) {
				    //break if we have found a transition that simulates with the best result possible
				    return true;
				}
			    }
			    return false;
			});

		    if (max_value >= min_value) {
			break;
		    }
		}

		/*if(min_value > std::numeric_limits<int>::lowest() && max_value == std::numeric_limits<int>::lowest()) {
		  std::cout << lts->name(t) << " does not simulate "
		  << lts->name(s) << " because of "
		  << lts->name(trs.src) << " => "
		  << lts->name(trs.target) << " ("
		  << tr_s_label << ")"; // << std::endl;
		  std::cout << "  Simulates? "
		  << q_simulates(t, trs.target);
		  std::cout << "  domnoop? "
		  << label_dominance.may_dominated_by_noop(tr_s_label, lts_id) << "   " << endl;
		  // label_dominance.dump(trs.label);
		  // for (auto trt : lts->get_transitions(t)) {
		  // std::cout << "Tried with: "
		  // 	  << lts->name(trt.src) << " => "
		  // 	  << lts->name(trt.target) << " ("
		  // 	  << trt.label << ")" << " label dom: "
		  // 	  << label_dominance.q_dominates(trt.label,
		  // 					 trs.label, lts_id)
		  // 	  << " target sim "
		  // 	  << q_simulates(trt.target, trs.target)
		  // 	  << std::endl;
		  // }
		  }*/
		min_value  = std::min(min_value, max_value);
		if (min_value <= lower_bound) {
		    return true;
		}
	    }
	    return false;
	});
    assert(min_value < std::numeric_limits<int>::max());

    min_value = std::max(min_value, lower_bound);

    // } else {
    // 	min_value = lower_bound;
    // }

    assert(min_value <= previous_value);
    
    
    if(min_value < previous_value) {
	//cout << "Updating " << lts->get_names()[s] << " <= " << lts->get_names()[t]
	// << " with " << min_value << " before " << previous_value << endl;


	update_value(t, s,  min_value);
	return true;
    }
    return false;
}

template <typename T>
int NumericSimulationRelation<T>::update_pair (int lts_id, const LabelledTransitionSystem * lts,
					       const NumericLabelRelation<T> & label_dominance,
					       const TauDistances<T> & tau_distances,
					       int s, int t) {
    assert (s != t && may_simulate(t, s)) ;
    
    T lower_bound = tau_distances.minus_shortest_path(t,s);

    T previous_value = q_simulates(t, s);
    
    // cout << "prev: " << previous_value << endl;
		    
    assert(lower_bound <= previous_value);
    if(lower_bound == previous_value) {
	return false;
    } 

    T min_value = previous_value;
    // if (lts->is_goal(t) || !lts->is_goal(s)) {

    //Check if really t simulates s
    //for each transition s--l->s':
    // a) with noop t >= s' and l dominated by noop?
    // b) exist t--l'-->t', t' >= s' and l dominated by l'?
    lts->applyPostSrc(s, [&](const LTSTransition & trs) {
	    for(int tr_s_label : lts->get_labels(trs.label_group)) {
		T max_value = std::numeric_limits<int>::lowest();
		for (int t2 : tau_distances.states_reachable_from(t)) {
		    T tau_distance = tau_distances.minus_shortest_path(t, t2);
                    
		    max_value = max(max_value, 
				    compare_noop(lts_id, trs.target, tr_s_label, t2, tau_distance, label_dominance));

		    if (max_value >= min_value) {
			continue; // Go to next transition
		    }
			       
		    lts->applyPostSrc(t2,[&](const LTSTransition & trt) {
			    for(int tr_t_label : lts->get_labels(trt.label_group)) {
				max_value = max(max_value, compare_transitions(lts_id, trs.target, tr_s_label, trt.target, tr_t_label, tau_distance, label_dominance));
				if (max_value >= min_value) {
				    //break if we have found a transition that simulates with the best result possible
				    return true;
				}
			    }
			    return false;
			});

		    if (max_value >= min_value) {
			break;
		    }
		}

		/*if(min_value > std::numeric_limits<int>::lowest() && max_value == std::numeric_limits<int>::lowest()) {
		  std::cout << lts->name(t) << " does not simulate "
		  << lts->name(s) << " because of "
		  << lts->name(trs.src) << " => "
		  << lts->name(trs.target) << " ("
		  << tr_s_label << ")"; // << std::endl;
		  std::cout << "  Simulates? "
		  << q_simulates(t, trs.target);
		  std::cout << "  domnoop? "
		  << label_dominance.may_dominated_by_noop(tr_s_label, lts_id) << "   " << endl;
		  // label_dominance.dump(trs.label);
		  // for (auto trt : lts->get_transitions(t)) {
		  // std::cout << "Tried with: "
		  // 	  << lts->name(trt.src) << " => "
		  // 	  << lts->name(trt.target) << " ("
		  // 	  << trt.label << ")" << " label dom: "
		  // 	  << label_dominance.q_dominates(trt.label,
		  // 					 trs.label, lts_id)
		  // 	  << " target sim "
		  // 	  << q_simulates(trt.target, trs.target)
		  // 	  << std::endl;
		  // }
		  }*/
		min_value  = std::min(min_value, max_value);

                if(min_value < -truncate_value) {
                    min_value = lower_bound;
                    return true;
                } else if (min_value <= lower_bound) {
		    return true;
		}
	    }
	    return false;
	});
    assert(min_value < std::numeric_limits<int>::max());

    min_value = std::max(min_value, lower_bound);

    // } else {
    // 	min_value = lower_bound;
    // }

    assert(min_value <= previous_value);
    
    
    if(min_value < previous_value) {
	//cout << "Updating " << lts->get_names()[s] << " <= " << lts->get_names()[t]
	// << " with " << min_value << " before " << previous_value << endl;

	update_value(t, s,  min_value);
	return true;
    }
    return false;
}


template <typename T>
int NumericSimulationRelation<T>::update (int lts_id, const LabelledTransitionSystem * lts,
					  const NumericLabelRelation<T> & label_dominance,
					  int max_time) {
    if (cancelled){
	cancel_simulation_computation(lts_id, lts); //check that tau-labels have not changed
	return 0;
    }
    
    assert(tau_labels);


    const auto & tau_distances = tau_labels->get_tau_distances(lts_id);
   
    int new_tau_distances_id = tau_distances.get_id();
    if(new_tau_distances_id != tau_distances_id) { //recompute_goal_respecting
	tau_distances_id = new_tau_distances_id;
	for (int s = 0; s < lts->size(); s++) {
            for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
		if (!lts->is_goal(t) && lts->is_goal(s)) {
		    if (tau_distances.get_goal_distance(t) == std::numeric_limits<int>::max()) {
			update_value(t, s, std::numeric_limits<int>::lowest());
		    } else {
			update_value(t, s, min(relation[s][t], -tau_distances.get_goal_distance(t)));
		    }
		}

	    }
	}
    }

  

    Timer timer;

    int num_iterations = 0;
    bool changes = true;    
    while (changes) {
	num_iterations ++;	
        changes = false;
        for (int s = 0; s < lts->size(); s++) {
            for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
		if (timer() > max_time) {
		    cout << "Computation of numeric simulation on LTS " << lts_id
			 << " with " << lts->size()
			 << " states cancelled after " << timer() << " seconds." << endl;
		    
		    cancel_simulation_computation(lts_id, lts);
		    return num_iterations;
		}

		// cout << "s: " << s << " t: " << t << endl;
		if (s != t  && may_simulate(t, s)) { // && !is_relation_stable[s][t]
		    changes |= update_pair (lts_id, lts, label_dominance, tau_distances, s, t);
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


template <typename T>
bool NumericSimulationRelation<T>::pruned(const State & state) const {
    return abs->get_abstract_state(state) == -1;
}

template <typename T>
int NumericSimulationRelation<T>::get_abstract_state_id (const vector<int> & s) const {
    return  abs->get_abstract_state(s);
}

template <typename T>
int NumericSimulationRelation<T>::get_abstract_state_id (const State & s) const {
    return  abs->get_abstract_state(s);
}

template <typename T>
T NumericSimulationRelation<T>::q_simulates (const State & t, const State & s) const{
    int tid = abs->get_abstract_state(t);
    int sid = abs->get_abstract_state(s);
    if(sid == -1 || tid == -1) {
	return std::numeric_limits<int>::lowest(); 
    }
    return q_simulates(tid, sid);
}

template <typename T> 
T NumericSimulationRelation<T>::q_simulates (const vector<int> & t, const vector<int> & s) const{
    int tid = abs->get_abstract_state(t);
    int sid = abs->get_abstract_state(s);
    if(sid == -1 || tid == -1) {
	return std::numeric_limits<int>::lowest(); 
    }
    return q_simulates(tid, sid);
}


template <typename T>
void NumericSimulationRelation<T>::dump(const vector<string> & names) const{ 
    cout << "SIMREL:" << endl;
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
            if(may_simulate(j, i) && i != j){
		cout << names[i] << " <= " << names[j] << " (" << q_simulates(j, i) << ")" << endl;
            }
        }
    }
}

template <typename T>
void NumericSimulationRelation<T>::dump() const{ 
    cout << "SIMREL:" << endl;
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
	    cout << q_simulates(j, i)  << " ";
	}
	cout << endl;
    }
}



template <typename T>
bool NumericSimulationRelation<T>::has_dominance() const {
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
	    if(i == j) continue;
	    if (relation[i][j] > std::numeric_limits<int>::lowest()) {
		return true;
	    }
	}
    }

    return false;
}


template <typename T>
bool NumericSimulationRelation<T>::has_positive_dominance() const {
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
	    if(i == j) continue;
	    if (relation[i][j] >= 0) {
		return true;
	    }
	}
    }

    return false;
}

template <typename T>
void NumericSimulationRelation<T>::statistics() const {
    map<T, int> values;
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

template <typename T>
void NumericSimulationRelation<T>::precompute_absstate_bdds(SymVariables * vars_){
    vars = vars_;
    abs->getAbsStateBDDs(vars, abs_bdds);
}

template <typename T>
void NumericSimulationRelation<T>::precompute_bdds(bool dominating, bool quantified, bool use_ADD) {
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
	cerr << "use_ADD not supported yet" << endl;
	exit(0); 
	// if(dominating) {
	//     dominating_adds.resize(abs->size(), vars->getADD(0));
	//     may_dominating_bdds.resize(abs->size(), vars->zeroBDD());
	// }else {
	//     dominated_adds.resize(abs->size(), vars->getADD(0));
	//     may_dominated_bdds.resize(abs->size(), vars->zeroBDD());
	// }
	
	// for(int i = 0; i < relation.size(); ++i){
	//     for(int j = 0; j < relation.size(); ++j){
	// 	T value = may_simulate(i, j) ? q_simulates (i, j) : std::numeric_limits<int>::lowest(); 

	// 	ADD val = vars->getADD(value);
	// 	if(dominating) {
	// 	    may_dominating_bdds[j] += abs_bdds[i];
	// 	    dominating_adds[j] += abs_bdds[i].Add()*val;
	// 	} else { 
	// 	    may_dominated_bdds[i] +=  abs_bdds[j];
	// 	    dominated_adds[i] += abs_bdds[j].Add()*val;
	// 	}
	//     }
	// }
	
    } else {
	if(dominating) {
	    //cout << "Precomputing dominating_bdd_maps ... " << endl;

	    dominating_bdd_maps.resize(abs->size());
	    may_dominating_bdds.resize(abs->size(), vars->zeroBDD()); 
	}else {
	    //cout << "Precomputing dominated_bdd_maps ... " << endl;

	    dominated_bdd_maps.resize(abs->size());
	    may_dominated_bdds.resize(abs->size(), vars->zeroBDD());
	}
	
	for(int i = 0; i < relation.size(); ++i){
	    for(int j = 0; j < relation.size(); ++j){
		if(may_simulate(i, j)) {
		    T val = q_simulates (i, j);
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

template <typename T>
BDD NumericSimulationRelation<T>::getSimulatedBDD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return dominated_bdds[absstate];
}

template <typename T>
BDD NumericSimulationRelation<T>::getSimulatingBDD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return dominating_bdds[absstate];
}


template <typename T>
const map<T, BDD> & NumericSimulationRelation<T>::getSimulatedBDDMap(const State & state) const{
    assert(!dominated_bdd_maps.empty());
    int absstate = abs->get_abstract_state(state);
    assert(absstate != -1); 
    return dominated_bdd_maps[absstate];
}

template <typename T>
const map<T, BDD> & NumericSimulationRelation<T>::getSimulatingBDDMap(const State & state) const{
    assert(!dominated_bdd_maps.empty());
    int absstate = abs->get_abstract_state(state);
    assert(absstate == -1);
    return dominating_bdd_maps[absstate];
}

template <typename T>
BDD NumericSimulationRelation<T>::getMaySimulatedBDD(const State & state) const{
    assert(!may_dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return may_dominated_bdds[absstate];
}

template <typename T>
BDD NumericSimulationRelation<T>::getMaySimulatingBDD(const State & state) const{
    assert(!may_dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->zeroBDD();
    else return may_dominating_bdds[absstate];
}

template <typename T>
ADD NumericSimulationRelation<T>::getSimulatedADD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->getADD(std::numeric_limits<int>::lowest());
    else return dominated_adds[absstate];
}

template <typename T>
ADD NumericSimulationRelation<T>::getSimulatingADD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return vars->getADD(std::numeric_limits<int>::lowest());
    else return dominating_adds[absstate];
}


template class NumericSimulationRelation<int>; 
template class NumericSimulationRelation<IntEpsilon>; 
