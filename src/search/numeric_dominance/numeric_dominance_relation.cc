#include "numeric_dominance_relation.h"

#include "numeric_simulation_relation.h" 
#include "../merge_and_shrink/abstraction.h" 
#include "../merge_and_shrink/labels.h" 
#include "../search_progress.h" 
#include "../merge_and_shrink/labelled_transition_system.h" 

using namespace std;



template <typename T> 
void NumericDominanceRelation<T>::init (const std::vector<Abstraction *> & abstractions){ 
    simulations.clear();
    
    simulation_of_variable.resize(g_variable_domain.size(), 0);
    for (int i = 0; i < abstractions.size(); ++i){
	simulations.push_back(init_simulation(abstractions[i])); 
	for (int v : abstractions[i]->get_varset()) {
	    simulation_of_variable[v] = i;
	}
    }
    parent.resize (g_variable_domain.size());
    parent_ids.resize (simulations.size());
    succ.resize (g_variable_domain.size());
    succ_ids.resize (simulations.size());
    initial_state.resize (g_variable_domain.size());
    initial_state_ids.resize (simulations.size());
    values_initial_state_against_parent.resize (simulations.size());

    set_initial_state(g_initial_state());
}

template <typename T> 
std::unique_ptr<NumericSimulationRelation<T>> 
NumericDominanceRelation<T>::init_simulation (Abstraction * _abs){
    auto res = make_unique<NumericSimulationRelation<T>> (_abs, truncate_value, tau_labels);
    res->init_goal_respecting();
    return res;
}

template <typename T> 
bool NumericDominanceRelation<T>::pruned_state(const State &state) const {
    for(auto & sim : simulations) {
        if(sim->pruned(state)){
            return true;
        }
    }
    return false;
}


// int NumericDominanceRelation<T>::get_cost(const State &state) const{
//     int cost = 0;
//     for(auto & sim : simulations) {
// 	int new_cost = sim->get_cost(state);
// 	if (new_cost == -1) return -1;
// 	cost = max (cost, new_cost);
//     }
//     return cost;
// }



// bool NumericDominanceRelation<T>::parent_dominates_successor(const State & parent, 
// 							  const Operator *op) const {
    
//     for(const auto & sim : simulations) {
// 	int val = sim->q_simulates(t, s);
// 	// cout << val << endl;
// 	if(val == std::numeric_limits<int>::lowest()) {
// 	    return false;
// 	}
// 	if(val < 0) {
// 	    sum_negatives += val;
// 	} else {
// 	    max_positive = std::max(max_positive, val);
// 	}

//     }
// }

template <typename T> 
bool NumericDominanceRelation<T>::dominates(const State &t, const State & s, int g_diff) const {
    T total_value = 0;
    for(const auto & sim : simulations) {
	assert(sim);
	T val = sim->q_simulates(t, s);
	if(val == std::numeric_limits<int>::lowest()) {
	    return false;
	}
	total_value += val;
    }
    
    return total_value - g_diff >= 0;

}

template <typename T> 
bool NumericDominanceRelation<T>::strictly_dominates(const State &t, const State & s) const {
    return dominates(t, s, 0) && !dominates(s, t, 0);
}


template <typename T> 
bool NumericDominanceRelation<T>::strictly_dominates_initial_state(const State &t) const {
    vector<int> copy_t (g_variable_domain.size());
    for(int i = 0; i < g_variable_domain.size(); ++i) {
	copy_t[i] = t[i];
    }
    return dominates_parent(copy_t, initial_state, 0) && !dominates_parent(initial_state, copy_t, 0);
}

template <typename T> 
bool NumericDominanceRelation<T>::dominates_parent(const vector<int> & state, 
						   const vector<int> & parent, 
						   int action_cost) const {
    T total_value = 0;
    for(const auto & sim : simulations) {
	T val = sim->q_simulates(state, parent);
	if(val == std::numeric_limits<int>::lowest()) {
	    return false;
	}
	total_value += val;
    }
    
    return total_value - action_cost >= 0;
}

template <typename T> 
void NumericDominanceRelation<T>::precompute_bdds(SymVariables * vars, 
					       bool dominating, bool quantified, bool use_ADD){
    Timer t;
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_bdds(dominating, quantified, use_ADD);
    }
    cout << "Precomputed  BDDs: " << t() << endl;
}


template <typename T> 
BDD NumericDominanceRelation<T>::getDominatedBDD(SymVariables * vars, const State &state, 
						 bool trade_off_dominance) const{

    if(!trade_off_dominance) {
	BDD res = vars->oneBDD();
    
	try{
	    for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
		res *= (*it)->getSimulatedBDD(state);
	    }
	}catch(BDDError e){
	    return vars->zeroBDD();
	}
	return res;
    } else  {
	BDD res = vars->zeroBDD();
	try{
	    map<T, BDD> mapa = getDominatedBDDMap (vars, state, true);
	    for(auto & entry : mapa) {
		assert (entry.first >= 0);
		res += entry.second;
	    }
	}catch(BDDError e){
	}	
	return res;
    }


}



template <typename T> 
BDD NumericDominanceRelation<T>::getDominatingBDD(SymVariables * vars, const State &state ) const{
    BDD res = vars->oneBDD();
    try{
        for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
            res *= (*it)->getSimulatingBDD(state);
        }
    }catch(BDDError e){
        return vars->zeroBDD();
    }

    return res;
}



template<typename T>
map<T, BDD> NumericDominanceRelation<T>::getDominatedBDDMap(SymVariables * vars, 
							    const State &state, 
							    bool only_positive ) const{   
    map<T, BDD> res; 
    res[T(0)] = vars->oneBDD();

    T accumulated_value = total_max_value;
    for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
	accumulated_value -= (*it)->get_max_value();
	const auto & sim_bdd_map = (*it)->getSimulatedBDDMap(state);
	map<T, BDD> new_res; 
	
	for(const auto & entry : sim_bdd_map) {
	    for(const auto & entry2 : res) {
		try{
		    if(entry.first != std::numeric_limits<int>::lowest()) {
			T value = entry.first + entry2.first;
			
			if(only_positive && value + accumulated_value < 0) {
			    continue;
			}

			if (new_res.count(value)) {
			    new_res[value] += entry.second*entry2.second;
			} else {
			    new_res[value] = entry.second*entry2.second;
			} 
		    }
		}catch(BDDError e){
		}
	    }
	}
	new_res.swap(res);
    }

    return res;
}



template <typename T> bool NumericDominanceRelation<T>::
action_selection_pruning(const State & state, std::vector<const Operator *> & applicable_operators,
			 SearchProgress & search_progress, OperatorCost cost_type) const {
    for(int i = 0; i < g_variable_domain.size(); ++i) {
	parent[i] = state[i];
    }
    for(int i = 0; i < simulations.size(); ++i) {
	parent_ids[i] = simulations[i]->get_abstract_state_id(parent);
    } 
    succ = parent;
    for (auto op : applicable_operators) {
	for(const auto & prepost : op->get_pre_post()){
	    succ[prepost.var] = prepost.post;
	    relevant_simulations.insert(simulation_of_variable[prepost.var]);
	}

	T total_value = 0;
	bool may_simulate = true;
	for (int sim :  relevant_simulations) {
	    int succ_id = simulations[sim]->get_abstract_state_id(succ);
	    if(succ_id == -1) {
		may_simulate = false;
		break;
	    }
	    T val = simulations[sim]->q_simulates(succ_id, parent_ids[sim]);
	    if(val == std::numeric_limits<int>::lowest()) {
		may_simulate = false;
		break;
	    }
	    total_value += val;
	}
      	relevant_simulations.clear();

	//TODO: Use adjusted cost instead.
	if(may_simulate && total_value - get_adjusted_action_cost(*op, cost_type) >= 0) {
	    search_progress.inc_action_selection(applicable_operators.size() - 1);
	    applicable_operators.clear();
	    applicable_operators.push_back(op);
	    return true;
	}
	    
	for(const auto & prepost : op->get_pre_post()){
	    succ[prepost.var] = parent[prepost.var];
	}
    }

    return false;
}

template <typename T> void NumericDominanceRelation<T>::
prune_dominated_by_parent_or_initial_state (const State & state, 
							     std::vector<const Operator *> & applicable_operators,
					    SearchProgress & search_progress, bool parent_ids_stored,
					    bool compare_against_parent,
					    bool compare_against_initial_state, OperatorCost cost_type) const {

    if(!parent_ids_stored) {
	for(int i = 0; i < g_variable_domain.size(); ++i) {
	    succ[i] = state[i];
	}
	if(compare_against_parent) {
	    parent = succ;
	for(int i = 0; i < simulations.size(); ++i) {
	    parent_ids[i] = simulations[i]->get_abstract_state_id(parent);
	} 
    }
    }

    vector<int> ts_initial_state_does_not_simulate_parent;
    T initial_state_against_parent = 0;
    if (compare_against_initial_state) {
	for(int i = 0; i < simulations.size(); ++i) {
	    values_initial_state_against_parent[i] =
		simulations[i]->q_simulates(initial_state_ids[i], parent_ids[i]);
	    if(values_initial_state_against_parent[i] == std::numeric_limits<int>::lowest()) {
		ts_initial_state_does_not_simulate_parent.push_back(i);
	    } else {
		initial_state_against_parent += values_initial_state_against_parent[i];
	    }
	}
    }

    int detected_dead_ends = 0;
    int ops_before = applicable_operators.size();
    applicable_operators.erase(std::remove_if(applicable_operators.begin(), 
					      applicable_operators.end(),
					      [&](const Operator * op){
						  for(const auto & prepost : op->get_pre_post()){
						      succ[prepost.var] = prepost.post;
						      relevant_simulations.insert(simulation_of_variable[prepost.var]);
						  }

						  bool proved_prunable = false;

						  //Check dead_ends
						  for (int sim :  relevant_simulations) {
						      succ_ids[sim] = simulations[sim]->get_abstract_state_id(succ);
						      if(succ_ids[sim] == -1) {
							  detected_dead_ends ++;
							  proved_prunable = true;
						      } 
						  }

						  if(!proved_prunable && compare_against_parent) {
						  T total_value = 0;
						  bool may_simulate = true;
						  for (int sim :  relevant_simulations) {
							  T val = simulations[sim]->q_simulates(parent_ids[sim], succ_ids[sim]);

	
							  if(val == std::numeric_limits<int>::lowest()) {
							      may_simulate = false;
							      break; 
							  }
							  total_value += val;
						      } 
						      
						      proved_prunable = may_simulate && (total_value >= 0 || total_value + get_adjusted_action_cost(*op, cost_type) > 0);
						  }

						  if (!proved_prunable && compare_against_initial_state && ts_initial_state_does_not_simulate_parent.size() <= relevant_simulations.size()) {
						      bool all_not_simulated_change = false;
						      for (int sim_must_change :  ts_initial_state_does_not_simulate_parent) {
							  bool found = false;
							  for(int sim : relevant_simulations) {
							      if (sim_must_change == sim) {
								  found = true;
								  break;
							      }
							  }
							  if (!found) {
							      all_not_simulated_change = false;
							      break; //proved no
							  }
						      }

						      if(all_not_simulated_change) {
							  T total_value = initial_state_against_parent;
							  bool may_simulate = true;
							  for (int sim :  relevant_simulations) {					      
							      T val = simulations[sim]->q_simulates(initial_state_ids[sim], succ_ids[sim]);
							  if(val == std::numeric_limits<int>::lowest()) {
							      may_simulate = false;
								  break; 
							  }
							  total_value += val;
							      if(values_initial_state_against_parent[sim] != std::numeric_limits<int>::lowest()) {
								  total_value -= values_initial_state_against_parent[sim];
							      }
							  }
							  proved_prunable = may_simulate && (total_value >= 0 || total_value + get_adjusted_action_cost(*op, cost_type) > 0);
						      }
						  }

						  relevant_simulations.clear();
						  for(const auto & prepost : op->get_pre_post()){
						      succ[prepost.var] = parent[prepost.var];
						  }

						  //TODO: Use adjusted cost instead.
						  return proved_prunable;
					      }),
			       applicable_operators.end());	
    if(ops_before > applicable_operators.size()) {
	search_progress.inc_dead_ends(detected_dead_ends);
	search_progress.inc_pruned((ops_before - applicable_operators.size()) - detected_dead_ends);
	//cout << "Pruned "  << ops_before  -applicable_operators.size() << " out of " << ops_before << endl;
    }
}


template <typename T> void NumericDominanceRelation<T>::
 set_initial_state (const State & state) {
    for(int i = 0; i < g_variable_domain.size(); ++i) {
	initial_state[i] = state[i];
    }
    for(int i = 0; i < simulations.size(); ++i) {
	initial_state_ids[i] = simulations[i]->get_abstract_state_id(initial_state);
    }
}

    
template class NumericDominanceRelation<int>; 
template class NumericDominanceRelation<IntEpsilon>; 
