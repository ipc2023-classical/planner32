#include "numeric_dominance_relation.h"

#include "numeric_simulation_relation.h" 
#include "../merge_and_shrink/abstraction.h" 
#include "../merge_and_shrink/labels.h" 
#include "../merge_and_shrink/labelled_transition_system.h" 

using namespace std;


template <typename T> 
void NumericDominanceRelation<T>::init (const std::vector<Abstraction *> & abstractions){ 
    simulations.clear();
    for (auto abs : abstractions){
	simulations.push_back(init_simulation(abs)); 
    }
}

template <typename T> 
std::unique_ptr<NumericSimulationRelation<T>> 
NumericDominanceRelation<T>::init_simulation (Abstraction * _abs){
    auto res = make_unique<NumericSimulationRelation<T>> (_abs, truncate_value);
    res->init_goal_respecting();
    return std::move (res);
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
	T val = sim->q_simulates(t, s);
	if(val == std::numeric_limits<int>::lowest()) {
	    return false;
	}
	total_value += val;
    }
    
    return total_value - g_diff >= 0;

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
BDD NumericDominanceRelation<T>::getDominatedBDD(SymVariables * vars, const State &state ) const{
    BDD res = vars->oneBDD();
    try{
        for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
            res *= (*it)->getSimulatedBDD(state);
        }
    }catch(BDDError e){
        return vars->zeroBDD();
    }
    return res;
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
							      const State &state ) const{
    map<T, BDD> res; 
    res[T(0)] = vars->oneBDD();

    for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
	const auto & sim_bdd_map = (*it)->getSimulatedBDDMap(state);
	map<T, BDD> new_res; 

	for(const auto & entry : sim_bdd_map) {
	    for(const auto & entry2 : res) {
		try{
		    if(entry.first != std::numeric_limits<int>::lowest()) {
			T value = entry.first + entry2.first;
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


template class NumericDominanceRelation<int>; 
template class NumericDominanceRelation<IntEpsilon>; 
