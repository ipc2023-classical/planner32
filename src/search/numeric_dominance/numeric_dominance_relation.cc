#include "numeric_dominance_relation.h"

#include "numeric_simulation_relation.h" 
#include "../merge_and_shrink/abstraction.h" 
#include "../merge_and_shrink/labels.h" 
#include "../merge_and_shrink/labelled_transition_system.h" 

using namespace std;


void NumericDominanceRelation::init (const std::vector<Abstraction *> & abstractions){ 
    simulations.clear();
    for (auto abs : abstractions){
	simulations.push_back(init_simulation(abs)); 
    }
}

std::unique_ptr<NumericSimulationRelation> 
NumericDominanceRelation::init_simulation (Abstraction * _abs){
    auto res = make_unique<NumericSimulationRelation> (_abs, truncate_value);
    res->init_goal_respecting();
    return std::move (res);
}


bool NumericDominanceRelation::pruned_state(const State &state) const {
    for(auto & sim : simulations) {
        if(sim->pruned(state)){
            return true;
        }
    }
    return false;
}


// int NumericDominanceRelation::get_cost(const State &state) const{
//     int cost = 0;
//     for(auto & sim : simulations) {
// 	int new_cost = sim->get_cost(state);
// 	if (new_cost == -1) return -1;
// 	cost = max (cost, new_cost);
//     }
//     return cost;
// }



// bool NumericDominanceRelation::parent_dominates_successor(const State & parent, 
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
bool NumericDominanceRelation::dominates(const State &t, const State & s, int g_diff) const {

    int total_value = 0;
    
    // int sum_negatives = 0;
    // int max_positive = 0;

    for(const auto & sim : simulations) {
	int val = sim->q_simulates(t, s);
	// cout << val << endl;
	if(val == std::numeric_limits<int>::lowest()) {
	    return false;
	}
	total_value += val;
	// if(val < 0) {
	//     sum_negatives += val;
	// } else {
	//     max_positive = std::max(max_positive, val);
	// }
    }
    //int total_value = sum_negatives + max_positive;

    // if(total_value -g_diff > 0 ||
    //    (total_value == g_diff && g_diff > 0)) {
    // 	cout << "Find domination with: " << total_value << " action_cost: " << g_diff << endl;
    // 	for (int var = 0; var < g_variable_domain.size(); ++var) {
    // 	    if(s[var] != t[var]) {
    // 		cout <<  "   " << g_fact_names[var][s[var]] << " -> " << g_fact_names[var][t[var]] << endl;
    // 	    }
    // 	}	     
    // }
    // cout << "Prune? " << total_value << " " << g_diff << ": " <<
    // 	(total_value -g_diff > 0 || (total_value == g_diff && g_diff > 0)) << endl;
    
    return total_value - g_diff > 0 ||
	(total_value == g_diff && g_diff > 0);
}



bool NumericDominanceRelation::dominates_parent(const vector<int> & state, const vector<int> & parent, int action_cost) const {
    int sum_negatives = 0;
    int max_positive = 0;

    for(const auto & sim : simulations) {
	int val = sim->q_simulates(state, parent);

	if(val == std::numeric_limits<int>::lowest()) {
	    return false;
	}
	if(val < 0) {
	    sum_negatives += val;
	} else {
	    max_positive = std::max(max_positive, val);
	}
    }    
    int total_value = sum_negatives + max_positive;
    
    return total_value - action_cost > 0 || (total_value == action_cost && action_cost > 0);
}

