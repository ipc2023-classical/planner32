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

std::unique_ptr<NumericSimulationRelation> NumericDominanceRelation::init_simulation (Abstraction * _abs){
    auto res = make_unique<NumericSimulationRelation> (_abs);
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


bool NumericDominanceRelation::dominates(const State &t, const State & s, int g_diff) const {

    int sum_negatives = 0;
    int max_positive = 0;

    for(auto & sim : simulations) {
	int val = sim->simulates(t, s);
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
    return total_value -g_diff > 0 ||
	(total_value == g_diff && g_diff > 0);
}

