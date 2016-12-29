#include "numeric_dominance_relation.h"

#include "numeric_simulation_relation.h" 
#include "abstraction.h" 
#include "labels.h" 
#include "labelled_transition_system.h" 

using namespace std;


void NumericDominanceRelation::init (const std::vector<Abstraction *> & abstractions){ 
    simulations.clear();
    for (auto abs : abstractions){
	simulations.push_back(init_simulation(abs)); 
    }
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


bool NumericDominanceRelation::dominates(const State &t, const State & s) const {
    for(auto & sim : simulations) {
	if (!sim->simulates(t, s)) {
	    return false;
	}
    }
    return true;
}

