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


int NumericDominanceRelation::get_cost(const State &state) const{
    int cost = 0;
    for(auto & sim : simulations) {
	int new_cost = sim->get_cost(state);
	if (new_cost == -1) return -1;
	cost = max (cost, new_cost);
    }
    return cost;
}


bool NumericDominanceRelation::dominates(const State &t, const State & s) const {
    for(auto & sim : simulations) {
	if (!sim->simulates(t, s)) {
	    return false;
	}
    }
    return true;
}


//l does not dominate l2 anymore, check if this changes the simulation relation
template <typename LR> 
bool NumericDominanceRelationLR<LR>::propagate_label_domination(int lts_id, 
						       const LabelledTransitionSystem * lts,
						       int l, int l2, 
						       SimulationRelation & simrel) const {
    for (int s = 0; s < lts->size(); s++) {
	for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
	    if (s != t && simrel.simulates(t, s)) {
		//Check if really t simulates s //for each transition s--l->s':
		// a) with noop t >= s' and l dominated by noop?
		// b) exist t--l'-->t', t' >= s' and l dominated by l'?
		bool not_simulates_anymore = lts->applyPostSrc(s, [&](const LTSTransition & trs) {
			if(trs.label != l2) return false;

			if(simrel.simulates (t, trs.target) &&
			   label_dominance.dominated_by_noop(trs.label, lts_id)) {
			    //cout << "Dominated by noop!" << endl;
			    return false;
			}
			bool found =
                            lts->applyPostSrc(t,[&](const LTSTransition & trt) {
				    if (trt.label == l) return false;
				    if(label_dominance.dominates(trt.label, trs.label, lts_id) &&
				       simrel.simulates(trt.target, trs.target)) {
					return true;
				    }
				    return false;
				});
			
			return !found;
		    });

		if(not_simulates_anymore) return false;
	    }
	}
    }
    return true;    
}
