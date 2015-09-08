#include "simulation_simple.h"

#include <cstdlib>
#include <iostream>

#include "simulation_relation.h"
#include "labelled_transition_system.h"
#include "lts_complex.h"
#include "label_relation.h"
#include "abstraction.h"


using namespace std;


std::unique_ptr<SimulationRelation> ComputeSimulationRelationSimple::init_simulation (Abstraction * _abs){
    std::unique_ptr<SimulationRelation> res (new SimulationRelation(_abs));
    res->init_goal_respecting();
    return std::move (res);
}

std::unique_ptr<SimulationRelation> 
ComputeSimulationRelationSimple::init_simulation_incremental (CompositeAbstraction * _abs, 
			     const SimulationRelation & simrel_one, 
			     const SimulationRelation & simrel_two){
    
    std::unique_ptr<SimulationRelation> res (new SimulationRelation(_abs));
    res->init_incremental(_abs, simrel_one, simrel_two);
    return std::move (res);
}


/*
 * THIS IMPLEMENTATION IS VERY INNEFICIENT
 * ONLY TO BE USED AS A PROOF OF CONCEPT
 */
template<typename LTS>
void ComputeSimulationRelationSimple::update_sim(int lts_id, const LTS * lts,
						 const LabelRelation & label_dominance, 
						 SimulationRelation & simrel) {
    bool changes = true;
    while (changes) {
        //std::cout << "looping" << std::endl;
        changes = false;
        for (int s = 0; s < lts->size(); s++) {
            for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
		if (s != t && simrel.simulates(t, s) && !simrel.fixed_simulates(t, s)) {
                    //cout << "Checking states " << lts->name(s) << " and " << lts->name(t) << endl;
                    //Check if really t simulates s
                    //for each transition s--l->s':
                    // a) with noop t >= s' and l dominated by noop?
                    // b) exist t--l'-->t', t' >= s' and l dominated by l'?
                    lts->applyPostSrc(s, [&](const LTSTransition & trs) {
                        //cout << "Checking transition " << trs.label << " " << g_operators[trs.label].get_name() << " to " << trs.target << endl;
                            if(simrel.simulates (t, trs.target) &&
                                    label_dominance.dominated_by_noop(trs.label, lts_id)) {
                                //cout << "Dominated by noop!" << endl;
                                return false;
                            }
                            bool found =
                            lts->applyPostSrc(t,[&](const LTSTransition & trt) {
                                        if(label_dominance.dominates(trt.label, trs.label, lts_id) &&
                                                simrel.simulates(trt.target, trs.target)) {
                                            return true;
                                        }
                                        return false;
                                    });

                            if(!found) {
                                changes = true;
                                simrel.remove(t, s);
                                /*std::cout << lts->name(t) << " does not simulate "
                                 << lts->name(s) << " because of "
                                 << lts->name(trs.src) << " => "
                                 << lts->name(trs.target) << " ("
                                 << trs.label << ")"; // << std::endl;
                                 std::cout << "  Simulates? "
                                 << simulates(trs.src, trs.target);
                                 std::cout << "  domnoop? "
                                 << label_dominance.dominated_by_noop(
                                 trs.label, lts_id) << "   ";
                                 label_dominance.dump(trs.label);*/
                                /*for (auto trt : lts->get_transitions(t)) {
                                 std::cout << "Tried with: "
                                 << lts->name(trt.src) << " => "
                                 << lts->name(trt.target) << " ("
                                 << trt.label << ")" << " label dom: "
                                 << label_dominance.dominates(trt.label,
                                 trs.label, lts_id)
                                 << " target sim "
                                 << simulates(trt.target, trs.target)
                                 << std::endl;
                                 }*/
                                return true;
                            }
                            return false;
                        });
                }
            }
        }
    }
}


//l does not dominate l2 anymore, check if this changes the simulation relation
bool ComputeSimulationRelationSimple::propagate_label_domination(int lts_id, 
								 const LabelledTransitionSystem * lts,
								 const LabelRelation & label_dominance, 
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
