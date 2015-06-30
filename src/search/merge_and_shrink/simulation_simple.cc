#include "simulation_simple.h"

#include "labelled_transition_system.h"
#include "lts_complex.h"

#include "label_relation.h"

#include <cstdlib>
#include <iostream>
#include "abstraction.h"

using namespace std;

SimulationRelationSimple::SimulationRelationSimple(Abstraction * _abs,
        bool incremental) :
        SimulationRelation(_abs) {
    if (incremental) {
        CompositeAbstraction* comp_abs = (CompositeAbstraction*) _abs;
        const Abstraction* abs_one = comp_abs->get_component(0);
        //const vector <bool> & goal_states_one = abs_one->get_goal_states();
        const SimulationRelation* simrel_one = abs_one->get_simulation_relation();
        assert(simrel_one);
        int num_one = simrel_one->num_states();
        const Abstraction* abs_two = comp_abs->get_component(1);
        //const vector <bool> & goal_states_two = abs_two->get_goal_states();
        const SimulationRelation* simrel_two = abs_two->get_simulation_relation();
        assert(simrel_two);
        int num_two = simrel_two->num_states();

        int num_states = relation.size();
        fixed_relation.resize(num_states);
        for (int i = 0; i < num_states; i++) {
            fixed_relation[i].resize(num_states, false);
        }

        for (int i = 0; i < num_one; i++) {
            for (int j = i + 1; j < num_one; j++) {
                if (simrel_one->simulates(i, j)) {
                    for (int x = 0; x < num_two; x++) {
                        int ip = comp_abs->get_abstract_state(i, x);
                        for (int y = 0; y < num_two; y++) {
                            if (simrel_two->simulates(x, y)) {
                                int jp = comp_abs->get_abstract_state(j, y);
                                fixed_relation[ip][jp] = true;
                            }
                        }
                    }
                }
                if (simrel_one->simulates(j, i)) {
                    for (int x = 0; x < num_two; x++) {
                        int ip = comp_abs->get_abstract_state(i, x);
                        for (int y = 0; y < num_two; y++) {
                            if (simrel_two->simulates(x, y)) {
                                int jp = comp_abs->get_abstract_state(j, y);
                                fixed_relation[jp][ip] = true;
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
 * THIS IMPLEMENTATION IS VERY INNEFICIENT
 * ONLY TO BE USED AS A PROOF OF CONCEPT
 */
template<typename LTS>
void SimulationRelationSimple::update_sim(int lts_id, const LTS * lts,
        const LabelRelation & label_dominance) {
    bool changes = true;
    while (changes) {
        //std::cout << "looping" << std::endl;
        changes = false;
        for (int s = 0; s < lts->size(); s++) {
            for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
                if (s != t && simulates(t, s)
                        && (fixed_relation.empty() || !fixed_relation[t][s]/*PIET-edit: was (likely incorrect): !fixed_relation[s][t]*/)) {
                    //cout << "Checking states " << lts->name(s) << " and " << lts->name(t) << endl;
                    //Check if really t simulates s
                    //for each transition s--l->s':
                    // a) with noop t >= s' and l dominated by noop?
                    // b) exist t--l'-->t', t' >= s' and l dominated by l'?
                    lts->applyPostSrc(s, [&](const LTSTransition & trs) {
                        //cout << "Checking transition " << trs.label << " " << g_operators[trs.label].get_name() << " to " << trs.target << endl;
                            if(simulates (t, trs.target) &&
                                    label_dominance.dominated_by_noop(trs.label, lts_id)) {
                                //cout << "Dominated by noop!" << endl;
                                return false;
                            }
                            bool found =
                            lts->applyPostSrc(t,[&](const LTSTransition & trt) {
                                        if(label_dominance.dominates(trt.label, trs.label, lts_id) &&
                                                simulates(trt.target, trs.target)) {
                                            return true;
                                        }
                                        return false;
                                    });

                            if(!found) {
                                changes = true;
                                remove(t, s);
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
bool SimulationRelationSimple::propagate_label_domination(int lts_id, 
							  const LabelledTransitionSystem * lts,
							  const LabelRelation & label_dominance, 
							  int l, int l2) const {
    for (int s = 0; s < lts->size(); s++) {
	for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
	    if (s != t && simulates(t, s)) {
		//Check if really t simulates s //for each transition s--l->s':
		// a) with noop t >= s' and l dominated by noop?
		// b) exist t--l'-->t', t' >= s' and l dominated by l'?
		bool not_simulates_anymore = lts->applyPostSrc(s, [&](const LTSTransition & trs) {
			if(trs.label != l2) return false;

			if(simulates (t, trs.target) &&
			   label_dominance.dominated_by_noop(trs.label, lts_id)) {
			    //cout << "Dominated by noop!" << endl;
			    return false;
			}
			bool found =
                            lts->applyPostSrc(t,[&](const LTSTransition & trt) {
				    if (trt.label == l) return false;
				    if(label_dominance.dominates(trt.label, trs.label, lts_id) &&
				       simulates(trt.target, trs.target)) {
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
