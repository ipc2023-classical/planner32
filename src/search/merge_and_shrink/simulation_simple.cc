#include "simulation_simple.h"

#include "labelled_transition_system.h"
#include "lts_efficient.h"

#include "label_relation.h"

#include <cstdlib>
#include <iostream>
#include "abstraction.h"

using namespace std;

SimulationRelationSimple::SimulationRelationSimple(Abstraction * _abs) :
        SimulationRelation(_abs) {
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
                        && (fixed_relation.empty() || !fixed_relation[s][t])) {
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

