#include "variable_partition_dominance.h"

#include "abstraction.h"
#include "abstraction_builder.h"
#include "ld_simulation.h"
#include "label_relation.h"
#include "labelled_transition_system.h"

#include <map>

using namespace std;

LocalLabelRelation::LocalLabelRelation(const SimulationRelation & atomic_dominance,
                                       std::shared_ptr<std::vector<int>>  label_cost_,
                                       LabelledTransitionSystem & lts) : 
    dom_by_noop (label_cost_->size(), true),
    label_to_equivalence_group (lts.compute_label_equivalence_groups()) {

    for(int l : lts.get_relevant_labels()) {
        for(auto tr : lts.get_transitions_label(l)){
            if (!atomic_dominance.simulates(tr.target, tr.src)) {
                dom_by_noop[l] = false;
                break;
            }
        }
    }
}


OutsideLabelRelation::OutsideLabelRelation(const vector<LocalLabelRelation> & local_dominances,
                                           std::shared_ptr<std::vector<int> > label_cost_,
                                           const std::vector<int> & pattern) :
    dom_by_noop (local_dominances[0].dom_by_noop.size(), true),
    label_cost(label_cost_) {

    vector <int> outside_pattern;
    for(size_t i = 0; i < local_dominances.size(); ++i) {
        if (std::find(pattern.begin(), pattern.end(), i) == pattern.end()) {
            outside_pattern.push_back(i);
        }
    }

    for(int i : outside_pattern) {
        for (size_t j = 0; j < dom_by_noop.size(); ++j) {
            dom_by_noop[j] &= local_dominances[i].dom_by_noop[j];
        }
    }

    label_to_equivalence_group.reserve(label_cost->size());
    
    vector<int> eq_class(outside_pattern.size());
    std::map<vector<int>, int> equivalence_classes;
    int group = 0;
    for (int l = 0; l < label_cost->size(); ++l) {
        for(size_t i = 0; i < outside_pattern.size(); ++i) {
            eq_class[i] = local_dominances[outside_pattern[i]].label_to_equivalence_group[l];
        }

        auto it = equivalence_classes.find(eq_class);
        if(it == equivalence_classes.end()) {
            equivalence_classes[eq_class] = group;
            label_to_equivalence_group.push_back(group);
            group++;
        } else {
            label_to_equivalence_group.push_back(it->second);
        }
    }
}

bool OutsideLabelRelation::dominates(int label, int label2) const {
    assert (label >= 0);
    assert (label2 >= 0);
    assert (label < label_to_equivalence_group.size());
    assert (label < label_cost->size());
    assert (label2 < label_to_equivalence_group.size());
    assert (label2 < label_cost->size());
    
    return label_to_equivalence_group[label] == label_to_equivalence_group[label2] &&
        (*label_cost)[label] <= (*label_cost)[label2];
}




VariablePartitionDominance::VariablePartitionDominance(Labels * labels,
                                                       const OutsideLabelRelation & label_dominance, 
                                                       const std::vector<int> & pattern) :
    abstraction (new PDBAbstraction(labels, pattern)) {
    abstraction->normalize();
    abstraction->compute_distances();

    LabelMap label_map (labels);
    LabelledTransitionSystem * lts = abstraction->get_lts(label_map);
    assert(lts);
   
    simulation_relation = std::make_unique<SimulationRelation> (abstraction.get());
    simulation_relation->init_goal_respecting();
    bool changes = true;
    while (changes) {
        changes = false;
        for (int s = 0; s < lts->size(); s++) {
            for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
                assert (simulation_relation);
		if (s != t && simulation_relation->simulates(t, s)) {
                    lts->applyPostSrc(s, [&](const LTSTransition & trs) {
                            if(simulation_relation->simulates (t, trs.target) &&
                               label_dominance.dominated_by_noop(trs.label)) {
                                return false;
                            }
                            bool found =
                                lts->applyPostSrc(t,[&](const LTSTransition & trt) {
                                        if(label_dominance.dominates(trt.label, trs.label) &&
                                           simulation_relation->simulates(trt.target, trs.target)) {
                                            return true;
                                        }
                                        return false;
                                    });

                            if(!found) {
                                changes = true;
                                simulation_relation->remove(t, s);

                                return true;
                            }
                            return false;
                        });
                }
            }
        }
    }
    
    abstraction->release_memory();    
}


PDBAbstraction & VariablePartitionDominance::get_abstraction() {
    assert((bool)abstraction);
    return *abstraction;
}

// const std::vector<int> & VariablePartitionDominance::get_dominated_states(const State & state) const {
//     return simulation_relation->get_dominated_states(state);
// }

// const std::vector<int> & VariablePartitionDominance::get_dominated_states(int s) const {
//     return simulation_relation->get_dominated_states(s);
// }



std::unique_ptr<VariablePartitionDominance>
VariablePartitionDominanceFactory::compute_dominance (const std::vector<int> & pattern) {

    cout << "Compute dominance for pattern: ";
    for (int v : pattern) {
        cout << v << " ";
    }
    cout << endl;
    assert(!atomic_label_relations.empty());
    OutsideLabelRelation outside_dominance (atomic_label_relations, label_cost, pattern);

    cout << "OutsideLabelRelation computed" << endl;
    return std::make_unique<VariablePartitionDominance> (atomic_dominance->get_labels(),
                                                         outside_dominance, pattern);
}

VariablePartitionDominanceFactory::VariablePartitionDominanceFactory
( const Options &opts) : builder_atomic(new AbsBuilderAtomic(opts)) {

}



void VariablePartitionDominanceFactory::init (bool unit_cost, OperatorCost cost_type) {
    std::vector<std::unique_ptr<Abstraction> > tmp; //tmp should be atomic abstractions but this requires major changes to ld_simulation and abstraction_builder.
    builder_atomic->build_abstraction (unit_cost, cost_type, atomic_dominance, tmp);
    const auto & atomic_abstractions = atomic_dominance->get_abstractions();
    cout << atomic_abstractions.size() << " atomic abstractions" << endl; 
    atomic_label_relations.reserve(atomic_abstractions.size());
    atomic_dominance->compute_final_simulation(SimulationType::SIMPLE,
                                               LabelDominanceType::NOOP,
                                               100000, false, false, 
                                               false, 
                                               false, false, true, 100000);
    Labels * labels = atomic_dominance->get_labels();
    LabelMap label_map (labels);
    label_cost = make_shared<vector<int>> (label_map.get_num_labels());
    for(size_t i = 0; i < label_cost->size(); ++i) {
        int old_label_id = label_map.get_old_id(i);
        (*label_cost)[i] = labels->get_label_by_index(old_label_id)->get_cost();
    }


    for(size_t i = 0; i < atomic_abstractions.size(); ++i) {

        LabelledTransitionSystem * lts = atomic_abstractions[i]->get_lts(label_map);

        assert (label_cost);
        assert (!label_cost->empty());
        assert (i < atomic_dominance->get_dominance_relation().size());

        atomic_label_relations.push_back(LocalLabelRelation(atomic_dominance->get_dominance_relation()[i],
                                                            label_cost, *lts));

    }

    atomic_dominance->release_memory();

    builder_atomic.reset();
    cout << "VariablePartitionDominanceFactory::init done " << endl;
}


void VariablePartitionDominanceFactory::add_options_to_parser(OptionParser &parser) {
    AbsBuilderAtomic::add_options_to_parser(parser);
}

