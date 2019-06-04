#ifndef VARIABLE_PARTITION_DOMINANCE_H
#define VARIABLE_PARTITION_DOMINANCE_H

#include <vector>
#include <memory>
#include <cassert>

#include "../operator_cost.h"

class PDBAbstraction;
class Abstraction;
class SimulationRelation;
class LDSimulation;
class LabelRelation;
class DominanceRelation;
class LabelMap;
class OptionParser;
class Options;
class State;
class Labels;
class LabelledTransitionSystem;
class AbsBuilderAtomic;

class LocalLabelRelation;

class OutsideLabelRelation {
    std::vector<int> dom_by_noop;
    // Mapping from labels to groups so that they have the same transitions in all other TRs
    std::vector<int> label_to_equivalence_group;
    std::shared_ptr<std::vector<int> > label_cost;
public:

    OutsideLabelRelation(const std::vector<LocalLabelRelation> & local_dominances,
                         std::shared_ptr<std::vector<int> > label_cost,
                         const std::vector<int> & pattern);

    bool dominated_by_noop(int label) const {
        return dom_by_noop[label];
    }

    bool dominates(int label, int label2) const;
};

class LocalLabelRelation {
    std::vector<int> dom_by_noop;
    // Mapping from labels to groups so that they have the same transitions in this TR
    std::vector<int> label_to_equivalence_group;    
public:
    LocalLabelRelation(const SimulationRelation & atomic_dominance,
                       std::shared_ptr<std::vector<int> > label_cost,
                       LabelledTransitionSystem & lts);

    friend class OutsideLabelRelation;    
};


class VariablePartitionDominance {
    std::unique_ptr<PDBAbstraction> abstraction;
    std::unique_ptr<SimulationRelation> simulation_relation;
public:
    VariablePartitionDominance(Labels * labels,
                               const OutsideLabelRelation & label_relation, 
                               const std::vector<int> & pattern);
    
    /* const std::vector<int> & get_dominated_states(const State & state) const; */
    /* const std::vector<int> & get_dominated_states(int s) const; */
    
    PDBAbstraction & get_abstraction();
    
};

class VariablePartitionDominanceFactory {
    std::unique_ptr<AbsBuilderAtomic> builder_atomic;
    std::vector<LocalLabelRelation> atomic_label_relations;
    std::unique_ptr<LDSimulation> atomic_dominance;
    std::shared_ptr<std::vector<int> > label_cost;
public:
    VariablePartitionDominanceFactory (const Options &opts);
    
    void init(bool unit_cost, OperatorCost cost_type);
    
    virtual std::unique_ptr<VariablePartitionDominance>
        compute_dominance (const std::vector<int> & pattern);

    static void add_options_to_parser(OptionParser &parser);
};


#endif
