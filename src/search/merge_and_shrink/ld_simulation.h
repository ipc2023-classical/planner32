#ifndef MERGE_AND_SHRINK_LD_SIMULATION_H
#define MERGE_AND_SHRINK_LD_SIMULATION_H

#include "../operator_cost.h"

#include "label_relation.h"
#include "simulation_relation.h"
#include "factored_simulation.h"

#include "../option_parser.h"

#include <memory>
#include <vector>

class Labels;
class MergeStrategy;
class Abstraction;
class SymManager;
class LabelledTransitionSystem;
class LTSComplex;
class ComputeSimulationRelation;

// Label dominance simulation
class LDSimulation {  
protected:
    const bool skip_simulation;
    const bool nold_simulation;

    const bool apply_simulation_shrinking;
    const bool apply_subsumed_transitions_pruning;
    const bool apply_label_dominance_reduction;
    const bool prune_dead_operators;
    const bool forbid_lr;

    const bool complex_simulation;
    const bool complex_lts;

    const bool use_expensive_statistics;
    const int limit_absstates_merge;
    const int limit_transitions_merge;

    const bool use_mas;
    const int limit_seconds_mas; //Limit of seconds for building the abstraction
    std::unique_ptr<MergeStrategy> merge_strategy;
    const bool use_bisimulation;
    const bool intermediate_simulations;
    const bool incremental_simulations;


    /* Parameters for constructing a final abstraction after the simulation */
    const bool compute_final_abstraction; 
    ShrinkStrategy *const shrink_strategy;
    const bool shrink_after_merge;
    const bool original_merge; //Forces the ld simulation to use the
			       //original merge,


    //TODO: Use unique_ptr here
    std::unique_ptr<Labels> labels;
    std::vector<Abstraction *> abstractions;
    FactoredSimulation simulations;
    std::unique_ptr<Abstraction> final_abstraction;

    std::vector<int> useless_vars;
    std::vector<bool> dead_labels;

    std::unique_ptr<ComputeSimulationRelation> alg_compute_simulation; 

    void dump_options() const;

    void build_abstraction();
    void compute_ld_simulation(bool incremental_step = false);
    void build_factored_systems ();

    std::vector<std::vector<int> > get_variable_partition_greedy();

    /* void compute_ld_simulation_after_merge(std::vector<Abstraction *> & _abstractions, */
    /*         FactoredSimulation & _simulations, */
    /*         const std::pair<int, int> & next_systems); */

    /* static void compute_ld_simulation(Labels * _labels, std::vector<Abstraction *> & _abstractions, */
    /*         FactoredSimulation & _simulations, bool no_ld = false); */

    /* static void compute_ld_simulation_complex(Labels * _labels, */
    /*         std::vector<Abstraction *> & _abstractions, */
    /*         FactoredSimulation & _simulations, bool no_ld = false); */


    template<typename LTS>
	void compute_ld_simulation(std::vector<LTS *> & _ltss,
				   const LabelMap & labelMap, 
				   LabelRelation & label_dominance, 
				   bool incremental_step);

    // If lts_id = -1 (default), then prunes in all ltss. If lts_id > 0,
    // prunes transitions dominated in all in all LTS, but other
    // transitions are only checked for lts_id
    int prune_subsumed_transitions(const LabelMap & labelMap,
				   LabelRelation & label_dominance, 
				   const std::vector<LabelledTransitionSystem *> & ltss, 
				   int lts_id);


    void remove_dead_labels(std::vector<Abstraction *> & abstractions);

    int remove_useless_abstractions(std::vector<Abstraction *> & abstractions) ;

    Abstraction * complete_heuristic() const;
public:
    LDSimulation(const Options &options);
    LDSimulation(bool unit_cost, const Options &options, OperatorCost cost_type);
    virtual ~LDSimulation();

    void initialize();

    bool pruned_state(const State &state) const;
    int get_cost(const State &state) const;

    inline FactoredSimulation & get_simulations() {
        return simulations;
    }

    void getVariableOrdering(std::vector <int> & var_order);

    static void add_options_to_parser(OptionParser &parser);
};


#endif
