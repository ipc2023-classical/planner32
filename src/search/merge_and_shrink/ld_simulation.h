#ifndef MERGE_AND_SHRINK_LD_SIMULATION_H
#define MERGE_AND_SHRINK_LD_SIMULATION_H

#include "../operator_cost.h"

#include "simulation_relation.h"
#include "dominance_relation.h"

#include "../option_parser.h"

#include <memory>
#include <vector>

class Labels;
class MergeStrategy;
class Abstraction;
class SymManager;
class LabelledTransitionSystem;
class LTSComplex;
class NumericDominanceRelation;

enum class LabelDominanceType {
    NONE, NOOP, NORMAL
};
enum class SimulationType {
    NONE, SIMPLE, COMPLEX
};

std::ostream & operator<<(std::ostream &os, const LabelDominanceType & m);
extern const std::vector<std::string> LabelDominanceTypeValues;
std::ostream & operator<<(std::ostream &os, const SimulationType & m);
extern const std::vector<std::string> SimulationTypeValues;


// Label dominance simulation
class LDSimulation {  
    friend class AbsBuilderPDB;
    
protected:
    std::unique_ptr<Labels> labels;
    std::vector<Abstraction *> abstractions;
    std::unique_ptr<DominanceRelation> dominance_relation;
    //std::unique_ptr<Abstraction> final_abstraction;

    std::vector<int> useless_vars;
    std::vector<bool> dead_labels;

    std::unique_ptr<DominanceRelation> create_dominance_relation(SimulationType simulation_type, 
								 LabelDominanceType label_dominance_type, 
								 int switch_off_label_dominance);

    void compute_ld_simulation(SimulationType simulation_type, 
			       LabelDominanceType label_dominance_type, 
			       int switch_off_label_dominance, bool complex_lts,
			       bool apply_subsumed_transitions_pruning, 
			       bool apply_label_dominance_reduction, 
			       bool apply_simulation_shrinking, bool preserve_all_optimal_plans,
			       bool incremental_step = false);

    std::vector<std::vector<int> > get_variable_partition_greedy();

    template<typename LTS> void compute_ld_simulation(std::vector<LTS *> & _ltss,
						      const LabelMap & labelMap, 
						      bool incremental_step);

    // If lts_id = -1 (default), then prunes in all ltss. If lts_id > 0,
    // prunes transitions dominated in all in all LTS, but other
    // transitions are only checked for lts_id
    int prune_subsumed_transitions(const LabelMap & labelMap,
				   const std::vector<LabelledTransitionSystem *> & ltss, 
				   int lts_id, bool preserve_all_optimal_plans);


    void remove_dead_labels(std::vector<Abstraction *> & abstractions);

    int remove_useless_abstractions(std::vector<Abstraction *> & abstractions) ;

    int remove_useless_atomic_abstractions(std::vector<Abstraction *> & abstractions) const;

    double estimated_memory_MB (std::vector<Abstraction * > all_abstractions) const;


public:
    LDSimulation(bool unit_cost, const Options &opts, OperatorCost cost_type);
    virtual ~LDSimulation();

    void init_atomic_abstractions();
    void init_factored_systems(const std::vector<std::vector<int> > & partition_vars);


    void build_abstraction(MergeStrategy * merge_strategy, int limit_absstates_merge, 
			   int limit_transitions_merge, bool original_merge,
			   ShrinkStrategy * shrink_strategy, bool forbid_lr, 
			   int limit_seconds_mas, int limit_memory_kb_total,
			   bool intermediate_simulations, bool incremental_simulations,
			   SimulationType simulation_type, 
			   LabelDominanceType label_dominance_type, 
			   int switch_off_label_dominance,
			   bool complex_lts, 
			   bool apply_subsumed_transitions_pruning, bool apply_label_dominance_reduction, 
			   bool apply_simulation_shrinking, 
			   bool preserve_all_optimal_plans, bool expensive_statistics );


    void compute_final_simulation(SimulationType simulation_type, 
				  LabelDominanceType label_dominance_type, 
				  int switch_off_label_dominance, bool intermediate_simulations, bool complex_lts, 
				  bool apply_subsumed_transitions_pruning, 
				  bool apply_label_dominance_reduction, bool apply_simulation_shrinking, bool preserve_all_optimal_plans); 

    void complete_heuristic(MergeStrategy * merge_strategy, 
			    ShrinkStrategy * shrink_strategy,
			    bool shrink_after_merge, int limit_second_mas, int limit_memory_mas, 
			    bool prune_dead_operators, bool use_expensive_statistics, 
			    std::vector<std::unique_ptr<Abstraction>> & res ) const;

    void prune_dead_ops () const {
	prune_dead_ops(abstractions);
    }

    void prune_dead_ops (const std::vector<Abstraction*> & all_abstractions) const ;
    bool pruned_state(const State &state) const;
    int get_cost(const State &state) const;
    bool dominates(const State &t, const State &s) const {
	return dominance_relation->dominates(t, s);
    }

    inline DominanceRelation & get_dominance_relation() {
        return *dominance_relation;
    }

    inline bool has_dominance_relation() {
        return dominance_relation.get() != nullptr;
    }

    std::unique_ptr<NumericDominanceRelation> compute_numeric_dominance_relation() const;


    void getVariableOrdering(std::vector <int> & var_order);

    static void add_options_to_parser(OptionParser &parser);

    void release_memory(); 
};


#endif
