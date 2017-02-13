#ifndef MERGE_AND_SHRINK_ABSTRACTION_BUILDER_H
#define MERGE_AND_SHRINK_ABSTRACTION_BUILDER_H

#include "ld_simulation.h"

#include "../option_parser.h"
#include <memory> 

class Labels;
class Abstraction;
class MergeStrategy;
class ShrinkStrategy;

class AbstractionBuilder {   
    const Options opts;
protected: 
    const bool expensive_statistics;
    const bool dump;

    const int limit_seconds_total;
    const int limit_memory_kb_total; //Limit of seconds for building the abstraction

public: 
    AbstractionBuilder(const Options &opts);
    virtual ~AbstractionBuilder() = default;

    void  init_ldsim (bool unit_cost, OperatorCost cost_type, 
		      std::unique_ptr <LDSimulation> & ldSim) const; 
    
    virtual void build_abstraction (bool unit_cost, OperatorCost cost_type, 
				    std::unique_ptr<LDSimulation> & ldSim, 
				    std::vector<std::unique_ptr<Abstraction> > & abstractions) const = 0;

    virtual void dump_options() const = 0; 

    static void add_options_to_parser(OptionParser &parser);
};


class AbsBuilderPDB : public AbstractionBuilder {

    const int limit_absstates_merge;

 public: 
    AbsBuilderPDB(const Options &opts);
    virtual ~AbsBuilderPDB() = default;

    virtual void build_abstraction (bool unit_cost, OperatorCost cost_type,
				    std::unique_ptr<LDSimulation> & ldSim, 
				    std::vector<std::unique_ptr<Abstraction> > & abstractions) const;

    virtual void dump_options() const {
	std::cout << "AbsBuilderPDB" << std::endl;
    } 
};

class AbsBuilderAtomic : public AbstractionBuilder {

public: 
    AbsBuilderAtomic(const Options &opts); 
    virtual  ~AbsBuilderAtomic() = default;
    

    virtual void build_abstraction (bool unit_cost, OperatorCost cost_type,
				    std::unique_ptr<LDSimulation> & ldSim, 
				    std::vector<std::unique_ptr<Abstraction> > & abstractions) const;

    
    virtual void dump_options() const {
	std::cout << "AbsBuilderAtomic" << std::endl;
    }
};


class AbsBuilderMAS : public AbstractionBuilder {
    std::unique_ptr<MergeStrategy> merge_strategy;
    std::unique_ptr<ShrinkStrategy> shrink_strategy;
    const bool shrink_after_merge;
     
    const int limit_seconds_mas; //Limit of seconds for building the abstraction

    bool prune_dead_operators;
    bool store_original_operators;

    const bool restart;
    const int num_abstractions; 

public: 
    AbsBuilderMAS(const Options &opts); 
    virtual ~AbsBuilderMAS() = default;

    virtual void build_abstraction (bool unit_cost, OperatorCost cost_type,
				    std::unique_ptr<LDSimulation> & ldSim, 
				    std::vector<std::unique_ptr<Abstraction> > & abstractions) const;

    virtual void dump_options() const {
	std::cout << "AbsBuilderMAS" << std::endl;
    }
}; 



class AbsBuilderDefault : public AbstractionBuilder {
    std::unique_ptr<MergeStrategy> merge_strategy;
    const bool original_merge;
    const int limit_absstates_merge;
    const int min_limit_absstates_merge;
    const int limit_transitions_merge; 

    const int limit_absstates_shrink; 
     
    const int limit_seconds_mas; //Limit of seconds for building the abstraction

    const int num_abstractions;
    const int switch_off_label_dominance;

public: 
    AbsBuilderDefault(const Options &opts); 
    virtual ~AbsBuilderDefault() = default;

    virtual void build_abstraction (bool unit_cost, OperatorCost cost_type,
				    std::unique_ptr<LDSimulation> & ldSim, 
				    std::vector<std::unique_ptr<Abstraction> > & abstractions) const;

    virtual void dump_options() const {
	std::cout << "AbsBuilderDefault" << std::endl;
    }
}; 



class AbsBuilderMasSimulation : public AbstractionBuilder{

    const SimulationType simulation_type; 
    const LabelDominanceType label_dominance_type; 
    /* //Allows to switch off label dominance from normal to noop if the */
    /* //number of labels is greater than this parameter. */
    const int switch_off_label_dominance;

    const bool apply_simulation_shrinking;
    bool apply_subsumed_transitions_pruning;
    const bool apply_label_dominance_reduction;
    bool prune_dead_operators;
    bool store_original_operators;

    const bool complex_lts;
    /* const bool use_mas; */
    /* /\* Parameters for constructing a final abstraction after the simulation *\/ */
    /* const bool compute_final_abstraction;  */

    std::unique_ptr<MergeStrategy> merge_strategy;
    const bool original_merge;
    const int limit_absstates_merge;
    const int min_limit_absstates_merge;
    const int limit_transitions_merge;

    bool intermediate_simulations;
    bool incremental_simulations;
    const bool compute_final_simulation;

    const bool forbid_lr;

    std::unique_ptr<ShrinkStrategy> shrink_strategy;
    const bool shrink_after_merge;

    const int limit_seconds_mas; //Limit of seconds for building the abstraction

    void compute_ld_simulation(LDSimulation & ldSim,
			       bool incremental_step = false) const;

public: 
    AbsBuilderMasSimulation(const Options &opts);
    virtual ~AbsBuilderMasSimulation() = default;

    virtual void build_abstraction (bool unit_cost, OperatorCost cost_type,
				    std::unique_ptr<LDSimulation> & ldSim, 
				    std::vector<std::unique_ptr<Abstraction> > & abstractions) const; 

    virtual void dump_options() const;
};



class AbsBuilderComposite : public AbstractionBuilder {
    const std::vector<AbstractionBuilder *> strategies;
public: 
    AbsBuilderComposite(const Options &opts); 
    virtual ~AbsBuilderComposite() = default;

    virtual void build_abstraction (bool unit_cost, OperatorCost cost_type,
				    std::unique_ptr<LDSimulation> & ldSim, 
				    std::vector<std::unique_ptr<Abstraction> > & abstractions) const;

    virtual void dump_options() const {
	for (auto st : strategies) {
	    std::cout << "  ST1: ";
	    st->dump_options();
	}
    }
};

#endif
