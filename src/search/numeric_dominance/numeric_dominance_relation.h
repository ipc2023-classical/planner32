#ifndef NUMERIC_DOMINANCE_NUMERIC_DOMINANCE_RELATION_H
#define NUMERIC_DOMINANCE_NUMERIC_DOMINANCE_RELATION_H

#include <vector>
#include <memory>
#include "../merge_and_shrink/abstraction.h" 
#include "../merge_and_shrink/labels.h" 
#include "numeric_simulation_relation.h"
#include "numeric_label_relation.h"

class LabelledTransitionSystem;
class LTSComplex;
class Operator;
class SearchProgress;

/*
 * Class that represents the collection of simulation relations for a
 * factored LTS. Uses unique_ptr so that it owns the simulations and
 * it cannot be copied away.
 */
template <typename T>
class NumericDominanceRelation {

    //Auxiliar data-structures to perform successor pruning 
    mutable std::set<int> relevant_simulations;
    mutable std::vector<int> parent, parent_ids, succ;

    const int truncate_value;
    NumericLabelRelation<T> label_dominance;

protected:
    std::vector<std::unique_ptr<NumericSimulationRelation<T>> > simulations;
    std::vector<int> simulation_of_variable;
    T total_max_value; 

    std::unique_ptr<NumericSimulationRelation<T>> init_simulation (Abstraction * _abs);

    template<typename LTS>
	void compute_ld_simulation_template(std::vector<LTS *> & _ltss,
					    const LabelMap & labelMap,
					    bool dump) {
	assert(_ltss.size() == simulations.size());
	Timer t;
	int num_iterations = 0;
	int num_inner_iterations = 0;

	std::cout << "Compute numLDSim on " << _ltss.size() << " LTSs." << std::endl;
	for (auto lts : _ltss) {
	    std::cout << lts->size() << " states and " << lts->num_transitions() << " transitions" << std::endl;
	}

	
        label_dominance.init(_ltss, *this, labelMap);
	
	std::cout << "  Init numLDSim in " << t() << "s: " << std::flush;
	do{
	    num_iterations++;
	    //label_dominance.dump();
		for (int i = 0; i < simulations.size(); i++){
		    num_inner_iterations += simulations[i]->update(i, _ltss[i], label_dominance);
		    //_dominance_relation[i]->dump(_ltss[i]->get_names());
		}
	    std::cout << " " << t() << "s" << std::flush;
	}while(label_dominance.update(_ltss, *this));
	std::cout << std::endl << "Numeric LDSim computed " << t() << std::endl;
	std::cout << "Numeric LDSim outer iterations: " << num_iterations << std::endl;
	std::cout << "Numeric LDSim inner iterations: " << num_inner_iterations << std::endl;
	
	std::cout << "------" << std::endl; 
        for(int i = 0; i < _ltss.size(); i++){
	    simulations[i]->statistics();
	    std::cout << "------" << std::endl; 
	}

	if(dump) {
	    std::cout << "------" << std::endl; 
	    for(int i = 0; i < _ltss.size(); i++){ 
		/*     //_ltss[i]->get_abstraction()->dump_names(); */
		simulations[i]->dump(_ltss[i]->get_names()); 
		std::cout << "------" << std::endl; 
		/*     label_dominance.dump(_ltss[i], i);  */
	    } 
	}
	//label_dominance.dump_equivalent();
	//label_dominance.dump_dominance();
	//exit(0);
	//}


	total_max_value = 0;
	for (auto &  sim : simulations){
	    total_max_value += sim->compute_max_value();
	}
    }

public:

    NumericDominanceRelation(Labels * labels, 
			     int truncate_value_, 
			     bool compute_tau_labels_with_noop_dominance) : 
								truncate_value(truncate_value_),
								label_dominance(labels,
										compute_tau_labels_with_noop_dominance) 
    {}


								
    bool action_selection_pruning (const State & state, 
				   std::vector<const Operator *> & applicable_operators,
				   SearchProgress & search_progress) const;
    void prune_dominated_by_parent (const State & state, 
				   std::vector<const Operator *> & applicable_operators,
				   SearchProgress & search_progress, bool parent_ids_stored) const;

 
    //Methods to use the dominance relation 
    bool pruned_state(const State &state) const;
    //int get_cost(const State &state) const;

    //bool parent_dominates_successor(const State & parent, const Operator *op) const;
    bool dominates(const State &t, const State & s, int g_diff) const;

    bool dominates_parent(const std::vector<int> & state, const std::vector<int> & parent, int action_cost) const;

    void init (const std::vector<Abstraction *> & abstractions);
    
    void compute_ld_simulation (std::vector<LabelledTransitionSystem *> & _ltss, 
				const LabelMap & labelMap, bool dump) {
	compute_ld_simulation_template(_ltss, labelMap, dump);
    }

    
    //Methods to obtain the BDD representation for pruning
    void precompute_bdds(SymVariables * vars, bool dominating, 
			 bool quantified, bool use_ADD);

    //Methods to obtain the BDD representation for pruning
    BDD getDominatedBDD(SymVariables * vars, const State &state, bool trade_off_dominance) const;
    BDD getDominatingBDD(SymVariables * vars, const State &state) const;
    std::map<T, BDD> getDominatedBDDMap(SymVariables * vars, const State &state, bool only_positive) const;
    //std::map<T, BDD> getDominatingBDDMap(SymVariables * vars, const State &state) const;

    /* map<int, BDD> getBDDMap(SymVariables * vars, const State &state, bool dominating); */
    /* ADD getADD(SymVariables * vars, const State &state, bool dominating); */


    //Methods to access the underlying simulation relations
    const std::vector<std::unique_ptr<NumericSimulationRelation<T>> > & get_simulations () const{
	return simulations;
    }

    int size () const {
	return simulations.size();
    }
    
    NumericSimulationRelation<T> & operator[](int index) {
        return *(simulations[index]);
    }

    const NumericSimulationRelation<T> & operator[](int index) const {
        return *(simulations[index]);
    }    
};

#endif
