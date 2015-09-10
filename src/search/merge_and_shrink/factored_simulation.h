#ifndef MERGE_AND_SHRINK_FACTORED_SIMULATION_H
#define MERGE_AND_SHRINK_FACTORED_SIMULATION_H

#include <vector>
#include <memory>
#include "label_relation.h" 

#include "../sym/sym_variables.h"
#include "simulation_relation.h"
#include "compute_simulation.h"

/*
 * Class that represents the collection of simulation relations for a
 * factored LTS. Uses unique_ptr so that it owns the simulations and
 * it cannot be copied away.
 */

class FactoredSimulation {
    friend class ComputeSimulationRelation; //To allow initializing
					    //the simulations

    std::vector<std::unique_ptr<SimulationRelation> > simulations;
    LabelRelation label_dominance;
  public: 
FactoredSimulation(Labels * labels) : label_dominance(labels) 
    {}
   
    template<typename LTS>
	void compute_ld_simulation(std::vector<LTS *> & _ltss,
				   const LabelMap & labelMap, 
				   bool incremental_step, bool nold_simulation, 
				   ComputeSimulationRelation  & alg_compute_simulation) {
	Timer t;

	if(!nold_simulation){
	    label_dominance.init(_ltss, *this, labelMap);
	}else{
	    label_dominance.init_identity(_ltss.size(), labelMap);
	}

        label_dominance.init(_ltss, *this, labelMap);
						
	std::cout << "Label dominance initialized: " << _ltss.size() <<" " << t() << std::endl;
	do{
	    std::cout << "LDsimulation loop: ";
	    //label_dominance.dump();
	    if (incremental_step) {
		// Should be enough to just update the last (i.e., the new) simulation here.
		alg_compute_simulation.update(simulations.size() - 1,
					      _ltss.back(), label_dominance, *(simulations.back()));
	    } else {
		for (int i = 0; i < simulations.size(); i++){
		    alg_compute_simulation.update(i, _ltss[i], label_dominance, *(simulations[i]));
		    //_dominance_relation[i]->dump(_ltss[i]->get_names());
		}
	    }
	    std::cout << " took " << t() << "s" << std::endl;
	    //return; //PIET-edit: remove this for actual runs; just here for debugging the complex stuff
	}while(label_dominance.update(_ltss, *this));
	//for(int i = 0; i < _ltss.size(); i++){
	//_ltss[i]->dump();
	//	_dominance_relation[i]->dump(_ltss[i]->get_names());
	//}
	//label_dominance.dump_equivalent();
	//label_dominance.dump_dominance();
	//exit(0);
	//}
    }

    bool propagate_transition_pruning(int lts_id,
				      const std::vector<LabelledTransitionSystem *> & ltss, 
				      int src, int label_id, int target) const;

    int prune_subsumed_transitions(std::vector<Abstraction *> & abstractions, 
				   const LabelMap & labelMap,
				   const std::vector<LabelledTransitionSystem *> & ltss, 
				   int lts_id);


    EquivalenceRelation* get_equivalent_labels_relation(const LabelMap & labelMap, 
							std::set<int> & dangerous_LTSs) const {
	return label_dominance.get_equivalent_labels_relation(labelMap, dangerous_LTSs);
    }
    
    void remove_useless();

    void dump_statistics () const;
    //Statistics of the factored simulation 
    int num_equivalences() const;
    int num_simulations() const;   
    //Computes the probability of selecting a random pair s, s' such
    //that s simulates s'.
    double get_percentage_simulations(bool ignore_equivalences) const;
    //Computes the probability of selecting a random pair s, s' such that
    //s is equivalent to s'.
    double get_percentage_equivalences() const;
    double get_percentage_equal() const;


    //Methods to obtain the BDD representation for pruning
    void precompute_dominated_bdds(SymVariables * vars);
    void precompute_dominating_bdds(SymVariables * vars);

    BDD getSimulatedBDD(SymVariables * vars, const State &state) const;
    BDD getSimulatingBDD(SymVariables * vars, const State &state) const;
    BDD getIrrelevantStates(SymVariables * vars) const;



    //Methods to access the underlying simulation relations
    int size () const {
	return simulations.size();
    }
    
    SimulationRelation & operator[](int index) {
        return *(simulations[index]);
    }

    const SimulationRelation & operator[](int index) const {
        return *(simulations[index]);
    }

    std::vector<std::unique_ptr<SimulationRelation> >::iterator begin() {
	return simulations.begin();
    }

    std::vector<std::unique_ptr<SimulationRelation> >::iterator end() {
	return simulations.end();
    }

    std::vector<std::unique_ptr<SimulationRelation> >::iterator rbegin() {
	return simulations.begin();
    }

    std::vector<std::unique_ptr<SimulationRelation> >::iterator rend() {
	return simulations.end();
    }


    std::vector<std::unique_ptr<SimulationRelation> >::const_iterator begin() const {
	return simulations.begin();
    }

    std::vector<std::unique_ptr<SimulationRelation> >::const_iterator end() const {
	return simulations.end();
    }


    std::vector<std::unique_ptr<SimulationRelation> >::const_iterator rbegin() const{
	return simulations.begin();
    }

    std::vector<std::unique_ptr<SimulationRelation> >::const_iterator rend() const{
	return simulations.end();
    }

    SimulationRelation * back() {
	return simulations.back().get();
    }

    void clear(){
	std::vector<std::unique_ptr<SimulationRelation>>().swap(simulations);
    }
};
#endif
