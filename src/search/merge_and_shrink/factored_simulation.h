#ifndef MERGE_AND_SHRINK_FACTORED_SIMULATION_H
#define MERGE_AND_SHRINK_FACTORED_SIMULATION_H

#include <vector>
#include <memory>
#include "../sym/sym_variables.h"

class SimulationRelation;

/*
 * Class that represents the collection of simulation relations for a
 * factored LTS. Uses unique_ptr so that it owns the simulations and
 * it cannot be copied away.
 */

class FactoredSimulation {
    friend class ComputeSimulationRelation; //To allow initializing
					    //the simulations

    std::vector<std::unique_ptr<SimulationRelation> > simulations;
  public: 
    FactoredSimulation() {}

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
