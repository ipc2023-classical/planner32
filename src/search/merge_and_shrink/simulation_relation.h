#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_H

#include <vector>
#include <string>
#include <iostream>
#include "../sym/sym_variables.h"

class Labels;
class Abstraction;
class CompositeAbstraction;

// First implementation of a simulation relation. 
class SimulationRelation {
protected:
    Abstraction * abs;

    //By now we assume that the partition is unitary... we can improve
    //this later with EquivalenceRelation
    std::vector<std::vector<bool> > relation;
    //To compute intermediate simulations. ]
    //If fixed_relation is set, then we can skip checking it
    std::vector<std::vector<bool> > fixed_relation;


    //For each abstract state, we create a BDD that represents all the
    //abstract states dominated by it and dominating it
    std::vector<BDD> dominated_bdds, dominating_bdds;

    //BDDs of each abstract state
    std::vector<BDD> abs_bdds;
    BDD zeroBDD;

    //Vectors of states dominated/dominating by each state. Lazily
    // computed when needed.
    std::vector<std::vector<int>> dominated_states, dominating_states;

public:
    SimulationRelation(Abstraction * _abs);

    void init_goal_respecting (); 
    void init_identity (); 
    void init_incremental(CompositeAbstraction * _abs, 
			  const SimulationRelation & simrel_one, 
			  const SimulationRelation & simrel_two); 

    SimulationRelation(CompositeAbstraction * _abs, 
		       const SimulationRelation & simrel_one, 
		       const SimulationRelation & simrel_two);

    virtual ~SimulationRelation();

    void apply_shrinking_to_table(const std::vector<int> & abstraction_mapping);

    inline bool simulates (int s, int t) const {
        return relation[s][t];
    }

    inline bool similar (int s, int t) const {
        return relation[s][t] && relation[t][s];
    }

    inline bool fixed_simulates(int s, int t) const {
	return !fixed_relation.empty() && fixed_relation[s][t];
     }
    inline void remove (int s, int t) {
        relation[s][t] = false;
    }

    inline const std::vector<std::vector<bool> > & get_relation () {
        return relation;
    }

    int num_equivalences() const;
    int num_simulations(bool ignore_equivalences) const;
    int num_states() const { 
        return relation.size();
    }

    int num_different_states() const;

    void reset();
    void dump(const std::vector<std::string> & names) const;
    void dump() const;

    BDD getSimulatedBDD(const State & state) const;
    BDD getSimulatingBDD(const State & state) const;
    BDD getIrrelevantStates(SymVariables * vars);

    void precompute_absstate_bdds(SymVariables * vars);
    void precompute_dominated_bdds();
    void precompute_dominating_bdds();

    inline const Abstraction * get_abstraction() const {
	return abs;
    }

    const std::vector<BDD> & get_dominated_bdds ();
    const std::vector<BDD> & get_dominating_bdds ();
    const std::vector<BDD> & get_abs_bdds() const;
    const std::vector <int> & get_varset() const;
    bool pruned(const State & state) const;
    int get_cost(const State & state) const;
    int get_index (const State & state) const;
    const std::vector<int> & get_dominated_states(const State & state);
    const std::vector<int> & get_dominating_states(const State & state);

    //Computes the probability of selecting a random pair s, s' such
    //that s simulates s'.
    double get_percentage_simulations(bool ignore_equivalences) const;

    //Computes the probability of selecting a random pair s, s' such that
    //s is equivalent to s'.
    double get_percentage_equivalences() const;

    void shrink();

    void compute_list_dominated_states();

};

#endif
