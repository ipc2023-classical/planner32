#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_H

#include <vector>
#include <string>
#include <iostream>
#include "abstraction.h"
#include "label_relation.h"
#include "../sym/sym_variables.h"

class Labels;

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

public:
    SimulationRelation(Abstraction * _abs);

    virtual ~SimulationRelation(){}

    virtual void update(int lts_id, const LabelledTransitionSystem * lts,
            const LabelRelation & label_dominance) = 0;

    virtual void update(int lts_id, const LTSEfficient * lts,
            const LabelRelation & label_dominance) = 0;

    inline bool simulates (int s, int t) const {
        return relation[s][t];
    }

    inline bool similar (int s, int t) const {
        return relation[s][t] && relation[t][s];
    }


    inline void remove (int s, int t) {
        relation[s][t] = false;
    }

    inline const std::vector<std::vector<bool> > & get_relation () {
        return relation;
    }


    int num_equivalences() const;
    int num_simulations(bool ignore_equivalences) const;
    int num_states() const{
        return abs_bdds.size();
    }

    int num_different_states() const;

    void reset();
    void dump(const std::vector<std::string> & names) const;

    BDD getSimulatedBDD(const State & state) const;
    BDD getSimulatingBDD(const State & state) const;

    void precompute_absstate_bdds(SymVariables * vars);
    void precompute_dominated_bdds();
    void precompute_dominating_bdds();

    inline const std::vector<BDD> & get_dominated_bdds () {
        if(dominated_bdds.empty()) precompute_dominated_bdds();
        return dominated_bdds;
    }

    inline const std::vector<BDD> & get_dominating_bdds () {
        if(dominating_bdds.empty()) precompute_dominating_bdds();
        return dominating_bdds;
    }


    inline const std::vector<BDD> & get_abs_bdds() const{
        return abs_bdds;
    }

    inline const std::vector <int> & get_varset() const {
        return abs->get_varset();
    }

    inline bool pruned(const State & state) const {
        return abs->get_abstract_state(state) == -1;
    }

    inline int get_cost(const State & state) const {
        return abs->get_cost(state);
    }

    //Computes the probability of selecting a random pair s, s' such
    //that s simulates s'.
    double get_percentage_simulations(bool ignore_equivalences) const;

    //Computes the probability of selecting a random pair s, s' such that
    //s is equivalent to s'.
    double get_percentage_equivalences() const;

    void shrink();

};

#endif
