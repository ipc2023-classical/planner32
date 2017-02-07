#ifndef NUMERIC_DOMINANCE_NUMERIC_SIMULATION_RELATION_H
#define NUMERIC_DOMINANCE_NUMERIC_SIMULATION_RELATION_H

#include <vector>
#include <string>
#include <iostream>
#include "../sym/sym_variables.h"
#include "numeric_label_relation.h"

class Labels;
class Abstraction;
class CompositeAbstraction;
class LabelledTransitionSystem;
class LTSTransition; 
class NumericSimulationRelation {
protected:
    const  Abstraction * abs;
    const int truncate_value; 
    
    std::vector<int> tau_labels;
    std::vector<std::vector<int> > distances_with_tau;
    //List of states for which distances_with_tau is not infinity
    std::vector<std::vector<int> > reachable_with_tau;

    std::vector<std::vector<int> > relation;

    //BDDs of each abstract state
    std::vector<BDD> abs_bdds;
    SymVariables * vars; 

    //For each abstract state, we create a BDD/ADD/BDDMap that represents all the
    //abstract states dominated by it and dominating it
    std::vector<BDD> may_dominated_bdds, may_dominating_bdds;
    std::vector<BDD> dominated_bdds, dominating_bdds;
    std::vector<ADD> dominated_adds, dominating_adds;
    std::vector<std::map<int, BDD> > dominated_bdd_maps, dominating_bdd_maps;


    int compare_noop(int lts_id, const LTSTransition & trs, int t,
		     int tau_distance, const NumericLabelRelation & label_dominance) const;


    int compare_transitions(int lts_id, const LTSTransition & trs, const LTSTransition & trt, 
			    int tau_distance, const NumericLabelRelation & label_dominance) const;

public:
    NumericSimulationRelation(Abstraction * _abs, int truncate_value);
    
    void init_goal_respecting (); 
    int update (int lts_id, const LabelledTransitionSystem * lts,
		 const NumericLabelRelation & label_dominance);


    bool pruned(const State & state) const;
	
    int  q_simulates (const State & t, const State & s) const;
    int  q_simulates (const std::vector<int> & t, const std::vector<int> & s) const;

    inline bool simulates (int s, int t) const {
        return relation[s][t] >= 0;
    }

    inline bool may_simulate (int s, int t) const {
	assert(s < relation.size());
	assert(t < relation[s].size());
        return relation[s][t] > std::numeric_limits<int>::lowest();
    }

    inline int q_simulates (int s, int t) const {
	if(s >= relation.size()) {
	    std::cout << s << std::endl;
	    std::cout << relation.size() << std::endl;
	}
	assert(s < relation.size());
	assert(t < relation[s].size());
	assert(s!=t || relation[s][t] == 0);
        return relation[s][t];
    }

    inline bool similar (int s, int t) const {
        return simulates(s, t) && simulates(t, s);
    }

    inline void update_value (int s, int t, int value) {

	if(value < -truncate_value) {
	    //std::cout << value << " rounded to -infty: " << truncate_value << std::endl;
	    value = std::numeric_limits<int>::lowest();
	}
        relation[s][t] = value;
    }

    inline void remove (int s, int t) {
        relation[s][t] = std::numeric_limits<int>::lowest();
    }

    inline const std::vector<std::vector<int> > & get_relation () {
        return relation;
    }

    inline int minus_shortest_path_with_tau(int from, int to) {
	assert(from < distances_with_tau.size());
	assert(to < distances_with_tau[from].size());
	return distances_with_tau[from][to] == std::numeric_limits<int>::max() 
	    ? std::numeric_limits<int>::lowest() 
	    : - distances_with_tau[from][to];
    }

    void precompute_shortest_path_with_tau(const LabelledTransitionSystem * lts, int lts_id,
					   const NumericLabelRelation & label_dominance);

    void dump(const std::vector<std::string> & names) const;
    void dump() const;

    void statistics() const;

    void precompute_absstate_bdds(SymVariables * vars);
    void precompute_bdds(bool dominating, bool quantified, bool use_ADD);


    BDD getSimulatedBDD(const State & state) const;
    BDD getSimulatingBDD(const State & state) const;
    const std::map<int, BDD> & getSimulatedBDDMap(const State & state) const;
    const std::map<int, BDD> & getSimulatingBDDMap(const State & state) const;
    BDD getMaySimulatedBDD(const State & state) const;
    BDD getMaySimulatingBDD(const State & state) const;
    ADD getSimulatedADD(const State & state) const;
    ADD getSimulatingADD(const State & state) const;


};

#endif
