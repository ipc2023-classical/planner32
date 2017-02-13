#ifndef NUMERIC_DOMINANCE_NUMERIC_SIMULATION_RELATION_H
#define NUMERIC_DOMINANCE_NUMERIC_SIMULATION_RELATION_H

#include <vector>
#include <string>
#include <iostream>
#include "../sym/sym_variables.h"
#include "numeric_label_relation.h"
#include "int_epsilon.h"

class Labels;
class Abstraction;
class CompositeAbstraction;
class LabelledTransitionSystem;
class LTSTransition; 

template <typename T> 
class NumericSimulationRelation {
protected:
    const  Abstraction * abs;
    const int truncate_value; 
    
    std::vector<int> tau_labels;
    std::vector<std::vector<T> > distances_with_tau;

    //List of states for which distances_with_tau is not infinity
    std::vector<std::vector<int> > reachable_with_tau;

    std::vector<std::vector<T> > relation;
    T max_relation_value;

    //BDDs of each abstract state
    std::vector<BDD> abs_bdds;
    SymVariables * vars; 

    //For each abstract state, we create a BDD/ADD/BDDMap that represents all the
    //abstract states dominated by it and dominating it
    std::vector<BDD> may_dominated_bdds, may_dominating_bdds;
    std::vector<BDD> dominated_bdds, dominating_bdds;
    std::vector<ADD> dominated_adds, dominating_adds;
    std::vector<std::map<T, BDD> > dominated_bdd_maps, dominating_bdd_maps;

    int get_label_cost (int label) const;

    T compare_noop(int lts_id, const LTSTransition & trs, int t,
		     T tau_distance, 
		     const NumericLabelRelation<T> & label_dominance) const;


    T compare_transitions(int lts_id, const LTSTransition & trs, const LTSTransition & trt, 
			  T tau_distance, const NumericLabelRelation<T> & label_dominance) const;

public:
    NumericSimulationRelation(Abstraction * _abs, int truncate_value);
    
    void init_goal_respecting (); 
    int update (int lts_id, const LabelledTransitionSystem * lts,
		 const NumericLabelRelation<T> & label_dominance);


    bool pruned(const State & state) const;
	
    T  q_simulates (const State & t, const State & s) const;
    T  q_simulates (const std::vector<int> & t, const std::vector<int> & s) const;

    int get_abstract_state_id(const State & t) const; 
    int get_abstract_state_id(const std::vector<int> & t) const; 

    inline bool simulates (int s, int t) const {
        return relation[s][t] >= 0;
    }

    inline bool may_simulate (int s, int t) const {
	assert(s < relation.size());
	assert(t < relation[s].size());
        return relation[s][t] > std::numeric_limits<int>::lowest();
    }

    inline T q_simulates (int s, int t) const {
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

    inline void update_value (int s, int t, T value) {
	if(value < -truncate_value) {
	    //std::cout << value << " rounded to -infty: " << truncate_value << std::endl;
	    value = std::numeric_limits<int>::lowest();
	}
        relation[s][t] = value;
    }

    inline const std::vector<std::vector<T> > & get_relation () {
        return relation;
    }

    T compute_max_value() {
	max_relation_value = 0; 
	for(const auto & row : relation) {
	    for (T value : row) {
		max_relation_value = std::max(max_relation_value, value);
	    }
	}
	return max_relation_value;
    }

    T get_max_value () const {
	return max_relation_value; 
    }

    inline T minus_shortest_path_with_tau(int from, int to) {
	assert(from < distances_with_tau.size());
	assert(to < distances_with_tau[from].size());
	assert((from != to && distances_with_tau[from][to] > 0) ||
	       (from == to && distances_with_tau[from][to] == 0));
	return distances_with_tau[from][to] == std::numeric_limits<int>::max() 
	    ? std::numeric_limits<int>::lowest() 
	    : - distances_with_tau[from][to];
    }

    void precompute_shortest_path_with_tau(const LabelledTransitionSystem * lts, int lts_id,
					   const NumericLabelRelation<T> & label_dominance);

    void dump(const std::vector<std::string> & names) const;
    void dump() const;

    void statistics() const;

    void precompute_absstate_bdds(SymVariables * vars);
    void precompute_bdds(bool dominating, bool quantified, bool use_ADD);

    BDD getSimulatedBDD(const State & state) const;
    BDD getSimulatingBDD(const State & state) const;
    const std::map<T, BDD> & getSimulatedBDDMap(const State & state) const;
    const std::map<T, BDD> & getSimulatingBDDMap(const State & state) const;
    BDD getMaySimulatedBDD(const State & state) const;
    BDD getMaySimulatingBDD(const State & state) const;
    ADD getSimulatedADD(const State & state) const;
    ADD getSimulatingADD(const State & state) const;
};


#endif
