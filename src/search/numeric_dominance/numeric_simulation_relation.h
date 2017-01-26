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
    Abstraction * abs;
    
    std::vector<int> tau_labels;
    std::vector<std::vector<int> > distances_with_tau;
    //List of states for which distances_with_tau is not infinity
    std::vector<std::vector<int> > reachable_with_tau;

    std::vector<std::vector<int> > relation;


    int compare_noop(int lts_id, const LTSTransition & trs, int t,
		     int tau_distance, const NumericLabelRelation & label_dominance) const;


    int compare_transitions(int lts_id, const LTSTransition & trs, const LTSTransition & trt, 
			    int tau_distance, const NumericLabelRelation & label_dominance) const;

public:
    NumericSimulationRelation(Abstraction * _abs);
    
    void init_goal_respecting (); 
    void update (int lts_id, const LabelledTransitionSystem * lts,
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
	assert(s < relation.size());
	assert(t < relation[s].size());
	assert(s!=t || relation[s][t] == 0);
        return relation[s][t];
    }

    inline bool similar (int s, int t) const {
        return simulates(s, t) && simulates(t, s);
    }

    inline void update_value (int s, int t, int value) {
	if(value < -1000) {
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
};

#endif
