#ifndef MERGE_AND_SHRINK_NUMERIC_SIMULATION_RELATION_H
#define MERGE_AND_SHRINK_NUMERIC_SIMULATION_RELATION_H

#include <vector>
#include <string>
#include <iostream>
#include "../sym/sym_variables.h"
#include "numeric_label_relation.h"

class Labels;
class Abstraction;
class CompositeAbstraction;
class LabelledTransitionSystem;

class NumericSimulationRelation {
protected:
    Abstraction * abs;


    std::vector<int> tau_labels;
    std::vector<std::vector<int> > distances_with_tau;

    std::vector<std::vector<int> > relation;
public:
    NumericSimulationRelation(Abstraction * _abs);
    
    void init_goal_respecting (); 
    void update (int lts_id, const LabelledTransitionSystem * lts,
			const NumericLabelRelation & label_dominance);


    bool simulates (const State & t, const State & s) const;

    inline bool simulates (int s, int t) const {
        return relation[s][t] >= 0;
    }

    inline bool may_simulate (int s, int t) const {
        return relation[s][t] > std::numeric_limits<int>::lowest();
    }

    inline int q_simulates (int s, int t) const {
        return relation[s][t];
    }

    inline bool similar (int s, int t) const {
        return simulates(s, t) && simulates(t, s);
    }

    inline void update_value (int s, int t, int value) {
        relation[s][t] = value;
    }

    inline void remove (int s, int t) {
        relation[s][t] = std::numeric_limits<int>::lowest();
    }

    inline const std::vector<std::vector<int> > & get_relation () {
        return relation;
    }

    inline int minus_shortest_path_with_tau(int from, int to) {
	return distances_with_tau[from][to] == std::numeric_limits<int>::max() 
	    ? std::numeric_limits<int>::lowest() 
	    : - distances_with_tau[from][to];
    }

    void precompute_shortest_path_with_tau(const LabelledTransitionSystem * lts, 
					   const std::vector<int> & tau_labels);

    void dump(const std::vector<std::string> & names) const;
    void dump() const;
};

#endif
