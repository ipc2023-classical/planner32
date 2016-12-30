#ifndef NUMERIC_DOMINANCE_NUMERIC_LABEL_RELATION_H
#define NUMERIC_DOMINANCE_NUMERIC_LABEL_RELATION_H

#include "../merge_and_shrink/label_relation.h"

#include "numeric_label_relation.h"

#include <iostream>
#include <vector>
#include "../merge_and_shrink/labels.h"
#include "../merge_and_shrink/label.h"

class LabelledTransitionSystem;
class NumericSimulationRelation;
class NumericDominanceRelation;

/* 
 * Label relation represents the preorder relations on labels that
 * occur in a set of LTS
 */ 
class NumericLabelRelation {
    Labels * labels;
    int num_labels;

    // Summary matrix for each l1, l2 indicating whether l1 dominates
    // l2 in all (-2), in none (-1) or only in i (i)
    std::vector<std::vector<int> > dominates_in;
    std::vector<int> dominates_noop_in, dominated_by_noop_in;

    //For each lts, matrix indicating whether l1 simulates l2, noop
    //simulates l or l simulates noop
    std::vector<std::vector<std::vector<int > > > lqrel;
    std::vector<std::vector<int> > simulated_by_irrelevant;
    std::vector<std::vector<int> > simulates_irrelevant;

    std::vector<std::vector<int> > tau_labels;

    bool update(int i, const LabelledTransitionSystem * lts, 
		const NumericSimulationRelation & sim);

    //Returns true if l1 simulates l2 in lts
    inline bool may_simulate (int l1, int l2, int lts) const{
        return dominates_in[l1][l2] !=  DOMINATES_IN_NONE &&
	    (dominates_in[l1][l2] == DOMINATES_IN_ALL ||
	     dominates_in[l1][l2] != lts);
    }

    inline int get_lqrel(int l1, int l2, int lts) {
	return lqrel[l1][l2][lts];
    }
    inline bool set_lqrel (int l1, int l2, int lts, int value){	
	assert(value != std::numeric_limits<int>::lowest() + 1);
	assert(dominates_in[l1][l2] == DOMINATES_IN_ALL || dominates_in[l1][l2] != lts);
	/* std::cout << value << " " << lqrel[l1][l2][lts] << std::endl; */
	assert(value <= lqrel[l1][l2][lts]);
	if(value < lqrel[l1][l2][lts]) {
	    lqrel[l1][l2][lts] = value;
	    if (value == std::numeric_limits<int>::lowest()) {
		if(dominates_in[l1][l2] == DOMINATES_IN_ALL){
		    dominates_in[l1][l2] = lts;
		}else if(dominates_in[l1][l2] != lts){
		    dominates_in[l1][l2] = DOMINATES_IN_NONE;
		}
	    }
	    return true;
	}
	return false;
    }

    inline bool set_simulated_by_irrelevant(int l, int lts, int value){
        //Returns if there were changes in dominated_by_noop_in
	assert(value <= simulated_by_irrelevant[l][lts]);
	assert(value != std::numeric_limits<int>::lowest() + 1);
	if(value < simulated_by_irrelevant[l][lts]) {
	    simulated_by_irrelevant[l][lts] = value;
	    if (value == std::numeric_limits<int>::lowest()){
		if(dominated_by_noop_in[l] == DOMINATES_IN_ALL){
		    dominated_by_noop_in[l] = lts;
		}else if(dominated_by_noop_in[l] != lts){
		    dominated_by_noop_in[l] = DOMINATES_IN_NONE;
		}
	    }
	    return true;
	}
        return false;
    }

    inline bool set_simulates_irrelevant(int l, int lts, int value){
	//std::cout << "simulates irrelevant: " << g_operators[l].get_name() << " in " << g_fact_names[lts][0] << ": " << value << std::endl;
	assert(value != std::numeric_limits<int>::lowest() + 1);
        //Returns if there were changes in dominates_noop_in
	assert(value <= simulates_irrelevant[l][lts]);
	if(value < simulates_irrelevant[l][lts]) {
	    simulates_irrelevant[l][lts] = value;

	    if (value == std::numeric_limits<int>::lowest()){
		if(dominates_noop_in[l] == DOMINATES_IN_ALL){
		    dominates_noop_in[l] = lts;
		}else if(dominates_noop_in[l] != lts){
		    dominates_noop_in[l] = DOMINATES_IN_NONE;
		}
	    }
	    return true;
	}
        return false;
    }

    int mix_numbers(const std::vector<int> & values, int lts) const; 

public:
    NumericLabelRelation(Labels * labels);

    //Initializes label relation (only the first time, to reinitialize call reset instead)
    void init(const std::vector<LabelledTransitionSystem *> & lts, const NumericDominanceRelation & sim, const LabelMap & labelMap);

    bool update(const std::vector<LabelledTransitionSystem*> & lts, const NumericDominanceRelation & sim);

    inline int get_num_labels() const {
        return num_labels;
    }

    inline bool may_dominated_by_noop (int l, int lts) const {
        return dominated_by_noop_in[l] == DOMINATES_IN_ALL || dominated_by_noop_in[l] == lts;
    }

    inline bool may_dominate_noop_in (int l, int lts) const {
        return dominates_noop_in[l] == DOMINATES_IN_ALL || dominates_noop_in[l] == lts;
    }

    //Returns true if l dominates l2 in lts (simulates l2 in all j \neq lts)
    inline bool may_dominate (int l1, int l2, int lts) const{
        return dominates_in[l1][l2] == DOMINATES_IN_ALL || (dominates_in[l1][l2] == lts);
    }

    //Returns true if l1 simulates l2 in lts
    int q_dominates (int l1, int l2, int lts) const {
	if (may_dominate(l1, l2, lts)) {
	    int res = mix_numbers(lqrel[l1][l2], lts);
	    assert(res != std::numeric_limits<int>::lowest());
	    return res;
	}else {
	    return std::numeric_limits<int>::lowest();
	}
    }
    
    int q_dominates_noop (int l, int lts) const{
	if (may_dominate_noop_in(l, lts)) {
	    int res = mix_numbers(simulates_irrelevant[l], lts);
	    assert(res != std::numeric_limits<int>::lowest());
	    return res;
	} else {
	    return std::numeric_limits<int>::lowest();
	}
    }
    int q_dominated_by_noop (int l, int lts) const {
	if (may_dominated_by_noop(l, lts)) {
	    int res = mix_numbers(simulated_by_irrelevant[l], lts);
	    assert(res != std::numeric_limits<int>::lowest());
	    assert(res != std::numeric_limits<int>::max());
	    return res;
	} else {
	    return std::numeric_limits<int>::lowest();
	}
    }

    const std::vector<int> & get_tau_labels(int lts) const {
	return tau_labels[lts];
    }
};

#endif
