#ifndef NUMERIC_DOMINANCE_NUMERIC_LABEL_RELATION_H
#define NUMERIC_DOMINANCE_NUMERIC_LABEL_RELATION_H

#include "../merge_and_shrink/label_relation.h"

#include "numeric_label_relation.h"

#include <iostream>
#include <vector>
#include "../merge_and_shrink/labels.h"
#include "../merge_and_shrink/label.h"
#include "int_epsilon.h"

class LabelledTransitionSystem;


template <typename T> class NumericDominanceRelation;
template <typename T> class NumericSimulationRelation;

/* 
 * Label relation represents the preorder relations on labels that
 * occur in a set of LTS
 */ 
template <typename T>
class NumericLabelRelation {
    const int compute_tau_labels_with_noop_dominance; 

    Labels * labels;
    int num_labels;
    int num_ltss;

    // Summary matrix for each l1, l2 indicating whether l1 dominates
    // l2 in all (-2), in none (-1) or only in i (i)
    std::vector<std::vector<int> > dominates_in;
    std::vector<int> dominates_noop_in, dominated_by_noop_in;

    //For each lts, matrix indicating whether l1 simulates l2, noop
    //simulates l or l simulates noop
    std::vector<std::vector<int> > position_of_label; //position that label l takes on lts
    std::vector<std::vector<int> > irrelevant_labels_lts;
    std::vector<std::vector<std::vector<T> > > lqrel;
    std::vector<std::vector<T> > simulated_by_irrelevant;
    std::vector<std::vector<T> > simulates_irrelevant;

    std::vector<std::vector<int> > tau_labels;
//Vector to indicate whether tau_labels changed since last time
    mutable std::vector<bool> tau_labels_changed; 

    bool update(int i, const LabelledTransitionSystem * lts, 
		const NumericSimulationRelation<T> & sim);

    //Returns true if l1 simulates l2 in lts
    inline bool may_simulate (int l1, int l2, int lts) const{
        return dominates_in[l1][l2] !=  DOMINATES_IN_NONE &&
	    (dominates_in[l1][l2] == DOMINATES_IN_ALL ||
	     dominates_in[l1][l2] != lts);
    }

    inline T get_lqrel(int l1, int l2, int lts) const {
	int pos1 = position_of_label[lts][l1];
	int pos2 = position_of_label[lts][l2];
	if(pos1 >= 0) {
	    if(pos2 >= 0) {
		return lqrel[lts][pos1][pos2];
	    }else {
		return simulates_irrelevant[lts][pos1];
	    }	  
	}else {
	    if(pos2 != -1) {
		return simulated_by_irrelevant[lts][pos2];
	    }else {
		return 0; //Both are irrelevant
	    }    
	}
    }

    inline bool set_lqrel (int l1, int l2, int lts, T value){	
	assert(value != std::numeric_limits<int>::lowest() + 1);
	assert(dominates_in[l1][l2] == DOMINATES_IN_ALL || dominates_in[l1][l2] != lts);
	assert(position_of_label[lts][l1] >= 0 &&  position_of_label[lts][l2] >= 0);

	assert(value <= lqrel[lts][position_of_label[lts][l1]][position_of_label[lts][l2]]);
	if(value < lqrel[lts][position_of_label[lts][l1]][position_of_label[lts][l2]]) {
	    lqrel[lts][position_of_label[lts][l1]][position_of_label[lts][l2]] = value;
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

    inline T get_simulated_by_irrelevant(int l, int lts) const {
	if(position_of_label[lts][l] >= 0) {
	    return simulated_by_irrelevant[lts][position_of_label[lts][l]];
	} else {
	    return 0;
	}
    }

    inline bool set_simulated_by_irrelevant(int l, int lts, T value){
        //Returns if there were changes in dominated_by_noop_in
	
	assert(position_of_label[lts][l] >= 0);
	int pos = position_of_label[lts][l];

	assert(value <= simulated_by_irrelevant[lts][pos]);
	assert(value != std::numeric_limits<int>::lowest() + 1);
	if(value < simulated_by_irrelevant[lts][pos]) {
	    simulated_by_irrelevant[lts][pos] = value;
	    if (value == std::numeric_limits<int>::lowest()) {
		if(dominated_by_noop_in[l] == DOMINATES_IN_ALL){
		    dominated_by_noop_in[l] = lts;
		}else if(dominated_by_noop_in[l] != lts){
		    dominated_by_noop_in[l] = DOMINATES_IN_NONE;
		}
		for(int l1 : irrelevant_labels_lts[lts]) {
		    if(dominates_in[l1][l] == DOMINATES_IN_ALL){
			dominates_in[l1][l] = lts;
		    }else if(dominates_in[l1][l] != lts){
			dominates_in[l1][l] = DOMINATES_IN_NONE;
		    }
		}
	    }
	    return true;
	}
        return false;
    }


    inline T get_simulates_irrelevant(int l, int lts) const {
	if(position_of_label[lts][l] >= 0) {
	    return simulates_irrelevant[lts][position_of_label[lts][l]];
	} else {
	    return 0;
	}
    }

    inline bool set_simulates_irrelevant(int l, int lts, T value){
	//std::cout << "simulates irrelevant: " << g_operators[l].get_name() << " in " << g_fact_names[lts][0] << ": " << value << std::endl;
	assert(value != std::numeric_limits<int>::lowest() + 1);
	
	assert(position_of_label[lts][l] >= 0); 
	int pos = position_of_label[lts][l];

        //Returns if there were changes in dominates_noop_in
	assert(value <= simulates_irrelevant[lts][pos]);
	if(value < simulates_irrelevant[lts][pos]) {
	    simulates_irrelevant[lts][pos] = value;

	    if (value == std::numeric_limits<int>::lowest()){
		if(dominates_noop_in[l] == DOMINATES_IN_ALL){
		    dominates_noop_in[l] = lts;
		}else if(dominates_noop_in[l] != lts){
		    dominates_noop_in[l] = DOMINATES_IN_NONE;
		}
		for(int l2 : irrelevant_labels_lts[lts]) {
		    if(dominates_in[l][l2] == DOMINATES_IN_ALL){
			dominates_in[l][l2] = lts;
		    }else if(dominates_in[l][l2] != lts){
			dominates_in[l][l2] = DOMINATES_IN_NONE;
		    }
		}
	    }
	    return true;
	}
        return false;
    }

    /* int mix_numbers(const std::vector<int> & values,  */
    /* 		    const std::vector<int> & values_irrelevant_labels,  */
    /* 		    int lts) const;  */

public:
    NumericLabelRelation(Labels * labels, bool compute_tau_labels_with_noop_dominance);

    //Initializes label relation (only the first time, to reinitialize call reset instead)
    void init(const std::vector<LabelledTransitionSystem *> & lts, 
	      const NumericDominanceRelation<T> & sim, 
	      const LabelMap & labelMap);

    bool update(const std::vector<LabelledTransitionSystem*> & lts, 
		const NumericDominanceRelation<T> & sim);

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
    T q_dominates (int l1, int l2, int lts) const {
	if (may_dominate(l1, l2, lts)) {
	    T total_sum = 0;
	    
	    for(int lts_id = 0; lts_id < num_ltss; ++lts_id) {
		if(lts_id != lts) {
		    assert (get_lqrel(l1, l2, lts_id) != std::numeric_limits<int>::lowest());
		    total_sum += get_lqrel(l1, l2, lts_id);
		}
	    }
    
	    return total_sum;    
	}else {
	    assert(false);
	    return std::numeric_limits<int>::lowest();
	}
    }
    
    T q_dominates_noop (int l, int lts) const{
	if (may_dominate_noop_in(l, lts)) {
	    T total_sum = 0;
	    
	    for(int lts_id = 0; lts_id < num_ltss; ++lts_id) {
		if(lts_id != lts) {
		    assert (get_simulates_irrelevant(l, lts_id) != std::numeric_limits<int>::lowest());
		    total_sum += get_simulates_irrelevant(l, lts_id);
		}
	    }
	    return total_sum;    
	} else {
	    assert(false);
	    return std::numeric_limits<int>::lowest();
	}
    }

    T q_dominated_by_noop (int l, int lts) const {
	if (may_dominated_by_noop(l, lts)) {
	    T total_sum = 0;
	    
	    for(int lts_id = 0; lts_id < num_ltss; ++lts_id) {
		if(lts_id != lts) {
		    assert (get_simulated_by_irrelevant(l, lts_id) != std::numeric_limits<int>::lowest());
		    total_sum += get_simulated_by_irrelevant(l, lts_id);
		}
	    }
	    return total_sum;    

	} else {
	    assert(false);
	    return std::numeric_limits<int>::lowest();
	}
    }

    const std::vector<int> & get_tau_labels(int lts) const {
	return tau_labels[lts];
    }

    bool have_tau_labels_changed (int lts) const {
	bool result = tau_labels_changed[lts];
	tau_labels_changed[lts] = false;
	return result;
    }

    void dump(const LabelledTransitionSystem * lts, int lts_id) const;
};

#endif
