#ifndef MERGE_AND_SHRINK_ALTERNATIVE_LABEL_RELATION_H
#define MERGE_AND_SHRINK_ALTERNATIVE_LABEL_RELATION_H

#include "label_relation.h"
#include "../utilities.h"

#include <iostream>
#include <vector>
#include "labels.h"
#include "label.h"

class LabelledTransitionSystem;

/* 
 * Label relation represents the preorder relations on labels that
 * occur in a set of LTS
 */ 
class AlternativeLabelRelation {
    Labels * labels;
    int num_labels;
    int num_ltss;

    // Summary matrix for each l1, l2 indicating whether l1 dominates
    // l2 in all (-2), in none (-1) or only in i (i)
    std::vector<std::vector<int> > dominates_in;
    std::vector<int> dominates_noop_in, dominated_by_noop_in;

    //For each lts, matrix indicating whether l1 simulates l2, noop
    //simulates l or l simulates noop
    std::vector<int> cost_of_label;
    std::vector<std::vector<int> > position_of_label; //position that label l takes on lts
    std::vector<std::vector<int> > irrelevant_labels_lts;
    std::vector<std::vector<std::vector<bool> > > lrel;
    std::vector<std::vector<bool> > simulated_by_irrelevant;
    std::vector<std::vector<bool> > simulates_irrelevant;

    bool update(int i, const LabelledTransitionSystem * lts, 
		const SimulationRelation & sim);

    inline bool get_lrel(int l1, int l2, int lts) const {
	int pos1 = position_of_label[lts][l1];
	int pos2 = position_of_label[lts][l2];
	if(pos1 >= 0) {
	    if(pos2 >= 0) {
		return lrel[lts][pos1][pos2];
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

    inline bool set_not_simulates (int l1, int l2, int lts){	
	if(lrel[lts][position_of_label[lts][l1]][position_of_label[lts][l2]]) {
	    lrel[lts][position_of_label[lts][l1]][position_of_label[lts][l2]] = false;
	    if (!dominates_in.empty()) {
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

    inline bool set_not_simulated_by_irrelevant(int l, int lts){
        //Returns if there were changes in dominated_by_noop_in
	assert(position_of_label[lts][l] >= 0);
	int pos = position_of_label[lts][l];

	if(simulated_by_irrelevant[lts][pos]) {
	    simulated_by_irrelevant[lts][pos] = false;

            if(dominated_by_noop_in[l] == DOMINATES_IN_ALL){
                dominated_by_noop_in[l] = lts;
            }else if(dominated_by_noop_in[l] != lts){
                dominated_by_noop_in[l] = DOMINATES_IN_NONE;
            }
            
            if(!dominates_in.empty()) {
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

    
    inline bool set_not_simulates_irrelevant(int l, int lts){
	int pos = position_of_label[lts][l];
        //Returns if there were changes in dominates_noop_in
	if(simulates_irrelevant[lts][pos]) {
	    simulates_irrelevant[lts][pos] = false;

            if(dominates_noop_in[l] == DOMINATES_IN_ALL){
                dominates_noop_in[l] = lts;
            }else if(dominates_noop_in[l] != lts){
                dominates_noop_in[l] = DOMINATES_IN_NONE;
            }
            if(!dominates_in.empty()) {
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



    
    inline bool get_simulated_by_irrelevant(int l, int lts) const {
	if(position_of_label[lts][l] >= 0) {
	    return simulated_by_irrelevant[lts][position_of_label[lts][l]];
	} else {
	    return 0;
	}
    }

    inline bool get_simulates_irrelevant(int l, int lts) const {
	if(position_of_label[lts][l] >= 0) {
	    return simulates_irrelevant[lts][position_of_label[lts][l]];
	} else {
	    return 0;
	}
    }

public:
    AlternativeLabelRelation(Labels * labels);

    //Initializes label relation (only the first time, to reinitialize call reset instead)
    void init(const std::vector<LabelledTransitionSystem *> & lts, 
	      const DominanceRelation & sim, 
	      const LabelMap & labelMap);

    bool update(const std::vector<LabelledTransitionSystem*> & lts, 
		const DominanceRelation & sim);

    void init(const std::vector<LTSComplex *> & ,
              const DominanceRelation & ,
              const LabelMap & ){
        std::cout << "LTSComplex not implemented." << std::endl;
        exit(EXIT_CRITICAL_ERROR);
    }

    bool update(const std::vector<LTSComplex*> & ,
                const DominanceRelation & ){
       std::cout << "LTSComplex not implemented." << std::endl;
        exit(EXIT_CRITICAL_ERROR);
    }


    inline int get_num_labels() const {
        return num_labels;
    }

    inline bool dominated_by_noop (int l, int lts) const {
        return dominated_by_noop_in[l] == DOMINATES_IN_ALL || dominated_by_noop_in[l] == lts;
    }

    inline int get_dominated_by_noop_in (int l) const {
        return dominated_by_noop_in[l];
    }


    inline bool dominates_noop (int l, int lts) const {
        return dominates_noop_in[l] == DOMINATES_IN_ALL || dominates_noop_in[l] == lts;
    }

    //Returns true if l1 simulates l2 in lts
    inline bool simulates (int l1, int l2, int lts) const{
	if(dominates_in.empty()) {
	    return get_lrel(l1, l2, lts);
	}
        return dominates_in[l1][l2] !=  DOMINATES_IN_NONE &&
	    (dominates_in[l1][l2] == DOMINATES_IN_ALL ||
	     dominates_in[l1][l2] != lts);
    }

    //Returns true if l1 simulates l2 in lts
    bool dominates (int l1, int l2, int lts) const {
        if(dominates_in.empty()) {
	    for(int lts_id = 0; lts_id < num_ltss; ++lts_id) {
		if(lts_id != lts && !get_lrel(l1, l2, lts_id)) {
		    assert(num_ltss > 1);
		    return false;
		}
	    }
	    return true;
	}
	
	assert(num_ltss > 1 || dominates_in[l1][l2] == DOMINATES_IN_ALL || (dominates_in[l1][l2] == lts));
        return dominates_in[l1][l2] == DOMINATES_IN_ALL || (dominates_in[l1][l2] == lts);
    }
    

    int get_label_cost (int label) const {
	return cost_of_label[label];
    }

    void dump(const LabelledTransitionSystem * lts, int lts_id) const;

    
    bool propagate_transition_pruning(int , 
				      const std::vector<LabelledTransitionSystem *> & , 
				      const DominanceRelation & , 
				      int , int , int ) const{
        std::cout << "propagate_transition_pruning not implemented." << std::endl;
        exit(EXIT_CRITICAL_ERROR);
        return false;
    }


    void kill_label(int ) {
        std::cout << "kill_label not implemented." << std::endl;
        exit(EXIT_CRITICAL_ERROR);
    }

    std::vector<int> get_labels_dominated_in_all() const;

    EquivalenceRelation* get_equivalent_labels_relation(const LabelMap & , std::set<int> & ) const {
        std::cout << "get_equivalent_labels_relation." << std::endl;
        exit(EXIT_CRITICAL_ERROR);
    }


};

#endif
