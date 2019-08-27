#ifndef MERGE_AND_SHRINK_LABEL_RELATION_IDENTITY_H
#define MERGE_AND_SHRINK_LABEL_RELATION_IDENTITY_H

#include <iostream>
#include <vector>
#include "labels.h"
#include "label.h"
#include "label_relation.h" //TODO: For DOMINATED_IN_NONE

class EquivalenceRelation;
/* class LTSComplex; */
class LabelledTransitionSystem;
class SimulationRelation;
class DominanceRelation;

/* 
 * Label relation represents the preorder relations on labels that
 * occur in a set of LTS
 */ 
class LabelRelationIdentity {
    Labels * labels;
    int num_labels;

public:
    LabelRelationIdentity (Labels * labels);

    void init(const std::vector<LabelledTransitionSystem *> & /*lts*/,
	      const DominanceRelation & /*sim*/,
	      const LabelMap & /*labelMap*/){}

    /* void init(const std::vector<LTSComplex *> & /\*lts*\/, */
    /*           const DominanceRelation & /\*sim*\/, */
    /*           const LabelMap & /\*labelMap*\/) {} */

    void reset() {}
    bool update(const std::vector<LabelledTransitionSystem*> & /*lts*/,
		const DominanceRelation & /*sim*/) {return false;}
    /* bool update(const std::vector<LTSComplex*> & /\*lts*\/, */
    /*     	const DominanceRelation & /\*sim*\/) {return false;} */


    void dump() const {}
    void dump(int /*label*/) const {}
    void dump_equivalent() const {}
    void dump_dominance() const {}


    inline int get_num_labels() const {
        return num_labels;
    }

    inline int get_dominated_by_noop_in (int /*l*/) const {
        return DOMINATES_IN_NONE;
    }

    inline bool dominated_by_noop (int /*l*/, int /*lts*/) const {
	return false;
    }

    inline bool dominates (int l1, int l2, int /*lts*/) const{
	return l1 == l2;
    }

    std::vector<int>  get_labels_dominated_in_all() const { 
	return std::vector<int>(); 
    }

    void kill_label(int /*l*/) {}



    void prune_operators(){}

    //dangerousLTSs returns the set of LTSs where labels that dominate each other could not be included in the equivalence relation. 
    EquivalenceRelation* get_equivalent_labels_relation(const LabelMap & labelMap, std::set<int> &  dangerous_LTSs) const;

    bool propagate_transition_pruning(int lts_id, 
				      const std::vector<LabelledTransitionSystem *> & ltss, 
				      const DominanceRelation & simulations, 
				      int src, int l1, int target) const; 

};

    

#endif
