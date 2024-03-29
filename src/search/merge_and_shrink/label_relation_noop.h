#ifndef MERGE_AND_SHRINK_LABEL_RELATION_NOOP_H
#define MERGE_AND_SHRINK_LABEL_RELATION_NOOP_H


#include "label_relation.h"

/* 
 * Label relation represents the preorder relations on labels that
 * occur in a set of LTS
 */ 
class LabelRelationNoop {
    Labels * labels;
    int num_labels;
    //Indicates whether labels are dominated by noop or other irrelevant
    //variables in theta
    std::vector<std::vector<bool> > simulated_by_irrelevant;
    std::vector<std::vector<bool> > simulates_irrelevant;

    std::vector<int> dominated_by_noop_in;
    //std::vector<int> dominates_noop_in;


    bool update(int i, const LabelledTransitionSystem * lts, 
		const SimulationRelation & sim);
    /* bool update(int i, const LTSComplex  * lts,  */
    /*     	const SimulationRelation & sim); */

    void merge_systems(int system_one, int system_two);
    void merge_labels();

    //void update_after_merge(int i, int j, LTS & lts);

    inline bool set_not_simulated_by_irrelevant(int l, int lts){
        //Returns if there were changes in dominated_by_noop_in
        simulated_by_irrelevant[l][lts] = false;
        if(dominated_by_noop_in[l] == DOMINATES_IN_ALL){
            dominated_by_noop_in[l] = lts;
            return true;
        }else if(dominated_by_noop_in[l] != lts){
            dominated_by_noop_in[l] = DOMINATES_IN_NONE;
            return true;
        }
        return false;
    }

public:
    LabelRelationNoop(Labels * labels);

    //Initializes label relation (only the first time, to reinitialize call reset instead)
    void init(const std::vector<LabelledTransitionSystem *> & lts,
            const DominanceRelation & sim,
            const LabelMap & labelMap);

    //Initializes label relation (only the first time, to reinitialize call reset instead)
    /* void init(const std::vector<LTSComplex *> & lts, */
    /*         const DominanceRelation & sim, */
    /*         const LabelMap & labelMap); */

    /* //Initializes label relation (only the first time, to reinitialize call reset instead) */
    /* void init_identity(int num_lts, const LabelMap & labelMap); */


    void reset();
    bool update(const std::vector<LabelledTransitionSystem*> & lts,
            const DominanceRelation & sim);
    /* bool update(const std::vector<LTSComplex*> & lts, */
    /*         const DominanceRelation & sim); */


    void dump() const;
    void dump(int label) const;
    void dump_equivalent() const;
    void dump_dominance() const;


    inline int get_num_labels() const {
        return num_labels;
    }
    void prune_operators();

    std::vector<int> get_labels_dominated_in_all() const;

    inline int get_dominated_by_noop_in (int l) const {
        return dominated_by_noop_in[l];
    }

    inline bool dominated_by_noop (int l, int lts) const {
        return dominated_by_noop_in[l] == DOMINATES_IN_ALL ||
                dominated_by_noop_in[l] == lts;
//        return dominated_by_noop_in[l] != DOMINATES_IN_NONE &&
//                (dominated_by_noop_in[l] == DOMINATES_IN_ALL ||
//                        dominated_by_noop_in[l] == lts);
    }

    //Returns true if l dominates l2 in lts (simulates l2 in all j \neq lts)
    inline bool dominates (int l1, int l2, int /*lts*/) const{
        return l1 == l2;
    }

    //dangerousLTSs returns the set of LTSs where labels that dominate each other could not be included in the equivalence relation. 
    EquivalenceRelation* get_equivalent_labels_relation(const LabelMap & labelMap, std::set<int> &  dangerous_LTSs) const;

    bool propagate_transition_pruning(int lts_id, 
				      const std::vector<LabelledTransitionSystem *> & ltss, 
				      const DominanceRelation & simulations, 
				      int src, int l1, int target) const; 


    void kill_label(int l) {
	dominated_by_noop_in[l] = DOMINATES_IN_NONE;
 	fill(begin(simulated_by_irrelevant[l]), end(simulated_by_irrelevant[l]), false);
 	fill(begin(simulates_irrelevant[l]), end(simulates_irrelevant[l]), false);
    }

};

    

#endif
