#include "label_relation.h"

#include "labels.h"
#include "simulation_relation.h"
#include "dominance_relation.h"
#include "labelled_transition_system.h"
#include "lts_complex.h"
#include "../equivalence_relation.h"
#include "../globals.h"
#include "../utilities.h"
#include "../debug.h"

using namespace std;

NumericLabelRelation::NumericLabelRelation(Labels * _labels) : labels (_labels), num_labels(_labels->get_size()){

}

void NumericLabelRelation::init(const std::vector<LabelledTransitionSystem *> & lts,
				const NumericDominanceRelation & sim,
				const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();

    std::vector<std::vector<int> >().swap(lqrel);
    std::vector<std::vector<int> >().swap(simulates_irrelevant);
    std::vector<std::vector<int> >().swap(simulated_by_irrelevant);
    std::vector<std::vector<int> > ().swap(dominates_in);
    std::vector<int> ().swap(dominated_by_noop_in);
    std::vector<int> ().swap(dominates_noop_in);

    simulates_irrelevant.resize(num_labels);
    simulated_by_irrelevant.resize(num_labels);
    lqrel.resize(num_labels);

    dominates_in.resize(num_labels);
    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    dominates_noop_in.resize(num_labels, DOMINATES_IN_ALL);

    for(int i = 0; i < num_labels; i++){
        simulates_irrelevant[i].resize(lts.size(), 0);
        simulated_by_irrelevant[i].resize(lts.size(), 0);
	lqrel[i].resize(num_labels);
	for(int j = 0; j < num_labels; j++){
	    lqrel[i][j].resize(lts.size(), 0);
	}
    }
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
        dominates_in[l1].resize(num_labels, DOMINATES_IN_ALL);
    }
    DEBUG_MSG(cout << "Update label dominance: " << num_labels
	      << " labels " << lts.size() << " systems." << endl;);
    
    for (int i = 0; i < lts.size(); ++i){
        update(i, lts[i], sim[i]);
    }
}



bool NumericLabelRelation::update(const std::vector<LabelledTransitionSystem*> & lts,
				  const NumericDominanceRelation & sim){
    bool changes = false;
    
    for (int i = 0; i < lts.size(); ++i){
        changes |= update(i, lts[i], sim[i]);
    }
    
    return changes;
}

/* TODO: This version is inefficient. It could be improved by
 * iterating only the right transitions (see TODO inside the loop)
 */
bool NumericLabelRelation::update(int lts_i, const LabelledTransitionSystem * lts,  const NumericSimulationRelation & sim){
    bool changes = false;
    for(int l2 : lts->get_relevant_labels()) {
        for(int l1 : lts->get_relevant_labels()){
            if(l1 != l2 && may_simulate(l1, l2, lts_i)){
		int min_value = std::numeric_limits<int>::max();
                //std::cout << "Check " << l1 << " " << l2 << std::endl;
                //std::cout << "Num transitions: " << lts->get_transitions_label(l1).size()
                //		    << " " << lts->get_transitions_label(l2).size() << std::endl;
                //Check if it really simulates
                //For each transition s--l2-->t, and every label l1 that dominates
                //l2, exist s--l1-->t', t <= t'?
                for(const auto & tr : lts->get_transitions_label(l2)){
		    int max_value = std::numeric_limits<int>::lowest();
                    //TODO: for(auto tr2 : lts->get_transitions_for_label_src(l1, tr.src)){
                    for(const auto & tr2 : lts->get_transitions_label(l1)){
                        if(tr2.src == tr.src && sim.may_simulate(tr2.target, tr.target)){
			    max_value = std::max(max_value, sim.q_simulates(tr2.target, tr.target));
			    if(max_value >= min_value){
				break; //Stop checking this tr
			    }
                        }
                    }
		    min_value = std::min(min_value, max_value);
		    changes |= set_lqrel(l1, l2, lts_i, min_value);
		    if(min_value == std::numeric_limits<int>::lowest()) {
			break; //Stop checking trs of l1, l2
		    }
                }
		assert(min_value != std::numeric_limits<int>::max());
            }
        }

        //Is l2 simulated by irrelevant_labels in lts?
	int old_value = simulated_by_irrelevant[l2][lts_i];
	if(old_value != std::numeric_limits::lowest()) {
	    int min_value = std::numeric_limits::max();
	    for(auto tr : lts->get_transitions_label(l2)){
		min_value = std::min(min_value, sim.q_simulates(tr.src, tr.target));
		if(min_value == std::numeric_limits::lowest()) {
		    break;
		}
	    }

	    if (min_value < old_value) {
		changes |= set_simulated_by_irrelevant(l2, lts_i, min_value);
		for (int l : lts->get_irrelevant_labels()){
		    changes |= set_lqrel(l, l2, lts_i, min_value);
		}
		old_value = min_value;
	    }
        }

        //Does l2 simulates irrelevant_labels in lts?
	old_value = simulates_irrelevant[l2][lts_i];
	if(old_value != std::numeric_limits::lowest()) {
	    int min_value = std::numeric_limits::max();
	    for(int s = 0; s < lts->size(); s++){
		int max_value = std::numeric_limits<int>::lowest();
		for(const auto & tr : lts->get_transitions_label(l2)) {
		    if(tr.src == s) {
			max_value = std::max(max_value, sim.q_simulates(tr.target, tr.src));
			if(max_value > min_value) {
			    break;
			}
		    }
		}
		min_value = std::min(min_value, max_value);
		if(min_value < old_value) {
		    old_value = min_value;
		    changes |= set_simulates_irrelevant(l2, lts_i, min_value);
		    for (int l : lts->get_irrelevant_labels()){
			changes |= set_lqrel(l2, l, lts_i, min_value);
		    }
		}
	    }
	}
    }

    return changes;
}



int NumericLabelRelation::q_simulates (int l1, int l2, int lts) const {
    if(dominates_in[l1][l2] !=  DOMINATES_IN_NONE &&
       (dominates_in[l1][l2] == DOMINATES_IN_ALL ||
	dominates_in[l1][l2] != lts)) {
	int sum_negatives = 0;
	int max_positive = 0;
	    
	for(int lts_id = 0; lts_id < lqrel[l1][l2].size(); ++lts_id) {
	    if(lts_id == lts) continue;
	    int val = lqrel[l1][l2][lts_id];
	    if(val < 0) {
		sum_negatives += val;
	    } else {
		max_positive = std::max(max_positive, val);
	    }
	}
	return sum_negatives + max_positive;
    }
}

int NumericLabelRelation::q_dominated_by_noop (int l1, int lts) const {
    if(dominated_by_noop_in[l1][l2] !=  DOMINATES_IN_NONE &&
       (dominates_in[l1][l2] == DOMINATES_IN_ALL ||
	dominates_in[l1][l2] != lts)) {
	int sum_negatives = 0;
	int max_positive = 0;
	    
	for(int lts_id = 0; lts_id < lqrel[l1][l2].size(); ++lts_id) {
	    if(lts_id == lts) continue;
	    int val = lqrel[l1][l2][lts_id];
	    if(val < 0) {
		sum_negatives += val;
	    } else {
		max_positive = std::max(max_positive, val);
	    }
	}
	return sum_negatives + max_positive;
    }
}


int NumericLabelRelation::mix_numbers (const std::vector<int> & values, int lts) const {
    int sum_negatives = 0;
    int max_positive = 0;
	    
    for(int lts_id = 0; lts_id < values.size(); ++lts_id) {
	if(lts_id == lts) continue;
	int val = values[lts_id];
	if(val < 0) {
	    sum_negatives += val;
	} else {
	    max_positive = std::max(max_positive, val);
	}
    }
    return sum_negatives + max_positive;    
}


