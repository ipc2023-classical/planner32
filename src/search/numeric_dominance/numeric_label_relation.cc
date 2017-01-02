#include "../merge_and_shrink/labels.h"
#include "numeric_simulation_relation.h"
#include "numeric_dominance_relation.h"
#include "../merge_and_shrink/labelled_transition_system.h"
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

    std::vector<std::vector<std::vector<int> > >().swap(lqrel);
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

    for(int i = 0; i < num_labels; i++) {
        simulates_irrelevant[i].resize(lts.size(), std::numeric_limits<int>::max());
        simulated_by_irrelevant[i].resize(lts.size(), std::numeric_limits<int>::max());
	lqrel[i].resize(num_labels);
	for(int j = 0; j < num_labels; j++) {
	    lqrel[i][j].resize(lts.size(), i == j ? 0 : std::numeric_limits<int>::max());
	}
    }
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
        dominates_in[l1].resize(num_labels, DOMINATES_IN_ALL);
    }
    cout << "Update label dominance: " << num_labels
	      << " labels " << lts.size() << " systems." << endl;
    
    for (int i = 0; i < lts.size(); ++i){
        update(i, lts[i], sim[i]);
    }

    cout << "Computing tau labels for " << lts.size() << endl;
    tau_labels.resize(lts.size());
    for(int l = 0; l < num_labels; l++){
	if(dominates_noop_in[l] == DOMINATES_IN_ALL){
	    for (int lts_id = 0; lts_id < lts.size(); ++lts_id){
		tau_labels[lts_id].push_back(l);
	    }
	} else if (dominates_noop_in[l] >= 0) {
	    tau_labels[dominates_noop_in[l]].push_back(l);
	} 
    }
    for (int lts_id = 0; lts_id < lts.size(); ++lts_id) {
	cout << "Number of tau labels: " << tau_labels[lts_id].size() << endl;
    }
}



bool NumericLabelRelation::update(const std::vector<LabelledTransitionSystem*> & lts,
				  const NumericDominanceRelation & sim) {
    bool changes = false;
    
    for (int i = 0; i < lts.size(); ++i){
        changes |= update(i, lts[i], sim[i]);
    }

    for (int lts_id = 0; lts_id < lts.size(); ++lts_id) {
	tau_labels[lts_id].erase(remove_if(tau_labels[lts_id].begin(),
					   tau_labels[lts_id].end(), [&](int label) {
					       return dominates_noop_in[label] != DOMINATES_IN_ALL &&
						   dominates_noop_in[label] != lts_id;
					   }),
				 tau_labels[lts_id].end());
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
			    if(max_value >= min_value) {
				break; //Stop checking this tr
			    }
                        }
                    }
		    min_value = std::min(min_value, max_value);
		    if(min_value == std::numeric_limits<int>::lowest()) {
			break; //Stop checking trs of l1, l2
		    }
                }
		changes |= set_lqrel(l1, l2, lts_i, min_value);		    
		assert(min_value != std::numeric_limits<int>::max());
            }
        }

        //Is l2 simulated by irrelevant_labels in lts?
	int old_value = simulated_by_irrelevant[l2][lts_i];
	if(old_value != std::numeric_limits<int>::lowest()) {
	    int min_value = std::numeric_limits<int>::max();
	    for(auto tr : lts->get_transitions_label(l2)){
		min_value = std::min(min_value, sim.q_simulates(tr.src, tr.target));
		if(min_value == std::numeric_limits<int>::lowest()) {
		    break;
		}
	    }

	    assert(min_value != std::numeric_limits<int>::max());

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
	if(old_value != std::numeric_limits<int>::lowest()) {
	    int min_value = std::numeric_limits<int>::max();
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
	    }
	    assert(min_value != std::numeric_limits<int>::max());
	    if(min_value < old_value) {
		old_value = min_value;
		changes |= set_simulates_irrelevant(l2, lts_i, min_value);
		for (int l : lts->get_irrelevant_labels()){
		    changes |= set_lqrel(l2, l, lts_i, min_value);
		}
	    }
	}
    }

    for (int l : lts->get_irrelevant_labels()) {
	set_simulates_irrelevant(l, lts_i, 0);
	set_simulated_by_irrelevant(l, lts_i, 0);
    }

    return changes;
}

int NumericLabelRelation::mix_numbers (const std::vector<int> & values, int lts) const {
    int sum_negatives = 0;
    int max_positive = 0;
	    
    for(int lts_id = 0; lts_id < values.size(); ++lts_id) {
	if(lts_id == lts) continue;
	int val = values[lts_id];
	// assert(val != std::numeric_limits<int>::max());
	if(val < 0) {
	    sum_negatives += val;
	} else {
	    max_positive = std::max(max_positive, val);
	}
    }
    
    return sum_negatives; //+ max_positive;    
}


void NumericLabelRelation::dump(const LabelledTransitionSystem * lts, int lts_id) const {
    for(int l2 : lts->get_relevant_labels()) {
        for(int l1 : lts->get_relevant_labels()) {
	    if(l2 == l1) {
		continue;
	    }
	    if (may_dominate(l2, l1, lts_id)) {
		cout << g_operators[l1].get_name() << " <= " << 
		    g_operators[l2].get_name() << " with " 
		     << q_dominates(l2, l1, lts_id) << endl;
	    }
	}
    }
}
