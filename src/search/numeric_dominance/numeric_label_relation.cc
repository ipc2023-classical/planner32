#include "../merge_and_shrink/labels.h"
#include "numeric_simulation_relation.h"
#include "numeric_dominance_relation.h"
#include "../merge_and_shrink/labelled_transition_system.h"
#include "../globals.h"
#include "../utilities.h"
#include "../debug.h"

using namespace std;

template <typename T>
NumericLabelRelation<T>::NumericLabelRelation(Labels * _labels, bool compute_tau_labels_with_noop_dominance_, bool compute_tau_labels_as_self_loops_everywhere_) : compute_tau_labels_with_noop_dominance (compute_tau_labels_with_noop_dominance_),compute_tau_labels_as_self_loops_everywhere (compute_tau_labels_as_self_loops_everywhere_),  labels (_labels), num_labels(_labels->get_size()){

}

template <typename T> 
void NumericLabelRelation<T>::init(const std::vector<LabelledTransitionSystem *> & lts,
				const NumericDominanceRelation<T> & sim,
				const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();
    num_ltss = lts.size();

    cout << "Init label dominance: " << num_labels
	      << " labels " << lts.size() << " systems." << endl;


    std::vector<std::vector<int> >().swap(position_of_label);
    std::vector<std::vector<std::vector<T> > >().swap(lqrel);
    std::vector<std::vector<T> >().swap(simulates_irrelevant);
    std::vector<std::vector<T> >().swap(simulated_by_irrelevant);

    irrelevant_labels_lts.resize(lts.size());
    position_of_label.resize(lts.size());
    simulates_irrelevant.resize(lts.size());
    simulated_by_irrelevant.resize(lts.size());
    lqrel.resize(lts.size());

    for (int i = 0; i < num_ltss; ++i){
        position_of_label[i].resize(num_labels, -1);
	irrelevant_labels_lts[i] = lts[i]->get_irrelevant_labels();

	int num_relevant_labels = 0;
	for(int l = 0; l < num_labels; l++) {
	    if(lts[i]->is_relevant_label(l)){
		position_of_label[i][l] = num_relevant_labels++;
	    }
	}

	cout << "Relevant labels: " << num_relevant_labels << endl;
	// for(int l : lts[i]->get_relevant_labels()) {
	//     assert(position_of_label[i][l]  >= 0);
	// }

	simulates_irrelevant[i].resize(num_relevant_labels, std::numeric_limits<int>::max());
        simulated_by_irrelevant[i].resize(num_relevant_labels, std::numeric_limits<int>::max());
	lqrel[i].resize(num_relevant_labels);

	for(int j = 0; j < num_relevant_labels; j++) {
	    lqrel[i][j].resize(num_relevant_labels, std::numeric_limits<int>::max());
	    lqrel[i][j][j] = 0;
	}
    }

    cout << "Dominating." << endl;
    std::vector<std::vector<int> > ().swap(dominates_in);
    std::vector<int> ().swap(dominated_by_noop_in);
    std::vector<int> ().swap(dominates_noop_in);
    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    dominates_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    
    if(num_labels < 7500) { // If we have more than 5000 labels, there is not enough space. 
	dominates_in.resize(num_labels);
	for (int l1 = 0; l1 < dominates_in.size(); ++l1){
	    dominates_in[l1].resize(num_labels, DOMINATES_IN_ALL);
	}
    }
    cout << "Update label dominance: " << num_labels
	      << " labels " << lts.size() << " systems." << endl;
    
    for (int i = 0; i < num_ltss; ++i){
        update(i, lts[i], sim[i]);
    }
	   

    cout << "Computing tau labels for " << lts.size() << endl;
    tau_labels_changed.resize(lts.size(), true);
    tau_labels.resize(lts.size());
    if(compute_tau_labels_with_noop_dominance) {
	cout << "Compute tau labels with noop dominance" << endl;
	for(int l = 0; l < num_labels; l++){
	    if(dominates_noop_in[l] == DOMINATES_IN_ALL){
		for (int lts_id = 0; lts_id < num_ltss; ++lts_id){
		    tau_labels[lts_id].push_back(l);
		}
	    } else if (dominates_noop_in[l] >= 0) {
		tau_labels[dominates_noop_in[l]].push_back(l);
	    } 
	}
    } else if (compute_tau_labels_as_self_loops_everywhere) {
	int num_tau_labels = 0;
	for(int l = 0; l < num_labels; l++) {
	    int transition_system_relevant = -1;
	    for (int lts_id = 0; lts_id < num_ltss; ++lts_id){
		if(lts[lts_id]->is_relevant_label(l)) {
		    if(transition_system_relevant == -1) {
			transition_system_relevant = lts_id;
		    }else {
			transition_system_relevant = -2;
			break;
		    }
		}
	    }
	    //TODO: check why there are labels irrelevant everywhere
	    // assert(transition_system_relevant != -1);
	    if(transition_system_relevant >= 0) {
		num_tau_labels ++;
		tau_labels[transition_system_relevant].push_back(l);
	    }
	}
	cout << "Computed tau labels as self-loops everywhere: " << num_tau_labels  << " / " << num_labels << endl;

    } else {
     	cout << "No tau labels" << endl;
    }

    // for (int lts_id = 0; lts_id < lts.size(); ++lts_id) {
    // 	cout << "Number of tau labels: " << tau_labels[lts_id].size() << endl;
    // }
}


template <typename T> 
bool NumericLabelRelation<T>::update(const std::vector<LabelledTransitionSystem*> & lts,
				     const NumericDominanceRelation<T> & sim) {
    bool changes = false;
    
    for (int i = 0; i < lts.size(); ++i){
        changes |= update(i, lts[i], sim[i]);
    }

    if(compute_tau_labels_with_noop_dominance) {
	for (int lts_id = 0; lts_id < lts.size(); ++lts_id) {
	    tau_labels[lts_id].erase(remove_if(tau_labels[lts_id].begin(),
					   tau_labels[lts_id].end(), [&](int label) {
					       return dominates_noop_in[label] != DOMINATES_IN_ALL &&
						   dominates_noop_in[label] != lts_id;
					       }),
				     tau_labels[lts_id].end());
	}
    }				 
		
    return changes;
}

template <typename T> 
bool NumericLabelRelation<T>::update(int lts_i, const LabelledTransitionSystem * lts,  
				     const NumericSimulationRelation<T> & sim){
    bool changes = false;
    for(int l2 : lts->get_relevant_labels()) {
	assert(lts->is_relevant_label(l2));
        for(int l1 : lts->get_relevant_labels()) {
	    assert(lts->is_relevant_label(l1));
            if(l1 != l2 && may_simulate(l1, l2, lts_i)){
		T min_value = std::numeric_limits<int>::max();
                //std::cout << "Check " << l1 << " " << l2 << std::endl;
                //std::cout << "Num transitions: " << lts->get_transitions_label(l1).size()
                //		    << " " << lts->get_transitions_label(l2).size() << std::endl;
                //Check if it really simulates
                //For each transition s--l2-->t, and every label l1 that dominates
                //l2, exist s--l1-->t', t <= t'?

                for(const auto & tr : lts->get_transitions_label(l2)){
		    T max_value = std::numeric_limits<int>::lowest();
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
	T old_value = get_simulated_by_irrelevant(l2, lts_i);
	if(old_value != T(std::numeric_limits<int>::lowest())) {
	    T min_value = std::numeric_limits<int>::max();
	    for(auto tr : lts->get_transitions_label(l2)){
		min_value = std::min(min_value, sim.q_simulates(tr.src, tr.target));
		if(min_value == std::numeric_limits<int>::lowest()) {
		    break;
		}
	    }

	    assert(min_value != std::numeric_limits<int>::max());

	    if (min_value < old_value) {
		changes |= set_simulated_by_irrelevant(l2, lts_i, min_value);
		// for (int l : lts->get_irrelevant_labels()){
		//     changes |= set_lqrel(l, l2, lts_i, min_value);
		// }
		old_value = min_value;
	    }
        }

        //Does l2 simulates irrelevant_labels in lts?
	old_value = get_simulates_irrelevant(l2, lts_i);
	if(old_value != std::numeric_limits<int>::lowest()) {
	    T min_value = std::numeric_limits<int>::max();
	    for(int s = 0; s < lts->size(); s++){
		T max_value = std::numeric_limits<int>::lowest();
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
		// for (int l : lts->get_irrelevant_labels()){
		//     changes |= set_lqrel(l2, l, lts_i, min_value);
		// }
	    }
	}
    }

    // for (int l : lts->get_irrelevant_labels()) {
    // 	set_simulates_irrelevant(l, lts_i, 0);
    // 	set_simulated_by_irrelevant(l, lts_i, 0);
    // }

    return changes;
}


template <typename T> 
void NumericLabelRelation<T>::dump(const LabelledTransitionSystem * lts, int lts_id) const {
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
	if (may_dominated_by_noop (l2, lts_id))  {
	    cout <<  g_operators[l2].get_name() << " dominates noop: " << q_dominated_by_noop(l2, lts_id) << endl;

	}
	if (may_dominate_noop_in (l2, lts_id))  {
	    cout <<  g_operators[l2].get_name() << " dominated by noop: " << q_dominates_noop(l2, lts_id) << endl;
	}
    }

}



template class NumericLabelRelation<int>; 
template class NumericLabelRelation<IntEpsilon>; 
