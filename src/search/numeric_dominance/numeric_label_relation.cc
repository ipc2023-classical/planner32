#include "numeric_label_relation.h"
#include "numeric_simulation_relation.h"
#include "numeric_dominance_relation.h"
#include "../merge_and_shrink/labels.h"
#include "../merge_and_shrink/labelled_transition_system.h"
#include "../globals.h"
#include "../utilities.h"
#include "../debug.h"

using namespace std;

template <typename T>
NumericLabelRelation<T>::NumericLabelRelation(Labels * _labels, int num_labels_to_use_dominates_in_) :
    labels (_labels), num_labels(_labels->get_size()), num_labels_to_use_dominates_in(num_labels_to_use_dominates_in_) 
{
    
}

template <typename T> 
bool NumericLabelRelation<T>::update(int lts_i, const LabelledTransitionSystem * lts,  
				     const NumericSimulationRelation<T> & sim){
    bool changes = false;
    for(LabelGroup lg2 (0); lg2.group < lts->get_num_label_groups(); ++lg2) {      
	for(LabelGroup lg1(0); lg1.group < lts->get_num_label_groups(); ++lg1) {
            if(lg1 != lg2 && may_simulate(lg1, lg2, lts_i)){
		T min_value = std::numeric_limits<int>::max();
                //std::cout << "Check " << l1 << " " << l2 << std::endl;
                //std::cout << "Num transitions: " << lts->get_transitions_label(l1).size()
                //		    << " " << lts->get_transitions_label(l2).size() << std::endl;
                //Check if it really simulates
                //For each transition s--l2-->t, and every label l1 that dominates
                //l2, exist s--l1-->t', t <= t'?

                for(const auto & tr : lts->get_transitions_label_group(lg2)){
		    T max_value = std::numeric_limits<int>::lowest();
                    //TODO: for(auto tr2 : lts->get_transitions_for_label_src(l1, tr.src)){
                    for(const auto & tr2 : lts->get_transitions_label_group(lg1)){
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
		
		changes |= set_lqrel(lg1, lg2, lts_i, lts, min_value);		    
		assert(min_value != std::numeric_limits<int>::max());
            }
        }

        //Is l2 simulated by irrelevant_labels in lts?
	T old_value = get_simulated_by_irrelevant(lg2, lts_i);
	if(old_value != T(std::numeric_limits<int>::lowest())) {
	    T min_value = std::numeric_limits<int>::max();
	    for(const auto & tr : lts->get_transitions_label_group(lg2)){
		min_value = std::min(min_value, sim.q_simulates(tr.src, tr.target));
		if(min_value == std::numeric_limits<int>::lowest()) {
		    break;
		}
	    }

	    assert(min_value != std::numeric_limits<int>::max());

	    if (min_value < old_value) {
		changes |= set_simulated_by_irrelevant(lg2, lts_i, lts,  min_value);
		// for (int l : lts->get_irrelevant_labels()){
		//     changes |= set_lqrel(l, l2, lts_i, min_value);
		// }
		old_value = min_value;
	    }
        }

        //Does l2 simulates irrelevant_labels in lts?
	old_value = get_simulates_irrelevant(lg2, lts_i);
	if(old_value != std::numeric_limits<int>::lowest()) {
	    T min_value = std::numeric_limits<int>::max();
	    for(int s = 0; s < lts->size(); s++){
		T max_value = std::numeric_limits<int>::lowest();
		for(const auto & tr : lts->get_transitions_label_group(lg2)) {
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
		changes |= set_simulates_irrelevant(lg2, lts_i, lts, min_value);
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
    int count = 0;
    for(int l2 : lts->get_relevant_labels()) {
        for(int l1 : lts->get_relevant_labels()) {
	    if(l2 == l1) {
		continue;
	    }
	    if (may_dominate(l2, l1, lts_id)) {
		cout << g_operators[l1].get_name() << " <= " << 
		    g_operators[l2].get_name() << " with " 
		     << q_dominates(l2, l1, lts_id) << endl;
		count ++;
	    }
	}
	if (may_dominated_by_noop (l2, lts_id))  {
	    cout <<  g_operators[l2].get_name() << " dominates noop: " << q_dominated_by_noop(l2, lts_id) << endl;
	    count ++;
	}
	if (may_dominate_noop_in (l2, lts_id))  {
	    cout <<  g_operators[l2].get_name() << " dominated by noop: " << q_dominates_noop(l2, lts_id) << endl;
	    count ++;
	}
    }
    cout << "Numeric label relation " << lts_id << " total count: " << count << endl; 

}


template class NumericLabelRelation<int>; 
template class NumericLabelRelation<IntEpsilon>; 
