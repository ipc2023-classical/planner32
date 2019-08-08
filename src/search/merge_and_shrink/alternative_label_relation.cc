#include "alternative_label_relation.h"

#include "labels.h"
#include "labelled_transition_system.h"
#include "simulation_relation.h"
#include "dominance_relation.h"

#include "../globals.h"
#include "../utilities.h"
#include "../debug.h"

using namespace std;


AlternativeLabelRelation::AlternativeLabelRelation(Labels * _labels) :
    labels (_labels), num_labels(_labels->get_size()){

}

 
void AlternativeLabelRelation::init(const std::vector<LabelledTransitionSystem *> & lts,
				const DominanceRelation & sim,
				const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();
    num_ltss = lts.size();

    cout << "Init label dominance: " << num_labels
	      << " labels " << lts.size() << " systems." << endl;

    std::vector<int> ().swap(cost_of_label);
    std::vector<std::vector<int> >().swap(position_of_label);
    std::vector<std::vector<std::vector<bool> > >().swap(lrel);
    std::vector<std::vector<bool> >().swap(simulates_irrelevant);
    std::vector<std::vector<bool> >().swap(simulated_by_irrelevant);

    irrelevant_labels_lts.resize(lts.size());
    position_of_label.resize(lts.size());
    simulates_irrelevant.resize(lts.size());
    simulated_by_irrelevant.resize(lts.size());
    lrel.resize(lts.size());

    cost_of_label.resize(num_labels);
    for(int l = 0; l < num_labels; l++) {
	cost_of_label[l] = labelMap.get_cost(l);
    }

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

	simulates_irrelevant[i].resize(num_relevant_labels, true);
        simulated_by_irrelevant[i].resize(num_relevant_labels, true);
	lrel[i].resize(num_relevant_labels);

	for(int j = 0; j < num_relevant_labels; j++) {
	    lrel[i][j].resize(num_relevant_labels, true);
	}
    }

    cout << "Dominating." << endl;
    std::vector<std::vector<int> > ().swap(dominates_in);
    std::vector<int> ().swap(dominated_by_noop_in);
    std::vector<int> ().swap(dominates_noop_in);
    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    dominates_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    
    if(num_labels < 0) { // If we have more than 5000 labels, there is not enough space. 
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

    // for (int lts_id = 0; lts_id < lts.size(); ++lts_id) {
    // 	cout << "Number of tau labels: " << tau_labels[lts_id].size() << endl;
    // }
}


 
bool AlternativeLabelRelation::update(const std::vector<LabelledTransitionSystem*> & lts,
				     const DominanceRelation & sim) {
    bool changes = false;
    
    for (int i = 0; i < lts.size(); ++i){
        changes |= update(i, lts[i], sim[i]);
    }
		
    return changes;
}

 
bool AlternativeLabelRelation::update(int lts_i, const LabelledTransitionSystem * lts,  
                                      const SimulationRelation & sim){

    bool changes = false;
    //cout << "UPDATE " << lts_i << " " << lts->get_relevant_labels().size() << endl;
    for(int l2 : lts->get_relevant_labels()) {
        for(int l1 : lts->get_relevant_labels()){
            if(l1 != l2 && simulates(l1, l2, lts_i)){
                //std::cout << "Check " << l1 << " " << l2 << std::endl;
                //std::cout << "Num transitions: " << lts->get_transitions_label(l1).size()
                //		    << " " << lts->get_transitions_label(l2).size() << std::endl;
                //Check if it really simulates
                //For each transition s--l2-->t, and every label l1 that dominates
                //l2, exist s--l1-->t', t <= t'?
                for(const auto & tr : lts->get_transitions_label(l2)){
                    bool found = false;
                    //TODO: for(auto tr2 : lts->get_transitions_for_label_src(l1, tr.src)){
                    for(const auto & tr2 : lts->get_transitions_label(l1)){
                        if(tr2.src == tr.src &&
                                sim.simulates(tr2.target, tr.target)){
                            found = true;
                            break; //Stop checking this tr
                        }
                    }
                    if(!found){
                        //std::cout << "Not sim " << l1 << " " << l2 << " " << i << std::endl;
                        set_not_simulates(l1, l2, lts_i);
                        changes = true;
                        break; //Stop checking trs of l1
                    }
                }
            }
        }

        //Is l2 simulated by irrelevant_labels in lts?
        if (get_simulated_by_irrelevant(l2, lts_i)) {
            for(auto tr : lts->get_transitions_label(l2)){
                if(!sim.simulates(tr.src, tr.target)) {
                    changes |= set_not_simulated_by_irrelevant(l2, lts_i);
                }
            }
        }

        //Does l2 simulates irrelevant_labels in lts?
        if(get_simulates_irrelevant(l2, lts_i)){
            for(int s = 0; s < lts->size(); s++){
                bool found = false;
                for(const auto & tr : lts->get_transitions_label(l2)){
                    if(tr.src == s && sim.simulates(tr.target, tr.src)) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    changes |= set_not_simulates_irrelevant(l2, lts_i);
                }
            }
        }
    }

    return changes;
}


void AlternativeLabelRelation::dump(const LabelledTransitionSystem * , int ) const {

}

std::vector<int> AlternativeLabelRelation::get_labels_dominated_in_all() const{
    std::vector<int> labels_dominated_in_all;

    //TODO: Implement case where dominates_in is empty
    //cout << "We have " << num_labels << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
    for (int l = 0; l < dominates_in.size(); ++l){
	//cout << "Check: " << l << endl;
        //labels->get_label_by_index(l)->dump();
        if (dominated_by_noop_in[l] == DOMINATES_IN_ALL){
            labels_dominated_in_all.push_back(l);
            continue;
        }

        // PIET-edit: Here we do not remove either label if they dominate each other in all LTSs.
        // for (int l2 = 0; l2 < dominates_in.size(); ++l2){
        //     if (l2 != l && dominates_in[l2][l] == DOMINATES_IN_ALL &&
        //             !dominates_in[l][l2] == DOMINATES_IN_ALL){
        //         labels_dominated_in_all.push_back(l);
        //         break;
        //     }
        // }

        // PIET-edit: Here we remove one of the two labels dominating each other in all LTSs.
        for (int l2 = 0; l2 < dominates_in.size(); ++l2){
	    //cout << " with : " << l2 << " " << dominates_in[l2][l] << "  " << dominates_in[l][l2];
	    // if ((l == 118 && l2 == 279) || (l2 == 118 && l == 279)) {
	    // 	cout << "HERE: " << dominates_in[l2][l] << " " << dominates_in[l][l2] << endl;
	    // 	cout << (l2 < l) << " " << (dominates_in[l2][l] == DOMINATES_IN_ALL) << " " << (dominates_in[l][l2] != DOMINATES_IN_ALL) << endl;
	    // 	cout <<  (l2 > l) << " " << (dominates_in[l2][l] == DOMINATES_IN_ALL) << endl;
	    // }
	
	    // if ( dominates_in[l2][l] == DOMINATES_IN_ALL && 
	    // 	 (dominates_in[l][l2] != DOMINATES_IN_ALL || l2 < l)) {
            if ((l2 < l && dominates_in[l2][l] == DOMINATES_IN_ALL && 
	    	 dominates_in[l][l2] != DOMINATES_IN_ALL)
	    	|| (l2 > l && dominates_in[l2][l] == DOMINATES_IN_ALL)) {
		//cout << " yes" << endl;			       
                labels_dominated_in_all.push_back(l);
                break;
            }
	    //cout << " no" << endl;
        }


        // PIET-edit: If we take proper care, this both dominating each other cannot ever happen.
//        for (int l2 = 0; l2 < dominates_in.size(); ++l2) {
//            if (l != l2 && dominates_in[l][l2] == DOMINATES_IN_ALL && dominates_in[l2][l] == DOMINATES_IN_ALL) {
//                cerr << "Error: two labels dominating each other in all abstractions. This CANNOT happen!" << endl;
//                cerr << l << "; " << l2 << endl;
//                exit(1);
//            }
//            if (dominates_in[l2][l] == DOMINATES_IN_ALL) {
//                labels_dominated_in_all.push_back(l);
//                break;
//            }
//        }
    }

    return labels_dominated_in_all;
}

