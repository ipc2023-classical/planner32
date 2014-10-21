#include "label_relation.h"

#include "labels.h"
#include "simulation_relation.h"
#include "labelled_transition_system.h"

using namespace std;

LabelRelation::LabelRelation(Labels * _labels) : labels (_labels){}

void LabelRelation::dump() const {

  for (int l = 0; l < dominates_in.size(); ++l){
    if (labels->is_label_reduced(l)) cout << "reduced";
    if (l < 10){
      cout << "l" << l << ": ";  dump(l);
    }else{
      cout << "l" << l << ":";  dump(l);
    }

  }
}

void LabelRelation::dump(int label) const {
  if(dominated_by_noop_in[label] >= 0 && dominated_by_noop_in[label] <= 9){
    cout << " Dominated by noop: " << dominated_by_noop_in[label] << ", labels: ";
  }else{
    cout << " Dominated by noop:" << dominated_by_noop_in[label] << ", labels: ";
  } 

  for (int l2 = 0; l2 < dominates_in[label].size(); ++l2){
    if(dominates_in[l2][label] >= 0 && dominates_in[l2][label] <= 9) cout <<" ";
    cout << dominates_in[l2][label] << " ";
  }
  cout  << endl;
}

void LabelRelation::prune_operators(){
    //cout << "We have " << labels->get_size() << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
    for (int l = 0; l < dominates_in.size(); ++l){
	//labels->get_label_by_index(l)->dump();	
	if (dominated_by_noop_in[l]== DOMINATES_IN_ALL){
	    cout << g_operators[l].get_name() << " is dominated by noop " << endl;
	}

	for (int l2 = 0; l2 < dominates_in.size(); ++l2){
	    if (l2 != l && dominates_in[l2][l] == DOMINATES_IN_ALL){
		cout << g_operators[l].get_name() << " is dominated by " << g_operators[l2].get_name() << endl;
	    }
	}
    }
}

void LabelRelation::get_labels_dominated_in_all(std::vector<int> & labels_dominated_in_all){
    //cout << "We have " << labels->get_size() << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
    for (int l = 0; l < dominates_in.size(); ++l){
	//labels->get_label_by_index(l)->dump();	
	if (dominated_by_noop_in[l] == DOMINATES_IN_ALL){
	    labels_dominated_in_all.push_back(l);
	    continue;
	}

	for (int l2 = 0; l2 < dominates_in.size(); ++l2){
	    if (l2 != l && dominates_in[l2][l] == DOMINATES_IN_ALL && 
		!dominates_in[l][l2] == DOMINATES_IN_ALL){
		labels_dominated_in_all.push_back(l);
		break;
	    }
	}
    }
}


void LabelRelation::reset(){
    for (int i = 0; i < labels->get_size(); i++){
	for(int j = 0; j < simulates_irrelevant[i].size(); j++){
	    simulates_irrelevant[i][j] = true;
	    simulated_by_irrelevant[i][j] = true;
	}
    }
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
	dominated_by_noop_in[l1] = DOMINATES_IN_ALL;
	for (int l2 = 0; l2 < dominates_in[l1].size(); ++l2){
	    if(labels->get_label_by_index(l1)->get_cost() <=
	       labels->get_label_by_index(l2)->get_cost()){
		dominates_in[l1][l2] = DOMINATES_IN_ALL;
	    }
	}
    }
}

void LabelRelation::init(const std::vector<LabelledTransitionSystem *> & lts,
	  const std::vector<SimulationRelation*> & sim){
    simulates_irrelevant.resize(labels->get_size());
    simulated_by_irrelevant.resize(labels->get_size());
    for(int i = 0; i < labels->get_size(); i++){
	simulates_irrelevant[i].resize(lts.size(), true);
	simulated_by_irrelevant[i].resize(lts.size(), true);
    }
    
    dominates_in.resize(labels->get_size());
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
	dominated_by_noop_in.resize(labels->get_size(), DOMINATES_IN_ALL);
	dominates_in[l1].resize(labels->get_size(), DOMINATES_IN_ALL);
	for (int l2 = 0; l2 < dominates_in[l1].size(); ++l2){
	    if(labels->get_label_by_index(l1)->get_cost() > 
	       labels->get_label_by_index(l2)->get_cost()){
		dominates_in[l1][l2] = DOMINATES_IN_NONE;
	    }
	}
    }
    cout << "Update label dominance: " << labels->get_size()  
	 << " labels " << lts.size() << " systems." << endl;
    for (int i = 0; i < lts.size(); ++i){
	update(i, lts[i], sim[i]);
    }
}

bool LabelRelation::update(const std::vector<LabelledTransitionSystem*> & lts,
	    const std::vector<SimulationRelation*> & sim){
    bool changes = false;
    for (int i = 0; i < lts.size(); ++i){
	changes |= update(i, lts[i], sim[i]);
    }
    return changes;
}

/* TODO: REALLY INEFFICIENT VERSION THAT SHOULD BE IMPROVED */
bool LabelRelation::update(int i, const LabelledTransitionSystem * lts, 
			   const SimulationRelation * sim){
    bool changes = false;
    for(int l2 : lts->get_relevant_labels()) {
	for(int l1 : lts->get_relevant_labels()){ 
	    if(l1 != l2 && simulates(l1, l2, i)){
		//std::cout << "Check " << l1 << " " << l2 << std::endl;
		//std::cout << "Num transitions: " << lts->get_transitions_label(l1).size() 
		//		    << " " << lts->get_transitions_label(l2).size() << std::endl;
		//Check if it really simulates
		//For each transition s--l2-->t, and every label l1 that dominates
		//l2, exist s--l1-->t', t <= t'?
		for(auto & tr : lts->get_transitions_label(l2)){
		    bool found = false;
		    //TODO: for(auto tr2 : lts->get_transitions_for_label_src(l1, tr.src)){
		    for(auto & tr2 : lts->get_transitions_label(l1)){
			if(tr2.src == tr.src &&
			   sim->simulates(tr2.target, tr.target)){
			    found = true;
			    break; //Stop checking this tr
			}
		    }
		    if(!found){
			//std::cout << "Not sim " << l1 << " " << l2 << " " << i << std::endl;
			set_not_simulates(l1, l2, i);
			changes = true;
			break; //Stop checking trs of l1
		    }
		}
	    }
	}

	//Is l2 simulated by irrelevant_labels in lts?
	for(auto tr : lts->get_transitions_label(l2)){
	    if (simulated_by_irrelevant[l2][i] && 
		!sim->simulates(tr.src, tr.target)) {
		changes |= set_not_simulated_by_irrelevant(l2, i);
		for (int l : lts->get_irrelevant_labels()){
		    if(simulates(l, l2, i)){
			changes = true;
			set_not_simulates(l, l2, i);
		    }
		}
	    }
	}

	//Does l2 simulates irrelevant_labels in lts?
	if(simulates_irrelevant[l2][i]){
	    for(int s = 0; s < lts->size(); s++){
		bool found = false;
		for(auto tr : lts->get_transitions_label(l2)){
		    if(tr.src == s && sim->simulates(tr.target, tr.src)) {
			found = true;
			break;
		    }
		}
		if(!found) {
		    simulates_irrelevant[l2][i] = false;
		    for (int l : lts->get_irrelevant_labels()){
			if(simulates(l2, l, i)){
			    set_not_simulates(l2, l, i);
			    changes = true;
			}
		    }
		}
	    }
	}
    }

    return changes;
}
