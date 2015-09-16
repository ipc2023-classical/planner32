#include "label_relation_noop.h"

#include "labels.h"
#include "simulation_relation.h"
#include "dominance_relation.h"
#include "labelled_transition_system.h"
#include "lts_complex.h"
#include "../equivalence_relation.h"
#include "../globals.h"
#include "../utilities.h"

using namespace std;

LabelRelationNoop::LabelRelationNoop(Labels * _labels) : labels (_labels), 
        num_labels(_labels->get_size()){
}

void LabelRelationNoop::dump_equivalent() const {

}


void LabelRelationNoop::dump_dominance() const {
    
}

void LabelRelationNoop::dump() const {
}

void LabelRelationNoop::dump(int label) const {
    cout << "Dump l: " << label << "; ";
    if(dominated_by_noop_in[label] >= 0 && dominated_by_noop_in[label] <= 9){
        cout << " Dominated by noop: " << dominated_by_noop_in[label] << ", labels: ";
    }else{
        cout << " Dominated by noop:" << dominated_by_noop_in[label] << ", labels: ";
    }
}

void LabelRelationNoop::prune_operators(){
    //cout << "We have " << num_labels << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
    for (int l = 0; l < dominated_by_noop_in.size(); ++l){
        //labels->get_label_by_index(l)->dump();
        if (dominated_by_noop_in[l]== DOMINATES_IN_ALL){
            cout << g_operators[l].get_name() << " is dominated by noop " << endl;
        }
    }
}


// void check_transitivity() {
//     //Check transitivity
//     for (int l = 0; l < dominates_in.size(); ++l){
// 	for (int l2 = 0; l2 < dominates_in.size(); ++l2){
// 	    if(dominates_in[l2][l] == DOMINATES_IN_ALL) {
// 		for (int l3 = 0; l3 < dominates_in.size(); ++l3){
// 		    if (dominates_in[l][l3] == DOMINATES_IN_ALL && 
// 			dominates_in[l2][l3] != DOMINATES_IN_ALL ){
// 			cerr << "Assertion error: label dominance in all is not transitive" << endl;
// 			exit(0);
// 		    }
// 		    if (dominates_in[l][l3] != DOMINATES_IN_NONE && 
// 			dominates_in[l2][l3] != DOMINATES_IN_ALL && 
// 			dominates_in[l2][l3] != dominates_in[l][l3]){
// 			cerr << "Assertion error: label dominance in all is not transitive 2" << endl;
// 			exit(0);
// 		    }
// 		}		
// 	    }
// 	}
	
//     }

//     for (int l = 0; l < dominates_in.size(); ++l){
// 	if(dominated_by_noop_in[l] == DOMINATES_IN_ALL) {
// 	    for (int l3 = 0; l3 < dominates_in.size(); ++l3){
// 		if (dominates_in[l][l3] == DOMINATES_IN_ALL && 
// 		    dominated_by_noop_in[l3] != DOMINATES_IN_ALL ){
// 		    cerr << "Assertion error: label dominance in all is not transitive for noops" << endl;
// 		    exit(0);
// 		}
// 		if (dominates_in[l][l3] != DOMINATES_IN_NONE && 
// 		    dominated_by_noop_in[l3] != DOMINATES_IN_ALL && 
// 		    dominated_by_noop_in[l3] != dominates_in[l][l3]){
// 		    cerr << "Assertion error: label dominance in all is not transitive for noops" << endl;
// 		    exit(0);
// 		}
// 	    }		
// 	}
//     }
    
//     cerr << "Transitivity checked" << endl;
// }

std::vector<int> LabelRelationNoop::get_labels_dominated_in_all() const{
    std::vector<int> labels_dominated_in_all;
    for (int l = 0; l < dominated_by_noop_in.size(); ++l){
        if (dominated_by_noop_in[l] == DOMINATES_IN_ALL){
            labels_dominated_in_all.push_back(l);
            continue;
        }
    }

    return labels_dominated_in_all;
}

void LabelRelationNoop::reset(){
    cout << "Error: reser of label relation has been disabled" << endl;
    exit(0);
    // for (int i = 0; i < num_labels; i++){
    // 	for(int j = 0; j < simulates_irrelevant[i].size(); j++){
    // 	    simulates_irrelevant[i][j] = true;
    // 	    simulated_by_irrelevant[i][j] = true;
    // 	}
    // }
    // for (int l1 = 0; l1 < dominates_in.size(); ++l1){
    // 	dominated_by_noop_in[l1] = DOMINATES_IN_ALL;
    // 	for (int l2 = 0; l2 < dominates_in[l1].size(); ++l2){
    // 	    if(labels->get_label_by_index(l1)->get_cost() <=
    // 	       labels->get_label_by_index(l2)->get_cost()){
    // 		dominates_in[l1][l2] = DOMINATES_IN_ALL;
    // 	    }
    // 	}
    // }
}

void LabelRelationNoop::init(const std::vector<LabelledTransitionSystem *> & lts,
        const DominanceRelation & sim,
        const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();

    //TODO: Only do this in the incremental step (reset) 
    //TODO: Is there a better way to reinitialize these? 
    std::vector<std::vector<bool> >().swap(simulates_irrelevant);
    std::vector<std::vector<bool> >().swap(simulated_by_irrelevant);
    std::vector<int> ().swap(dominated_by_noop_in);

    simulates_irrelevant.resize(num_labels);
    simulated_by_irrelevant.resize(num_labels);
    for(int i = 0; i < num_labels; i++){
        simulates_irrelevant[i].resize(lts.size(), true);
        simulated_by_irrelevant[i].resize(lts.size(), true);
    }

    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_ALL);
        cout << "Update label dominance: " << num_labels
            << " labels " << lts.size() << " systems." << endl;

    for (int i = 0; i < lts.size(); ++i){
        update(i, lts[i], sim[i]);
    }
}

void LabelRelationNoop::init(const std::vector<LTSComplex *> & lts,
        const DominanceRelation & sim,
        const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();

    simulates_irrelevant.resize(num_labels);
    simulated_by_irrelevant.resize(num_labels);
    for(int i = 0; i < num_labels; i++){
        simulates_irrelevant[i].resize(lts.size(), true);
        simulated_by_irrelevant[i].resize(lts.size(), true);
    }

    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    
    cout << "Update label dominance: " << num_labels
            << " labels " << lts.size() << " systems." << endl;
    for (int i = 0; i < lts.size(); ++i){
        update(i, lts[i], sim[i]);
    }
}


bool LabelRelationNoop::update(const std::vector<LabelledTransitionSystem*> & lts,
        const DominanceRelation & sim){
    bool changes = false;
    
    for (int i = 0; i < lts.size(); ++i){
        changes |= update(i, lts[i], sim[i]);
    }
    
    return changes;
}

bool LabelRelationNoop::update(const std::vector<LTSComplex*> & lts,
        const DominanceRelation & sim){
    bool changes = false;
    for (int i = 0; i < lts.size(); ++i){
        changes |= update(i, lts[i], sim[i]);
    }
    return changes;
}

/* TODO: This version is inefficient. It could be improved by
 * iterating only the right transitions (see TODO inside the loop)
 */
bool LabelRelationNoop::update(int i, const LabelledTransitionSystem * lts, 
        const SimulationRelation & sim){
    bool changes = false;
    for(int l2 : lts->get_relevant_labels()) {
        //Is l2 simulated by irrelevant_labels in lts?
        for(auto tr : lts->get_transitions_label(l2)){
            if (simulated_by_irrelevant[l2][i] &&
                    !sim.simulates(tr.src, tr.target)) {
                changes |= set_not_simulated_by_irrelevant(l2, i);
            }
        }

        //Does l2 simulates irrelevant_labels in lts?
        if(simulates_irrelevant[l2][i]){
            for(int s = 0; s < lts->size(); s++){
                bool found = false;
                for(const auto & tr : lts->get_transitions_label(l2)){
                    if(tr.src == s && sim.simulates(tr.target, tr.src)) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    simulates_irrelevant[l2][i] = false;
                }
            }
        }
    }

    return changes;
}

bool LabelRelationNoop::update(int i, const LTSComplex * lts, 
        const SimulationRelation & sim){
    bool changes = false;
    for(int l2 : lts->get_relevant_labels()) {
        
        //Is l2 simulated by irrelevant_labels in lts?
        if (simulated_by_irrelevant[l2][i] &&
                lts->applyPost(l2,
                        //Exists s-l2-> y s.t. s is not simulated by t
                        [&](const LTSTransition & tr){
            return !sim.simulates(tr.src, tr.target);
        })){
            changes |= set_not_simulated_by_irrelevant(l2, i);
	}

        //Does l2 simulates irrelevant_labels in lts?
        if(simulates_irrelevant[l2][i]){
            for(int s = 0; s < lts->size(); s++){
                if(!lts->applyPost(l2,
                        [&](const LTSTransition & tr){
                    return tr.src == s && sim.simulates(tr.target, tr.src);
                })){
                    simulates_irrelevant[l2][i] = false;
                }
            }
        }
    }

    return changes;
}

EquivalenceRelation * 
LabelRelationNoop::get_equivalent_labels_relation(const LabelMap & labelMap, 
					      set<int> &  /*dangerous_LTSs*/) const {
    list<Block> rel;
    for (int l1 = 0; l1 < num_labels; l1++){
	Block eq;
	eq.insert(labelMap.get_old_id(l1));
	rel.push_back(eq);
    }

    return new EquivalenceRelation(rel.size(), rel);
}


/* Returns true if we succeeded in propagating the effects of pruning a transition in lts i. */
bool LabelRelationNoop::propagate_transition_pruning(int lts_id, 
						 const vector<LabelledTransitionSystem *> & ltss, 
						 const DominanceRelation & simulations,
						 int src, int l1, int target) const {
    LabelledTransitionSystem * lts = ltss[lts_id];
    const SimulationRelation & sim = simulations[lts_id]; 

    if(simulates_irrelevant[l1][lts_id]) {
	//For each transition from src, check if anything has changed
	bool found = lts->applyPostSrc(src, [&](const LTSTransition & tr){
		if (l1 == tr.label) { //Same label
		    if(tr.target == target) return false;
		    if(sim.simulates(tr.target, tr.src)){
			//There is another transition with the same label which simulates noop
			return true;
		    }
		}
		return false;
	});
	if(!found) {
	    return false;
	}
    }

    //TODO: This should be moved somewhere else, but it is convinient to place it here. 
    lts->kill_transition (src, l1, target);
    return true;
}
