#include "label_relation.h"

#include "labels.h"
#include "simulation_relation.h"
#include "labelled_transition_system.h"
#include "lts_efficient.h"
#include "../equivalence_relation.h"
#include "../globals.h"

using namespace std;

LabelRelation::LabelRelation(Labels * _labels) : labels (_labels), 
        num_labels(_labels->get_size()){

}

void LabelRelation::dump_equivalent() const {
    vector<bool> redundant(dominates_in.size(), false);
    int num_redundant = 0;
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
        for (int l2 = l1+1; l2 < dominates_in.size(); ++l2){
            if(!redundant[l2] && dominates_in[l1][l2] != DOMINATES_IN_NONE &&
                    dominates_in[l2][l1] == dominates_in[l1][l2]){
                redundant[l2] =true;
                num_redundant ++;
                cout << l1 << " equivalent to " << l2 << " in " << dominates_in[l1][l2] << endl;
            }
        }
    }
    cout << "Redundant labels: " << num_redundant << endl;
}


void LabelRelation::dump_dominance() const {
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
        for (int l2 = 0; l2 < dominates_in.size(); ++l2){
            if(dominates_in[l1][l2] != DOMINATES_IN_NONE &&
                    dominates_in[l2][l1] != dominates_in[l1][l2]){
                cout << l1 << " dominates " << l2 << " in " << dominates_in[l1][l2] << endl;
                cout << g_operators[l1].get_name() << " dominates " << g_operators[l2].get_name() << endl;
            }
        }
    }
}

void LabelRelation::dump() const {
    for (int l = 0; l < dominates_in.size(); ++l){
        //if (labels->is_label_reduced(l)) cout << "reduced";
        if (l < 10){
            cout << "l" << l << ": ";  dump(l);
        }else{
            cout << "l" << l << ":";  dump(l);
        }

    }
}

void LabelRelation::dump(int label) const {
    cout << "Dump l: " << label << "; ";
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
    //cout << "We have " << num_labels << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
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
    //cout << "We have " << num_labels << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
    for (int l = 0; l < dominates_in.size(); ++l){
        //labels->get_label_by_index(l)->dump();
        if (dominated_by_noop_in[l] == DOMINATES_IN_ALL){
            labels_dominated_in_all.push_back(l);
            continue;
        }

        // PIET-edit: Here we do not remove either label if they dominate each other in all LTSs.
        for (int l2 = 0; l2 < dominates_in.size(); ++l2){
            if (l2 != l && dominates_in[l2][l] == DOMINATES_IN_ALL &&
                    !dominates_in[l][l2] == DOMINATES_IN_ALL){
                labels_dominated_in_all.push_back(l);
                break;
            }
        }
        // PIET-edit: Here we remove one of the two labels dominating each other in all LTSs.
//        for (int l2 = 0; l2 < dominates_in.size(); ++l2){
//            if ((l2 < l && dominates_in[l2][l] == DOMINATES_IN_ALL && !dominates_in[l][l2] == DOMINATES_IN_ALL)
//                    || (l2 > l && dominates_in[l2][l] == DOMINATES_IN_ALL)) {
//                labels_dominated_in_all.push_back(l);
//                break;
//            }
//        }
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
}

void LabelRelation::reset(){
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

void LabelRelation::init(const std::vector<LabelledTransitionSystem *> & lts,
        const std::vector<SimulationRelation*> & sim,
        const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();
    simulates_irrelevant.resize(num_labels);
    simulated_by_irrelevant.resize(num_labels);
    for(int i = 0; i < num_labels; i++){
        simulates_irrelevant[i].resize(lts.size(), true);
        simulated_by_irrelevant[i].resize(lts.size(), true);
    }

    dominates_in.resize(num_labels);
    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
        int old_l1 = labelMap.get_old_id(l1);
        dominates_in[l1].resize(num_labels, DOMINATES_IN_ALL);
        for (int l2 = 0; l2 < dominates_in[l1].size(); ++l2){
            int old_l2 = labelMap.get_old_id(l2);
            if(labels->get_label_by_index(old_l1)->get_cost() >
            labels->get_label_by_index(old_l2)->get_cost()){
                dominates_in[l1][l2] = DOMINATES_IN_NONE;
            }
        }
    }
    cout << "Update label dominance: " << num_labels
            << " labels " << lts.size() << " systems." << endl;


    for (int i = 0; i < lts.size(); ++i){
        update(i, lts[i], sim[i]);
    }

}

void LabelRelation::init(const std::vector<LTSEfficient *> & lts,
        const std::vector<SimulationRelation*> & sim,
        const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();

    simulates_irrelevant.resize(num_labels);
    simulated_by_irrelevant.resize(num_labels);
    for(int i = 0; i < num_labels; i++){
        simulates_irrelevant[i].resize(lts.size(), true);
        simulated_by_irrelevant[i].resize(lts.size(), true);
    }

    dominates_in.resize(num_labels);
    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_ALL);
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
        int old_l1 = labelMap.get_old_id(l1);
        dominates_in[l1].resize(num_labels, DOMINATES_IN_ALL);
        for (int l2 = 0; l2 < dominates_in[l1].size(); ++l2){
            int old_l2 = labelMap.get_old_id(l2);

            if(labels->get_label_by_index(old_l1)->get_cost() >
            labels->get_label_by_index(old_l2)->get_cost()){
                dominates_in[l1][l2] = DOMINATES_IN_NONE;
            }
        }
    }

    cout << "Update label dominance: " << num_labels
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

bool LabelRelation::update(const std::vector<LTSEfficient*> & lts,
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
                for(const auto & tr : lts->get_transitions_label(l2)){
                    bool found = false;
                    //TODO: for(auto tr2 : lts->get_transitions_for_label_src(l1, tr.src)){
                    for(const auto & tr2 : lts->get_transitions_label(l1)){
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
                for(const auto & tr : lts->get_transitions_label(l2)){
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

bool LabelRelation::update(int i, const LTSEfficient * lts, 
        const SimulationRelation * sim){
    bool changes = false;
    for(int l2 : lts->get_relevant_labels()) {
        for(int l1 : lts->get_relevant_labels()){
            if(l1 != l2 && simulates(l1, l2, i)){
                //std::cout << "Check " << l1 << " " << l2 << std::endl;
                //std::cout << "Num transitions: " << lts->get_transitions_label(l1).size()
                //		    << " " << lts->get_transitions_label(l2).size() << std::endl;
                //Check if it really simulates
                //For each transition s--l2-->t, and every label l1 that simulates
                //l2, exist s--l1-->t', t <= t'?
                lts->applyPost(l2,
                        [&](const LTSTransition & tr){
                    if(!lts->applyPost(l1, tr.src,
                            [&](const LTSTransition & tr2){
                        return sim->simulates(tr2.target, tr.target);
                    })){
                        //std::cout << "Not sim " << l1 << " " << l2 << " " << i << std::endl;
                        set_not_simulates(l1, l2, i);
                        changes = true;
                        return true;
                    }
                    return false;
                });
            }
        }

        //Is l2 simulated by irrelevant_labels in lts?
        if (simulated_by_irrelevant[l2][i] &&
                lts->applyPost(l2,
                        //Exists s-l2-> y s.t. s is not simulated by t
                        [&](const LTSTransition & tr){
            return !sim->simulates(tr.src, tr.target);
        })){
            changes |= set_not_simulated_by_irrelevant(l2, i);
            for (int l : lts->get_irrelevant_labels()){
                if(simulates(l, l2, i)){
                    changes = true;
                    set_not_simulates(l, l2, i);
                }
            }
        }

        //Does l2 simulates irrelevant_labels in lts?
        if(simulates_irrelevant[l2][i]){
            for(int s = 0; s < lts->size(); s++){
                if(!lts->applyPost(l2,
                        [&](const LTSTransition & tr){
                    return tr.src == s && sim->simulates(tr.target, tr.src);
                })){
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





void LabelRelation::init_identity(int num_lts, const LabelMap & labelMap){
    num_labels = labelMap.get_num_labels();
    simulates_irrelevant.resize(num_labels);
    simulated_by_irrelevant.resize(num_labels);
    for(int i = 0; i < num_labels; i++){
        simulates_irrelevant[i].resize(num_lts, false);
        simulated_by_irrelevant[i].resize(num_lts, false);
    }

    dominates_in.resize(num_labels);
    dominated_by_noop_in.resize(num_labels, DOMINATES_IN_NONE);
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
        dominates_in[l1].resize(num_labels, DOMINATES_IN_NONE);
        dominates_in[l1][l1] = DOMINATES_IN_ALL;
    }
}

EquivalenceRelation * LabelRelation::get_equivalent_labels_relation(const LabelMap & labelMap) const {
    list<Block> rel;
    vector<bool> captured_labels(num_labels, false);
    for (int l1 = 0; l1 < num_labels; l1++){
        if (captured_labels[l1])
            continue;
        Block eq;
        eq.insert(labelMap.get_old_id(l1));
        captured_labels[l1] = true;
        for (int l2 = l1 + 1; l2 < num_labels; l2++) {
            if (captured_labels[l2])
                // was already marked as equivalent to some other label, so cannot be equivalent to l1
                continue;
            if (dominates_in[l1][l2] != DOMINATES_IN_NONE && dominates_in[l2][l1] != DOMINATES_IN_NONE &&
                    (dominates_in[l1][l2] == DOMINATES_IN_ALL ||
                            dominates_in[l2][l1] == DOMINATES_IN_ALL ||
                            dominates_in[l1][l2] == dominates_in[l2][l1])) {
                eq.insert(labelMap.get_old_id(l2));
                //cout << l2 << " eq " << l1 << endl;
                captured_labels[l2] = true;
            }
        }
        rel.push_back(eq);
    }
    return new EquivalenceRelation(rel.size(), rel);
}
