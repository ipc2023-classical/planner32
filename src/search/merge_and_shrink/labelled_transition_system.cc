#include "labelled_transition_system.h"

#include "abstraction.h"
#include "labels.h"

using namespace std;

LabelledTransitionSystem::LabelledTransitionSystem (Abstraction * _abs, const LabelMap & labelMap) : 
    abs(_abs),  num_states(_abs->size()),goal_states(_abs->get_goal_states()){
    
    int num_labels = labelMap.get_num_labels();
    const vector<bool> &was_rel_label = abs->get_relevant_labels();

    for(int i = 0; i < num_states; i++){
	name_states.push_back(abs->description(i));
    }

    transitions_src.resize(abs->size());
    transitions_label.resize(num_labels);
 
    for (int label_no = 0; label_no < num_labels; label_no++) {
	int old_label = labelMap.get_old_id(label_no);
	if(was_rel_label[old_label]){
	    const vector<AbstractTransition> &abs_tr = abs->get_transitions_for_label(old_label);

	    if (!abs_tr.empty()) {
		relevant_labels.push_back(label_no);
		for (int j = 0; j < abs_tr.size(); j++) {
		    LTSTransition t (abs_tr[j].src, abs_tr[j].target, label_no);
		    transitions.push_back(t);
		    transitions_src[t.src].push_back(t);
		    transitions_label[t.label].push_back(t);
		}
	    } else {
		//dead_label
	    }
	}else{
	    irrelevant_labels.push_back(label_no);
	    /*for(int i = 0; i < num_states; i++){
	      LTSTransition t (i, i, label_no);
	      transitions.push_back(t);
	      transitions_src[t.src].push_back(t);
	      transitions_label[t.label].push_back(t);
	      }*/
	}
    }
}

void LabelledTransitionSystem::kill_transition(int src, int label, int target) {
    LTSTransition t(src, target, label);
    kill_from_vector(t, transitions);
    kill_from_vector(t, transitions_src[src]);
    kill_from_vector(t, transitions_label[label]);
}


void LabelledTransitionSystem::kill_label(int l) {
    cout << "KILL " << l << endl;
    std::vector <LTSTransition>().swap(transitions_label[l]);

    transitions.erase(std::remove_if(begin(transitions),
				     end(transitions),
				     [&](LTSTransition & t){
					 return t.label == l;
				     }), end(transitions));
    
    for (auto & trs : transitions_src){
    	trs.erase(std::remove_if(begin(trs),
    				 end(trs),
    				 [&](LTSTransition & t){
    				     return t.label == l;
    				 }), end(trs));
    }


    relevant_labels.erase(remove(begin(relevant_labels), end(relevant_labels), l), end(relevant_labels));
    irrelevant_labels.erase(remove(begin(irrelevant_labels), end(irrelevant_labels), l), end(irrelevant_labels));   
}


vector<int> LabelledTransitionSystem::compute_label_equivalence_groups() {
    //Ensure that transitions from label are sorted
    for(int l : relevant_labels) {
        std::sort(transitions_label[l].begin(),
                  transitions_label[l].end(),
                  [](const LTSTransition & tr, const LTSTransition & tr2) -> bool {
                      return tr.src < tr2.src || (tr.src == tr2.src && tr.target < tr2.target);
                  });
    }
    
    //All labels go to group 0 by default (irrelevant labels group)
    vector<int> label_to_equivalence_group (transitions_label.size(), 0);
    std::vector<int> representative_label_from_group(1);
    for(int l : relevant_labels) {
        int candidate_group = 1;
        for (; candidate_group < representative_label_from_group.size(); ++candidate_group) {
            // Pick a label from the candidate group to check if both labels have the same transitions
            int l2 = representative_label_from_group[candidate_group];

            if (transitions_label[l].size() == transitions_label[l2].size() &&
                std::equal(transitions_label[l].begin(),
                           transitions_label[l].end(),
                           transitions_label[l2].begin(),
                           [](const LTSTransition & tr, const LTSTransition & tr2) {
                               return tr.src == tr2.src && tr.target == tr2.target;
                           })) {
                //We have found an equivalent label group
                break;
            }          
        }

        label_to_equivalence_group[l] = candidate_group;
        if (candidate_group == representative_label_from_group.size()) {
            representative_label_from_group.push_back(l);
        }
    }
    
    return label_to_equivalence_group;
}




