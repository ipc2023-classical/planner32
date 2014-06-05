#include "labelled_transition_system.h"

#include "abstraction.h"

using namespace std;

LabelledTransitionSystem::LabelledTransitionSystem (Abstraction * abs) : 
  num_states(abs->size()),
  goal_states(abs->get_goal_states()){

  const vector<bool> is_rel_label (abs->get_relevant_labels());

  for(int i = 0; i < num_states; i++){
    name_states.push_back(abs->description(i));
  }
  transitions_src.resize(abs->size());
  transitions_label.resize(abs->get_num_labels());
 
  for (int label_no = 0; label_no < abs->get_num_labels(); label_no++) {
    if(is_rel_label[label_no]){
      relevant_labels.push_back(label_no);
      const vector<AbstractTransition> &abs_tr = abs->get_transitions_for_label(label_no);
      for (int j = 0; j < abs_tr.size(); j++) {
	LTSTransition t (abs_tr[j].src, abs_tr[j].target, label_no);
	transitions.push_back(t);
	transitions_src[t.src].push_back(t);
	transitions_label[t.label].push_back(t);
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


 
