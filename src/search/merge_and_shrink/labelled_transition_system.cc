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

    label_group_of_label.resize(num_labels, LabelGroup(-1));
    label_groups.reserve(num_labels);

    transitions_src.resize(abs->size());
    transitions_label_group.reserve(num_labels);

    for (int label_no = 0; label_no < num_labels; label_no++) {
	int old_label = labelMap.get_old_id(label_no);
	if(was_rel_label[old_label]){
	    const vector<AbstractTransition> &abs_tr = abs->get_transitions_for_label(old_label);

	    if (!abs_tr.empty()) {
		vector<TSTransition> transitions_label;
		relevant_labels.push_back(label_no);
		for (int j = 0; j < abs_tr.size(); j++) {
		    transitions_label.push_back(TSTransition(abs_tr[j].src, abs_tr[j].target));
		}
        	std::sort(transitions_label.begin(), transitions_label.end());

		for(int g = 0; g < transitions_label_group.size(); ++g){
		    if(transitions_label_group[g] == transitions_label) {
			assert(g < label_groups.size());
			label_groups[g].push_back(label_no);
			label_group_of_label[label_no] = LabelGroup(g);
			break;
		    }
		}
		if (label_group_of_label[label_no].dead()) {
		    LabelGroup new_group (transitions_label_group.size());
		    for (const TSTransition & tr : transitions_label) {
			transitions.push_back(LTSTransition (tr.src, tr.target, new_group));
			transitions_src[tr.src].push_back(LTSTransition (tr.src, tr.target, new_group));
		    }
		    transitions_label_group.push_back(move(transitions_label));
		    label_groups.push_back(std::vector<int>());
		    label_groups[new_group.group].push_back(label_no);
		    label_group_of_label[label_no] = new_group;
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

#ifndef NDEBUG
    for (int label_no = 0; label_no < num_labels; label_no++) {
	is_relevant_label(label_no);
    }
#endif
}

void LabelledTransitionSystem::kill_transition(int src, int label, int target) {
    auto group = label_group_of_label[label];

    if(label_groups[group.group].size() == 1) {
	LTSTransition t(src, target, LabelGroup(group));
	kill_from_vector(t, transitions);
	kill_from_vector(t, transitions_src[src]);
	kill_from_vector(TSTransition(src, target), transitions_label_group[group.group]);
    } else {
	LabelGroup new_group (transitions_label_group.size());
	transitions_label_group.push_back(transitions_label_group[group.group]);
	kill_from_vector(TSTransition(src, target), transitions_label_group[new_group.group]);
	for (const auto & t :  transitions_label_group[new_group.group]) {
	    transitions.push_back(LTSTransition(t.src, t.target, new_group));
	    transitions_src[t.src].push_back(LTSTransition(t.src, t.target, new_group));
	}
	label_groups[group.group].erase(remove(begin(label_groups[group.group]),
					 end(label_groups[group.group]), label),
				  end(label_groups[group.group]));
	label_groups[new_group.group].push_back(label);
	label_group_of_label[label] = new_group;
    }
}


void LabelledTransitionSystem::kill_label(int l) {
    cout << "KILL " << l << endl;
    auto group = label_group_of_label[l];
    if (group.dead()) {
	irrelevant_labels.erase(remove(begin(irrelevant_labels), end(irrelevant_labels), l), end(irrelevant_labels));
    } else {
	label_group_of_label[l] = LabelGroup(-1);

	relevant_labels.erase(remove(begin(relevant_labels), end(relevant_labels), l), end(relevant_labels));
	label_groups[group.group].erase(remove(begin(label_groups[group.group]), end(label_groups[group.group]), l), end(label_groups[group.group]));
	if(label_groups[group.group].empty()) {
	    //Kill group
	    std::vector <TSTransition>().swap(transitions_label_group[l]);

	    transitions.erase(std::remove_if(begin(transitions),
					     end(transitions),
					     [&](LTSTransition & t){
						 return t.label_group == group;
					     }), end(transitions));

	    for (auto & trs : transitions_src){
		trs.erase(std::remove_if(begin(trs),
					 end(trs),
					 [&](LTSTransition & t){
					     return t.label_group == group;
					 }), end(trs));
	    }

	}
    }
}



void LabelledTransitionSystem::dump() const {
    for (int s = 0; s < size(); s++) {
	applyPostSrc(s, [&](const LTSTransition & trs) {
		cout << trs.src << " -> " << trs.target << " (" << trs.label_group.group << ":";
		for(int tr_s_label : get_labels(trs.label_group)) {
		    cout << " " << tr_s_label;
		}
		cout << ")\n";
		return false;
	    });
    }

}

bool LabelledTransitionSystem::is_self_loop_everywhere_label(int label) const {
    if (!is_relevant_label(label)) return true;
    const auto & trs = get_transitions_label(label);
    if (trs.size() < num_states) return false;

    // This assumes that there is no repeated transition
    int num_self_loops = 0;
    for (const auto & tr : trs) {
        if (tr.src == tr.target) {
            num_self_loops ++;
        }
    }

    assert(num_self_loops <= num_states);
    return num_self_loops == num_states;
}
