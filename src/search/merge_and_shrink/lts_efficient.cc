#include "lts_efficient.h"

#include "abstraction.h"

using namespace std;

LTSEfficient::LTSEfficient (Abstraction * _abs) : 
    abs(_abs),  num_states(_abs->size()),goal_states(_abs->get_goal_states()){

    const vector<bool> & is_rel_label = abs->get_relevant_labels();

    for(int i = 0; i < num_states; i++){
	name_states.push_back(abs->description(i));
    }

    for (int label_no = 0; label_no < abs->get_num_labels(); label_no++) {
	if(is_rel_label[label_no]){
	    relevant_labels.push_back(label_no);
	    const vector<AbstractTransition> &abs_tr = abs->get_transitions_for_label(label_no);
	    for (int j = 0; j < abs_tr.size(); j++) {
		LTSTransitionEfficient t (abs_tr[j].src, abs_tr[j].target, label_no);
		transitionsPre.push_back(t);
		transitionsPost.push_back(t);
	    }
	}else{
	    irrelevant_labels.push_back(label_no);
	}
    }

    std::sort(begin(transitionsPre), 
	      end(transitionsPre), 
	      [] (const LTSTransitionEfficient & t, const LTSTransitionEfficient & t2){
		  return t.target < t2.target || (t.target == t2.target && t.label < t2.label)
		      || (t.target == t2.target && t.label == t2.label && t.src < t2.src);
	      });
    std::sort(begin(transitionsPost), 
	      end(transitionsPost), 
	      [] (const LTSTransitionEfficient & t, const LTSTransitionEfficient & t2){
		  return t.src < t2.src || (t.src == t2.src && t.label < t2.label)
		      || (t.src == t2.src && t.label == t2.label && t.target < t2.target);
	      });

    //cout << "Set Transitions pre: " << endl;
    set_sl(transitionsPre, qaPre, qaPre_map, [](const LTSTransitionEfficient & t){    
	    return t.target;
	});
    //cout << "Set Transitions post: " << endl;
    set_sl(transitionsPost, qaPost, qaPost_map, [](const LTSTransitionEfficient & t){    
	    return t.src;
	});
    
    // cout << "Generated lts efficient. Pre: " << endl;
    // for(const auto & t : transitionsPre){
    // 	cout << t << endl;
    // }
    // cout << "Post: " << endl;
    // for(const auto & t : transitionsPost){
    // 	cout << t << endl;
    // }
}

void LTSEfficient::set_sl(vector <LTSTransitionEfficient> & transitions,
			  vector <Qa> & qa, std::map<int, std::map<int, int> > & qaMap,
			  function<int (const LTSTransitionEfficient &)> fget) {
    int s, l, index = 0, qaindex= 0;
    for(int i = 0; i < transitions.size(); i++){
	//cout << "Introducing t: " << transitions[i] << endl;
	//transitions[i].sl = qaindex;
	if(i==0 || s != fget(transitions[i]) || l != transitions[i].label){
	    if(i > 0){
		qaMap[l][s] = qaindex;
		//cout << "Added qa: " << s << " " << l << " qa: " << qaindex << endl;
		Qa q {qaindex++, s, l, index, i-1};
		qa.push_back(q);
		index = i;
		//transitions[i].sl = qaindex;
	    }
	    s = fget(transitions[i]);
	    l = transitions[i].label;
	}
    }
    if(transitions.size() > 0){
	int lastIndex = transitions.size()-1;
	qaMap[l][s] = qaindex;
	//cout << "Added qa: " << s << " " << l << " qa: " << qaindex << endl;
	Qa q {qaindex, s, l, index, lastIndex};
	qa.push_back(q);
    }
}


void LTSEfficient::dump_names() const {
    cout << "LTS Names: "; 
    for (const auto & n : name_states) cout << " " << n; 
    cout << endl;
}


std::ostream & operator << (std::ostream& o , const LTSTransitionEfficient & t){
    return o  << t.src << " -- " << t.label << " --> " << " " << t.target; //<< " (" << t.sl << ")";
}


