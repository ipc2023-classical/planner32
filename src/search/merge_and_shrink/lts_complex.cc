#include "lts_complex.h"

#include "abstraction.h"
#include "labels.h"

using namespace std;

LTSComplex::LTSComplex (Abstraction * _abs, const LabelMap & labelMap) : 
    abs(_abs),  num_states(_abs->size()),goal_states(_abs->get_goal_states()){

    int num_labels = labelMap.get_num_labels();
    const vector<bool> is_rel_label(num_labels, false);    
    const vector<bool> &was_rel_label = abs->get_relevant_labels();

    for(int i = 0; i < num_states; i++){
	name_states.push_back(abs->description(i));
    }

    for (int label_no = 0; label_no < labelMap.get_num_labels(); label_no++) {
	int old_label = labelMap.get_old_id(label_no);
	if(was_rel_label[old_label]){
	    relevant_labels.push_back(label_no);
	    const vector<AbstractTransition> &abs_tr = abs->get_transitions_for_label(old_label);
	    for (int j = 0; j < abs_tr.size(); j++) {
		LTSTransition t (abs_tr[j].src, abs_tr[j].target, label_no);
		transitionsPre.push_back(t);
		transitionsPost.push_back(t);
	    }
	}else{
	    irrelevant_labels.push_back(label_no);
	}
    }

    std::sort(begin(transitionsPre), 
	      end(transitionsPre), 
	      [] (const LTSTransition & t, const LTSTransition & t2){
		  return t.target < t2.target || (t.target == t2.target && t.label < t2.label)
		      || (t.target == t2.target && t.label == t2.label && t.src < t2.src);
	      });
    std::sort(begin(transitionsPost), 
	      end(transitionsPost), 
	      [] (const LTSTransition & t, const LTSTransition & t2){
		  return t.src < t2.src || (t.src == t2.src && t.label < t2.label)
		      || (t.src == t2.src && t.label == t2.label && t.target < t2.target);
	      });

    qaPre_map.resize(num_labels);
    qaPost_map.resize(num_labels);
    rangePostSrc.resize(num_states);
    rangePreTarget.resize(num_states);
    //cout << "Set Transitions pre: " << endl;
    set_sl(transitionsPre, qaPre, qaPre_map, rangePreTarget, [](const LTSTransition & t){    
	    return t.target;
	});
    //cout << "Set Transitions post: " << endl;
    set_sl(transitionsPost, qaPost, qaPost_map, rangePostSrc, [](const LTSTransition & t){    
	    return t.src;
	});
    
    // cout << "Generated lts complex. Pre: " << endl;
    // for(const auto & t : transitionsPre){
    // 	cout << t << endl;
    // }
    // cout << "Post: " << endl;
    // for(const auto & t : transitionsPost){
    // 	cout << t << endl;
    // }
}

void LTSComplex::set_sl(vector <LTSTransition> & transitions,
			  vector <Qa> & qa, std::vector<std::map<int, int> > & qaMap,
			  vector<pair<int, int> > & rangeStates,
			  function<int (const LTSTransition &)> fget) {
    int s = 0, l = 0; /* set to 0 only to avoid warning */
    int index = 0, qaindex= 0;
    int indexS = 0;
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
		if(s != fget(transitions[i])){
		    rangeStates[s] = pair<int, int>(indexS, i-1);
		    indexS = i;
		}
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

void LTSComplex::dump_names() const {
    cout << "LTS Names: "; 
    for (const auto & n : name_states) cout << " " << n; 
    cout << endl;
}

std::ostream & operator << (std::ostream& o , const LTSTransition & t){
    return o  << t.src << " -- " << t.label << " --> " << " " << t.target; //<< " (" << t.sl << ")";
}
