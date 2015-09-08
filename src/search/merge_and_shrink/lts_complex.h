#ifndef LTS_COMPLEX_H
#define LTS_COMPLEX_H

#include "labelled_transition_system.h"//For getting LTSTransition

#include <functional>

#include <vector>
#include <map>
#include <string>

typedef int AbstractStateRef;
class Abstraction;
class SimulationRelation;
class LabelMap;

class Qa {
public:
    int index;
    int state;
    int label;
    int b,e; //range of transitions with this Qa

    bool markRemove; //Mark to indicate whether it has been inserted in the remove list

    Qa(int i, int s, int l, int b_, int e_) :
        index(i), state(s), label(l), b(b_), e(e_),  markRemove(false) {}

    Qa(const Qa & o) = default;
    Qa(): index(0), state(0), label(0), b(0), e(0),  markRemove(false){}

    int size() const {
        return e - b + 1;
    }
};

//Alvaro: Class added to implement the simple simulation
class LTSComplex {
    Abstraction * abs;

    //Duplicated from abstraction
    int num_states;
    std::vector <bool> goal_states;
    AbstractStateRef init_state;
    std::vector<int> relevant_labels;
    std::vector<int> irrelevant_labels;
    std::vector <std::string> name_states;

    //List of transitions sorted (src, label) and (target, label)
    std::vector <LTSTransition> transitionsPost, transitionsPre;
    std::vector<Qa> qaPre, qaPost;
    std::vector<std::map <int, int> > qaPre_map, qaPost_map;
    std::vector<std::pair<int, int> > rangePostSrc, rangePreTarget;
    //Map [label][state] -> qa

    void set_sl(std::vector <LTSTransition> & transitions,
            std::vector<Qa> & qa, std::vector<std::map<int, int> > & qaMap,
		std::vector<std::pair<int, int> > & rangeStates,
		std::function<int (const LTSTransition &)> fget);

public:
    LTSComplex (Abstraction * abs, const LabelMap & labelMap);
    ~LTSComplex(){}

    const std::vector<bool> & get_goal_states() const {
        return goal_states;
    }

    inline int size() const{
        return num_states;
    }

    int num_transitions() const{
        return transitionsPost.size();
    }

    const std::vector <LTSTransition> & get_transitions_post() const {
        return transitionsPost;
    }

    const std::vector <LTSTransition> & get_transitions_pre() const {
        return transitionsPre;
    }

    bool hasQaPre(int label, int target) const {
        return qaPre_map [label].count(target);
        /* auto it = qaPre_map.find(label); */
        /* if (it != qaPre_map.end()) { */
        /*     return (*it).second.count(target); */
        /* } */
        /* return false; */
        //return qaPre_map.count(label) && qaPre_map.at(label).count(target);
    }

    bool hasQaPost(int label, int src) const {
        return qaPost_map [label].count(src);
        /* std::map<int, std::map<int, int> >::const_iterator it = qaPost_map.find(label); */
        /* if (it != qaPost_map.end()) { */
        /*     return (*it).second.count(src); */
        /* } */
        /* return false; */
        //return qaPost_map.count(label) && qaPost_map.at(label).count(src);
    }

    const Qa & get_qa_pre(int label, int target) const {
	return qaPre[qaPre_map[label].at(target)];
        //return qaPre[qaPre_map.at(label).at(target)];
    }

    const Qa & get_qa_post(int label, int src) const {
	return qaPost[qaPost_map[label].at(src)];
        //return qaPost[qaPost_map.at(label).at(src)];
    }

    const std::vector<Qa> & get_qa_post() const {
        return qaPost;
    }

    const std::vector<Qa> & get_qa_pre() const {
        return qaPre;
    }

    int get_pos_qa_pre(int label, int target) const {
	auto it = qaPre_map[label].find(target);
	if (it != qaPre_map[label].end())
	    return (*it).second;
	return -1;
	
        /* std::map<int, std::map<int, int> >::const_iterator it = qaPre_map.find(label); */
        /* if (it != qaPre_map.end()) { */
        /*     std::map<int, int>::const_iterator it2 = (*it).second.find(target); */
        /*     if (it2 != (*it).second.end()) */
        /*         return (*it2).second; */
        /* } */
        /* return -1; */
        //if(!hasQaPre(label, target)) return -1;
        //return qaPre_map.at(label).at(target);
    }

    int get_pos_qa_post(int label, int src) const {
	auto it = qaPost_map[label].find(src);
	if (it != qaPost_map[label].end())
	    return (*it).second;
	return -1;

        /* std::map<int, std::map<int, int> >::const_iterator it = qaPost_map.find(label); */
        /* if (it != qaPost_map.end()) { */
        /*     std::map<int, int>::const_iterator it2 = (*it).second.find(src); */
        /*     if (it2 != (*it).second.end()) */
        /*         return (*it2).second; */
        /* } */
        /* return -1; */
        //if(!hasQaPost(label, src)) return -1;
        //return qaPost_map.at(label).at(src);
    }

    /* //Given s-l> t and l' check whether exists s-l'> t', t' >= t */
    /* bool check (int label, const LTSTransition & tr,  */
    /* 	      const SimulationRelation * sim) const { */
    /*     if(qaPost_map.count(label) && qaPost_map.at(label).count(tr.src)){ */
    /* 	  const Qa & qa = qaPost_map.at(label).at(tr.src); */
    /* 	  for(int i = qa.b; i <= qa.e; ++i){ */
    /* 	      if(sim->simulates(transitionsPost[i].target, tr.target)){ */
    /* 		  return true; */
    /* 	      } */
    /* 	  } */
    /*     } */
    /*     return false; */
    /* } */


    //For each transition labelled with l, apply a function. If returns true, applies a break
    bool applyPost(int label, int src,
		   std::function<bool(const LTSTransition & tr)> && f) const {
	auto it = qaPost_map[label].find(src);
	if (it != qaPost_map[label].end()) {
	    const Qa & qa = qaPost[(*it).second];
	    for(int i = qa.b; i <= qa.e; ++i){
		if(f(transitionsPost[i])) return true;

	    }
	}
        return false;
    }


    //For each transition labelled with l, apply a function. If returns true, applies a break
    bool applyPost(int label,
            std::function<bool(const LTSTransition & tr)> && f) const {
	for (const auto & qai : qaPost_map[label]) {
	    const Qa & qa = qaPost[qai.second];
	    for(int i = qa.b; i <= qa.e; ++i){
		if(f(transitionsPost[i])) return true;
	    }
	}
        return false;
    }

    //For each transition labelled with l, apply a function. If returns true, applies a break
    /* bool applyPostSrc(int src, */
    /* 		      std::function<bool(const LTSTransition & tr)> && f) const { */
    /*     for(auto & pm : qaPost_map){ */
    /*         auto it = pm.find(src); */
    /*         if (it != pm.end()) { */
    /*             const Qa & qa = qaPost[(*it).second]; */
    /*             for(int i = qa.b; i <= qa.e; ++i){ */
    /*                 if(f(transitionsPost[i])) return true; */
    /*             } */
    /*         } */
    /*     } */
    /*     return false; */
    /* } */

    /* //For each transition labelled with l, apply a function. If returns true, applies a break */
    /* bool applyPreTarget(int target, */
    /*         std::function<bool(const LTSTransition & tr)> && f) const { */
    /*     for(auto & pm : qaPre_map){ */
    /*         auto it = pm.find(target); */
    /*         if (it != pm.end()) { */
    /*             const Qa & qa = qaPre[(*it).second]; */
    /*             for(int i = qa.b; i <= qa.e; ++i){ */
    /*                 if(f(transitionsPre[i])) return true; */
    /*             } */
    /*         } */
    /*     } */
    /*     return false; */
    /* } */

    bool applyPostSrc(int src,
    		      std::function<bool(const LTSTransition & tr)> && f) const {
	const std::pair<int, int> &p  = rangePostSrc[src];
	for(int i = p.first; i <= p.second; ++i){
	    if(f(transitionsPost[i])) return true;
	}
        return false;
    }

    bool applyPreTarget(int target,
    		      std::function<bool(const LTSTransition & tr)> && f) const {
	const std::pair<int, int> &p  = rangePreTarget[target];
	for(int i = p.first; i <= p.second; ++i){
	    if(f(transitionsPost[i])) return true;
	}
        return false;
    }



    void dump_names() const;

    const std::vector<std::string> & get_names () const {
        return name_states;
    }

    const std::string & name (int s) const {
        return name_states[s];
    }

    const std::vector<int> & get_irrelevant_labels() const {
        return irrelevant_labels;
    }

    const std::vector<int> & get_relevant_labels() const {
        return relevant_labels;
    }

    inline const Abstraction * get_abstraction() const {
        return abs;
    }

};

#endif
