#ifndef MERGE_AND_SHRINK_LABELLED_TRANSITION_SYSTEM_H
#define MERGE_AND_SHRINK_LABELLED_TRANSITION_SYSTEM_H

#include <functional>
#include <vector>
#include <string>
#include <algorithm>    // std::find

typedef int AbstractStateRef;
class Abstraction;
class LabelMap;

class LTSTransition {
 public:
  AbstractStateRef src, target;
  int label;
 LTSTransition(AbstractStateRef _src, AbstractStateRef _target, int _label) : 
  src(_src), target(_target), label(_label) {
  }

 LTSTransition(const LTSTransition & t) : 
  src(t.src), target(t.target), label(t.label) {
  }

  bool operator==(const LTSTransition &other) const {
    return src == other.src && target == other.target && label == other.label;
  }

  bool operator!=(const LTSTransition &other) const {
    return !(*this == other);
  }

  bool operator<(const LTSTransition &other) const {
    return src < other.src || (src == other.src && target < other.target) 
      || (src == other.src && target == other.target && label < other.label);
  }

  bool operator>=(const LTSTransition &other) const {
    return !(*this < other);
  }
};

//Alvaro: Class added to implement the simple simulation
class LabelledTransitionSystem {
  Abstraction * abs;

  //Duplicated from abstraction
  int num_states;
  std::vector <bool> goal_states;
  AbstractStateRef init_state;
  std::vector<int> relevant_labels;
  std::vector<int> irrelevant_labels;

  std::vector <std::string> name_states;
  std::vector <LTSTransition> transitions;
  std::vector<std::vector <LTSTransition> > transitions_src;
  std::vector<std::vector <LTSTransition> > transitions_label;


  inline void kill_from_vector(const LTSTransition & t, std::vector <LTSTransition> & v) {
      auto it = std::find(begin(v), end(v), t);
      *it = v.back();
      v.pop_back();
  }

 public:
  LabelledTransitionSystem (Abstraction * abs, const LabelMap & labelMap);
  ~LabelledTransitionSystem(){}

  const std::vector<bool> & get_goal_states() const {
    return goal_states;
  }

  inline int size() const{
    return num_states;
  }

  int num_transitions() const{
      return transitions.size();
  }

  const std::vector<LTSTransition> & get_transitions() const{
    return transitions;
  }

  const std::vector<LTSTransition> & get_transitions(int sfrom) const {
    return transitions_src[sfrom];
  }

  const std::vector<LTSTransition> & get_transitions_label(int label) const {
    return transitions_label[label];
  }

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

  inline Abstraction * get_abstraction()  {
    return abs;
  }

  void kill_transition(int src, int label, int target); 

  //For each transition labelled with l, applya a function. If returns true, applies a break
  bool applyPostSrc(int from,
		    std::function<bool(const LTSTransition & tr)> && f) const {
      for(const auto & tr : transitions_src[from]){
	  if(f(tr)) return true;
      }
      return false;
  }
};

#endif
