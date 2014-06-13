#ifndef MERGE_AND_SHRINK_LABELLED_TRANSITION_SYSTEM_H
#define MERGE_AND_SHRINK_LABELLED_TRANSITION_SYSTEM_H


#include <vector>
#include <string>



class Abstraction;


class LTSTransition {
 public:
  int src, target, label;
 LTSTransition(int _src, int _target, int _label) : 
  src(_src), target(_target), label(_label) {
  }

 LTSTransition(const LTSTransition & t) : 
  src(t.src), target(t.target), label(t.label) {
  }

};

//Alvaro: Class added to implement the simple simulation
class LabelledTransitionSystem {
  Abstraction * abs;

  int num_states;
  std::vector <bool> goal_states;
  std::vector<int> relevant_labels;
  std::vector<int> irrelevant_labels;

  std::vector <std::string> name_states;
  std::vector <LTSTransition> transitions;
  std::vector<std::vector <LTSTransition> > transitions_src;
  std::vector<std::vector <LTSTransition> > transitions_label;

 public:
  LabelledTransitionSystem (Abstraction * abs);
  ~LabelledTransitionSystem(){}

  const std::vector<bool> & get_goal_states() const {
    return goal_states;
  }

  inline int size() const{
    return num_states;
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

  inline bool is_relevant_label(int l) const {
    return relevant_labels [l];
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
