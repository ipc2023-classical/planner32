
#ifndef MERGE_AND_SHRINK_LABELLED_TRANSITION_SYSTEM_H
#define MERGE_AND_SHRINK_LABELLED_TRANSITION_SYSTEM_H

#include <functional>
#include <vector>
#include <string>
#include <algorithm>    // std::find
#include <cassert>

typedef int AbstractStateRef;
class Abstraction;
class LabelMap;

struct LabelGroup {
    int group;

    explicit LabelGroup(int g) : group(g) {
    }
    
    LabelGroup(const LabelGroup & other) = default;

    bool operator==(const LabelGroup &other) const {
	return group == other.group;
    }

    bool operator!=(const LabelGroup &other) const {
	return !(*this == other);
    }

    LabelGroup & operator++() {
	group++;
	return *this;
    }

    bool operator<(const  LabelGroup &other) const {
	return group < other.group;
    }

    bool dead() const {
	return group == -1;
    }
};
class LTSTransition {
 public:
  AbstractStateRef src, target;
  LabelGroup label_group;
 LTSTransition(AbstractStateRef _src, AbstractStateRef _target, LabelGroup _label) : 
  src(_src), target(_target), label_group(_label) {
  }

  LTSTransition(const LTSTransition & t) = default;

  bool operator==(const LTSTransition &other) const {
    return src == other.src && target == other.target && label_group == other.label_group;
  }

  bool operator!=(const LTSTransition &other) const {
    return !(*this == other);
  }

  bool operator<(const LTSTransition &other) const {
    return src < other.src || (src == other.src && target < other.target) 
      || (src == other.src && target == other.target && label_group < other.label_group);
  }

  bool operator>=(const LTSTransition &other) const {
    return !(*this < other);
  }
};

struct TSTransition {
    AbstractStateRef src, target;
TSTransition(AbstractStateRef _src, AbstractStateRef _target) : 
    src(_src), target(_target) {
    }

    TSTransition(const TSTransition & t) = default;

    bool operator==(const TSTransition &other) const {
	return src == other.src && target == other.target;
    }

    bool operator!=(const TSTransition &other) const {
	return !(*this == other);
    }

    bool operator<(const TSTransition &other) const {
	return src < other.src || (src == other.src && target < other.target) ;
    }

    bool operator>=(const TSTransition &other) const {
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

  std::vector<std::vector<int> > label_groups;
  std::vector<LabelGroup> label_group_of_label;
  std::vector <std::string> name_states;
  std::vector <LTSTransition> transitions;
  std::vector<std::vector <LTSTransition> > transitions_src;
  std::vector<std::vector <TSTransition> > transitions_label_group;


  template <typename T> 
  inline void kill_from_vector(const T & t, std::vector <T> & v) {
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

  bool is_goal(int state) const {
      return goal_states[state];
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

  const std::vector<TSTransition> & get_transitions_label(int label) const {
    return transitions_label_group[label_group_of_label[label].group];
  }
  const std::vector<TSTransition> & get_transitions_label_group(LabelGroup label_group) const {
    return transitions_label_group[label_group.group];
  }

  const std::vector<std::string> & get_names () const {
    return name_states;
  }

  const std::string & name (int s) const {
    return name_states[s];
  }

  int is_relevant_label(int label) const {
#ifndef NDEBUG
      bool relevant1 = std::find(relevant_labels.begin(), relevant_labels.end(), label) != relevant_labels.end();
      bool relevant2 = std::find(irrelevant_labels.begin(), irrelevant_labels.end(), label) == irrelevant_labels.end();
      bool relevant3 = !label_group_of_label[label].dead() &&
	  !transitions_label_group[label_group_of_label[label].group].empty();
      assert (relevant1 == relevant2);
      assert (relevant3 == relevant2);
#endif
      return std::find(relevant_labels.begin(), relevant_labels.end(), label) != relevant_labels.end();
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

  void kill_label(int l);

  void kill_transition(int src, int label, int target); 

  //For each transition labelled with l, applya a function. If returns true, applies a break
  bool applyPostSrc(int from,
		    std::function<bool(const LTSTransition & tr)> && f) const {
      for(const auto & tr : transitions_src[from]){
	  if(f(tr)) return true;
      }
      return false;
  }
  const std::vector<int> & get_labels(LabelGroup label_group) const {
      return label_groups[label_group.group];
  }

  LabelGroup get_group_label(int label) const {
      return LabelGroup(label_group_of_label[label]);
  }

  int get_num_label_groups() const {
      return label_groups.size();	
  }

  const std::vector<LabelGroup> & get_group_of_label() const {
      return label_group_of_label;
  }

  void dump() const; 
};

#endif
