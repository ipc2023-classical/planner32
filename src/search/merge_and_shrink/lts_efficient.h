#ifndef LTS_EFFICIENT_H
#define LTS_EFFICIENT_H

#include <functional>

#include <vector>
#include <map>
#include <string>

typedef int AbstractStateRef;
class Abstraction;
class SimulationRelation;

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

class LTSTransitionEfficient {
 public:
  AbstractStateRef src, target;
  int label;
  //int sl; //Summary of src, label pair

 LTSTransitionEfficient(AbstractStateRef _src, AbstractStateRef _target, int _label) : 
  src(_src), target(_target), label(_label)/*, sl(0)*/ {
  }

  LTSTransitionEfficient(const LTSTransitionEfficient & t) = default;

  bool operator==(const LTSTransitionEfficient &other) const {
    return src == other.src && target == other.target && label == other.label;
  }

  bool operator!=(const LTSTransitionEfficient &other) const {
    return !(*this == other);
  }

  bool operator<(const LTSTransitionEfficient &other) const {
    return src < other.src || (src == other.src && target < other.target) 
      || (src == other.src && target == other.target && label < other.label);
  }

  bool operator>=(const LTSTransitionEfficient &other) const {
    return !(*this < other);
  }

  friend std::ostream & operator << (std::ostream& o , const LTSTransitionEfficient & t);
};

//Alvaro: Class added to implement the simple simulation
class LTSEfficient {
  Abstraction * abs;

  //Duplicated from abstraction
  int num_states;
  std::vector <bool> goal_states;
  AbstractStateRef init_state;
  std::vector<int> relevant_labels;
  std::vector<int> irrelevant_labels;
  std::vector <std::string> name_states;

  //List of transitions sorted (src, label) and (target, label)
  std::vector <LTSTransitionEfficient> transitionsPost, transitionsPre;
  std::vector<Qa> qaPre, qaPost;
  std::map<int, std::map <int, int> > qaPre_map, qaPost_map;
  //Map [label][state] -> qa

  void set_sl(std::vector <LTSTransitionEfficient> & transitions, 
	      std::vector<Qa> & qa, std::map<int, std::map<int, int> > & qaMap,
	      std::function<int (const LTSTransitionEfficient &)> fget);

 public:
  LTSEfficient (Abstraction * abs);
  ~LTSEfficient(){}

  const std::vector<bool> & get_goal_states() const {
    return goal_states;
  }

  inline int size() const{
    return num_states;
  }

  int num_transitions() const{
      return transitionsPost.size();
  }

  const std::vector <LTSTransitionEfficient> & get_transitions_post() const {
      return transitionsPost;
  }

  const std::vector <LTSTransitionEfficient> & get_transitions_pre() const {
      return transitionsPre;
  }

  bool hasQaPre(int label, int target) const {
      return qaPre_map.count(label) && qaPre_map.at(label).count(target);
  }

  bool hasQaPost(int label, int src) const {
      return qaPost_map.count(label) && qaPost_map.at(label).count(src);
  }

  const Qa & get_qa_pre(int label, int target) const {
      return qaPre[qaPre_map.at(label).at(target)];
  }

  const Qa & get_qa_post(int label, int src) const {
      return qaPost[qaPost_map.at(label).at(src)];
  }

  const std::vector<Qa> & get_qa_post() const {
      return qaPost;
  }

  const std::vector<Qa> & get_qa_pre() const {
      return qaPre;
  }

  int get_pos_qa_pre(int label, int target) const {
      if(!hasQaPre(label, target)) return -1;
      return qaPre_map.at(label).at(target);
  }

  int get_pos_qa_post(int label, int src) const {
      if(!hasQaPost(label, src)) return -1;
      return qaPost_map.at(label).at(src);
  }

  /* //Given s-l> t and l' check whether exists s-l'> t', t' >= t */
  /* bool check (int label, const LTSTransitionEfficient & tr,  */
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
			       std::function<bool(const LTSTransitionEfficient & tr)> && f) const {
      if(hasQaPost(label, src)){
	  const Qa & qa = get_qa_post(label, src);
	  for(int i = qa.b; i <= qa.e; ++i){
	      if(f(transitionsPost[i])) return true;
	  
	  }
      }
      return false;
  }


  //For each transition labelled with l, apply a function. If returns true, applies a break
  bool applyPost(int label, 
		 std::function<bool(const LTSTransitionEfficient & tr)> && f) const {
      if(qaPost_map.count(label)){
	  for (auto & qai : qaPost_map.at(label)){
	      const Qa & qa = qaPost[qai.second];
	      for(int i = qa.b; i <= qa.e; ++i){
		  if(f(transitionsPost[i])) return true;
	      }
	  }
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

  inline Abstraction * get_abstraction()  {
      return abs;
  }
 
};

#endif
