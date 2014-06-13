#ifndef MERGE_AND_SHRINK_LABEL_RELATION_H
#define MERGE_AND_SHRINK_LABEL_RELATION_H

#include <iostream>
#include <vector>
#include "labels.h"
#include "label.h"

const int DOMINATES_IN_ALL = -2;
const int DOMINATES_IN_NONE = -1;

/* 
 * Label relation represents the preorder relations on labels that
 * occur in a set of LTS
 */ 
class LabelRelation {
  Labels * labels;
 public:

  //For each lts, matrix indicating whether l1 simulates l2
  //std::vector<std::vector<std::vector<bool > > > dominates;
  
  //Matrix for each l1, l2 indicating whether l1 dominates l2 in all
  //(-2), in none (-1) or only in i (i)
  std::vector<std::vector<int> > dominates_in;

  //Indicates whether labels are dominated by noop or other irrelevant
  //variables in theta
  std::vector<std::vector<bool> > simulated_by_irrelevant;
  std::vector<std::vector<bool> > simulates_irrelevant;

  std::vector<int> dominated_by_noop_in;

  LabelRelation(Labels * labels);

  template <typename LTS, typename SimRel> 
    void init(const std::vector<LTS*> & lts,
	      const std::vector<SimRel*> & sim){
    simulates_irrelevant.resize(labels->get_size());
    simulated_by_irrelevant.resize(labels->get_size());
    for(int i = 0; i < labels->get_size(); i++){
      simulates_irrelevant[i].resize(lts.size(), true);
      simulated_by_irrelevant[i].resize(lts.size(), true);
    }
    

    dominates_in.resize(labels->get_size());
    for (int l1 = 0; l1 < dominates_in.size(); ++l1){
      dominated_by_noop_in.resize(labels->get_size(), DOMINATES_IN_ALL);
      dominates_in[l1].resize(labels->get_size(), DOMINATES_IN_ALL);
      for (int l2 = 0; l2 < dominates_in[l1].size(); ++l2){
	if(labels->get_label_by_index(l1)->get_cost() > 
	   labels->get_label_by_index(l2)->get_cost()){
	  dominates_in[l1][l2] = DOMINATES_IN_NONE;
	}
      }
    }

    for (int i = 0; i < lts.size(); ++i){
      update(i, lts[i], sim[i]);
    }

  }

  template <typename LTS, typename SimRel> 
    bool update(const std::vector<LTS*> & lts,
		const std::vector<SimRel*> & sim){
    bool changes = false;
    for (int i = 0; i < lts.size(); ++i){
      changes |= update(i, lts[i], sim[i]);
    }
    return changes;
  }

  template <typename LTS, typename SimRel>
    bool update(int i, const LTS * lts, const SimRel * sim){
    bool changes = false;
    for(int l2 : lts->get_relevant_labels()) {
      for(int l1 : lts->get_relevant_labels()){ 
	if(l1 != l2 && simulates(l1, l2, i)){
	  //std::cout << "Check " << l1 << " " << l2 << std::endl;
	  //std::cout << "Num transitions: " << lts->get_transitions_label(l1).size() 
	  //		    << " " << lts->get_transitions_label(l2).size() << std::endl;
	  //Check if it really simulates
	  //For each transition s--l2-->t, and evey label l1 that dominates
	  //l2, exist s--l1-->t', t <= t'?
	  for(auto tr : lts->get_transitions_label(l2)){
	    bool found = false;
	    //TODO: for(auto tr2 : lts->get_transitions_for_label_src(l1, tr.src)){
	    for(auto tr2 : lts->get_transitions_label(l1)){
	      if(tr2.src == tr.src &&
		 sim->simulates(tr2.target, tr.target)){
		found = true;
		break; //Stop checking this tr
	      }
	    }
	    if(!found){
	      //std::cout << "Not sim " << l1 << " " << l2 << " " << i << std::endl;
	      set_not_simulates(l1, l2, i);
	      changes = true;
	      break; //Stop checking trs of l1
	    }
	  }
	}
      }

      //Is l2 simulated by irrelevant_labels in lts?
      for(auto tr : lts->get_transitions_label(l2)){
	if (simulated_by_irrelevant[l2][i] && 
	    !sim->simulates(tr.src, tr.target)) {
	  changes |= set_not_simulated_by_irrelevant(l2, i);
	  for (int l : lts->get_irrelevant_labels()){
	    if(simulates(l, l2, i)){
	      changes = true;
	      set_not_simulates(l, l2, i);
	    }
	  }
	}
      }
      //Does l2 simulates irrelevant_labels in lts?
      if(simulates_irrelevant[l2][i]){
	for(int s = 0; s < lts->size(); s++){
	  bool found = false;
	  for(auto tr : lts->get_transitions_label(l2)){
	    if(tr.src == s && sim->simulates(tr.target, tr.src)) {
	      found = true;
	      break;
	    }
	  }
	  if(!found) {
	    simulates_irrelevant[l2][i] = false;
	    for (int l : lts->get_irrelevant_labels()){
	      if(simulates(l2, l, i)){
		set_not_simulates(l2, l, i);
		changes = true;
	      }
	    }
	  }
	}
      } 
    }

    return changes;
  }
  
  //void update_after_merge(int i, int j, LTS & lts);
  void dump() const;
  void dump(int label) const;

  //Returns true if l dominates l2 in lts (simulates l2 in all j \neq lts)
  inline bool dominates (int l1, int l2, int lts){
    return dominates_in[l1][l2] != DOMINATES_IN_NONE &&
      (dominates_in[l1][l2] == DOMINATES_IN_ALL || 
       dominates_in[l1][l2] == lts);
  }

  //Returns true if l1 simulates l2 in lts
  inline bool simulates (int l1, int l2, int lts){
    return dominates_in[l1][l2] !=  DOMINATES_IN_NONE &&
      (dominates_in[l1][l2] == DOMINATES_IN_ALL || 
       dominates_in[l1][l2] != lts);
  }

  inline void set_not_simulates (int l1, int l2, int lts){
    //std::cout << "Not simulates: " << l1 << " to " << l2 << " in " << lts << std::endl;
    if(dominates_in[l1][l2] == DOMINATES_IN_ALL){
      dominates_in[l1][l2] = lts;
    }else if(dominates_in[l1][l2] != lts){
      dominates_in[l1][l2] = DOMINATES_IN_NONE;
    }else{
      std::cerr << "ERROR: RECOMPUTING INNECESSARILY" << std::endl;
      std::exit(0);
    }
  }

  inline bool set_not_simulated_by_irrelevant(int l, int lts){
    //Returns if there were changes in dominated_by_noop_in
    simulated_by_irrelevant[l][lts] = false;
    if(dominated_by_noop_in[l] == DOMINATES_IN_ALL){
      dominated_by_noop_in[l] = lts;
      return true;
    }else if(dominated_by_noop_in[l] != lts){
      dominated_by_noop_in[l] = DOMINATES_IN_NONE;
      return true;
    }
    return false;
  }

  inline bool dominated_by_noop (int l, int lts){
    return dominated_by_noop_in[l] != DOMINATES_IN_NONE &&
      (dominated_by_noop_in[l] == DOMINATES_IN_ALL || 
       dominated_by_noop_in[l] == lts);
  }

};

#endif
