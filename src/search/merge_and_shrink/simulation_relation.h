#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_H

#include <vector>
#include <string>
#include <iostream>
#include "abstraction.h"
#include "label_relation.h"
#include "../sym/sym_variables.h"

class Labels;

// First implementation of a simulation relation. 
class SimulationRelation{
  const Abstraction * abs;

  //By now we assume that the partition is unitary... we can improve
  //this later with EquivalenceRelation
  std::vector<std::vector<bool> > relation;

  //For each abstract state, we create a BDD that represents all the
  //abstract states dominated by it and dominating it
  std::vector<BDD> dominated_bdds, dominating_bdds; 

  //BDDs of each abstract state
  std::vector<BDD> abs_bdds;
  BDD zeroBDD;

 public:
  SimulationRelation(const Abstraction * _abs, int num_states,
		     const std::vector<bool> & goal_states);

  inline bool simulates (int s, int t) const {
    return relation[s][t];
  }

  inline bool similar (int s, int t) const {
    return relation[s][t] && relation[t][s];
  }


  inline void remove (int s, int t) {
    relation[s][t] = false;
  }

  int num_equivalences() const;
  int num_simulations() const;
  int num_states() const{
    return abs_bdds.size();
  }

  int num_different_states() const;

  /*
   * THIS IMPLEMENTATION IS VERY INNEFICIENT
   * ONLY TO BE USED AS A PROOF OF CONCEPT
   */
  //Template class method that gets as input a list of LTSs and
  // returns a list of simulation relations over them
  template<class LTS>
    static void compute_label_dominance_simulation(const std::vector<LTS *> & lts_list, 
						   Labels * labels, 
						   std::vector<SimulationRelation *> & res){
  //Initialization phase
  for(auto & lts : lts_list){
    //Create initial goal-respecting relation
    res.push_back(new SimulationRelation(lts->get_abstraction(),
					 lts->size(), lts->get_goal_states()));
  }
  LabelRelation label_dominance(labels);
  label_dominance.init(lts_list, res);

  do{
    std::cout << "Loop" << std::endl;
    //label_dominance.dump();
    for (int i = 0; i < lts_list.size(); i++){
      res[i]->update(i, lts_list[i], label_dominance);
      //res[i]->dump(lts_list[i]->get_names());
    }
  }while(label_dominance.update(lts_list, res));
  //label_dominance.dump();
}

  template<class LTS>  
    void update(int lts_id, const LTS * lts, LabelRelation & label_dominance){
    bool changes = true;
    while(changes){
      changes = false;
      for (int s = 0; s < lts->size(); s++){
	for (int t = 0; t < lts->size(); t++){ //for each pair of states t, s
	  if(s != t && simulates(t, s)){
	    //Check if really t simulates s
	    //for each transition s--l->s':
	    // a) with noop t >= s' and l dominated by noop?
	    // b) exist t--l'-->t', t' >= s' and l dominated by l'?
	    for (auto trs : lts->get_transitions(s)){
	      if(simulates (t, trs.target) && 
		 label_dominance.dominated_by_noop(trs.label, lts_id)){
		continue;
	      }
	      bool found = false;
	      for (auto trt : lts->get_transitions(t)){
		if(label_dominance.dominates(trt.label, trs.label, lts_id) && 
		   simulates(trt.target, trs.target)){		
		  found = true;
		  break;
		}
	      }
	      if(!found){
		changes = true;
		remove(t, s);
		/*std::cout << lts->name(t) << " does not simulate " <<  lts->name(s) 
			  << " because of " <<
		  lts->name(trs.src)  << " => " << lts->name(trs.target) << " (" << trs.label << ")";// << std::endl;
		std::cout << "  Simulates? "<<simulates (trs.src, trs.target);
		std::cout << "  domnoop? "<<label_dominance.dominated_by_noop(trs.label, lts_id) << "   ";
		label_dominance.dump(trs.label);
		for (auto trt : lts->get_transitions(t)){
		  std::cout << "Tried with: " << lts->name(trt.src)  << " => " << lts->name(trt.target) << " (" << trt.label << ")" << " label dom: "<< label_dominance.dominates(trt.label, trs.label, lts_id) << " target sim " << simulates(trt.target, trs.target) << std::endl;
		}*/

		break;
	      }
	    }
	  }
	}
      }
    }
  }
  void dump(const std::vector<std::string> & names) const;
 
 BDD getSimulatedBDD(const State & state) const;
 BDD getSimulatingBDD(const State & state) const;

 void precompute_absstate_bdds(SymVariables * vars);
 void precompute_dominated_bdds();
 void precompute_dominating_bdds();

 inline const std::vector<BDD> & get_dominated_bdds () {
   if(dominated_bdds.empty()) precompute_dominated_bdds();
   return dominated_bdds;
 }

 inline const std::vector<BDD> & get_dominating_bdds () {
   if(dominating_bdds.empty()) precompute_dominating_bdds();
   return dominating_bdds;
 }


 inline const std::vector<BDD> & get_abs_bdds() const{
   return abs_bdds;
 }

 inline const std::vector <int> & get_varset() const {
   return abs->get_varset();
 }

 inline bool pruned(const State & state) const {
   return abs->get_abstract_state(state) == -1;
 }
 
};








#endif
