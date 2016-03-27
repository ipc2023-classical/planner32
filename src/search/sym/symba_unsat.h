#ifndef SYM_SYMBA_UNSAT_H
#define SYM_SYMBA_UNSAT_H

#include "../search_engine.h"
#include "sym_controller.h"
#include "sym_hnode.h"
#include "sym_bdexp.h"
#include "sym_enums.h"

#include <set>
#include <map>
#include <memory>

class SymBreadthFirstSearch;
class UCTNode;

class SymBAUnsat : public SearchEngine, public SymController{
  Dir searchDir; //Direction of search in the original state space
  Dir abstractDir; //Direction of search in the abstract state space

//Common parameters to every hierarchy policy
  const SymParamsMgr mgrParams; 
  const SymParamsSearch searchParams; //Parameters to perform the abstract searches
  const double phTime, phMemory;

//Maximum time and nodes to perform the whole? step? relaxation process 
  const int maxRelaxTime, maxRelaxNodes;
 
//How to compute the TRs of the abstract state space.
  const AbsTRsStrategy absTRsStrategy;
 
//Parameters to decide the relaxation 
  const bool perimeterPDBs;  //Initializes explorations with the one being relaxed.
  const double ratioRelaxTime, ratioRelaxNodes; 

//Whether the ph should use mutexes
  const bool use_mutex_in_abstraction;

  const double shouldAbstractRatio;
  const int maxNumAbstractions;

  //Constant for UCT formula
  double UCT_C;
  UCTRewardType rewardType;

  bool add_abstract_to_ongoing_searches;

  int numAbstractions;
  // List of hierarchy policies to derive new abstractions
  //std::vector <SymPH *> phs;
  //std::unique_ptr<UCTTree> ph;

  std::vector<std::unique_ptr<UCTNode> > nodes;
  std::map<std::set<int>, UCTNode *> nodesByPattern;

  bool askHeuristic();  

  //Statistics; 
  double time_step_abstract, time_step_original, 
      time_select_exploration, time_notify_mutex, 
      time_init;

  std::vector<SymBreadthFirstSearch *> ongoing_searches;

  SymBreadthFirstSearch * selectExploration();

  std::vector <BDD> dead_end_fw, dead_end_bw;

  void insertDeadEnds(BDD bdd, bool isFW);

  //std::pair<UCTNode *, bool> relax(std::vector<UCTNode *> & uct_trace);
  UCTNode * relax(UCTNode * node,  bool fw, std::vector<UCTNode *> & uct_trace); 

  void notifyFinishedAbstractSearch(SymBreadthFirstSearch * currentSearch, double time_spent,
				    const std::vector<UCTNode *> & uct_trace); 

  void notifyFinishedAbstractSearch(SymBreadthFirstSearch * currentSearch) {
      return notifyFinishedAbstractSearch(currentSearch, 0, std::vector<UCTNode *> ());
  }

  double computeReward (const BDD & bdd, double time_spent) const; 

  bool chooseDirection() const;
 public:

 UCTNode * getUCTNode (const std::set<int> & pattern);

  SymBAUnsat(const Options &opts);
  virtual ~SymBAUnsat(){}

  virtual void initialize();
  virtual int step();

  virtual void print_options() const;
  virtual void statistics() const;

  virtual bool proves_task_unsolvable() const {
      return true;
  }

  UCTNode * getRoot () const {
      return nodes[0].get();
  }

};



#endif
