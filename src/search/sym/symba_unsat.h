#ifndef SYM_SYMBA_UNSAT_H
#define SYM_SYMBA_UNSAT_H

#include "sym_engine.h"
#include <set>

class SymBreadthFirstSearch;

class SymBAUnsat : public SymEngine{
  int currentPH;

  //Statistics; 
  double time_step_abstract, time_step_original, time_select_exploration;

  std::vector<std::unique_ptr<SymBreadthFirstSearch> > ongoing_searches;

  SymBreadthFirstSearch * selectExploration();

  std::vector <BDD> dead_end_fw, dead_end_bw;

  void insertDeadEnds(BDD bdd, bool isFW) {
      if(isFW){
	  dead_end_fw.push_back(bdd);
      } else {
	  dead_end_bw.push_back(bdd);
      }
      
  }

 public:
  SymBAUnsat(const Options &opts);
  virtual ~SymBAUnsat(){}

  virtual void initialize();
  virtual int step();

  virtual void print_options() const;
  virtual void statistics() const;

  virtual bool proves_task_unsolvable() const {
      return true;
  }
};



#endif
