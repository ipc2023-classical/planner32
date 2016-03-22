#ifndef SYM_BREADTH_FIRST_SEARCH_H
#define SYM_BREADTH_FIRST_SEARCH_H

#include <vector>
#include <iostream>

#include "sym_bucket.h"
#include "sym_exploration.h"

class SymBreadthFirstSearch : public SymExploration  {  
  Bucket open;   // States in open 
  Bucket closed; // States in closed
  SymStepCostEstimation estimation;

  void filterDuplicates(Bucket & bucket) {
      for (BDD & bdd : bucket)
	  for(const BDD & c : closed)
	      bdd *=  !c;
  }
 public:
  SymBreadthFirstSearch(const SymParamsSearch & params);
  SymBreadthFirstSearch(const SymBreadthFirstSearch & ) = delete;
  SymBreadthFirstSearch(SymBreadthFirstSearch &&) = default;
  SymBreadthFirstSearch& operator=(const SymBreadthFirstSearch& ) = delete;
  SymBreadthFirstSearch& operator=(SymBreadthFirstSearch &&) = default;
  ~SymBreadthFirstSearch() {}


  bool init(SymManager * manager, bool forward);

  BDD pop();

  inline bool finished() const {
    return open.empty(); 
  }

  virtual bool stepImage(int maxTime, int maxNodes);

  virtual long nextStepTime() const;
  virtual long nextStepNodes() const;
  virtual long nextStepNodesResult() const;

  virtual bool isUseful(double ratio) const;
  virtual bool isSearchableWithNodes(int maxNodes) const; 
  void violated(TruncatedReason reason , double time, int maxTime, int maxNodes);

};
#endif 
