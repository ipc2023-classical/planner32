#ifndef SYM_BREADTH_FIRST_SEARCH_H
#define SYM_BREADTH_FIRST_SEARCH_H

#include <vector>
#include <iostream>

#include "sym_bucket.h"
#include "sym_exploration.h"

class SymBreadthFirstSearch : public SymExploration  {  
    Bucket open;   // States in open 
    //Bucket closed; // States in closed
    BDD closedTotal;
    
    SymStepCostEstimation estimation;

    SymBreadthFirstSearch * parent;

    void filterDuplicates(BDD & bdd) {
	/* for(const BDD & c : closed) */
	bdd *=  !closedTotal;
    }
    
    void filterDuplicates(Bucket & bucket) {
	for (BDD & bdd : bucket)
	    filterDuplicates(bdd);
    }

    void close (const BDD & bdd) {
	closedTotal += bdd;
	//closed.push_back(bdd);
	//mgr->mergeBucket(closed);
    }
 public:
  SymBreadthFirstSearch(const SymParamsSearch & params);
  SymBreadthFirstSearch(const SymBreadthFirstSearch & ) = delete;
  SymBreadthFirstSearch(SymBreadthFirstSearch &&) = default;
  SymBreadthFirstSearch& operator=(const SymBreadthFirstSearch& ) = delete;
  SymBreadthFirstSearch& operator=(SymBreadthFirstSearch &&) = default;
  ~SymBreadthFirstSearch() {}

  bool init(SymManager * manager, bool forward);

  bool init(SymBreadthFirstSearch * other, SymManager * manager, 
	    int maxRelaxTime, int maxRelaxNodes);


  BDD pop();

  virtual bool finished() const {
    return open.empty(); 
  }
  
  BDD getUnreachableStates() const ;

  bool foundSolution () const {
      if (parent && parent->foundSolution()) return true;

      BDD target = (fw ? mgr->getGoal() : mgr->getInitialState());
      return !((closedTotal*target).IsZero()); 
      /* for (auto & bdd : closed) { */
      /* 	  if (!((bdd*target).IsZero())) { */
      /* 	      return true; */
      /* 	  } */
      /* } */
  }

  bool isBetter(const SymBreadthFirstSearch & other) const{
      return nextStepTime() < other.nextStepTime();
  }

  virtual void getHeuristic(std::vector<ADD> & /*heuristics*/,
			    std::vector <int> & /*maxHeuristicValues*/) const {
      
  }

  virtual bool stepImage(int maxTime, int maxNodes);

  virtual long nextStepTime() const;
  virtual long nextStepNodes() const;
  virtual long nextStepNodesResult() const;

  virtual bool isUseful(double ratio) const;
  virtual bool isSearchableWithNodes(int maxNodes) const; 
  void violated(TruncatedReason reason , double time, int maxTime, int maxNodes);

  void notifyMutexes (const BDD & bdd);

};
#endif 
