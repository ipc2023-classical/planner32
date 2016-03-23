#ifndef SYM_EXPLORATION_H
#define SYM_EXPLORATION_H

#include "../debug.h" 

#include "sym_manager.h"
#include "sym_estimate.h"
#include "sym_util.h"
#include <vector>
#include <map>
#include <memory>

class SymManager;

//We use this enumerate to know why the current operation was truncated
enum class TruncatedReason {
  FILTER_MUTEX, MERGE_BUCKET, MERGE_BUCKET_COST, IMAGE_ZERO, IMAGE_COST
};
std::ostream & operator<<(std::ostream &os, const TruncatedReason & dir);


class SymExpStatistics {
public:
    double image_time, image_time_failed;
    double time_heuristic_evaluation;
    int num_steps_succeeded; 
    double step_time;

    SymExpStatistics() :
    image_time (0), 
	image_time_failed  (0), time_heuristic_evaluation(0), 
	num_steps_succeeded  (0), step_time(0) {  }
	

    void add_image_time(double t) {
	image_time += t;
	num_steps_succeeded += 1;
    }

    void add_image_time_failed(double t) {
	image_time += t;
	image_time_failed += t;
	num_steps_succeeded += 1;
    }
};

class SymExploration  { 
protected: 
  //Attributes that characterize the search:
  SymManager * mgr;            //Symbolic manager to perform bdd operations
  SymParamsSearch p;
  bool fw; //Direction of the search. true=forward, false=backward 

  SymExpStatistics stats;

 public: 
  SymExploration (const SymParamsSearch & params);

  inline bool isFW() const{
    return fw;
  }

  inline bool isAbstracted() const{
    return mgr->getAbstraction() != nullptr &&
      mgr->getAbstraction()->isAbstracted();
  }

  inline bool isOriginal() const{
      return mgr->getAbstraction() == nullptr ||
	  !mgr->getAbstraction()->isAbstracted();
  }

  SymAbstraction * getAbstraction() const{
    return mgr->getAbstraction();
    
  }

  inline bool isUseful() const {
      return isUseful(p.ratioUseful);
  }

  inline bool isSearchable() const {
      return isSearchableWithNodes(p.maxStepNodes);
  }

  inline bool isSearchableAfterRelax(int num_relaxations) const {
    double maxNodes = p.maxStepNodes;
    if(num_relaxations){
      maxNodes *= pow(p.ratioAfterRelax, num_relaxations);
    }
    return isSearchableWithNodes((int)maxNodes);
  }

  bool step(){
    return stepImage(p.getAllotedTime(nextStepTime()), 
		     p.getAllotedNodes(nextStepNodesResult()));
  }

  virtual bool stepImage(int maxTime, int maxNodes) = 0;

  virtual long nextStepTime() const = 0;
  virtual long nextStepNodes() const = 0;
  virtual long nextStepNodesResult() const = 0;

  void statistics() const;


  virtual bool isUseful(double ratio) const = 0;
  virtual bool isSearchableWithNodes(int maxNodes) const = 0;   

};
#endif // SYMBOLIC_EXPLORATION

