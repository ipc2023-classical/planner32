#ifndef SYM_ASTAR_H
#define SYM_ASTAR_H

#include "../debug.h" 
#include "sym_exploration.h"
#include "sym_heuristic.h"
#include "sym_manager.h"
#include "sym_bucket.h"
#include "sym_astar_open.h"
#include "sym_astar_closed.h"
#include "sym_estimate.h"
#include "sym_util.h"
#include <vector>
#include <map>
#include <memory>

/*
 * This class allows to perform a BDD search.  It is designed to
 * mantain the current state in the search.  We consider four
 * different points at which we may truncate the search:
 * pop(), filter_mutex(), expand_zero(), expand_cost()
 * We mantain 3 BDDs to know the current state: Stmp, S and Szero.
 * Briefly:
 * 1) if Sfilter, Szero and S are empty => pop() => Szero.
 * 2) else if Stfilter => filter_mutex() => Szero
 * 3) else if Szero => expand_zero => S (passing by Sfilter)
 * 4) else (S must have something) => expand_cost()
 * 
 * Zero cost operators have been expanded iff !S.IsZero() && Szero.IsZero()
 */
class SymController;
class SymBDExp;

class SymAstar : public SymExploration  {  
    friend class SymAstarOpen;
    SymAstar * parent; //Parent of the search
    SymBDExp * bdExp;

  //Current state of the search:
  SymAstarOpen open_list;

  Bucket Sfilter;   //current g-bucket without duplicates and h-classified (still not filtered mutexes)
  Bucket Smerge;    // bucket before applying merge
  Bucket Szero;     // bucket to expand 0-cost transitions
  Bucket S;         // bucket to expand cost transitions

  //bucket to store temporary image results in expand_zero() and expand_cost()
  //For each BDD in Szero or S, stores a map with pairs <cost, resImage>
  std::vector<std::map<int, Bucket>> Simg;

  std::unique_ptr<SymAstarClosed> closed;    // Closed list 
  int f, g;            // f and g value of current bucket (S, Szero, Sfilter and Smerge)

  //acceptedValues: f, g pairs that must be checked because there are
  //some heuristics that may make some states to have them 
  std::map<int, std::set<int>> acceptedValues; 

  //To seek the next f,g bucket, we have two different sets: hValues.
  //hValues contains all the possible h-values returned by the heuristics. 
  //std::set<int> hValues;   //Possible h-values 
  SymAstarClosed * perfectHeuristic;                  //The perfect heuristic for this state space
  std::vector <std::shared_ptr<SymHeuristic> > heuristics;  //List of non-perfect heuristics
  std::set<int> hValuesExplicit;
  
  SymStepCostEstimation estimationCost, estimationZero;//Time/nodes estimated
  // NOTE: This was used to estimate the time and nodes needed to
  //perform a step in case that the next bucket is still not prepared.
  //Now, we always prepare the next bucket and when that fails no
  //estimation is needed (the exploration is deemed as not searchable
  //and is worse than any other exploration which has its next bucket
  //to expand ready)
  //SymStepCostEstimation estimationDisjCost, estimationDisjZero;
  bool lastStepCost; //If the last step was a cost step (to know if we are in estimationDisjCost or Zero
  SymController * engine; //Access to the bound and notification of new solutions

  SymExpStatistics stats;

  const bool storeBDDsToDisk; 

  bool bucketReady() const {
    /*cout << "bucket ready " << !(Szero.empty() && S.empty() && 
      Sfilter.empty() && Smerge.empty()) << std::endl;*/
    return !(Szero.empty() && S.empty() && Sfilter.empty() && Smerge.empty());
  }

  inline bool expansionReady() const {
    return Sfilter.empty() && Smerge.empty() && 
      !(Szero.empty() && S.empty());
  }
  virtual bool initialization() const{
    return g==0 && lastStepCost;
  }


  /*
   * Check generated or closed states with other frontiers.  In the
   * original state space we obtain a solution (maybe suboptimal if
   * the states are not closed). 
   */
  void checkCutOriginal(Bucket & bucket, int g);

  void closeStates(Bucket & bucket, int g);

  /*Get the next set from the open list and update g and f.
    Remove duplicate and spurious states. */
  void pop();

  bool prepareBucket(/*int maxTime, int maxNodes, bool afterPop*/);

  
  /* Apply 0-cost operators over Szero. */
  /* Puts the result on Szero. */
  /* Includes the result on S or open (depends on the heuristic value).*/  
  bool expand_zero(int maxTime, int maxNodes);
  
  /* Apply cost-operators over S. */
  /* Insert S on closed. */
  /* Insert successors on open.  */
  bool expand_cost(int maxTime, int maxNodes);

  // Returns the subset with h_value h
  BDD compute_heuristic(const BDD & from, int fVal, int hVal, bool store_eval); 

  void computeEstimation(bool prepare);

  //void debug_pop();
  
  //////////////////////////////////////////////////////////////////////////////
 public:
  SymAstar(SymController * eng, const SymParamsSearch & params);
  SymAstar(const SymAstar & ) = delete;
  SymAstar(SymAstar &&) = default;
  SymAstar& operator=(const SymAstar& ) = delete;
  SymAstar& operator=(SymAstar &&) = default;
  ~SymAstar() {}


  virtual bool finished() const {
      return open_list.empty() && !bucketReady(); 
  }

  const SymAstarOpen & getOpen() const {
    return open_list;
  }

  virtual bool stepImage(int maxTime, int maxNodes);


  bool init(SymBDExp * exp, SymManager * manager, bool fw); //Init forward or backward search

  //Initialize another search process by reutilizing information of this search
  //calls to 5 methods are needed.
  //1) init(), prepares the data of the other exploration.
  void init(SymBDExp * exp, SymAstar * other);
  //2) init2() reopens closed states in other frontier and initializes g, f
  //Should be called right after init is executed on both frontiers.
  void init2(SymAstar * opposite);
  //Then, relaxFrontier only relaxes the first bucket to expand. 
  //The caller should check if expansion is feasible and useful
  //Finally, all the open list is relaxed to the new abstract state space
  bool relaxFrontier(SymManager * manager, int maxTime, int maxNodes);
  bool relax_open(int maxTime, int maxNodes){
      return open_list.relax(maxTime, maxNodes);
  }
  void relaxClosed();

  void addHeuristic(std::shared_ptr<SymHeuristic> heuristic);
  void setPerfectHeuristic(SymAstarClosed * h);

  //Adds a new heuristic to evaluate States
  void setChild(SymAstar * child){
      closed->addChild(child->getClosed());
  }


  virtual void getHeuristic(std::vector<ADD> & heuristics,
			    std::vector <int> & maxHeuristicValues) const {
      closed->getHeuristic(heuristics, maxHeuristicValues);
  }

  // void getUsefulExplorations(set <SymAstar *> & explorations, double minRatioUseful);

  //double computeRatioUseful(SymHeuristic * h) const;

  //void notifyH(SymAstarClosed * heur, int value, bool isNotClosed);
  //void notifyF(SymAstarClosed * heur, int value);
  void notifyPrunedBy(int fVal, int gVal);
  void notify(const Bucket & bucket, int fNotClosed = 0); //May prune  
  void notifyNotClosed(int fValue, int hValue);

  void getPossiblyUsefulExplorations(std::vector <SymAstar *> & potentialExps){
    perfectHeuristic->getUsefulExps(potentialExps);
  }

  bool isBetter(const SymAstar & other) const;

  SymAstar * getOpposite() const{
    if(perfectHeuristic)
      return perfectHeuristic->getExploration();
    else 
      return nullptr;
  }

  virtual bool isSearchableWithNodes(int maxNodes) const;

  using SymExploration::isUseful;

  bool isUseful(const std::vector<BDD> & evalStates, 
		       const std::vector<BDD> & newFrontier, 
		       double ratio) const;

  virtual bool isUseful(double ratio) const {
    return !isAbstracted() || closed->isUseful(ratio);
  }

  /* inline bool isOtherUseful(Bucket & closedAbstract,  */
  /* 			    double ratio) const { */
  /*     assert (!isAbstracted()); // We should only call this method */
  /* 				//over original state space searches */
  /*     double rUseful = ratioUseful(closedAbstract); */
  /*     return rUseful > 0 && rUseful >= ratio  ; */
  /* } */

  double ratioUseful(Bucket & bucket) const;

  // Pointer to the closed list Used to set as heuristic of other explorations.
  inline SymAstarClosed * getClosed() const{
    return closed.get();
  }

  inline SymController * getEngine() const {
      return engine;
  }

  inline const SymManager * get_mgr() const{
    return mgr;
  }

  inline SymAstar * getParent() const{
    return parent;
  }

  inline Bucket getSfilter() const{
    return Sfilter;
  }

  inline Bucket getSmerge() const{
    return Smerge;
  }

  inline Bucket getSzero() const{
    return Szero;
  }

  inline Bucket getS() const{
    return S;
  }

  inline void getFrontier(Bucket & res) const {
      for(const BDD & bdd : Sfilter){
	  res.push_back(bdd);
      }
      for(const BDD & bdd : Smerge){
	  res.push_back(bdd);
      }
      for(const BDD & bdd : S){
	  res.push_back(bdd);
      }
      for(const BDD & bdd : Szero){
	  res.push_back(bdd);
      }
  }

  inline SymBDExp * getBDExp() const{
    return bdExp;
  }

  inline int getF() const{
    return f;
  }

  inline int getG() const{
    return g;
  }

  inline int getH() const{
    return f-g;
  }

  inline BDD getClosedTotal(){
    return closed->getClosed();
  }

  inline BDD notClosed(){
    return closed->notClosed();
  }

  void desactivate(){
    closed->desactivate();
  }

  void filterDuplicates(Bucket & bucket) {
      //For each BDD in the bucket, get states with f
      for(auto & bdd : bucket){
	  bdd *= closed->notClosed();
	  DEBUG_MSG(std::cout << ", duplicates: " << bdd.nodeCount(););
	  if(perfectHeuristic && 
	     perfectHeuristic->getFNotClosed() == std::numeric_limits<int>::max()){
	      bdd *= perfectHeuristic->getClosed();
	      DEBUG_MSG(std::cout << ", dead ends: " << bdd.nodeCount(););
	  }
      }
  }


  void filterHeuristic (Bucket & bucket, int fVal, int hVal, 
			Bucket & res, bool store_eval = true) {
      Timer time_h;
      for(int i = 0; i < bucket.size(); ++i){
	  BDD bddH = bucket[i];
	
	  //We left in bucket all the states that have been pruned
	  bucket[i] = compute_heuristic(bucket[i], fVal, hVal, store_eval);
	
	  if (!bucket[i].IsZero()) {
	      //bddH contains all the extracted states (those that fit fVal and hVal)
	      DEBUG_MSG(std::cout << "Pruning thanks to the heuristic: " << bddH.nodeCount(););
	      bddH -= bucket[i];
	      DEBUG_MSG(std::cout << " => " << bddH.nodeCount() << std::endl;);	
	  }
	
	  DEBUG_MSG(std::cout << ", h="<< hVal << ", extracted: " << bddH.nodeCount() 
		    << ", left: " << bucket[i].nodeCount() << std::endl;);
	  if(!bddH.IsZero()){
	      res.push_back(bddH);
	  }
      }
      removeZero(bucket);

      stats.time_heuristic_evaluation += time_h();
  }

  //Do not accept any larger value on this diagonal Only
  //applicable on the original state space because of the
  //usage of nipping.
  bool rejectLargerG (int f, int g) const {
      return !isAbstracted() && (perfectHeuristic->getHNotClosed() > f - g ||
				 perfectHeuristic->getFNotClosed() > f);
  }

  bool acceptFG(int f, int g) const {
      assert (f >= 0);
      assert (g >= 0);
      assert (f >= g);
      return  hValuesExplicit.count(f- g) ||
	  (acceptedValues.count(f) && acceptedValues.at(f).count(g)) ||
	  perfectHeuristic->accept(f, f - g);
  }

  std::pair<int, int> getAcceptedUpperBound() {
      std::pair<int, int> upper_bound {std::numeric_limits<int>::max(), 
	      open_list.minG()};
      const auto & candidates = acceptedValues.upper_bound(f);
      acceptedValues.erase(begin(acceptedValues), candidates);
      if(candidates != end(acceptedValues) && !candidates->second.empty()){
	  upper_bound = {candidates->first, *(candidates->second.begin())};
      }
      return upper_bound;
  }

  virtual long nextStepTime() const;
  virtual long nextStepNodes() const;
  virtual long nextStepNodesResult() const;

  //Returns the nodes that have been expanded by the algorithm (closed without the current frontier)
  BDD getExpanded() const;
  void getNotExpanded(Bucket & res) const;

  /* void write(const std::string & file) const; */
  /* void init(SymBDExp * exp, SymManager * manager,  const std::string & file); */

  void filterMutex (Bucket & bucket) {
      mgr->filterMutex(bucket, fw, initialization());
  }
 private: 
  
  virtual void printFrontier() const;

  int frontierNodes() const{
    if(!Szero.empty()){
      return nodeCount(Szero);
    }else if (!S.empty()){
      return nodeCount(S);
    }else{
      return nodeCount(Sfilter) + nodeCount(Smerge);
    }
  }

  double frontierStates() const{
    if(!Szero.empty()){
      return mgr->stateCount(Szero);
    }else if (!S.empty()){
      return mgr->stateCount(S);
    }else{
      return mgr->stateCount(Sfilter) + mgr->stateCount(Smerge);
    }
  }


  int frontierBuckets() const{
    if(!Szero.empty()){
      return Szero.size();
    }else if (!S.empty()){
      return S.size();
    }else{
      return Sfilter.size() + Smerge.size();
    }
  }

  void violated(TruncatedReason reason , double time, 
		int maxTime, int maxNodes);

  friend std::ostream & operator<<(std::ostream &os, const SymAstar & bdexp);
  
};
#endif

