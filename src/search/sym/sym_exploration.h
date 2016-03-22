#ifndef SYM_EXPLORATION_H
#define SYM_EXPLORATION_H

#include "../debug.h" 
#include "sym_heuristic.h"
#include "sym_manager.h"
#include "sym_bucket.h"
#include "sym_open.h"
#include "sym_closed.h"
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
    friend class SymOpen;

  //Attributes that characterize the search:
  SymBDExp * bdExp;
  SymManager * mgr;            //Symbolic manager to perform bdd operations
  SymParamsSearch p;
  bool fw; //Direction of the search. true=forward, false=backward
  SymExploration * parent; //Parent of the search 

  //Current state of the search:
  SymOpen open_list;

  Bucket Sfilter;   //current g-bucket without duplicates and h-classified (still not filtered mutexes)
  Bucket Smerge;    // bucket before applying merge
  Bucket Szero;     // bucket to expand 0-cost transitions
  Bucket S;         // bucket to expand cost transitions

  //bucket to store temporary image results in expand_zero() and expand_cost()
  //For each BDD in Szero or S, stores a map with pairs <cost, resImage>
  std::vector<std::map<int, Bucket>> Simg;

  std::unique_ptr<SymClosed> closed;    // Closed list 
  int f, g;            // f and g value of current bucket (S, Szero, Sfilter and Smerge)

  //acceptedValues: f, g pairs that must be checked because there are
  //some heuristics that may make some states to have them 
  std::map<int, std::set<int>> acceptedValues; 

  //To seek the next f,g bucket, we have two different sets: hValues.
  //hValues contains all the possible h-values returned by the heuristics. 
  //std::set<int> hValues;   //Possible h-values 
  SymClosed * perfectHeuristic;                  //The perfect heuristic for this state space
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

  bool bucketReady() const {
    /*cout << "bucket ready " << !(Szero.empty() && S.empty() && 
      Sfilter.empty() && Smerge.empty()) << std::endl;*/
    return !(Szero.empty() && S.empty() && Sfilter.empty() && Smerge.empty());
  }

  inline bool expansionReady() const {
    return Sfilter.empty() && Smerge.empty() && 
      !(Szero.empty() && S.empty());
  }
  inline bool initialization() const{
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
  SymExploration(SymController * eng, const SymParamsSearch & params);
  SymExploration(const SymExploration & ) = delete;
  SymExploration(SymExploration &&) = default;
  SymExploration& operator=(const SymExploration& ) = delete;
  SymExploration& operator=(SymExploration &&) = default;
  ~SymExploration() {}


  inline bool finished() const {
      return open_list.empty() && !bucketReady(); 
  }

  const SymOpen & getOpen() const {
    return open_list;
  }

  bool stepImage(){
    return stepImage(p.getAllotedTime(nextStepTime()), 
		     p.getAllotedNodes(nextStepNodesResult()));
  }
  bool stepImage(int maxTime, int maxNodes);


  bool init(SymBDExp * exp, SymManager * manager, bool fw); //Init forward or backward search

  //Initialize another search process by reutilizing information of this search
  //calls to 5 methods are needed.
  //1) init(), prepares the data of the other exploration.
  void init(SymBDExp * exp, SymExploration * other);
  //2) init2() reopens closed states in other frontier and initializes g, f
  //Should be called right after init is executed on both frontiers.
  void init2(SymExploration * opposite);
  //Then, relaxFrontier only relaxes the first bucket to expand. 
  //The caller should check if expansion is feasible and useful
  //Finally, all the open list is relaxed to the new abstract state space
  bool relaxFrontier(SymManager * manager, int maxTime, int maxNodes);
  bool relax_open(int maxTime, int maxNodes){
      return open_list.relax(maxTime, maxNodes);
  }
  void relaxClosed();

  void addHeuristic(std::shared_ptr<SymHeuristic> heuristic);
  void setPerfectHeuristic(SymClosed * h);

  //Adds a new heuristic to evaluate States
  void setChild(SymExploration * child){
      closed->addChild(child->getClosed());
  }

  // void getUsefulExplorations(set <SymExploration *> & explorations, double minRatioUseful);

  //double computeRatioUseful(SymHeuristic * h) const;

  //void notifyH(SymClosed * heur, int value, bool isNotClosed);
  //void notifyF(SymClosed * heur, int value);
  void notifyPrunedBy(int fVal, int gVal);
  void notify(const Bucket & bucket, int fNotClosed = 0); //May prune  
  void notifyNotClosed(int fValue, int hValue);

  void getPossiblyUsefulExplorations(std::vector <SymExploration *> & potentialExps){
    perfectHeuristic->getUsefulExps(potentialExps);
  }

  bool isBetter(const SymExploration & other) const;

  SymExploration * getOpposite() const{
    if(perfectHeuristic)
      return perfectHeuristic->getExploration();
    else 
      return nullptr;
  }

  inline bool isSearchable() const{
      return isSearchableWithNodes(p.maxStepNodes);
  }

  inline bool isSearchableAfterRelax(int num_relaxations) const{
    double maxNodes = p.maxStepNodes;
    if(num_relaxations){
      maxNodes *= pow(p.ratioAfterRelax, num_relaxations);
    }
    return isSearchableWithNodes((int)maxNodes);
  }

  bool isSearchableWithNodes(int maxNodes) const;


  bool isUseful(const std::vector<BDD> & evalStates, 
		       const std::vector<BDD> & newFrontier, 
		       double ratio) const;

  inline bool isUseful() const {
    return isUseful(p.ratioUseful);
  }

  inline bool isUseful(double ratio) const {
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
  inline SymClosed * getClosed() const{
    return closed.get();
  }

  inline const SymManager * get_mgr() const{
    return mgr;
  }

  inline SymExploration * getParent() const{
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

  inline int getF() const{
    return f;
  }

  inline int getG() const{
    return g;
  }

  inline int getH() const{
    return f-g;
  }

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

  inline SymBDExp * getBDExp() const{
    return bdExp;
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

  void filterMutex (Bucket & bucket) {
      mgr->filterMutexBucket(bucket, fw, initialization(), 
			     p.max_pop_time, p.max_pop_nodes);
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

  void mergeBucket(Bucket & bucket) {
      mergeBucket(bucket, p.max_pop_time, p.max_pop_nodes);
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

  long nextStepTime() const;
  long nextStepNodes() const;
  long nextStepNodesResult() const;

  //Returns the nodes that have been expanded by the algorithm (closed without the current frontier)
  BDD getExpanded() const;
  void getNotExpanded(Bucket & res) const;

  void write(const std::string & file) const;
  void init(SymBDExp * exp, SymManager * manager,  const std::string & file);


  inline SymController * getEngine() const{
    return engine;
  }


  void statistics() const;
  
 private: 

  double stateCount(const Bucket & bucket) const {
      double sum = 0;
      for(const BDD & bdd : bucket){
	  sum += mgr->getVars()->numStates(bdd);
      }
      return sum;
  }


  void shrinkBucket(Bucket & bucket, int maxNodes);

  //  void addCountStates(SymHeuristic * h, const Bucket & bucket, double & possible, double & total) const;    

  //BDD closedByParents(const BDD & bdd, SymHeuristic * heur) const;
  //BDD notClosedByParents(const BDD & bdd, SymHeuristic * heur) const;

  //Extract states with h-value from list, removing duplicates
  /* void extract_states(Bucket & list, int fVal, int hVal, */
  /* 		      Bucket & res, bool duplicates);  */

  /* bool extract_states(Bucket & list, const Bucket & pruned, Bucket & res) const;  */


  bool mergeBucket(Bucket & bucket, int maxTime, int maxNodes){
      auto mergeBDDs = [] (BDD bdd, BDD bdd2, int maxNodes){
	  return bdd.Or(bdd2, maxNodes);
      };
      merge(mgr->getVars(), bucket, mergeBDDs, maxTime, std::min(maxNodes, p.max_disj_nodes));
      removeZero(bucket); //Be sure that we do not contain only the zero BDD
    
      return maxNodes >= p.max_disj_nodes || bucket.size() <= 1;
  }


  //Extract states without h-value from list using sym heuristic 
  /* void extract_states(Bucket & list, int fVal, int hVal, */
  /* 		      Bucket & res, SymClosed * heur);  */



  /* bool stepFDiagonal(int maxTime, int maxNodes); */
  /* bool explore(int maxTime, int maxNodes, */
  /* 	       int maxStepTime, BDD target = BDD(), */
  /* 	       int fTarget = 0, int gTarget = 0); */

  /* bool stepExpand(int maxTime, int maxNodes); */
  
  void printFrontier() const;

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
      return stateCount(Szero);
    }else if (!S.empty()){
      return stateCount(S);
    }else{
      return stateCount(Sfilter) + stateCount(Smerge);
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

  friend std::ostream & operator<<(std::ostream &os, const SymExploration & bdexp);
  
};
#endif // SYMBOLIC_EXPLORATION

