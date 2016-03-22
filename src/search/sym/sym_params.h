#ifndef SYM_PARAMS_H
#define SYM_PARAMS_H

#include "sym_enums.h" 

//TODO: move variable_order_finder somewhere else
#include "../merge_and_shrink/variable_order_finder.h" 
class OptionParser;
class Options;

//Parameters are divided in two classes. SymParamsMgr for mgr params
//and SymParamsSearch for parameters that control the characteristics
//of a search. 

class SymParamsMgr{

 public:
    //TODO: Rename merge strategy type to VarOrder

    // Var order used by the algorithm. By default: MERGE_LINEAR_LEVEL
    // so that we follow the preprocessor order (which should be
    // optimized)
    VariableOrderType variable_ordering; 
  //Parameters to initialize the CUDD manager
  long cudd_init_nodes; //Number of initial nodes
  long cudd_init_cache_size; //Initial cache size
  long cudd_init_available_memory; //Maximum available memory (bytes)

  //Parameters to generate the TRs
  int max_tr_size;
  int max_tr_time;

  //Parameters to generate the mutex BDDs
  MutexType mutex_type;
  int max_mutex_size;
  int max_mutex_time; 

    SymParamsMgr();
  SymParamsMgr(const Options & opts);
  static void add_options_to_parser(OptionParser &parser);
  static void add_options_to_parser_simulation(OptionParser &parser);
  void print_options() const ;
};

class SymParamsSearch{
 public:
  //By default max<int>. If lower, it allows for skip disjunction if the size of S is greater.
  int max_disj_nodes;

  //Parameters for sym_step_cost_estimation
  double min_estimation_time;         // Dont estimate if time is lower
  double penalty_time_estimation_sum; // violated_time = sum + time*mult
  double penalty_time_estimation_mult;
  double penalty_nodes_estimation_sum;// violated_nodes = sum + nodes*mult
  double penalty_nodes_estimation_mult;

  //Time and memory bounds for auxiliary operations
  int max_pop_nodes;
  int max_pop_time;

  //Parameters to control isUseful() and isSearchable()
  int maxStepTime, maxStepNodes;
  double ratioUseful;

  // Parameters to decide the alloted time for a step
  // alloted = max(minAlloted, estimated*multAlloted)
  int minAllotedTime, minAllotedNodes;   // min alloted time and nodes to a step
  int maxAllotedTime, maxAllotedNodes;   // min alloted time and nodes to a step
  double ratioAllotedTime, ratioAllotedNodes; // factor to multiply the estimation

  double ratioAfterRelax;

    bool non_stop;
  SymParamsSearch(const Options & opts);

  static void add_options_to_parser(OptionParser &parser, int maxStepTime, int maxStepNodes);
  //Parameters with default values for hierarchy policies
  static void add_options_to_parser_abstractions(OptionParser &parser, int maxStepTime, int maxStepNodes);
  void print_options() const ;  

  inline double getAllotedTime(double estimatedTime) const{
    return std::min(maxAllotedTime,
		    std::max<int>(estimatedTime*ratioAllotedTime, minAllotedTime));
  }

  inline double getAllotedNodes(double estimatedNodes) const{
    return std::min(maxAllotedNodes, std::max<int>(estimatedNodes*ratioAllotedNodes, minAllotedNodes));
  }

  void inheritParentParams(const SymParamsSearch & other){
    maxStepTime = std::min(maxStepTime, other.maxStepTime);
    maxStepNodes = std::min(maxStepNodes, other.maxStepNodes);
    }

    inline bool get_non_stop () const {
	return non_stop;
  }
};

#endif
