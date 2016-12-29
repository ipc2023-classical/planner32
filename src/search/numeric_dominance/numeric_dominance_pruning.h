#ifndef NUMERIC_DOMINANCE_NUMERIC_DOMINANCE_PRUNING_H
#define NUMERIC_DOMINANCE_NUMERIC_DOMINANCE_PRUNING_H

#include "../prune_heuristic.h"
#include "../sym/sym_variables.h"
#include "../sym/sym_params.h"

#include "numeric_dominance_relation.h"

class LDSimulation;
class AbstractionBuilder;
class SymVariables;
class SymManager;
class Abstraction;

enum class PruningType {Expansion, Generation, Parent, None};
std::ostream & operator<<(std::ostream &os, const PruningType & m);
extern const std::vector<std::string> PruningTypeValues;

class NumericDominancePruning : public PruneHeuristic {  
 protected:
  //Parameters to control the pruning
  const SymParamsMgr mgrParams; //Parameters for SymManager configuration.

  bool initialized;
  const bool remove_spurious_dominated_states;
  const bool insert_dominated;
  const PruningType pruning_type;

  /*
   * Three parameters help to decide whether to apply dominance
   * pruning or not. Dominance pruning is used until
   * min_insertions_desactivation are performed. At that moment, if
   * the ratio pruned/checked is lower than min_desactivation_ratio
   * the pruning is desactivated. If not, the pruning remains
   * activated until the planner finishes.
   */
  const int min_insertions_desactivation;
  const double min_desactivation_ratio;
  
  std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  
  std::unique_ptr<SymManager> mgr;    //The symbolic manager to handle mutex BDDs

  std::unique_ptr<AbstractionBuilder> abstractionBuilder;
  std::unique_ptr<LDSimulation> ldSimulation;
  std::unique_ptr<NumericDominanceRelation> numeric_dominance_relation;
  std::vector<std::unique_ptr<Abstraction> > abstractions;

  bool all_desactivated;
  bool activation_checked;

  int states_inserted; //Count the number of states inserted
  int states_checked; //Count the number of states inserted
  int states_pruned; //Count the number of states pruned
  int deadends_pruned; //Count the number of dead ends detected

  void dump_options() const;

  /* Methods to help concrete classes */
  BDD getBDDToInsert(const State &state);

  //Methods to keep dominated states in explicit search
  //Check: returns true if a better or equal state is known
  virtual bool check (const State & state, int g) = 0;
  virtual void insert (const State & state, int g) = 0;

  inline bool is_activated() {
      if(!activation_checked && states_inserted > min_insertions_desactivation){
	  activation_checked = true;
	  all_desactivated = states_pruned == 0 || 
	      states_pruned < states_checked*min_desactivation_ratio;
	  std::cout << "Simulation pruning " << (all_desactivated ? "desactivated: " : "activated: ")
		    << states_pruned << " pruned " << states_checked << " checked " << 
	      states_inserted << " inserted " << deadends_pruned << " deadends " << std::endl;
      }
      return !all_desactivated;
  }

 public:
  virtual void initialize();

  //Methods for pruning explicit search
  virtual bool prune_generation(const State &state, int g, const State &parent, int action_cost);
  virtual bool prune_expansion (const State &state, int g);

  virtual bool is_dead_end(const State &state);

  virtual int compute_heuristic(const State &state);
  NumericDominancePruning(const Options &opts);
  virtual ~NumericDominancePruning();

  virtual void print_statistics();

  virtual bool proves_task_unsolvable() const {
      return true;
  }
};

class NumericDominancePruningBDDMap : public NumericDominancePruning {
    std::map<int, BDD> closed;
public:
    NumericDominancePruningBDDMap (const Options &opts) : 
    NumericDominancePruning(opts)
    {}
    virtual ~NumericDominancePruningBDDMap () = default;

    //Methods to keep dominated states in explicit search
    virtual bool check (const State & state, int g);
    virtual void insert (const State & state, int g);
};


#endif
