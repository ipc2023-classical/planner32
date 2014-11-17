#ifndef MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H
#define MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H

#include "../prune_heuristic.h"
#include "../sym/sym_variables.h"
#include "../sym/sym_params.h"

class LDSimulation;
class SymVariables;
class SymManager;

enum class PruningDD {BDD_MAP, ADD, BDD, BDD_MAP_DISJ};
std::ostream & operator<<(std::ostream &os, const PruningDD & m);
extern const std::vector<std::string> PruningDDValues;

enum class PruningType {Expansion, Generation};
std::ostream & operator<<(std::ostream &os, const PruningType & m);
extern const std::vector<std::string> PruningTypeValues;

class SimulationHeuristic : public PruneHeuristic {  
 protected:
  //Parameters to control the pruning
  const SymParamsMgr mgrParams; //Parameters for SymManager configuration.

  bool initialized;
  const bool remove_spurious_dominated_states;
  const bool insert_dominated;
  const PruningType pruning_type;
  bool print_desactivation;

  /*
   * Three parameters help to decide whether to apply dominance
   * pruning or not. Dominance pruning is used until min_expansions
   * are performed. Afterwards, pruning is only used if the ratio
   * pruned/inserted is greater than min_pruning_ratio. If the ratio
   * pruned/inserted is lower than min_insert_ratio, we use the states
   * already inserted to prune but avoid inserting more states.
   */
  const double min_pruning_ratio, min_insert_ratio, min_deadend_ratio;
  const int min_insertions, min_deadends;

  std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  
  std::unique_ptr<SymManager> mgr;    //The symbolic manager to handle mutex BDDs

  std::unique_ptr<LDSimulation> ldSimulation;
  int states_inserted; //Count the number of states inserted
  int states_pruned; //Count the number of states pruned
  int deadends_pruned; //Count the number of dead ends detected

  void dump_options() const;

  /* Methods to help concrete classes */
  BDD getBDDToInsert(const State &state);

  //Methods to keep dominated states in explicit search
  //Check: returns true if a better or equal state is known
  virtual bool check (const State & state, int g) = 0;
  virtual void insert (const State & state, int g) = 0;

  /* //Methods to keep dominated states in symbolic search */
  /* virtual BDD check (const BDD & bdd, int g) = 0; */
  /* virtual void insert (const BDD & bdd, int g) = 0; */

  //void build_abstraction();

  inline bool insert_is_activated() const {
      return states_inserted < min_insertions || 
	  states_pruned >= states_inserted*min_insert_ratio;
  }
  inline bool prune_is_activated() const {
      return states_inserted < min_insertions || 
	  states_pruned >= states_inserted*min_pruning_ratio;
  }

  inline bool deadend_is_activated() const {
      return states_inserted < min_deadends || 
	  deadends_pruned >= states_inserted*min_deadend_ratio
	  || prune_is_activated();
  }

 public:
  virtual void initialize();

  //Methods for pruning explicit search
  virtual bool prune_generation(const State &state, int g);
  virtual bool prune_expansion (const State &state, int g);

  virtual bool is_dead_end(const State &state);

  virtual int compute_heuristic(const State &state);
  SimulationHeuristic(const Options &opts);
  virtual ~SimulationHeuristic();
};

class SimulationHeuristicBDDMap : public SimulationHeuristic {
  std::map<int, BDD> closed;
 public:
  SimulationHeuristicBDDMap (const Options &opts) : 
  SimulationHeuristic(opts)
  {}
  virtual ~SimulationHeuristicBDDMap (){}

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g);
  virtual void insert (const State & state, int g);
};

class SimulationHeuristicBDDMapDisj : public SimulationHeuristic {
    std::map<int, std::vector<BDD> > closed;
 public:
  SimulationHeuristicBDDMapDisj (const Options &opts) : 
  SimulationHeuristic(opts)
  {}
  virtual ~SimulationHeuristicBDDMapDisj (){}

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g);
  virtual void insert (const State & state, int g);
};


class SimulationHeuristicBDD : public SimulationHeuristic {
  BDD closed, closed_inserted;
  bool initialized;
 public:
  SimulationHeuristicBDD (const Options &opts) : 
  SimulationHeuristic(opts), initialized(false)
  {}
  virtual ~SimulationHeuristicBDD (){}

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g);
  virtual void insert (const State & state, int g);
};

/* class SimulationHeuristicADD : public SimulationHeuristic { */
/*   ADD closed; */
/*  public: */
/*   SimulationHeuristicADD (const Options &opts,  */
/* 			  bool _insert_dominated) :  */
/*   SimulationHeuristic(opts, _insert_dominated) */
/*   {} */
/*   virtual ~SimulationHeuristicADD (){} */

/*     virtual bool check (const State & state, int g); */
/*     virtual void insert (const State & state, int g); */
/* }; */

#endif
