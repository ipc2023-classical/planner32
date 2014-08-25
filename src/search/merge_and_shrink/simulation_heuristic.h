#ifndef MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H
#define MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H

#include "../prune_heuristic.h"
#include "../sym/sym_variables.h"
#include "../sym/sym_params.h"

class LDSimulation;
class SymVariables;
class SymManager;

enum class PruningDD {BDD_MAP, ADD, BDD};
std::ostream & operator<<(std::ostream &os, const PruningDD & m);
extern const std::vector<std::string> PruningDDValues;

enum class PruningType {Expansion, Generation};
std::ostream & operator<<(std::ostream &os, const PruningType & m);
extern const std::vector<std::string> PruningTypeValues;

class SimulationHeuristic : public PruneHeuristic {  
 protected:
  //Parameters to control the pruning
  const SymParamsMgr mgrParams; //Parameters for SymManager configuration.
  const bool remove_spurious_dominated_states;
  const bool insert_dominated;
  const PruningType pruning_type;

  std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  
  std::unique_ptr<SymManager> mgr;    //The symbolic manager to handle mutex BDDs

  std::unique_ptr<LDSimulation> ldSimulation;
  int states_inserted; //Count the number of states inserted.

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
  virtual int compute_heuristic(const State &state);
 public:
    virtual void initialize();

    //Methods for pruning explicit search
    virtual bool prune_generation(const State &state, int g);
    virtual bool prune_expansion (const State &state, int g);

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

class SimulationHeuristicBDD : public SimulationHeuristic {
  BDD closed, closed_inserted;
 public:
  SimulationHeuristicBDD (const Options &opts) : 
  SimulationHeuristic(opts), closed(vars->zeroBDD()), closed_inserted(vars->zeroBDD())
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
