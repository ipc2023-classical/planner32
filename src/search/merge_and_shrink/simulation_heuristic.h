#ifndef MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H
#define MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H

#include "../prune_heuristic.h"
#include "../sym/sym_variables.h"
#include "../sym/sym_params.h"

class MergeStrategy;
class Abstraction;
class Labels;
class SimulationRelation;
class SymVariables;
class SymManager;
class SymTransition;

enum class PruningType {BDD_MAP, ADD, BDD};
std::ostream & operator<<(std::ostream &os, const PruningType & m);
extern const std::vector<std::string> PruningTypeValues;

class SimulationHeuristic : public PruneHeuristic {  
 protected:
  const bool use_expensive_statistics;
  
  //Parameters to control the pruning
  const SymParamsMgr mgrParams; //Parameters for SymManager configuration.
  const bool remove_spurious_dominated_states;
  const bool insert_dominated;

  //Parameters to control the simulation
  const int limit_absstates_merge;
  MergeStrategy *const merge_strategy;
  const bool use_bisimulation;

  std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  
  std::unique_ptr<SymManager> mgr;    //The symbolic manager to handle mutex BDDs

  //TODO: Should labels be deleted after the simulation algorithm? If
  //yes, declare there. If no, consider using unique_ptr
  Labels *labels;

  std::unique_ptr<SymTransition> tr; //TR that computes dominated states
  int states_inserted; //Count the number of states inserted.

  void dump_options() const;

  int num_equivalences() const;
  int num_simulations() const;

  std::vector<Abstraction *> abstractions;
  std::vector<SimulationRelation *> simulations;


  /* Methods to help concrete classes */
  BDD getBDDToInsert(const State &state);

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g) = 0;
  virtual void insert (const State & state, int g) = 0;

  /* //Methods to keep dominated states in symbolic search */
  /* virtual BDD check (const BDD & bdd, int g) = 0; */
  /* virtual void insert (const BDD & bdd, int g) = 0; */

  void build_abstraction();
  virtual int compute_heuristic(const State &state);
 public:
    virtual void initialize(bool explicit_search);

    //Methods for pruning explicit search
    virtual bool prune_generation(const State &state, int g);
    virtual bool prune_expansion (const State &state, int g);

    virtual SymTransition * getTR(SymVariables * vars);

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
