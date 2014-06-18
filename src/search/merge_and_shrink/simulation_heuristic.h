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
class SymTransition;

enum class PruningType {BDD_DOMINATED, ADD_DOMINATED, BDD_DOMINATED_BLIND, 
    BDD_CLOSED, ADD_CLOSED, BDD_CLOSED_BLIND};
std::ostream & operator<<(std::ostream &os, const PruningType & m);
extern const std::vector<std::string> PruningTypeValues;

class SimulationHeuristic : public PruneHeuristic {  

  MergeStrategy *const merge_strategy;
  const bool use_expensive_statistics;

 protected:
  std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  

  bool insert_dominated;
  int limit_absstates_merge;
  SymParamsMgr mgrParams; //Parameters for SymManager configuration.
  Labels *labels;

  std::unique_ptr<SymTransition> tr; //TR that computes dominated states

  void dump_options() const;

  int num_equivalences() const;
  int num_simulations() const;

  std::vector<Abstraction *> abstractions;
  std::vector<SimulationRelation *> simulations;

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

    SimulationHeuristic(const Options &opts, bool _insert_dominated);
    virtual ~SimulationHeuristic();
};

class SimulationHeuristicBDDMap : public SimulationHeuristic {
  std::map<int, BDD> closed;
 public:
  SimulationHeuristicBDDMap (const Options &opts, 
			     bool _insert_dominated) : 
  SimulationHeuristic(opts, _insert_dominated)
  {}
  virtual ~SimulationHeuristicBDDMap (){}

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g);
  virtual void insert (const State & state, int g);
};

class SimulationHeuristicBDD : public SimulationHeuristic {
  BDD closed;
 public:
  SimulationHeuristicBDD (const Options &opts, 
			  bool _insert_dominated) : 
  SimulationHeuristic(opts, _insert_dominated), closed(vars->zeroBDD())
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
