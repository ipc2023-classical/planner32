#ifndef MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H
#define MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H

#include "../heuristic.h"
#include "../sym/sym_variables.h"
#include "../sym/sym_params.h"
class Abstraction;
class Labels;
class SimulationRelation;
class SymVariables;

class SimulationHeuristic : public Heuristic {
  std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  
  SymParamsMgr mgrParams; //Parameters for SymManager configuration.
  Labels *labels;
  BDD closed;
  //std::map<int, BDD> closed;
  std::vector<Abstraction *> abstractions;
  std::vector<SimulationRelation *> simulations;    
  
  void dump_options() const;
 protected:
    virtual void initialize();
    virtual int compute_heuristic(const State &state);
 public:
    SimulationHeuristic(const Options &opts);
    ~SimulationHeuristic();
};

#endif
