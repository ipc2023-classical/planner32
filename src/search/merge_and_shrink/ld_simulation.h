#ifndef MERGE_AND_SHRINK_LD_SIMULATION_H
#define MERGE_AND_SHRINK_LD_SIMULATION_H

#include "../operator_cost.h"

#include "label_relation.h"
#include "../sym/sym_variables.h"
#include "../option_parser.h"

#include <memory>
#include <vector>

class Labels;
class MergeStrategy;
class Abstraction;
class SimulationRelation;
class SymManager;
class LabelledTransitionSystem;

// Label dominance simulation
class LDSimulation {  
 protected:
  const bool use_expensive_statistics;
  const int limit_absstates_merge;
  const int limit_seconds_mas; //Limit of seconds for building the abstraction

  std::unique_ptr<MergeStrategy> merge_strategy;
  const bool use_bisimulation;
  const bool intermediate_simulations;

  //TODO: Use unique_ptr here
  std::unique_ptr<Labels> labels;
  std::vector<Abstraction *> abstractions;
  std::vector<SimulationRelation *> simulations;

  // std::unique_ptr<LabelRelation>  label_dominance;
  //  std::vector<LabelledTransitionSystem *> ltss;

  void dump_options() const;
  int num_equivalences() const;
  int num_simulations() const;

  void build_abstraction();
  void compute_ld_simulation();

  void compute_ld_simulation_after_merge(std::vector<Abstraction *> & _abstractions, 
					 std::vector<SimulationRelation *> & _simulations, 
					 const std::pair<int, int> & next_systems);

  static void compute_ld_simulation(Labels * _labels, std::vector<LabelledTransitionSystem *> & _ltss,
				    std::vector<SimulationRelation *> & _simulations);

  static void compute_ld_simulation(Labels * _labels, std::vector<Abstraction *> & _abstractions, 
				    std::vector<SimulationRelation *> & _simulations);


  static int prune_irrelevant_transitions(std::vector<LabelledTransitionSystem *> & _ltss, 
					   std::vector<SimulationRelation *> & _simulations, 
					   LabelRelation & label_dominance);
 public:
    LDSimulation(const Options &options);
    LDSimulation(bool unit_cost, const Options &options, OperatorCost cost_type);
    virtual ~LDSimulation();

    void initialize();
    void precompute_dominated_bdds(SymVariables * vars);
    void precompute_dominating_bdds(SymVariables * vars);

    bool pruned_state(const State &state) const;
    int get_cost(const State &state) const;

    BDD getSimulatedBDD(SymVariables * vars, const State &state) const;
    BDD getSimulatingBDD(SymVariables * vars, const State &state) const;

    inline std::vector<SimulationRelation *> get_simulations(){
      return simulations;
    }
    
    //Computes the probability of selecting a random pair s, s' such
    //that s simulates s'.
    double get_percentage_simulations(bool ignore_equivalences) const;

    //Computes the probability of selecting a random pair s, s' such that
    //s is equivalent to s'.
    double get_percentage_equivalences() const;

    double get_percentage_equal() const;


    static void add_options_to_parser(OptionParser &parser);
};


#endif
