#ifndef SYM_PRUNE_HEURISTIC_H
#define SYM_PRUNE_HEURISTIC_H
#include "sym_variables.h"
#include "../option_parser.h"
#include <memory>
class SymManager;
class LDSimulation;
class SymTransition;

class SymPruneHeuristic {
  const bool prune_irrelevant; 
  const bool dominance_pruning; 
    
  std::unique_ptr<LDSimulation> ldSimulation;
  std::unique_ptr<SymTransition> tr; //TR that computes dominated states
  
 public:
  virtual void initialize(SymManager * mgr);
  BDD simulatedBy(const BDD & bdd);

  SymPruneHeuristic(const Options &opts);

  ~SymPruneHeuristic();

  bool use_dominance_pruning() const{
      return dominance_pruning;
  }
};

#endif
