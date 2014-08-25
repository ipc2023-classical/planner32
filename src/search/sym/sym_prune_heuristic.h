#ifndef SYM_PRUNE_HEURISTIC_H
#define SYM_PRUNE_HEURISTIC_H
#include "sym_variables.h"
#include "../option_parser.h"
#include <memory>
class SymManager;
class LDSimulation;
class SymTransition;

class SymPruneHeuristic {
  std::unique_ptr<LDSimulation> ldSimulation;
  std::unique_ptr<SymTransition> tr; //TR that computes dominated states
  
 public:
  virtual void initialize(SymManager * mgr);
  BDD simulatedBy(const BDD & bdd);

  SymPruneHeuristic(const Options &opts);
  ~SymPruneHeuristic();
};

#endif
