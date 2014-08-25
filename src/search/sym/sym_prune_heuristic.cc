#include "sym_prune_heuristic.h"

#include "../merge_and_shrink/ld_simulation.h"
#include "sym_transition.h"

using namespace std;

SymPruneHeuristic::SymPruneHeuristic(const Options &opts) : 
  ldSimulation(new LDSimulation(opts)){
  
}

SymPruneHeuristic::~SymPruneHeuristic(){}

void SymPruneHeuristic::initialize(SymManager * mgr) {
  if(!tr){
    tr = unique_ptr<SymTransition>(new SymTransition(mgr,
						     ldSimulation->get_simulations()));
  }
}

BDD SymPruneHeuristic::simulatedBy(const BDD & bdd) {
  return tr->image(bdd);
}

