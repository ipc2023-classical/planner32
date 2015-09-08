#include "sym_prune_heuristic.h"

#include "../merge_and_shrink/ld_simulation.h"
#include "sym_transition.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "sym_manager.h"

using namespace std;

SymPruneHeuristic::SymPruneHeuristic(const Options &opts) :
    prune_irrelevant(opts.get<bool>("prune_irrelevant")), 
    dominance_pruning(opts.get<bool>("dominance_pruning")), 
    ldSimulation(new LDSimulation(opts)){
    
    ldSimulation->initialize();
}

SymPruneHeuristic::~SymPruneHeuristic(){}

void SymPruneHeuristic::initialize(SymManager * mgr) {
    cout << "Initialize sym prune heuristic" << endl;
  if(!tr && dominance_pruning){
    tr = unique_ptr<SymTransition>(new SymTransition(mgr,
						     ldSimulation->get_simulations()));
  }

  if(prune_irrelevant){
      cout << "Computing irrelevant states BDD " << g_timer() << endl;
      BDD irrelevantStates = ldSimulation->get_simulations().getIrrelevantStates(mgr->getVars());
      cout << "Irrelevant states BDD: " << irrelevantStates.nodeCount() << " " << g_timer() << endl;
      //Prune irrelevant states in both directions
      mgr->addDeadEndStates(true, irrelevantStates);
      mgr->addDeadEndStates(false, irrelevantStates);
  }
}

BDD SymPruneHeuristic::simulatedBy(const BDD & bdd) {
    if(dominance_pruning)
	return tr->image(bdd);
    else
	return bdd;
}

static SymPruneHeuristic *_parse(OptionParser &parser) {
    parser.document_synopsis("Simulation prune heuristic", "");
    
    parser.add_option<bool>("prune_irrelevant",
            "Activate removing irrelevant states from the search",
            "false");

    parser.add_option<bool>("dominance_pruning",
            "Activate dominance pruning",
            "false");

    LDSimulation::add_options_to_parser(parser);

    Options opts = parser.parse();

    if (parser.dry_run()) {
        return 0;
    } else {
        return new SymPruneHeuristic (opts);
    }
}

static Plugin <SymPruneHeuristic> _plugin("simulation", _parse);
