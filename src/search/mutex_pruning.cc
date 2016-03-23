#include "mutex_pruning.h"
#include "merge_and_shrink/opt_order.h"
#include "globals.h"
#include "sym/sym_variables.h"
#include "sym/sym_manager.h"
#include "causal_graph.h"
#include "option_parser.h"
#include "plugin.h"

using namespace std;

MutexPruning::MutexPruning(const Options &opts)
    : PruneHeuristic(opts), mgrParams(opts){
}

void MutexPruning::initialize() {
    if(!vars){
	vars = unique_ptr<SymVariables> (new SymVariables());
        vector <int> var_order;
	InfluenceGraph::compute_gamer_ordering(var_order);
        vars->init(var_order, mgrParams);
	//mgr = unique_ptr<SymManager>(new SymManager (vars.get(), nullptr, mgrParams));
	SymManager mgr  (vars.get(), nullptr, mgrParams);
	mutex_bdds = mgr.getNotMutexBDDs(true);
    }
}

static PruneHeuristic *_parse(OptionParser &parser) {
    parser.document_synopsis("Mutex pruning", "");

    Heuristic::add_options_to_parser(parser);
    SymParamsMgr::add_options_to_parser_simulation(parser);

    Options opts = parser.parse();

    if (parser.dry_run()) {
        return nullptr;
    }
    
    return new MutexPruning(opts);
}


static Plugin<PruneHeuristic> _plugin("mutex", _parse);

bool MutexPruning::is_dead_end(const State &state){
    for (BDD & bdd : mutex_bdds) {
	if(vars->isIn(state, bdd)) {
	    return true;
	}
    }
    return false;
}


void MutexPruning::print_statistics()
 {
 }
