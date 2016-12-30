#include "numeric_dominance_pruning.h"

#include "../merge_and_shrink/abstraction.h"
#include "../merge_and_shrink/labels.h"
#include "../merge_and_shrink/labelled_transition_system.h"
#include "../merge_and_shrink/abstraction_builder.h"

#include "numeric_simulation_relation.h"

#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../state.h"
#include "../timer.h"

#include "../sym/sym_transition.h"
#include "../sym/sym_manager.h"

#include <cassert>
#include <vector>
#include <limits>

using namespace std;

NumericDominancePruning::NumericDominancePruning(const Options &opts)
: PruneHeuristic(opts), 
  mgrParams(opts), initialized(false),
  remove_spurious_dominated_states(opts.get<bool>("remove_spurious")),
  insert_dominated(opts.get<bool>("insert_dominated")),
  pruning_type(NumericPruningType(opts.get_enum("pruning_type"))),
  min_insertions_desactivation(opts.get<int>("min_insertions")),
  min_desactivation_ratio(opts.get<double>("min_desactivation_ratio")),
  vars(new SymVariables()), abstractionBuilder(opts.get<AbstractionBuilder *>("abs")),
  all_desactivated(false), activation_checked(false), 
  states_inserted(0), states_checked(0), states_pruned(0), deadends_pruned(0) {
}

NumericDominancePruning::~NumericDominancePruning() {
}

void NumericDominancePruning::dump_options() const {
    cout << "Type pruning: " << pruning_type << endl;
}

void NumericDominancePruning::initialize() {
    if(!initialized){
	dump_options();
        initialized = true;
	abstractionBuilder->build_abstraction(is_unit_cost_problem() ||
					      cost_type == OperatorCost::ZERO,
					      cost_type, ldSimulation, abstractions);
	cout << "LDSimulation finished" << endl;

	if(pruning_type != NumericPruningType::None) {
	    numeric_dominance_relation = ldSimulation->compute_numeric_dominance_relation();
	}

	ldSimulation->release_memory();

	// if(pruning_type != NumericPruningType::None && pruning_type != NumericPruningType::Parent) {
	//     vector <int> var_order;
	//     ldSimulation->getVariableOrdering(var_order);

	//     vars->init(var_order, mgrParams);
	//     if(remove_spurious_dominated_states){
	// 	mgr = unique_ptr<SymManager> (new SymManager(vars.get(), nullptr, mgrParams, cost_type));
	// 	mgr->init();
	//     }
	//     mgrParams.print_options();

	//     if(insert_dominated){
	// 	ldSimulation->get_dominance_relation().precompute_dominated_bdds(vars.get());
		
	//     }else{
	// 	ldSimulation->get_dominance_relation().precompute_dominating_bdds(vars.get());
	//     }
	// }
        cout << "Completed preprocessing: " << g_timer() << endl;
    }
}

bool NumericDominancePruning::is_dead_end(const State &state) {
    if(ldSimulation && ldSimulation->has_dominance_relation() &&
       /*is_activated() && */ldSimulation->get_dominance_relation().pruned_state(state)){
        deadends_pruned ++;
        return true;
    }
    
    for (auto & abs : abstractions) {
	if(abs->is_dead_end(state)) {
	    return true;
	}
    }
    return false;
}

int NumericDominancePruning::compute_heuristic(const State &state) {
    int cost = (ldSimulation && ldSimulation->has_dominance_relation()) ? 
	ldSimulation->get_dominance_relation().get_cost(state) : 0;
    if (cost == -1)
        return DEAD_END;

    for (auto & abs : abstractions) {
	if(abs->is_dead_end(state)) {
	    return DEAD_END;
	}
	cost = max (cost, abs->get_cost(state));
    }

    return cost;
}

void NumericDominancePruning::prune_applicable_operators(const State & state, int /*g*/,
					std::vector<const Operator *> & applicable_operators) {

    if(pruning_type == NumericPruningType::Successor ||
       pruning_type == NumericPruningType::ParentSuccessor) {
	vector<int> parent (g_variable_domain.size());
	for(int i = 0; i < g_variable_domain.size(); ++i) {
	    parent[i] = state[i];
	}
	vector<int> succ (parent);
	for (auto op : applicable_operators) {
	    for(const auto & prepost : op->get_pre_post()){
		succ[prepost.var] = prepost.post;
	    }

	    //TODO: Use adjusted cost instead.
	    if(numeric_dominance_relation->dominates_parent(succ, parent, op->get_cost())) {
		//cout << (applicable_operators.size() - 1) << " operators pruned because I have " << op->get_name() << endl;
		applicable_operators.clear();
		applicable_operators.push_back(op);
		return;
	    }
	    
	    for(const auto & prepost : op->get_pre_post()){
		succ[prepost.var] = parent[prepost.var];
	    }
	}
    }    
}

bool NumericDominancePruning::prune_generation(const State &state, int /*g*/, const State &parent, int action_cost ) {
    if(pruning_type == NumericPruningType::None) return false;
    if(!is_activated()) return false;
    
    // if(states_inserted%1000 == 0){
    // 	cout << "Deadend is still activated. "
    // 	     << "  States inserted: " << states_inserted
    // 	     << "  Min insertions: " << min_insertions << " " << min_deadends
    // 	     << "  States pruned: " << states_pruned
    // 	     << "  Limit pruned: " << states_inserted*min_pruning_ratio
    // 	     << "  Deadends pruned: " << deadends_pruned
    // 	     << "  Limit deadends: " << states_inserted*min_deadend_ratio
    // 	     << "  Prune activated: " << prune_is_activated()
    // 	     << "  Deadend activated: " << deadend_is_activated()
    // 	     <<  endl;
    // }

    if(numeric_dominance_relation->pruned_state(state)){
        return true;
    } 

    if(pruning_type == NumericPruningType::Parent || pruning_type == NumericPruningType::ParentSuccessor) {
	return numeric_dominance_relation->dominates(parent, state, action_cost);
    }

    //a) Check if state is in a BDD with g.closed <= g
    // states_checked ++;
    // if (check(state, g)){
    //     states_pruned ++;
    //     return true;
    // }

    // if(numeric_dominance_relation->dominates(parent, state)) {
    // 	cerr << "Fatal error." << endl;
    // 	cout << "Parent: "; parent.dump_pddl();
    // 	cout << "State: "; state.dump_pddl();

    //     insert(parent, g);

    // 	if (check(state, g)){
    // 	    return true;
    // 	}



    // 	exit(EXIT_CRITICAL_ERROR);
    // }


    // //b) Insert state and other states dominated by it
    // if(pruning_type == NumericPruningType::Generation){
    //     insert(state, g);
    //     states_inserted ++;
    // }
    return false;
}

bool NumericDominancePruning::prune_expansion (const State &/*state*/, int /*g*/){
    if(pruning_type == NumericPruningType::None ||  pruning_type == NumericPruningType::Parent
       || pruning_type == NumericPruningType::ParentSuccessor
       || pruning_type == NumericPruningType::Successor) return false;
    // if(!is_activated()) return false;
    
    // //a) Check if state is in a BDD with g.closed <= g
    // states_checked++;
    // if(check(state, g)){
    //     states_pruned ++;
    //     return true;
    // }
    // //b) Insert state and other states dominated by it
    // if(pruning_type == NumericPruningType::Expansion){
    //     insert(state, g);
    //     states_inserted ++;
    // }
    return false;
}

BDD NumericDominancePruning::getBDDToInsert(const State &state){
    // if(insert_dominated){
    //     BDD res = numeric_dominance_relation.getSimulatedBDD(vars.get(), state);
    //     if(remove_spurious_dominated_states){
    //         res = mgr->filter_mutex(res, true, 1000000, true);
    //         res = mgr->filter_mutex(res, false, 1000000, true);
    //     }
    //     if(pruning_type == NumericPruningType::Generation){
    //         res -= vars->getStateBDD(state); //Remove the state
    //     }else if (vars->numStates(res) == 1) {
    //         //Small optimization: If we have a single state, not include it
    //         return vars->zeroBDD();
    //     }
    //     return res;
    // }else{
         return vars->getStateBDD(state);
    // }
}

static PruneHeuristic *_parse(OptionParser &parser) {
    parser.document_synopsis("Simulation heuristic", "");
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional_effects", "supported (but see note)");
    parser.document_language_support("axioms", "not supported");
    parser.document_property("admissible", "yes");
    parser.document_property("consistent", "yes");
    parser.document_property("safe", "yes");
    parser.document_property("preferred operators", "no");
    parser.document_note(
            "Note",
            "Conditional effects are supported directly. Note, however, that "
            "for tasks that are not factored (in the sense of the JACM 2014 "
            "merge-and-shrink paper), the atomic abstractions on which "
            "merge-and-shrink heuristics are based are nondeterministic, "
            "which can lead to poor heuristics even when no shrinking is "
            "performed.");

    parser.add_option<bool>("remove_spurious",
            "If activated, remove spurious states from the sets of simulated/simulating states",
            "true");

    Heuristic::add_options_to_parser(parser);
    SymParamsMgr::add_options_to_parser_simulation(parser);

    parser.add_option<AbstractionBuilder *>(
            "abs",
            "abstraction builder",
            "");
    //LDSimulation::add_options_to_parser(parser);

    parser.add_enum_option
    ("pruning_type", NumericPruningTypeValues,
            "Implementation of the simulation pruning: "
            "Expansion: prunes states when they are simulated by an expanded state"
            "Generation: prunes states when they are simulated by a generated state ",
            "expansion");

    parser.add_option<bool>("insert_dominated",
            "Whether we store the set of dominated states (default) or just the set of closed.",
            "true");

    parser.add_option<double>("min_desactivation_ratio",
            "Ratio of pruned/checked needed to continue pruning the search.",
            "0.0");

    parser.add_option<int>("min_insertions",
            "States are inserted and pruning until this limit. Afterwards, depends on the ratios",
            "1000");

    Options opts = parser.parse();

    if (parser.dry_run()) {
        return 0;
    } else {
	return new NumericDominancePruningBDDMap (opts);
    }
}


static Heuristic *_parse_h(OptionParser &parser) {
    return static_cast<Heuristic *> (_parse(parser));
}


static Plugin<PruneHeuristic> _plugin("num_simulation", _parse);
static Plugin<Heuristic> _plugin_h("num_simulation", _parse_h);


void NumericDominancePruningBDDMap::insert (const State & state, int g){
    //Timer t;
    BDD res = getBDDToInsert(state);
    //time_bdd += t();
    if (!closed.count(g)){
        closed[g] = res;
    }else{
        closed[g] += res;
    }
    /*time_insert += t();
    if(states_inserted % 1000 == 0){
	cout << "SimulationClosed: ";
	for (auto & entry : closed){
	    cout << " " << entry.first << "("<<entry.second.nodeCount()<< "),";
	}
	cout << " after " << states_inserted << endl;
	cout << time_bdd << ", " << time_insert << " and " << time_check << " of " << g_timer() << endl;
	}*/

}

bool NumericDominancePruningBDDMap::check (const State & /* state*/, int /* g*/){
    //Timer t;
    // //CODE TO TEST Simulation TR
    // BDD dominatedBDD = vars->oneBDD();
    // for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
    //   dominatedBDD *= (*it)->getSimulatedBDD(state);
    // }

    // BDD sBDD = vars->getStateBDD(state);
    // BDD trBDD = tr->image(sBDD);
    // if (trBDD != dominatedBDD) {
    //   cerr << "ERROR: the TR does not work " << endl;
    //   exit(0);
    // }

    // //END CODE TO TEST
    // if(insert_dominated){
    //     auto sb = vars->getBinaryDescription(state);
    //     for(auto & entry : closed){
    //         if(entry.first > g) break;
    //         if(!(entry.second.Eval(sb).IsZero())){
    //             return true;
    //         }
    //     }
    // }else{
    //     BDD simulatingBDD = numeric_dominance_relation.getSimulatingBDD(vars.get(), state);
    //     for(auto & entry : closed){
    //         if(entry.first > g) break;
    //         if(!((entry.second*simulatingBDD).IsZero())){
    //             return true;
    //         }
    //     }
    // }
    //time_check += t();

    return false;
}

void NumericDominancePruning::print_statistics()
 {
     if(mgr){
	 cout << "Dominance BDD nodes: " << mgr->totalNodes() << endl;
	 cout << "Dominance BDD memory: " << mgr->totalMemory() << endl;
     }
 }


std::ostream & operator<<(std::ostream &os, const NumericPruningType & pt){
    switch(pt){
    case NumericPruningType::Expansion: return os << "expansion";
    case NumericPruningType::Generation: return os << "generation";
    case NumericPruningType::Parent: return os << "parent";

    case NumericPruningType::Successor: return os << "successor";
    case NumericPruningType::ParentSuccessor: return os << "parent_successor";
    case NumericPruningType::None: return os << "none";
    default:
        std::cerr << "Name of NumericPruningTypeStrategy not known";
        exit(-1);
    }
}

const std::vector<std::string> NumericPruningTypeValues {
    "expansion", "generation", "parent", "successor", "parent_successor", "none"
};
