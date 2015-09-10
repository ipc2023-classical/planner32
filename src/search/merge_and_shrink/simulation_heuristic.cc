#include "simulation_heuristic.h"

#include "abstraction.h"
#include "labels.h"
#include "merge_strategy.h"

#include "simulation_relation.h"
#include "labelled_transition_system.h"
#include "ld_simulation.h"
#include "shrink_bisimulation.h"

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

SimulationHeuristic::SimulationHeuristic(const Options &opts)
: PruneHeuristic(opts), 
  mgrParams(opts), initialized(false),
  remove_spurious_dominated_states(opts.get<bool>("remove_spurious")),
  insert_dominated(opts.get<bool>("insert_dominated")),
  pruning_type(PruningType(opts.get_enum("pruning_type"))),
  min_insertions_desactivation(opts.get<int>("min_insertions")),
  min_desactivation_ratio(opts.get<double>("min_desactivation_ratio")),
  vars(new SymVariables()), ldSimulation(new LDSimulation(opts)),
  all_desactivated(false), activation_checked(false), 
  states_inserted(0), states_checked(0), states_pruned(0), deadends_pruned(0) {
}

SimulationHeuristic::~SimulationHeuristic() {
}

void SimulationHeuristic::dump_options() const {
    cout << "Type pruning: " << pruning_type;
}

void SimulationHeuristic::initialize() {
    if(!initialized){
        initialized = true;
        ldSimulation->initialize();

        vector <int> var_order;
        ldSimulation->getVariableOrdering(var_order);

        vars->init(var_order, mgrParams);
        if(remove_spurious_dominated_states){
            mgr = unique_ptr<SymManager> (new SymManager(vars.get(), nullptr, mgrParams));
            mgr->init();
        }
        mgrParams.print_options();

        if(insert_dominated){
            ldSimulation->get_dominance_relation().precompute_dominated_bdds(vars.get());
        }else{
            ldSimulation->get_dominance_relation().precompute_dominating_bdds(vars.get());
        }
        cout << "Completed preprocessing: " << g_timer() << endl;
    }
}

bool SimulationHeuristic::is_dead_end(const State &state) {
    if(is_activated() && ldSimulation->pruned_state(state)){
        deadends_pruned ++;
        return true;
    }
    return false;
}

int SimulationHeuristic::compute_heuristic(const State &state) {
    int cost = ldSimulation->get_cost(state);
    if (cost == -1)
        return DEAD_END;
    return cost;
}

bool SimulationHeuristic::prune_generation(const State &state, int g) {
    if(pruning_type == PruningType::None) return false;
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

    if(ldSimulation->pruned_state(state)){
        return true;
    }
    //a) Check if state is in a BDD with g.closed <= g
    states_checked ++;
    if (check(state, g)){
        states_pruned ++;
        return true;
    }

    //b) Insert state and other states dominated by it
    if(pruning_type == PruningType::Generation){
        insert(state, g);
        states_inserted ++;
    }
    return false;
}

bool SimulationHeuristic::prune_expansion (const State &state, int g){
    if(pruning_type == PruningType::None) return false;
    if(!is_activated()) return false;
    
    //a) Check if state is in a BDD with g.closed <= g
    states_checked++;
    if(check(state, g)){
        states_pruned ++;
        return true;
    }
    //b) Insert state and other states dominated by it
    if(pruning_type == PruningType::Expansion){
        insert(state, g);
        states_inserted ++;
    }
    return false;
}

BDD SimulationHeuristic::getBDDToInsert(const State &state){
    if(insert_dominated){
        BDD res = ldSimulation->get_dominance_relation().getSimulatedBDD(vars.get(), state);
        if(remove_spurious_dominated_states){
            res = mgr->filter_mutex(res, true, 1000000, true);
            res = mgr->filter_mutex(res, false, 1000000, true);
        }
        if(pruning_type == PruningType::Generation){
            res -= vars->getStateBDD(state); //Remove the state
        }else if (vars->numStates(res) == 1) {
            //Small optimization: If we have a single state, not include it
            return vars->zeroBDD();
        }
        return res;
    }else{
        return vars->getStateBDD(state);
    }
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
    LDSimulation::add_options_to_parser(parser);

    parser.add_enum_option
    ("pruning_type", PruningTypeValues,
            "Implementation of the simulation pruning: "
            "Expansion: prunes states when they are simulated by an expanded state"
            "Generation: prunes states when they are simulated by a generated state ",
            "expansion");

    parser.add_enum_option
    ("pruning_dd", PruningDDValues,
            "Implementation data structure of the simulation pruning: "
            "BDD_MAP: (default) inserts in a map of BDDs all the dominated states "
            "ADD: inserts in an ADD all the dominated states "
            "BDD: inserts in a BDD all the dominated states (does"
            "not consider the g-value so it is not safe to use it with A* search)"
            "BDD_MAP_DISJ: allows a disjunctive partitioning of BDDs in BDD_MAP "
            "SKYLINE_BDD_MAP, SKYLINE_BDD: inserts the real states instead of dominated states"
            "BDD_MAP");

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
        switch(PruningDD(opts.get_enum("pruning_dd"))){
        case PruningDD::BDD_MAP: return new SimulationHeuristicBDDMap (opts);
        //case PruningType::ADD_DOMINATED: return new SimulationHeuristicADD (opts, true);
        case PruningDD::BDD: return new SimulationHeuristicBDD (opts);
        case PruningDD::BDD_MAP_DISJ: return new SimulationHeuristicBDDMapDisj (opts);
        case PruningDD::SKYLINE_BDD_MAP: return new SimulationHeuristicSkylineBDDMap (opts);
        case PruningDD::SKYLINE_BDD: return new SimulationHeuristicSkylineBDD (opts);
        default:
            std::cerr << "Name of PruningTypeStrategy not known";
            exit(-1);
        }

        return nullptr;
    }
}


static Heuristic *_parse_h(OptionParser &parser) {
    return static_cast<Heuristic *> (_parse(parser));
}


static Plugin<PruneHeuristic> _plugin("simulation", _parse);
static Plugin<Heuristic> _plugin_h("simulation", _parse_h);


std::ostream & operator<<(std::ostream &os, const PruningDD & pt){
    switch(pt){
    case PruningDD::BDD_MAP: return os << "BDD map";
    case PruningDD::ADD: return os << "ADD";
    case PruningDD::BDD: return os << "BDD";
    case PruningDD::BDD_MAP_DISJ: return os << "BDDmapDisj";
    case PruningDD::SKYLINE_BDD_MAP: return os << "SkylineBDDmap";
    case PruningDD::SKYLINE_BDD: return os << "SkylineBDD";
    default:
        std::cerr << "Name of PruningTypeStrategy not known";
        exit(-1);
    }
}

std::ostream & operator<<(std::ostream &os, const PruningType & pt){
    switch(pt){
    case PruningType::Expansion: return os << "expansion";
    case PruningType::Generation: return os << "generation";
    case PruningType::None: return os << "none";
    default:
        std::cerr << "Name of PruningTypeStrategy not known";
        exit(-1);
    }
}


const std::vector<std::string> PruningDDValues {
    "BDD_MAP", "ADD", "BDD", "BDD_MAP_DISJ", "SKYLINE_BDD_MAP", "SKYLINE_BDD"
};

const std::vector<std::string> PruningTypeValues {
    "expansion", "generation", "none"
};


float time_insert = 0, time_check = 0, time_bdd = 0;
void SimulationHeuristicBDDMap::insert (const State & state, int g){
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

bool SimulationHeuristicBDDMap::check (const State & state, int g){
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
    if(insert_dominated){
        auto sb = vars->getBinaryDescription(state);
        for(auto & entry : closed){
            if(entry.first > g) break;
            if(!(entry.second.Eval(sb).IsZero())){
                return true;
            }
        }
    }else{
        BDD simulatingBDD = ldSimulation->get_dominance_relation().getSimulatingBDD(vars.get(), state);
        for(auto & entry : closed){
            if(entry.first > g) break;
            if(!((entry.second*simulatingBDD).IsZero())){
                return true;
            }
        }
    }
    //time_check += t();

    return false;
}
void SimulationHeuristicBDD::insert (const State & state, int /*g*/){
    if(!initialized){
        closed = vars->zeroBDD();
        closed_inserted = vars->zeroBDD();
        initialized=true;
    }

    closed += getBDDToInsert(state);

    /*if(use_expensive_statistics && states_inserted % 1000 == 0){
      cout << "SimulationClosed: " << vars->numStates(closed) << "   " << closed.nodeCount() 
      << " after " << states_inserted << " using " <<closed_inserted.nodeCount() << endl;
      }*/
}

bool SimulationHeuristicBDD::check (const State & state, int /*g*/){
    if(!initialized){
        closed = vars->zeroBDD();
        closed_inserted = vars->zeroBDD();
        initialized=true;
    }

    if(insert_dominated){
        auto sb = vars->getBinaryDescription(state);
        return !(closed.Eval(sb).IsZero());
    } else{
        BDD simulatingBDD = ldSimulation->get_dominance_relation().getSimulatingBDD(vars.get(), state);
        return !((closed*simulatingBDD).IsZero());
    }
}


// void SimulationHeuristicBDD::insert (const BDD & bdd, int /*g*/){
//   if(insert_dominated){
//     closed += tr->image(bdd);
//   }else{
//     closed += bdd;
//   }
// }

// BDD SimulationHeuristicBDD::check (const BDD & bdd, int /*g*/){
//   if(insert_dominated){
//     return bdd*!closed;
//   }else{
//     cerr << "Pruning with not insert dominated not supported yet in symbolic search" << endl;
//     exit(0);
//   }
// }

// void SimulationHeuristicBDDMap::insert (const BDD & /*bdd*/, int /*g*/){
//   cerr << "Pruning with BDDMap not supported yet in symbolic search" << endl;
//   exit(0);
// }

// BDD SimulationHeuristicBDDMap::check (const BDD & /*bdd*/, int /*g*/){
//   cerr << "Pruning with BDDMap not supported yet in symbolic search" << endl;
//   exit(0);
// }



void SimulationHeuristicBDDMapDisj::insert (const State & state, int g){
    Timer t;
    BDD res = getBDDToInsert(state);
    time_bdd += t();

    if(!closed[g].empty() && closed[g].back().nodeCount() <= 1000){
        closed[g][closed[g].size() - 1] += res;
    }else{
        closed[g].push_back(res);
    }

    time_insert += t();
    if(states_inserted % 1000 == 0){
        cout << "SimulationClosed: ";
        for (auto & entry : closed){
            cout << " " << entry.first << "("<<entry.second.size() << "),";
        }
        cout << " after " << states_inserted << endl;
        cout << time_bdd << ", " << time_insert << " and " << time_check << " of " << g_timer() << endl;
    }
}

bool SimulationHeuristicBDDMapDisj::check (const State & state, int g){
    Timer t;

    if(insert_dominated){
        auto sb = vars->getBinaryDescription(state);
        for(auto & entry : closed){
            if(entry.first > g) break;
            for(auto & bdd : entry.second) {
                if(!(bdd.Eval(sb).IsZero())){
                    return true;
                }
            }
        }
    }else{
        BDD simulatingBDD = ldSimulation->get_dominance_relation().getSimulatingBDD(vars.get(), state);
        for(auto & entry : closed){
            if(entry.first > g) break;
            for(auto & bdd : entry.second) {
                if(!((bdd*simulatingBDD).IsZero())){
                    return true;
                }
            }
        }
    }
    time_check += t();

    return false;
}



bool SimulationHeuristicSkylineBDDMap::check (const State & state, int g){
    auto & sims = ldSimulation->get_dominance_relation();

    if(closed.empty()){
        closed.resize(sims.size());
	for(int part = 0; part < sims.size(); part++) {
	    closed[part].resize(sims[part].num_states());
	}
    }

    for (int gval : g_values){
	if(gval > g){
	    return false;
	}
        
	BDD res = vars->oneBDD();

	for(int part = 0; !res.IsZero() && part < sims.size(); part++) {
	    if(insert_dominated){
		int val = sims[part].get_index(state);
		if (closed[part][val].count(gval)){		    
		    res *= closed[part][val][gval];
		}else{
		    return false;
		}
	    } else {
		BDD dom = vars->zeroBDD(); 
		for (int valp : sims[part].get_dominating_states (state)){
		    if (closed[part][valp].count(gval)){
			dom += closed[part][valp][gval];
		    }
		}
		res *= dom;
	    }
	}

	if (!res.IsZero()){
	    return true;
	}
    }
    return false;
}

void SimulationHeuristicSkylineBDDMap::insert (const State & state, int g){
    g_values.insert(g);
    auto & sims = ldSimulation->get_dominance_relation();

    if(closed.empty()){
        closed.resize(sims.size());
	for(int part = 0; part < sims.size(); part++) {
	    closed[part].resize(sims[part].num_states());
	}
    }

    BDD bdd = vars->getStateBDD(state);
    
    for(int part = 0; part < sims.size(); part++){
	if (insert_dominated){
	    for (int valp : sims[part].get_dominated_states (state)){
		if (!closed[part][valp].count(g)){
		    closed[part][valp][g] = bdd;
		}else{
		    closed[part][valp][g] += bdd;
		}
	    }
	} else {
	    int val = sims[part].get_index(state);
	    if (!closed[part][val].count(g)){
		closed[part][val][g] = bdd;
	    }else{
		closed[part][val][g] += bdd;
	    }
	}
    }
}






bool SimulationHeuristicSkylineBDD::check (const State & state, int /*g*/){
    auto & sims = ldSimulation->get_dominance_relation();

    if(closed.empty()){
        closed.resize(sims.size());
	for(int part = 0; part < sims.size(); part++) {
	    closed[part].resize(sims[part].num_states(), vars->zeroBDD());
	}
    }
      
    BDD res = vars->oneBDD();

    for(int part = 0; !res.IsZero() && part < sims.size(); part++) {
	if(insert_dominated){		     
	    int val = sims[part].get_index(state);
	    res *= closed[part][val];
	} else {
	    BDD dom = vars->zeroBDD();
	    for (int valp : sims[part].get_dominating_states (state)){
		dom += closed[part][valp];
	    }
	    res *= dom;
	}
    }

    return !res.IsZero();
}

    void SimulationHeuristicSkylineBDD::insert (const State & state, int /*g*/) {
    auto & sims = ldSimulation->get_dominance_relation ();

    if(closed.empty()){
        closed.resize(sims.size());
	for(int part = 0; part < sims.size(); part++) {
	    closed[part].resize(sims[part].num_states(), vars->zeroBDD());
	}
    }

    BDD bdd = vars->getStateBDD(state);
    for(int part = 0; part < sims.size(); part++){
	if (insert_dominated){
	    for (int valp : sims[part].get_dominated_states (state)){
		closed[part][valp] += bdd;
	    }
	} else {
	    int val = sims[part].get_index(state);
	    closed[part][val] += bdd;
	}
    }
}
