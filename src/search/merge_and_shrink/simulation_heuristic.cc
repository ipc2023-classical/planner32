#include "simulation_heuristic.h"

#include "abstraction.h"
#include "labels.h"
#include "merge_strategy.h"

#include "simulation_relation.h"
#include "labelled_transition_system.h"

#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../state.h"
#include "../timer.h"

#include "../sym/sym_transition.h"

#include <cassert>
#include <vector>
using namespace std;

SimulationHeuristic::SimulationHeuristic(const Options &opts, bool _insert_dominated)
  : PruneHeuristic(opts), merge_strategy(opts.get<MergeStrategy *>("merge_strategy")), 
    use_expensive_statistics(opts.get<bool>("expensive_statistics")),
    vars(new SymVariables()),  insert_dominated(_insert_dominated), 
    limit_absstates_merge(opts.get<int>("limit_merge")),
    mgrParams(opts) {
    labels = new Labels(is_unit_cost_problem(), opts, cost_type);
    vector <int> var_order; 
    for(int i = 0; i < g_variable_name.size(); i++){
      var_order.push_back(i);
    }
    vars->init(var_order, mgrParams);
    mgrParams.print_options();
}

SimulationHeuristic::~SimulationHeuristic() {
    delete merge_strategy;
    delete labels;
    for(auto abs : abstractions){
      delete abs;
    }

    for(auto sim : simulations){
      delete sim;
    }
    cout << "Expensive statistics: "
         << (use_expensive_statistics ? "enabled" : "disabled") << endl;

    if (use_expensive_statistics) {
        string dashes(79, '=');
        cerr << dashes << endl
             << ("WARNING! You have enabled extra statistics for "
            "merge-and-shrink heuristics.\n"
            "These statistics require a lot of time and memory.\n"
            "When last tested (around revision 3011), enabling the "
            "extra statistics\nincreased heuristic generation time by "
            "76%. This figure may be significantly\nworse with more "
            "recent code or for particular domains and instances.\n"
            "You have been warned. Don't use this for benchmarking!")
        << endl << dashes << endl;
    }
}

void SimulationHeuristic::dump_options() const {
    merge_strategy->dump_options();
    labels->dump_options();
    cout << "Type prunning: ";

}
void SimulationHeuristic::build_abstraction() {
    // TODO: We're leaking memory here in various ways. Fix this.
    //       Don't forget that build_atomic_abstractions also
    //       allocates memory.

    // vector of all abstractions. entries with 0 have been merged.
    vector<Abstraction *> all_abstractions;
    all_abstractions.reserve(g_variable_domain.size() * 2 - 1);
    Abstraction::build_atomic_abstractions(all_abstractions, labels);

    // cout << "Shrinking atomic abstractions..." << endl;
    // for (size_t i = 0; i < all_abstractions.size(); ++i) {
    //     all_abstractions[i]->compute_distances();
    //     if (!all_abstractions[i]->is_solvable())
    //         return all_abstractions[i];
    //     shrink_strategy->shrink_atomic(*all_abstractions[i]);
    // }

    cout << "Merging abstractions..." << endl;

    while (!merge_strategy->done()) {
      pair<int, int> next_systems = merge_strategy->get_next(all_abstractions, 
							     limit_absstates_merge);
      int system_one = next_systems.first;
      if(system_one == -1){
	break; //No pairs to be merged under the limit
      }
        Abstraction *abstraction = all_abstractions[system_one];
        assert(abstraction);
        int system_two = next_systems.second;
        assert(system_one != system_two);
        Abstraction *other_abstraction = all_abstractions[system_two];
        assert(other_abstraction);

        // Note: we do not reduce labels several times for the same abstraction
        // bool reduced_labels = false;
        // if (shrink_strategy->reduce_labels_before_shrinking()) {
        //     labels->reduce(make_pair(system_one, system_two), all_abstractions);
        //     reduced_labels = true;
        //     abstraction->normalize();
        //     other_abstraction->normalize();
        //     abstraction->statistics(use_expensive_statistics);
        //     other_abstraction->statistics(use_expensive_statistics);
        // }

        // distances need to be computed before shrinking
        abstraction->compute_distances();
        other_abstraction->compute_distances();
        // if (!abstraction->is_solvable())
        //     return abstraction;
        // if (!other_abstraction->is_solvable())
        //     return other_abstraction;

        //shrink_strategy->shrink_before_merge(*abstraction, *other_abstraction);
        abstraction->statistics(use_expensive_statistics);
        other_abstraction->statistics(use_expensive_statistics);

        // if (!reduced_labels) {
        //     labels->reduce(make_pair(system_one, system_two), all_abstractions);
        // }
        abstraction->normalize();
        other_abstraction->normalize();
        // if (!reduced_labels) {
        //     // only print statistics if we just possibly reduced labels
        //     other_abstraction->statistics(use_expensive_statistics);
        //     abstraction->statistics(use_expensive_statistics);
        // }

        Abstraction *new_abstraction = new CompositeAbstraction(labels,
                                                                abstraction,
                                                                other_abstraction);

        abstraction->release_memory();
        other_abstraction->release_memory();

        new_abstraction->statistics(use_expensive_statistics);

        all_abstractions[system_one] = 0;
        all_abstractions[system_two] = 0;
        all_abstractions.push_back(new_abstraction);
    }

    for (size_t i = 0; i < all_abstractions.size(); ++i) {
        if (all_abstractions[i]) {
	  abstractions.push_back(all_abstractions[i]);
	  all_abstractions[i]->compute_distances();
	  all_abstractions[i]->statistics(use_expensive_statistics);
	  //all_abstractions[i]->release_memory();
        }
    }
}


void SimulationHeuristic::initialize(bool explicit_checker) {
    Timer timer;
    cout << "Initializing simulation heuristic..." << endl;
    dump_options();
    verify_no_axioms();
 
    //build_abstraction();
    Abstraction::build_atomic_abstractions(abstractions, labels);

    cout << "Building LTS" << endl;
    vector<LabelledTransitionSystem *> lts;
    for (auto a : abstractions){
      lts.push_back(new LabelledTransitionSystem(a));
    }
    
    cout << "Computing simulation..." << endl;
    SimulationRelation::compute_label_dominance_simulation(lts, labels, simulations);

    for(int i = 0; i < simulations.size(); i++){ 
      simulations[i]->dump(lts[i]->get_names()); 
    } 

    for (auto l : lts){
      delete l;
    }

    if(explicit_checker){
      for (auto sim : simulations){
	sim->precompute_absstate_bdds(vars.get());
	if(insert_dominated){
	  sim->precompute_dominated_bdds();
	}else{
	  sim->precompute_dominating_bdds();
	}
      }
    }

    cout << "Done initializing simulation heuristic [" << timer << "]"
         << endl;
}

SymTransition * SimulationHeuristic::getTR(SymVariables * _vars){
  if (!tr){
    for (auto sim : simulations){
      sim->precompute_absstate_bdds(_vars);
    }
    tr = unique_ptr<SymTransition> {new SymTransition(_vars, simulations)};
    cout << "Simulation TR: " << *tr << endl;
  }
  return tr.get();
}

int SimulationHeuristic::compute_heuristic(const State & /*state*/) {
  return 0;
}

bool SimulationHeuristic::prune_generation(const State &state, int g){
  for(auto sim : simulations) {
    if(sim->pruned(state)){
      return true;
    }
  }
  return check(state, g);
}

bool SimulationHeuristic::prune_expansion (const State &state, int g){
  //a) Check if state is in a BDD with g.closed <= g
  if(check(state, g)){
    return true;
  }
  //b) Insert state and other states dominated by it
  insert(state, g);
  return false;
}

// BDD SimulationHeuristic::prune_generation(const BDD &bdd, int g){
//   return check(bdd, g);
// }

// BDD SimulationHeuristic::prune_expansion (const BDD &bdd, int g){
//   //a) Get subset of not dominated states with g.closed <= g
//   BDD res = check(bdd, g);
//   //b) Insert state and other states dominated by it
//   insert(res, g);
//   return res;
// }


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

    parser.add_option<MergeStrategy *>(
        "merge_strategy",
        "merge strategy; choose between merge_linear and merge_dfp",
        "merge_linear");


    vector<string> label_reduction_method;
    label_reduction_method.push_back("NONE");
    label_reduction_method.push_back("OLD");
    label_reduction_method.push_back("TWO_ABSTRACTIONS");
    label_reduction_method.push_back("ALL_ABSTRACTIONS");
    label_reduction_method.push_back("ALL_ABSTRACTIONS_WITH_FIXPOINT");
    parser.add_enum_option("label_reduction_method", label_reduction_method,
                           "label reduction method: "
                           "none: no label reduction will be performed "
                           "old: emulate the label reduction as desribed in the "
                           "IJCAI 2011 paper by Nissim, Hoffmann and Helmert."
                           "two_abstractions: compute the 'combinable relation' "
                           "for labels only for the two abstractions that will "
                           "be merged next and reduce labels."
                           "all_abstractions: compute the 'combinable relation' "
                           "for labels once for every abstraction and reduce "
                           "labels."
                           "all_abstractions_with_fixpoint: keep computing the "
                           "'combinable relation' for labels iteratively for all "
                           "abstractions until no more labels can be reduced.",
                           "ALL_ABSTRACTIONS_WITH_FIXPOINT");

    vector<string> label_reduction_system_order;
    label_reduction_system_order.push_back("REGULAR");
    label_reduction_system_order.push_back("REVERSE");
    label_reduction_system_order.push_back("RANDOM");
    parser.add_enum_option("label_reduction_system_order", label_reduction_system_order,
                           "order of transition systems for the label reduction methods "
                           "that iterate over the set of all abstractions. only useful "
                           "for the choices all_abstractions and all_abstractions_with_fixpoint "
                           "for the option label_reduction_method.", "RANDOM");
    parser.add_option<bool>("expensive_statistics",
                            "show statistics on \"unique unlabeled edges\" (WARNING: "
                            "these are *very* slow, i.e. too expensive to show by default "
                            "(in terms of time and memory). When this is used, the planner "
                            "prints a big warning on stderr with information on the performance impact. "
                            "Don't use when benchmarking!)",
                            "false");

    parser.add_option<bool>("expensive_statistics",
                            "show statistics on \"unique unlabeled edges\" (WARNING: "
                            "these are *very* slow, i.e. too expensive to show by default "
                            "(in terms of time and memory). When this is used, the planner "
                            "prints a big warning on stderr with information on the performance impact. "
                            "Don't use when benchmarking!)",
                            "false");


    parser.add_option<int>("limit_merge",
			    "limit on the number of abstract states after the merge"
			    "By default: 1, does not perform any merge",
                            "1");


    Heuristic::add_options_to_parser(parser);
    SymParamsMgr::add_options_to_parser(parser);

    parser.add_enum_option
      ("pruning_type", PruningTypeValues,
       "Implementation of the simulation pruning: "
       "BDD_DOMINATED: (default) inserts in a map of BDDs all the dominated states "
       "ADD_DOMINATED: inserts in an ADD all the dominated states "
       "BDD_DOMINATED_BLIND: inserts in a BDD all the dominated states (does"
       "not consider the g-value so it is not safe to use it with A* search)"
       "BDD_CLOSED: inserts in a map of BDDs only the reached states "
       "ADD_CLOSED: inserts in an ADD only the reached states "
       "BDD_CLOSED_BLIND: inserts in a BDD only the reached states (does"
       "not consider the g-value so it is not safe to use it with A* search)", 
       "BDD_DOMINATED");

    
    Options opts = parser.parse();

    if (parser.dry_run()) {
        return 0;
    } else {
      
      switch(PruningType(opts.get_enum("pruning_type"))){
      case PruningType::BDD_DOMINATED: return new SimulationHeuristicBDDMap (opts, true);
	//case PruningType::ADD_DOMINATED: return new SimulationHeuristicADD (opts, true);
      case PruningType::BDD_DOMINATED_BLIND: return new SimulationHeuristicBDD (opts, true);
      case PruningType::BDD_CLOSED: return new SimulationHeuristicBDDMap (opts, false);
	//case PruningType::ADD_CLOSED:  return new SimulationHeuristicADD (opts, false);
      case PruningType::BDD_CLOSED_BLIND: return new SimulationHeuristicBDD (opts, false);
      default:
	std::cerr << "Name of PruningTypeStrategy not known";
	exit(-1);  
      }
      return nullptr;
    }
}

static Plugin<PruneHeuristic> _plugin("simulation", _parse);



std::ostream & operator<<(std::ostream &os, const PruningType & pt){
  switch(pt){
  case PruningType::BDD_DOMINATED: return os << "BDD dominated";
  case PruningType::ADD_DOMINATED: return os << "ADD dominated";
  case PruningType::BDD_DOMINATED_BLIND: return os << "BDD dominated blind";
  case PruningType::BDD_CLOSED: return os << "BDD closed";
  case PruningType::ADD_CLOSED: return os << "ADD closed";
  case PruningType::BDD_CLOSED_BLIND: return os << "BDD closed blind";
  default:
    std::cerr << "Name of PruningTypeStrategy not known";
    exit(-1);
  }
}

const std::vector<std::string> PruningTypeValues {
  "BDD_DOMINATED",   "ADD_DOMINATED",   "BDD_DOMINATED_BLIND",   
  "BDD_CLOSED",   "ADD_CLOSED",   "BDD_CLOSED_BLIND"
};


void SimulationHeuristicBDDMap::insert (const State & state, int g){
  BDD res = vars->oneBDD();
  if(insert_dominated){
    for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
      res *= (*it)->getSimulatedBDD(state); 
    }
  }else{
    res = vars->getStateBDD(state); 
  }

  if (!closed.count(g)){
    closed[g] = res;
  }else{
    closed[g] += res;
  }
}

bool SimulationHeuristicBDDMap::check (const State & state, int g){

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
    for(auto entry : closed){
      if(entry.first > g) break;
      if(!(entry.second.Eval(sb).IsZero())){
	return true;
      }
    }
  }else{
    BDD simulatingBDD = vars->oneBDD();
    for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
      simulatingBDD *= (*it)->getSimulatingBDD(state); 
    }
    for(auto entry : closed){
      if(entry.first > g) break;
      if(!((entry.second*simulatingBDD).IsZero())){
	return true;
      }
    }
  }
  return false;
}

void SimulationHeuristicBDD::insert (const State & state, int /*g*/){
  if(insert_dominated){
    BDD res = vars->oneBDD();
    for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
      res *= (*it)->getSimulatedBDD(state); 
    }
    closed += res;
  }else{
    closed += vars->getStateBDD(state);
  }

}

bool SimulationHeuristicBDD::check (const State & state, int /*g*/){
  if(insert_dominated){
    auto sb = vars->getBinaryDescription(state);
    return !(closed.Eval(sb).IsZero());
  } else{
    BDD simulatingBDD = vars->oneBDD();
    for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
      simulatingBDD *= (*it)->getSimulatingBDD(state); 
    }
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
