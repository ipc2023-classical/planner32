#include "simulation_heuristic.h"

#include "abstraction.h"
#include "labels.h"

#include "simulation_relation.h"
#include "labelled_transition_system.h"

#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../state.h"
#include "../timer.h"

#include <cassert>
#include <vector>
using namespace std;

SimulationHeuristic::SimulationHeuristic(const Options &opts)
    : Heuristic(opts) {
    labels = new Labels(is_unit_cost_problem(), opts, cost_type);
}

SimulationHeuristic::~SimulationHeuristic() {
    delete labels;
    for(auto abs : abstractions){
      delete abs;
    }

    for(auto sim : simulations){
      delete sim;
    }
}

void SimulationHeuristic::dump_options() const {
    labels->dump_options();
}


void SimulationHeuristic::initialize() {
    Timer timer;
    cout << "Initializing simulation heuristic..." << endl;
    dump_options();
    verify_no_axioms();
 
    Abstraction::build_atomic_abstractions(abstractions, labels);

    cout << "Building LTS" << endl;
    vector<LabelledTransitionSystem *> lts;
    for (auto a : abstractions){
      lts.push_back(new LabelledTransitionSystem(a));
    }
    
    cout << "Computing simulation..." << endl;
    SimulationRelation::compute_label_dominance_simulation(lts, labels, simulations);

    cout << "Done initializing simulation heuristic [" << timer << "]"
         << endl;

    for(int i = 0; i < simulations.size(); i++){ 
      simulations[i]->dump(lts[i]->get_names()); 
    } 

    for (auto l : lts){
      delete l;
    }

    exit(0);
}


int SimulationHeuristic::compute_heuristic(const State &/*state*/) {
  return 0;
}


static Heuristic *_parse(OptionParser &parser) {
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

    Heuristic::add_options_to_parser(parser);
    Options opts = parser.parse();

    if (parser.dry_run()) {
        return 0;
    } else {
        SimulationHeuristic *result = new SimulationHeuristic(opts);
        return result;
    }
}

static Plugin<Heuristic> _plugin("simulation", _parse);
