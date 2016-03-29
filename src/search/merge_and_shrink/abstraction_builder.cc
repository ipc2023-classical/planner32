#include "abstraction_builder.h" 

#include "labels.h"
#include "abstraction.h"
#include "dominance_relation.h"

#include "shrink_strategy.h"
#include "merge_strategy.h"
#include "variable_partition_finder.h"

#include "../option_parser.h"
#include "../plugin.h"

using namespace std;

void AbstractionBuilder::init_ldsim (bool unit_cost, OperatorCost cost_type, 
				     std::unique_ptr <LDSimulation> & ldSim) const {
    ldSim.reset (new LDSimulation(unit_cost, opts, cost_type));
}

void AbsBuilderComposite::build_abstraction (bool unit_cost, OperatorCost cost_type,
					     std::unique_ptr<LDSimulation> & ldSim, 
					     std::vector<std::unique_ptr<Abstraction> > & abstractions) const {
    for (auto st : strategies) {
	st->build_abstraction(unit_cost, cost_type, ldSim, abstractions);
    }
}


void AbsBuilderPDB::build_abstraction (bool unit_cost, OperatorCost cost_type,
				       std::unique_ptr<LDSimulation> & ldSim, 
				       std::vector<std::unique_ptr<Abstraction> > & /*abstractions*/) const {
    if (ldSim) {
	cerr << "Error: AbsBuilderPDB can only be used to initialize the abstractions" << endl;
	exit(EXIT_INPUT_ERROR);
    }

    init_ldsim(unit_cost, cost_type, ldSim);

    VariablePartitionGreedy v(limit_absstates_merge);
    ldSim->init_factored_systems(v.get_partition());
}


void AbsBuilderAtomic::build_abstraction (bool unit_cost, OperatorCost cost_type, 
					  std::unique_ptr<LDSimulation> & ldSim, 
					  std::vector<std::unique_ptr<Abstraction> > & /*abstractions*/) const {
    init_ldsim(unit_cost, cost_type, ldSim);

    if(!ldSim) {
	ldSim->init_atomic_abstractions();
    }

}

void AbsBuilderMasSimulation::build_abstraction (bool unit_cost, OperatorCost cost_type,
						 std::unique_ptr<LDSimulation> & ldSim, 
						 std::vector<std::unique_ptr<Abstraction> > & /*abstractions*/) const {
    if(!ldSim) {
	init_ldsim(unit_cost, cost_type, ldSim);
    }

    ldSim->build_abstraction(merge_strategy.get(), limit_absstates_merge, 
			     limit_transitions_merge, original_merge,
			     shrink_strategy.get(), forbid_lr, 
			     limit_seconds_mas, 
			     intermediate_simulations, incremental_simulations, 
			     simulation_type, 
			     label_dominance_type, 
			     switch_off_label_dominance, 
			     complex_lts, 
			     apply_subsumed_transitions_pruning, apply_label_dominance_reduction, 
			     apply_simulation_shrinking, expensive_statistics);


    if(compute_final_simulation)
	ldSim->compute_final_simulation(simulation_type, 
					label_dominance_type, 
					switch_off_label_dominance, 
					intermediate_simulations, complex_lts, 
					apply_subsumed_transitions_pruning, 
					apply_label_dominance_reduction, apply_simulation_shrinking,
					prune_dead_operators); 
}


void AbsBuilderMAS::build_abstraction (bool unit_cost, OperatorCost cost_type,
						 std::unique_ptr<LDSimulation> & ldSim, 
						 std::vector<std::unique_ptr<Abstraction> > & abstractions) const {
    if(!ldSim) {
	init_ldsim(unit_cost, cost_type, ldSim);
    }

    ldSim->complete_heuristic(merge_strategy.get(), shrink_strategy.get(), shrink_after_merge, 
			      limit_seconds_mas, expensive_statistics,abstractions);

} 



void AbsBuilderMasSimulation::dump_options() const {
    cout << "AbsBuilderMasSimulation" << endl;
    merge_strategy->dump_options();
    if (shrink_strategy) shrink_strategy->dump_options();
    else cout << " no shrinking" << endl;

    cout << "Expensive statistics: "
            << (expensive_statistics ? "enabled" : "disabled") << endl;

    if (expensive_statistics) {
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

AbstractionBuilder::AbstractionBuilder(const Options &opts_)  : 
    opts (opts_), expensive_statistics (opts.get<bool>("expensive_statistics"))  {     
}

AbsBuilderAtomic::AbsBuilderAtomic(const Options &opts) : 
    AbstractionBuilder (opts) {
}

AbsBuilderPDB::AbsBuilderPDB(const Options &opts)  : 
    AbstractionBuilder (opts), 
    limit_absstates_merge (opts.get<int>("limit_absstates_merge")) {     
}

AbsBuilderComposite::AbsBuilderComposite(const Options &opts) : 
    AbstractionBuilder (opts), 
    strategies(opts.get_list<AbstractionBuilder *>("strategies")) {
}


AbsBuilderMasSimulation::AbsBuilderMasSimulation(const Options &opts) : 
    AbstractionBuilder(opts), 
    simulation_type(SimulationType(opts.get_enum("simulation_type"))),
    label_dominance_type(LabelDominanceType(opts.get_enum("label_dominance_type"))),
    switch_off_label_dominance(opts.get<int>("switch_off_label_dominance")), 
    apply_simulation_shrinking(opts.get<bool>("apply_simulation_shrinking")),
    apply_subsumed_transitions_pruning(opts.get<bool>("apply_subsumed_transitions_pruning")),
    apply_label_dominance_reduction(opts.get<bool>("apply_label_dominance_reduction")),
    prune_dead_operators(opts.get<bool>("prune_dead_operators")),
    complex_lts(opts.get<bool>("complex_lts")),	
    merge_strategy(opts.get<MergeStrategy *>("merge_strategy")), 
    original_merge(opts.get<bool>("original_merge")),				   
    limit_absstates_merge(opts.get<int>("limit_merge")),
    limit_transitions_merge(opts.get<int>("limit_transitions_merge")), 
    intermediate_simulations(opts.get<bool>("intermediate_simulations")),
    incremental_simulations(opts.get<bool>("incremental_simulations")),
    compute_final_simulation(opts.get<bool>("compute_final_simulation")),
    forbid_lr(opts.get<bool>("forbid_lr")),
    shrink_strategy(opts.get<ShrinkStrategy *>("shrink_strategy")),
    shrink_after_merge(opts.get<bool>("shrink_after_merge")), 
    limit_seconds_mas(opts.get<int>("limit_seconds")) {
    if (incremental_simulations && !intermediate_simulations) {
        cerr << "Error: To use incremental calculation of simulations, intermediate simulations must be used!" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }

    Abstraction::store_original_operators = opts.get<bool>("store_original_operators");

    if (!prune_dead_operators && Abstraction::store_original_operators) {
        cerr << "Error: Why do you want to store operators if you don't prune them?" << endl;
        exit(1);
    }
}



AbsBuilderMAS::AbsBuilderMAS(const Options &opts) : 
    AbstractionBuilder(opts), 
    merge_strategy(opts.get<MergeStrategy *>("merge_strategy")), 
    shrink_strategy(opts.get<ShrinkStrategy *>("shrink_strategy")),
    shrink_after_merge(opts.get<bool>("shrink_after_merge")), 
    limit_seconds_mas(opts.get<int>("limit_seconds")) {
}


static AbstractionBuilder *_parse_composite(OptionParser &parser) {
    AbstractionBuilder::add_options_to_parser(parser);

    parser.add_list_option<AbstractionBuilder *>("strategies");

    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    opts.verify_list_non_empty<AbstractionBuilder *>("strategies");

    if (!parser.dry_run())
        return new AbsBuilderComposite(opts);
    else
        return 0;
}

static Plugin<AbstractionBuilder> _plugin_composite("builder_composite", _parse_composite);


static AbstractionBuilder *_parse_pdb(OptionParser &parser) {

    AbstractionBuilder::add_options_to_parser(parser);

    parser.add_option<int>("limit_absstates_merge",
			   "maximum number of states",
			   "10000");

    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    if (!parser.dry_run())
        return new AbsBuilderPDB(opts);
    else
        return 0;
}

static Plugin<AbstractionBuilder> _plugin_pdb("builder_pdb", _parse_pdb);




static AbstractionBuilder *_parse_atomic(OptionParser &parser) {

    AbstractionBuilder::add_options_to_parser(parser);

    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    if (!parser.dry_run())
        return new AbsBuilderAtomic(opts);
    else
        return 0;
}

static Plugin<AbstractionBuilder> _plugin_atomic("builder_atomic", _parse_atomic);



void AbstractionBuilder::add_options_to_parser(OptionParser &parser) {
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
            "all_abstractions_with_fixpoaint: keep computing the "
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
}



static AbstractionBuilder *_parse_massim(OptionParser &parser) {

    AbstractionBuilder::add_options_to_parser(parser);

    parser.add_option<int>("limit_merge",
            "limit on the number of abstract states after the merge"
            "By default: 1, does not perform any merge",
            "1");

    parser.add_option<int>("limit_transitions_merge",
            "limit on the number of transitions after the merge"
            "By default: 0: no limit at all",
            "0");

    parser.add_option<int>("limit_seconds",
            "limit the number of seconds for building the merge and shrink abstractions"
            "By default: 600 ",
            "600");


    parser.add_option<bool>("use_bisimulation",
            "If activated, use bisimulation to shrink abstractions before computing the simulation",
            "true");

    parser.add_option<bool>("intermediate_simulations",
            "Compute intermediate simulations and use them for shrinking",
            "false");

    parser.add_option<bool>("compute_final_simulation",
            "Compute intermediate simulations and use them for shrinking",
            "true");

    parser.add_option<bool>("incremental_simulations",
            "Compute incremental simulations and use them for shrinking",
            "false");

    parser.add_option<bool>("use_mas",
            "Use MaS to derive the factoring (or the factored strategy)",
            "true");

    parser.add_option<MergeStrategy *>(
            "merge_strategy",
            "merge strategy; choose between merge_linear and merge_dfp",
            "merge_linear");

    parser.add_option<bool>("complex_lts",
            "Use the complex method for LTS representation",
            "false");

    parser.add_option<bool>("apply_simulation_shrinking",
            "Perform simulation shrinking",
            "false");

    parser.add_option<bool>("apply_subsumed_transitions_pruning",
            "Perform pruning of subsumed transitions, based on simulation shrinking. Note: can only be used if simulation shrinking is applied!",
            "false");

    parser.add_option<bool>("apply_label_dominance_reduction",
            "Perform label reduction based on found label dominances",
            "false");

    parser.add_option<bool>("prune_dead_operators",
            "Prune all operators that are dead in some abstraction. Note: not yet implemented; so far, only the number of dead operators is returned!",
            "false");


    parser.add_option<bool>("forbid_lr",
            "Disable lr from the first part",
            "false");

    parser.add_option<bool>("store_original_operators",
            "Store the original operators for each transition in an abstraction",
            "false");

    parser.add_option<bool>("shrink_after_merge",
                            "If true, performs the shrinking after merge instead of before",
                            "false");

    parser.add_option<bool>("original_merge",
                            "Whether it continues merging variables after the next recommended merge has exceeded size",
                            "false");

    parser.add_option<ShrinkStrategy *>("shrink_strategy",
					"shrink strategy; ", 
					"none");

    parser.add_enum_option("simulation_type", SimulationTypeValues ,
			   "type of simulation implementation: NONE, SIMPLE or COMPLEX .", "SIMPLE" );


    parser.add_enum_option("label_dominance_type", LabelDominanceTypeValues ,
			   "type of simulation implementation: NONE, NOOP or NORMAL.", "NORMAL" );  


    parser.add_option<int>("switch_off_label_dominance",
            "disables label dominance if there are too many labels"
            "By default: 1000, to avoid memory errors",
            "1000");

    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    if (!parser.dry_run())
        return new AbsBuilderMasSimulation(opts);
    else
        return 0;
}


static Plugin<AbstractionBuilder> _plugin_massim("builder_massim", _parse_massim);


static AbstractionBuilder *_parse_mas(OptionParser &parser) {

    AbstractionBuilder::add_options_to_parser(parser);

    parser.add_option<int>("limit_seconds",
            "limit the number of seconds for building the merge and shrink abstractions"
            "By default: 600 ",
            "600");

    parser.add_option<MergeStrategy *>(
            "merge_strategy",
            "merge strategy; choose between merge_linear and merge_dfp",
            "merge_linear");

    parser.add_option<bool>("shrink_after_merge",
                            "If true, performs the shrinking after merge instead of before",
                            "false");

    parser.add_option<ShrinkStrategy *>("shrink_strategy",
					"shrink strategy; ", 
					"none");

    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    if (!parser.dry_run())
        return new AbsBuilderMAS(opts);
    else
        return 0;
}


static Plugin<AbstractionBuilder> _plugin_mas("builder_mas", _parse_mas);
