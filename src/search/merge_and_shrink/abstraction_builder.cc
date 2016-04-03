#include "abstraction_builder.h" 

#include "labels.h"
#include "abstraction.h"
#include "dominance_relation.h"

#include "shrink_strategy.h"
#include "shrink_bisimulation.h"
#include "shrink_composite.h"
#include "shrink_own_labels.h"



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
    Abstraction::store_original_operators = store_original_operators;

    if(!ldSim) {
	init_ldsim(unit_cost, cost_type, ldSim);
    }

    int remaining_time = max(0, min<int>(limit_seconds_mas, limit_seconds_total - g_timer()));


    ldSim->build_abstraction(merge_strategy.get(), limit_absstates_merge, 
			     limit_transitions_merge, original_merge,
			     shrink_strategy.get(), forbid_lr, 
			     remaining_time, limit_memory_kb_total, 
			     intermediate_simulations, incremental_simulations, 
			     simulation_type, 
			     label_dominance_type, 
			     switch_off_label_dominance, 
			     complex_lts, 
			     apply_subsumed_transitions_pruning, apply_label_dominance_reduction, 
			     apply_simulation_shrinking, /*preserve_all_optimal_plans*/ false,
				 expensive_statistics);


    if(compute_final_simulation)
	ldSim->compute_final_simulation(simulation_type, 
					label_dominance_type, 
					    switch_off_label_dominance, 
					    intermediate_simulations, complex_lts, 
					    apply_subsumed_transitions_pruning, 
					    apply_label_dominance_reduction, apply_simulation_shrinking,
					    /*preserve_all_optimal_plans*/false); 

    if(prune_dead_operators) ldSim->prune_dead_ops();

}


void AbsBuilderMAS::build_abstraction (bool unit_cost, OperatorCost cost_type,
				       std::unique_ptr<LDSimulation> & ldSim, 
				       std::vector<std::unique_ptr<Abstraction> > & abstractions) const {
    Abstraction::store_original_operators = store_original_operators;
    
    for(int i = 0; i < num_abstractions ; ++i) {
	int remaining_time = max(0, min<int>(limit_seconds_mas, limit_seconds_total - g_timer()));
	if(remaining_time <= 0) break;

	if(restart) {
	    ldSim->release_memory();
	    std::unique_ptr<LDSimulation> tmpldSim;
	    init_ldsim(unit_cost, cost_type, tmpldSim);
	    tmpldSim->init_atomic_abstractions();

	    tmpldSim->complete_heuristic(merge_strategy.get(),  shrink_strategy.get(),shrink_after_merge, 
					 remaining_time, limit_memory_kb_total, prune_dead_operators, expensive_statistics, abstractions);

	} else {
	    if(!ldSim) {
		init_ldsim(unit_cost, cost_type, ldSim);
		ldSim->init_atomic_abstractions();
	    }

	    ldSim->complete_heuristic(merge_strategy.get(),  shrink_strategy.get(), shrink_after_merge, 
				      remaining_time, limit_memory_kb_total, prune_dead_operators,  expensive_statistics,abstractions);
	
	}
    }
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
    opts (opts_), expensive_statistics (opts.get<bool>("expensive_statistics")), 
    limit_seconds_total(opts.get<int>("limit_seconds_total")), 
    limit_memory_kb_total(opts.get<int>("limit_memory_kb")) {     
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
    store_original_operators(opts.get<bool>("store_original_operators")),
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

    if (opts.get<bool>("incremental_pruning")){
	apply_subsumed_transitions_pruning=true;
	prune_dead_operators=true; 
	store_original_operators = true;
	intermediate_simulations=true;  
	incremental_simulations=true;
    }

    if (incremental_simulations && !intermediate_simulations) {
        cerr << "Error: To use incremental calculation of simulations, intermediate simulations must be used!" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }

    if (!prune_dead_operators && Abstraction::store_original_operators) {
        cerr << "Error: Why do you want to store operators if you don't prune them?" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }
}



AbsBuilderMAS::AbsBuilderMAS(const Options &opts) : 
    AbstractionBuilder(opts), 
    merge_strategy(opts.get<MergeStrategy *>("merge_strategy")), 
    shrink_strategy(opts.get<ShrinkStrategy *>("shrink_strategy")),
    shrink_after_merge(opts.get<bool>("shrink_after_merge")),
    limit_seconds_mas(opts.get<int>("limit_seconds")), 
    prune_dead_operators(opts.get<bool>("prune_dead_operators")),
    store_original_operators(opts.get<bool>("store_original_operators")),
    restart(opts.get<bool>("restart")), 
    num_abstractions(opts.get<int>("num_abstractions")) {
}



AbsBuilderDefault::AbsBuilderDefault(const Options &opts) : 
    AbstractionBuilder(opts), 
    merge_strategy(opts.get<MergeStrategy *>("merge_strategy")), original_merge(opts.get<bool>("original_merge")),				   
    limit_absstates_merge(opts.get<int>("limit_merge")),
    limit_transitions_merge(opts.get<int>("limit_transitions_merge")), 
    limit_absstates_shrink(opts.get<int>("limit_shrink")),

    limit_seconds_mas(opts.get<int>("limit_seconds")), 
    num_abstractions(opts.get<int>("num_abstractions")), 
    switch_off_label_dominance(opts.get<int>("switch_off_label_dominance")) {
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

    parser.add_option<int>("label_reduction_max_time",
			   "limit the number of seconds for label reduction",
			   "60");

    parser.add_option<bool>("expensive_statistics",
			    "show statistics on \"unique unlabeled edges\" (WARNING: "
			    "these are *very* slow, i.e. too expensive to show by default "
			    "(in terms of time and memory). When this is used, the planner "
			    "prints a big warning on stderr with information on the performance impact. "
			    "Don't use when benchmarking!)",
			    "false");


    parser.add_option<int>("limit_seconds_total",
			   "limit the number of seconds for building the merge and shrink abstractions"
			   "By default: 1400, reserving ~100 seconds for the preprocessor and ~300 for search  ",
			   "1400");

    parser.add_option<int>("limit_memory_kb",
			   "limit the memory for building the merge and shrink abstractions"
			   "By default:  ",
			   "2000000");


}



static AbstractionBuilder *_parse_massim(OptionParser &parser) {

    AbstractionBuilder::add_options_to_parser(parser);

    parser.add_option<int>("limit_seconds",
			   "limit the number of seconds for each iteration"
			   "By default: 300 ",
			   "300");

    parser.add_option<int>("limit_merge",
			   "limit on the number of abstract states after the merge"
			   "By default: 1, does not perform any merge",
			   "1");


    parser.add_option<bool>("original_merge",
                            "Whether it continues merging variables after the next recommended merge has exceeded size",
                            "false");

    parser.add_option<int>("limit_transitions_merge",
			   "limit on the number of transitions after the merge"
			   "By default: 0: no limit at all",
			   "0");


    parser.add_option<bool>("intermediate_simulations",
			    "Compute intermediate simulations and use them for shrinking",
			    "false");

    parser.add_option<bool>("compute_final_simulation",
			    "Compute intermediate simulations and use them for shrinking",
			    "true");

    parser.add_option<bool>("incremental_simulations",
			    "Compute incremental simulations and use them for shrinking",
			    "false");

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
			    "true");


    parser.add_option<bool>("forbid_lr",
			    "Disable lr from the first part",
			    "false");

    parser.add_option<bool>("store_original_operators",
			    "Store the original operators for each transition in an abstraction",
			    "false");

    parser.add_option<bool>("shrink_after_merge",
                            "If true, performs the shrinking after merge instead of before",
                            "false");

    parser.add_option<bool>("incremental_pruning",
                            "Sets to true apply_subsumed_transitions_pruning, prune_dead_operators, store_original_operators, intermediate_simulations, and incremental_simulations",
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
			   "200");

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
			   "limit the number of seconds for each iteration"
			   "By default: 300 ",
			   "300");

    parser.add_option<MergeStrategy *>(
	"merge_strategy",
	"merge strategy; choose between merge_linear and merge_dfp",
	"merge_dfp");

    parser.add_option<bool>("shrink_after_merge",
                            "If true, performs the shrinking after merge instead of before",
                            "false");

    parser.add_option<ShrinkStrategy *>("shrink_strategy",
					"shrink strategy; ", 
					"none");

    parser.add_option<bool>("restart",
                            "If true, starts from atomic abstraction heuristics",
                            "false");


    parser.add_option<int>("num_abstractions",
                            "how many abstractions should be generated",
                            "1");


    parser.add_option<bool>("store_original_operators",
			    "Store the original operators for each transition in an abstraction",
			    "false");

    parser.add_option<bool>("prune_dead_operators",
			    "Prune all operators that are dead in some abstraction. Note: not yet implemented; so far, only the number of dead operators is returned!",
			    "true");



    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    if (!parser.dry_run())
        return new AbsBuilderMAS(opts);
    else
        return 0;
}


static Plugin<AbstractionBuilder> _plugin_mas("builder_mas", _parse_mas);





static AbstractionBuilder *_parse_default(OptionParser &parser) {

    AbstractionBuilder::add_options_to_parser(parser);

    parser.add_option<int>("limit_seconds",
			   "limit the number of seconds for each iteration"
			   "By default: 300 ",
			   "300");

    parser.add_option<MergeStrategy *>(
	"merge_strategy",
	"merge strategy; choose between merge_linear and merge_dfp",
	"merge_dfp");


    parser.add_option<int>("num_abstractions",
                            "how many abstractions should be generated",
                            "1");


    parser.add_option<int>("limit_merge",
			   "limit on the number of abstract states after the merge"
			   "By default: 100000",
			   "100000");

    parser.add_option<int>("limit_shrink",
			   "limit on the number of abstract states for shrinking",
			   "100000");


    parser.add_option<bool>("original_merge",
                            "Whether it continues merging variables after the next recommended merge has exceeded size",
                            "false");

    parser.add_option<int>("limit_transitions_merge",
			   "limit on the number of transitions after the merge", 
			    "100000");

    parser.add_option<int>("switch_off_label_dominance",
			   "disables label dominance if there are too many labels"
			   "By default: 1000, to avoid memory errors",
			   "200");



    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    if (!parser.dry_run())
        return new AbsBuilderDefault(opts);
    else
        return 0;
}


static Plugin<AbstractionBuilder> _plugin_default("builder", _parse_default);







void AbsBuilderDefault::build_abstraction (bool unit_cost, OperatorCost cost_type,
					   std::unique_ptr<LDSimulation> & ldSim_, 
					   std::vector<std::unique_ptr<Abstraction> > & abstractions) const {

    Abstraction::store_original_operators = true;

    bool preserve_all_optimal_plans = ldSim_.get() != nullptr;
    LDSimulation * ldSim = nullptr;
    std::unique_ptr<LDSimulation> tmpldSim;

    if(!ldSim) {
	init_ldsim(unit_cost, cost_type, ldSim_);
	ldSim = ldSim_.get();
    } else {
	ldSim_->release_memory();
	init_ldsim(unit_cost, cost_type, tmpldSim);
	ldSim = tmpldSim.get(); 
    }


    int remaining_time = max(0, min<int>(limit_seconds_mas, limit_seconds_total - g_timer()));


   // 1) Incremental simulations without shrinking or label reduction

    cout << "1) Incremental simulations without shrinking or label reduction. Max states: "
	 << limit_absstates_merge << " transitions: " << limit_transitions_merge  << endl;

    ldSim->build_abstraction(merge_strategy.get(), limit_absstates_merge, 
			     limit_transitions_merge, /*original_merge*/ true,
			     nullptr, /*forbid_lr*/ true, 
			     remaining_time, limit_memory_kb_total, 
			     /*intermediate_simulations */true, /*incremental_simulations*/ true, 
			     SimulationType::SIMPLE, 
			     LabelDominanceType::NORMAL, 
			     switch_off_label_dominance, 
			     /*complex_lts*/ false, 
			     /*apply_subsumed_transitions_pruning*/ true, /*apply_label_dominance_reduction*/false, 
			     /*apply_simulation_shrinking*/false, preserve_all_optimal_plans, expensive_statistics);


    ldSim->compute_final_simulation(SimulationType::SIMPLE, 
				    LabelDominanceType::NORMAL,
				    switch_off_label_dominance, 
				    true, false, 
				    true, 
				    false, false, preserve_all_optimal_plans); 

    ldSim->prune_dead_ops();

    // 2) Incremental simulations with shrinking and label reduction
    cout << "2) Incremental simulations with shrinking and label reduction " << endl;
    unique_ptr<ShrinkStrategy> bisim (ShrinkBisimulation::create_default(true));

    remaining_time = max(0, min<int>(limit_seconds_mas, limit_seconds_total - g_timer()));

    ldSim->build_abstraction(merge_strategy.get(), limit_absstates_merge, 
			     limit_transitions_merge, original_merge,
			     bisim.get(), /*forbid_lr*/ false, 
			     remaining_time, limit_memory_kb_total, 
			     /*intermediate_simulations */true, /*incremental_simulations*/ true, 
			     SimulationType::SIMPLE, 
			     LabelDominanceType::NORMAL, 
			     switch_off_label_dominance, 
			     /*complex_lts*/ false, 
			     /*apply_subsumed_transitions_pruning*/ true, /*apply_label_dominance_reduction*/false, 
			     /*apply_simulation_shrinking*/false, preserve_all_optimal_plans, expensive_statistics);

    ldSim->compute_final_simulation(SimulationType::SIMPLE, 
				    LabelDominanceType::NORMAL,
				    switch_off_label_dominance, 
				    true, false, 
				    true, 
				    false, false, preserve_all_optimal_plans); 
    
    ldSim->prune_dead_ops();

    cout << "3) Complete abstractions" << endl;

    // 3) Complete abstractions 
    unique_ptr<ShrinkStrategy> bisim_limit (limit_absstates_shrink == 0 ? 
					    ShrinkBisimulation::create_default(true) :
					    ShrinkBisimulation::create_default(true, limit_absstates_shrink));

    unique_ptr<ShrinkStrategy> own (ShrinkOwnLabels::create_default());
    std::vector<ShrinkStrategy *> strategies; 
    strategies.push_back(own.get());
    strategies.push_back(bisim_limit.get());
    unique_ptr<ShrinkStrategy> shrink_combined (ShrinkComposite::create_default(strategies));
    


    Abstraction::store_original_operators = false;
    
    for(int i = 0; i < num_abstractions ; ++i) {
	int remaining_time = max(0, min<int>(limit_seconds_mas, limit_seconds_total - g_timer()));
	if(remaining_time <= 0) break;

	ldSim->complete_heuristic(merge_strategy.get(), shrink_combined.get(),  /*shrink_after_merge*/ false, 
				  remaining_time, limit_memory_kb_total, /*prune_dead_operators*/true,  
				  expensive_statistics, abstractions);
	
    }

    ldSim->prune_dead_ops();
	

}

