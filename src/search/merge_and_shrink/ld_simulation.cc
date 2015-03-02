#include "ld_simulation.h"

#include "abstraction.h"
#include "shrink_bisimulation.h"
#include "simulation_simple.h"
#include "simulation_identity.h"
#include "simulation_efficient.h"
#include "simulation_efficient_nold.h"
#include "merge_strategy.h"
#include "labelled_transition_system.h"
#include "opt_order.h"
#include "../globals.h"
#include "../causal_graph.h"
#include "label_reducer.h"
#include <boost/dynamic_bitset.hpp>

#include "variable_partition_finder.h"

using namespace std;

LDSimulation::LDSimulation(bool unit_cost, const Options &opts, OperatorCost cost_type) : 
		                          skip_simulation(opts.get<bool>("skip_simulation")),
		                          nold_simulation(opts.get<bool>("nold_simulation")),
		                          apply_simulation_shrinking(opts.get<bool>("apply_simulation_shrinking")),
		                          apply_subsumed_transitions_pruning(opts.get<bool>("apply_subsumed_transitions_pruning")),
		                          apply_label_dominance_reduction(opts.get<bool>("apply_label_dominance_reduction")),
		                          prune_dead_operators(opts.get<bool>("prune_dead_operators")),
					  forbid_lr(opts.get<bool>("forbid_lr")),	                          
                                          efficient_simulation(opts.get<bool>("efficient_simulation")),
		                          efficient_lts(opts.get<bool>("efficient_lts")),
		                          use_expensive_statistics(opts.get<bool>("expensive_statistics")),
		                          limit_absstates_merge(opts.get<int>("limit_merge")),
		                          limit_transitions_merge(opts.get<int>("limit_transitions_merge")),
		                          use_mas(opts.get<bool>("use_mas")),
		                          limit_seconds_mas(opts.get<int>("limit_seconds")),
		                          merge_strategy(opts.get<MergeStrategy *>("merge_strategy")),
		                          use_bisimulation(opts.get<bool>("use_bisimulation")),
		                          intermediate_simulations(opts.get<bool>("intermediate_simulations")),
		                          incremental_simulations(opts.get<bool>("incremental_simulations")),
		                          compute_final_abstraction(opts.get<bool>("compute_final_abstraction")),
    shrink_strategy(opts.get<ShrinkStrategy *>("shrink_strategy")),
    shrink_after_merge(opts.get<bool>("shrink_after_merge")) , 
    labels (new Labels(unit_cost, opts, cost_type)) //TODO: c++14::make_unique
{
    /*if (apply_subsumed_transitions_pruning && (! apply_simulation_shrinking && ! intermediate_simulations)) {
        cerr << "Error: can only apply pruning of subsumed transitions if simulation shrinking (either at the end or in an intermediate fashion) is used!" << endl;
        exit(1);
	}*/
    /*if (apply_subsumed_transitions_pruning && ! labels->applies_perfect_label_reduction()) {
        cerr << "Error: can only apply pruning of subsumed transitions if perfect label reduction is applied" << endl;
        exit(1);
    }*/
    if (incremental_simulations && !intermediate_simulations) {
        cerr << "Error: To use incremental calculation of simulations, intermediate simulations must be used!" << endl;
        exit(1);
    }
    if (incremental_simulations && efficient_simulation) {
        cerr << "Error: Support for incremental calculation of simulations not yet implemented in (supposedly) efficient simulation relation!" << endl;
        exit(1);
    }

    Abstraction::store_original_operators = opts.get<bool>("store_original_operators");

    if (!prune_dead_operators && Abstraction::store_original_operators) {
        cerr << "Error: Why do you want to store operators if you don't prune them?" << endl;
        exit(1);
    }

}

LDSimulation::LDSimulation(const Options &opts) : 
    		                        skip_simulation(opts.get<bool>("skip_simulation")),
    		                        nold_simulation(opts.get<bool>("nold_simulation")),
    		                        apply_simulation_shrinking(opts.get<bool>("apply_simulation_shrinking")),
    		                        apply_subsumed_transitions_pruning(opts.get<bool>("apply_subsumed_transitions_pruning")),
                                    apply_label_dominance_reduction(opts.get<bool>("apply_label_dominance_reduction")),
                                    prune_dead_operators(opts.get<bool>("prune_dead_operators")),
    		                     forbid_lr(opts.get<bool>("forbid_lr")),	                          
		                        efficient_simulation(opts.get<bool>("efficient_simulation")),
    		                        efficient_lts(opts.get<bool>("efficient_lts")),
    		                        use_expensive_statistics(opts.get<bool>("expensive_statistics")),
    		                        limit_absstates_merge(opts.get<int>("limit_merge")),
    		                        limit_transitions_merge(opts.get<int>("limit_transitions_merge")),
    		                        use_mas(opts.get<bool>("use_mas")),
    		                        limit_seconds_mas(opts.get<int>("limit_seconds")),
    		                        merge_strategy(opts.get<MergeStrategy *>("merge_strategy")),
    		                        use_bisimulation(opts.get<bool>("use_bisimulation")),
    		                        intermediate_simulations(opts.get<bool>("intermediate_simulations")),
		incremental_simulations(opts.get<bool>("incremental_simulations")), 
    compute_final_abstraction(opts.get<bool>("compute_final_abstraction")), 
    shrink_strategy(opts.get<ShrinkStrategy *>("shrink_strategy")),
		shrink_after_merge(opts.get<bool>("shrink_after_merge")){
    /*if (apply_subsumed_transitions_pruning && (! apply_simulation_shrinking && ! intermediate_simulations)) {
        cerr << "Error: can only apply pruning of subsumed transitions if simulation shrinking (either at the end or in an intermediate fashion) is used!" << endl;
        exit(1);
	}*/ 
   if (incremental_simulations && !intermediate_simulations) {
        cerr << "Error: To use incremental calculation of simulations, intermediate simulations must be used!" << endl;
        exit(1);
    }
    if (incremental_simulations && efficient_simulation) {
        cerr << "Error: Support for incremental calculation of simulations not yet implemented in (supposedly) efficient simulation relation!" << endl;
        exit(1);
    }

    //TODO: Copied from heuristic.cc. Move to the PlanningTask class
    //that someday will be added to FD.
    OperatorCost cost_type = OperatorCost::NORMAL;
    bool is_unit_cost = true;
    for (size_t i = 0; i < g_operators.size(); ++i) {
        if (get_adjusted_action_cost(g_operators[i], cost_type) != 1) {
            is_unit_cost = false;
            break;
        }
    }
    labels = unique_ptr<Labels> (new Labels(is_unit_cost, opts, cost_type)); //TODO: c++14::make_unique
    /*if (apply_subsumed_transitions_pruning && ! labels->applies_perfect_label_reduction()) {
        cerr << "Error: can only apply pruning of subsumed transitions if perfect label reduction is applied" << endl;
        exit(1);
    }*/
    Abstraction::store_original_operators = opts.get<bool>("store_original_operators");
    if (!prune_dead_operators && Abstraction::store_original_operators) {
        cerr << "Error: Why do you want to store operators if you don't prune them?" << endl;
        exit(1);
    }
}

LDSimulation::~LDSimulation(){
    for(auto abs : abstractions){
        delete abs;
    }
    for(auto sim : simulations){
        delete sim;
    }
}

void LDSimulation::build_factored_systems() {
    VariablePartitionGreedy v(limit_absstates_merge);
    v.get_partition();
    v.dump();

    for (auto factor : v.get_partition()) {
        PDBAbstraction * abs_factor = new PDBAbstraction(labels.get(),
                factor);
        abstractions.push_back(abs_factor);
        abs_factor->normalize();
        abs_factor->compute_distances();
        abs_factor->statistics(use_expensive_statistics);
    }
}

void LDSimulation::remove_dead_labels(vector<Abstraction *> & abstractions){
    vector<int> new_dead_labels; 
    dead_labels.resize(labels->get_size(), false);
    for (auto abs : abstractions) {
	if(abs) abs->get_dead_labels(dead_labels, new_dead_labels);
    }
    
    if(new_dead_labels.size() > 0){
	for(auto l : new_dead_labels){
	    for (auto abs : abstractions) {
		if(abs) abs->prune_transitions_dominated_label_all(l);
	    }
	}

	for (auto abs : abstractions) {
	    if(abs) abs->reset_lts();
	}
    }
}

int LDSimulation::remove_useless_abstractions(vector<Abstraction *> & abstractions, 
					      vector<SimulationRelation *> & simulations) {
    remove_dead_labels(abstractions);
    int removed_abstractions = 0;
    for(int i =0; i < abstractions.size(); ++i){
	if(abstractions[i] && abstractions[i]->is_useless()){
	    useless_vars.insert(end(useless_vars), begin(abstractions[i]->get_varset()), 
				end(abstractions[i]->get_varset()));
	    if (i < simulations.size() && simulations[i]){
		delete simulations[i];
		simulations[i] = nullptr;
	    }
	    abstractions[i]->release_memory();
	    delete abstractions[i];
	    abstractions[i] = nullptr;
	    
	    removed_abstractions++;
	}
    }
    return removed_abstractions;			    
}

// Just main loop copied from merge_and_shrink heuristic 
Abstraction * LDSimulation::complete_heuristic() const {
    cout << "Complete heuristic Initialized with " << abstractions.size() << " abstractions" << endl;
    //Insert atomic abstractions in the first g_variable_domain
    //variables, filling with nullptr. Then composite abstractions
    vector<Abstraction *> all_abstractions(g_variable_domain.size(), nullptr);
    int index_atomic = 0;
    for (auto a : abstractions) {
	if (a->get_varset().size() == 1) {
	    all_abstractions[index_atomic++] = a;
	}else{
	    all_abstractions.push_back(a);
	}
    }
    merge_strategy->set_remaining_merges(abstractions.size() - 1);
    if(abstractions.size() > 1){
	labels->reduce(make_pair(0, 1), all_abstractions); 
// With the reduction methods we use here, this should just apply label reduction on all abstractions
    }
    while (!merge_strategy->done()) {
        pair<int, int> next_systems = merge_strategy->get_next(all_abstractions);
        int system_one = next_systems.first;
        int system_two = next_systems.second;
	cout << " NEXT SYSTEMS: " << system_one <<  " " << system_two << endl;
        assert(system_one != system_two);


        Abstraction *abstraction = all_abstractions[system_one];
        assert(abstraction);
        Abstraction *other_abstraction = all_abstractions[system_two];
        assert(other_abstraction);


	if(!shrink_after_merge){
	    // Note: we do not reduce labels several times for the same abstraction
	    bool reduced_labels = false;
	    if (shrink_strategy->reduce_labels_before_shrinking()) {
		labels->reduce(make_pair(system_one, system_two), all_abstractions);
		reduced_labels = true;
		abstraction->normalize();
		other_abstraction->normalize();
		abstraction->statistics(use_expensive_statistics);
		other_abstraction->statistics(use_expensive_statistics);
	    }

	    // distances need to be computed before shrinking
	    abstraction->compute_distances();
	    other_abstraction->compute_distances();
	    if (!abstraction->is_solvable())
		return abstraction;
	    if (!other_abstraction->is_solvable())
		return other_abstraction;

	    shrink_strategy->shrink_before_merge(*abstraction, *other_abstraction);
	    // TODO: Make shrink_before_merge return a pair of bools
	    //       that tells us whether they have actually changed,
	    //       and use that to decide whether to dump statistics?
	    // (The old code would print statistics on abstraction iff it was
	    // shrunk. This is not so easy any more since this method is not
	    // in control, and the shrink strategy doesn't know whether we want
	    // expensive statistics. As a temporary aid, we just print the
	    // statistics always now, whether or not we shrunk.)
	    abstraction->statistics(use_expensive_statistics);
	    other_abstraction->statistics(use_expensive_statistics);

	    if (!reduced_labels) {
		labels->reduce(make_pair(system_one, system_two), all_abstractions);
	    }
	    abstraction->normalize();
	    other_abstraction->normalize();
	    if (!reduced_labels) {
		// only print statistics if we just possibly reduced labels
		other_abstraction->statistics(use_expensive_statistics);
		abstraction->statistics(use_expensive_statistics);
	    }
	}else{
	    abstraction->normalize();
	    other_abstraction->normalize();
	}

        Abstraction *new_abstraction = new CompositeAbstraction(labels.get(),
                                                                abstraction,
                                                                other_abstraction);

        abstraction->release_memory();
        other_abstraction->release_memory();

        new_abstraction->statistics(use_expensive_statistics);

        all_abstractions[system_one] = 0;
        all_abstractions[system_two] = 0;
        all_abstractions.push_back(new_abstraction);

	/* this can help pruning unreachable/irrelevant states before starting on label reduction
	 * problem before: label reduction ran out of memory if unreachable/irrelevant states not
	 * pruned beforehand (at least in some instances)
	 * possible downside: Too many transitions here
	 */
	new_abstraction->compute_distances();
	if (!new_abstraction->is_solvable())
	    return new_abstraction;

	if(shrink_after_merge){
            labels->reduce(make_pair(all_abstractions.size() - 1, 
				     all_abstractions.size() - 1), all_abstractions);
	    new_abstraction->normalize();
	    shrink_strategy->shrink(*new_abstraction, numeric_limits<int>::max(), true);
	    assert (new_abstraction->is_solvable());
	}
    }

    assert(all_abstractions.size() == g_variable_domain.size() * 2 - 1);
    Abstraction *res_abstraction = 0;
    for (size_t i = 0; i < all_abstractions.size(); ++i) {
        if (all_abstractions[i]) {
            if (res_abstraction) {
                cerr << "Found more than one remaining abstraction!" << endl;
                exit_with(EXIT_CRITICAL_ERROR);
            }
            res_abstraction = all_abstractions[i];
            assert(i == all_abstractions.size() - 1);
        }
    }

    res_abstraction->compute_distances();
    if (!res_abstraction->is_solvable())
        return res_abstraction;

    res_abstraction->statistics(use_expensive_statistics);
    res_abstraction->release_memory();

    return res_abstraction;
    
}

void LDSimulation::build_abstraction() {
    // TODO: We're leaking memory here in various ways. Fix this.
    //       Don't forget that build_atomic_abstractions also
    //       allocates memory.
    Timer t;
    // vector of all abstractions. entries with 0 have been merged.
    vector<Abstraction *> all_abstractions;
    all_abstractions.reserve(g_variable_domain.size() * 2 - 1);
    vector<SimulationRelation *> all_simulations;
    if (incremental_simulations)
        all_simulations.reserve(g_variable_domain.size() * 2 - 1);
    Abstraction::build_atomic_abstractions(all_abstractions, labels.get());

    // compute initial simulations, based on atomic abstractions

    unique_ptr<ShrinkStrategy> shrink_strategy;
    if(use_bisimulation){
        shrink_strategy = unique_ptr<ShrinkStrategy>(ShrinkBisimulation::create_default());
    }

    if (intermediate_simulations) {
	if(!forbid_lr){
	    cout << "Reduce labels: " << labels->get_size() << " t: " << t() << endl;
	    labels->reduce(make_pair(0, 1), all_abstractions); // With the reduction methods we use here, this should just apply label reduction on all abstractions
	    cout << "Normalize: " << t() << endl;
	    for (auto abs : all_abstractions) {
		if(abs){
		    abs->normalize();
		    abs->statistics(use_expensive_statistics);
		}
	    }
	}
        cout << "Simulation-shrinking atomic abstractions..." << endl;
        for (size_t i = 0; i < all_abstractions.size(); ++i) {
            if (all_abstractions[i]) {
                abstractions.push_back(all_abstractions[i]);
            }
        }
        // initialize simulations
        compute_ld_simulation();
        if (incremental_simulations) {
            for (int i = 0; i < simulations.size(); i++) {
                all_simulations.push_back(simulations[i]);
            }
        }
    } else if (use_bisimulation) {
	if(!forbid_lr){
	    cout << "Reduce labels: " << labels->get_size() << " t: " << t() << endl;
	    labels->reduce(make_pair(0, 1), all_abstractions); // With the reduction methods we use here, this should just apply label reduction on all abstractions
	    cout << "Normalize: " << t() << endl;
	    for (auto abs : all_abstractions) {
		if(abs){
		    abs->normalize();
		    abs->statistics(use_expensive_statistics);
		}
	    }
	}
        // do not use bisimulation shrinking on atomic abstractions if simulations are used
        cout << "Bisimulation-shrinking atomic abstractions..." << endl;
        for (size_t i = 0; i < all_abstractions.size(); ++i) {
	    if(all_abstractions[i]){
		all_abstractions[i]->compute_distances();
		// if (!all_abstractions[i]->is_solvable())
		//  return all_abstractions[i];
		shrink_strategy->shrink_atomic(*all_abstractions[i]);
	    }
        }
    }

    int remaining_abstractions = all_abstractions.size();
    remaining_abstractions -= remove_useless_abstractions(all_abstractions, all_simulations);

    cout << "Merging abstractions..." << endl;

    while (!merge_strategy->done() && t() <= limit_seconds_mas && remaining_abstractions > 1) {

        pair<int, int> next_systems = merge_strategy->get_next(all_abstractions,
                limit_absstates_merge,
                limit_transitions_merge);
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

        cout << "Merge: " << t() << endl;

        //TODO: UPDATE SIMULATION WHEN DOING INCREMENTAL COMPUTATION
        Abstraction *new_abstraction = new CompositeAbstraction(labels.get(),
                abstraction,
                other_abstraction);

        abstraction->release_memory();
        other_abstraction->release_memory();

	remaining_abstractions --;
        new_abstraction->statistics(use_expensive_statistics);

        all_abstractions[system_one] = 0;
        all_abstractions[system_two] = 0;
        all_abstractions.push_back(new_abstraction);

        // Note: we do not reduce labels several times for the same abstraction
        bool reduced_labels = false;
        if (shrink_strategy && shrink_strategy->reduce_labels_before_shrinking()) {
	    remove_dead_labels(all_abstractions);
	    if(!forbid_lr){
		cout << "Reduce labels: " << labels->get_size() << " t: " << t() << endl;
		//labels->reduce(make_pair(system_one, system_two), all_abstractions);
		if(remaining_abstractions == 1){
		    labels->reduce_to_cost();
		}else{
		    labels->reduce(make_pair(0, 1), all_abstractions); // With the reduction methods we use here, this should just apply label reduction on all abstractions
		}
		reduced_labels = true;
	    }
            cout << "Normalize: " << t() << endl;
            /*abstraction->normalize();
            other_abstraction->normalize();
            abstraction->statistics(use_expensive_statistics);
            other_abstraction->statistics(use_expensive_statistics);*/
            new_abstraction->normalize();
            new_abstraction->statistics(use_expensive_statistics);
        }

        cout << "Compute distances: " << t() << endl;
        // distances need to be computed before shrinking
        //abstraction->compute_distances();
        //other_abstraction->compute_distances();
        new_abstraction->compute_distances();
        // if (!abstraction->is_solvable())
        //     return abstraction;
        // if (!other_abstraction->is_solvable())
        //     return other_abstraction;

        if(shrink_strategy){
            /* PIET-edit: Well, we actually want to have bisim shrinking AFTER merge */
            cout << "Shrink: " << t() << endl;
            //shrink_strategy->shrink_before_merge(*abstraction, *other_abstraction);
            shrink_strategy->shrink(*new_abstraction, new_abstraction->size(), true); /* force bisimulation shrinking */
            //abstraction->statistics(use_expensive_statistics);
            //other_abstraction->statistics(use_expensive_statistics);
            new_abstraction->statistics(use_expensive_statistics);
        }

        if ((shrink_strategy || intermediate_simulations || apply_subsumed_transitions_pruning) && !reduced_labels) {
	    remove_dead_labels(all_abstractions);
	    if(!forbid_lr){	
		//labels->reduce(make_pair(system_one, system_two), all_abstractions);
		if(remaining_abstractions == 1){
		    labels->reduce_to_cost();
		}else{
		    labels->reduce(make_pair(0, 1), all_abstractions); // With the reduction methods we use here, this should just apply label reduction on all abstractions
		}
	    }
	    for(auto a : all_abstractions) if (a) a->normalize(); 
	    
	}else{	    
	    cout << "Normalize: " << t() << endl;
	    //abstraction->normalize();
	    //other_abstraction->normalize();
	    new_abstraction->normalize();
	}
        if (!reduced_labels) {
            // only print statistics if we just possibly reduced labels
            //other_abstraction->statistics(use_expensive_statistics);
            //abstraction->statistics(use_expensive_statistics);
            new_abstraction->statistics(use_expensive_statistics);
        }

        cout << "Next it: " << t() << endl;
        if(intermediate_simulations){
            /* PIET-edit: What to do in case of incremental calculations? My guess is: Same as with abstractions.
             * That is, set the simulation relations corresponding to the just merged abstractions to null pointers
             * and in compute_ld_simulation add a single new simulation relation at the end of the vector.
             */
            abstractions.clear();
            for (size_t i = 0; i < all_abstractions.size(); ++i) {
                if (all_abstractions[i]) {
                    abstractions.push_back(all_abstractions[i]);
                }
            }
            vector<SimulationRelation *>().swap(simulations);
            if (incremental_simulations) {
                for (size_t i = 0; i < all_simulations.size(); ++i) {
                    if (all_abstractions[i]) {
                        // There should be one simulation for each abstraction (plus the two simulations for the just merged abstractions)
                        simulations.push_back(all_simulations[i]);
                    }
                }
                compute_ld_simulation(true);
                // Now we can also free the memory of the simulations corresponding to the two just merged abstractions
                delete all_simulations[system_one];
                all_simulations[system_one] = 0;
                delete all_simulations[system_two];
                all_simulations[system_two] = 0;
                all_simulations.push_back(simulations.back());
            } else {
                compute_ld_simulation();
            }
        }

	remaining_abstractions -= remove_useless_abstractions(all_abstractions, all_simulations);
    }

    if (intermediate_simulations) {
        abstractions.clear();
    }
    for (size_t i = 0; i < all_abstractions.size(); ++i) {
//        if (intermediate_simulations) {
//            abstractions.clear();
//        }
        if (all_abstractions[i]) {
            all_abstractions[i]->compute_distances();
            all_abstractions[i]->statistics(use_expensive_statistics);
            abstractions.push_back(all_abstractions[i]);
            //all_abstractions[i]->release_memory();
        }
    }

    cout << "Partition: " << endl;
    for (auto a : abstractions){
        const auto & varset = a->get_varset();
        int size = 1;
        for(int v : varset){
            cout << " " << v << "(" << g_fact_names[v][0] << ")";
            size *= g_variable_domain[v];
        }
        cout << " (" << size << ")" << endl;
    }
    cout << endl;
}

void LDSimulation::compute_ld_simulation(bool incremental_step) {

    LabelMap labelMap (labels.get());

    cout << "Building LTSs and Simulation Relations" << endl;
    vector<LabelledTransitionSystem *> ltss_simple;
    vector<LTSEfficient *> ltss_efficient;
    for (auto a : abstractions){
        a->compute_distances();
        int lts_size, lts_trs;
        if(efficient_lts){
            ltss_efficient.push_back(a->get_lts_efficient(labelMap));
            lts_size= ltss_efficient.back()->size();
            lts_trs= ltss_efficient.back()->num_transitions();
        }else{
            ltss_simple.push_back(a->get_lts(labelMap));
            lts_size= ltss_simple.back()->size();
            lts_trs= ltss_simple.back()->num_transitions();
        }
        cout << "LTS built: " << lts_size << " states " 
                << lts_trs << " transitions "
                << labelMap.get_num_labels() << " num_labels"  << endl;

        if (!incremental_step) {
            if (efficient_simulation) {
                if (nold_simulation) {
                    simulations.push_back(
                            new SimulationRelationEfficientNoLD(a));
                } else {
                    simulations.push_back(new SimulationRelationEfficient(a));
                }
            } else {
                //Create initial goal-respecting relation
                simulations.push_back(new SimulationRelationSimple(a));
            }
        }
    }
    if (incremental_step) {
        // We only allow for incremental calculation with the simple simulation relation so far.
        // Enough to add the new simulation relation to the end of the vector.
        simulations.push_back(new SimulationRelationSimple(abstractions.back()));
    }

    LabelRelation label_dominance (labels.get());
    if(efficient_lts){
        compute_ld_simulation(ltss_efficient, labelMap, label_dominance, incremental_step);
    }else{
        compute_ld_simulation(ltss_simple, labelMap, label_dominance, incremental_step);
    }

    //for(int i = 0; i < simulations.size(); i++)
        //simulations[i]->dump(ltss_simple[i]->get_names());

    if (apply_subsumed_transitions_pruning) {
	Timer t;
	int lts_id = incremental_step ? simulations.size() -1 : -1;

        int num_pruned_trs = prune_subsumed_transitions(labelMap, label_dominance, ltss_simple, lts_id/*TODO: Hack lts_efficient will not work ever */);

        //if(num_pruned_trs){
        std::cout << num_pruned_trs << " transitions in the LTSs were pruned. " << t() << std::endl;
        //_labels->prune_irrelevant_labels();
    }
   
    if (apply_label_dominance_reduction) {
	set<int> dangerous_LTSs;
	//labels->reduce(make_pair(0, 1), abstractions);

	labels->reduce(labelMap, label_dominance, dangerous_LTSs);
	cout << "Labels reduced. Dangerous for: " << dangerous_LTSs.size() << endl;
	//for (auto v : dangerous_LTSs) cout << v << " ";
	//cout << endl;

	for (auto abs : abstractions) {
	    // normalize here is necessary, as otherwise
	    // compute_distances might remove more transitions than it
	    // should (e.g., in nomystery-opt11:p06)
	    abs->normalize();
	}
	if (apply_simulation_shrinking) {
	    if (incremental_step) {
		// Should be enough to just shrink the new abstraction (using the new simulation relation).
		if(!dangerous_LTSs.count(simulations.size() - 1)){
		    simulations.back()->shrink();
		}
	    } else {
		//cout << "Shrink all" << endl;
		for (int i = 0; i < simulations.size(); i++) {
		    if(!dangerous_LTSs.count(i)){
			simulations[i]->shrink();
		    }
		}
	    }       
	}

    }

    remove_dead_labels(abstractions);

    for (auto abs : abstractions) {
        // normalize here is necessary, as otherwise compute_distances might remove more transitions than it should (e.g., in nomystery-opt11:p06)
        abs->normalize();
        abs->compute_distances();
    }

    for (auto abs : abstractions) {
        abs->normalize();
        std::cout << "Final LTS: " << abs->size() << " states " << abs->total_transitions() << " transitions " << abs->get_num_nonreduced_labels() << " / " << abs->get_num_labels() << " labels still alive" << std::endl;
    }
}


//Compute a simulation from 0
// void LDSimulation::compute_ld_simulation(Labels * _labels, 
// 					 vector<Abstraction *> & _abstractions, 
// 					 vector<SimulationRelation *> & _simulations, bool no_ld) {

//     LabelMap labelMap (_labels);

//     cout << "Building LTSs and Simulation Relations" << endl;
//     vector<LabelledTransitionSystem *> ltss;
//     for (auto a : _abstractions){
//         ltss.push_back(a->get_lts(_labels));
//         cout << "LTS built: " << ltss.back()->size() << " states " << ltss.back()->num_transitions() << " transitions " << _labels->get_size() << " num_labels"  << endl;
//         //Create initial goal-respecting relation
//         _simulations.push_back(new SimulationRelationSimple(a));
//     }


//     compute_ld_simulation(_labels, ltss, _simulations, labelMap, no_ld);

//     for(int i = 0; i < _simulations.size(); i++)
//         _simulations[i]->dump(ltss[i]->get_names());
// }

//template void LDSimulation::compute_ld_simulation<LabelledTransitionSystem>
//(Labels * _labels,
//        std::vector<LabelledTransitionSystem *> & _ltss,
//        std::vector<SimulationRelation *> & _simulations,
//        const LabelMap & labelMap, bool no_ld);
//
//template void LDSimulation::compute_ld_simulation <LTSEfficient>
//(Labels * _labels,
//        std::vector<LTSEfficient *> & _ltss,
//        std::vector<SimulationRelation *> & _simulations,
//        const LabelMap & labelMap, bool no_ld);

template<typename LTS>
void LDSimulation::compute_ld_simulation(std::vector<LTS *> & _ltss,
        const LabelMap & labelMap, LabelRelation & label_dominance,
        bool incremental_step) {
    Timer t;

    if(!nold_simulation){
        label_dominance.init(_ltss, simulations, labelMap);
    }else{
        label_dominance.init_identity(_ltss.size(), labelMap);
    }
    std::cout << "Label dominance initialized: " << t() << std::endl;
    do{
        std::cout << "LDsimulation loop: ";
        //label_dominance.dump();
        if (incremental_step) {
            // Should be enough to just update the last (i.e., the new) simulation here.
            simulations.back()->update(simulations.size() - 1, _ltss.back(), label_dominance);
        } else {
            for (int i = 0; i < simulations.size(); i++){
                simulations[i]->update(i, _ltss[i], label_dominance);
                //_simulations[i]->dump(_ltss[i]->get_names());
            }
        }
        std::cout << " took " << t() << "s" << std::endl;
        //return; //PIET-edit: remove this for actual runs; just here for debugging the efficient stuff
    }while(label_dominance.update(_ltss, simulations));
    //for(int i = 0; i < _ltss.size(); i++){
    //_ltss[i]->dump();
    //	_simulations[i]->dump(_ltss[i]->get_names());
    //}
    //label_dominance.dump_equivalent();
    //label_dominance.dump_dominance();
    //exit(0);
    //}
}

// void LDSimulation::compute_ld_simulation_efficient(Labels * _labels, 
//         vector<Abstraction *> & _abstractions,
//         vector<SimulationRelation *> & _simulations, bool no_ld) {
//     cout << "Building LTSs and Simulation Relations" << endl;
//     vector<LTSEfficient *> ltss;
//     LabelMap labelMap(_labels);
//     for (auto a : _abstractions){
//         ltss.push_back(a->get_lts_efficient(labelMap));
//         cout << "LTS built: " << ltss.back()->size() << " states " << ltss.back()->num_transitions() << " transitions " << _labels->get_size() << " num_labels"  << endl;
//         //Create initial goal-respecting relation
// 	if(no_ld){
// 	    _simulations.push_back(new SimulationRelationEfficientNoLD(a));
// 	}else{
// 	    _simulations.push_back(new SimulationRelationEfficient(a));
// 	}
//     }

//     LabelRelation TMPlabel (_labels);
//     if(!no_ld)
// 	TMPlabel.init(ltss, _simulations, labelMap);
//     Timer t;
//     compute_ld_simulation(_labels, ltss, _simulations, labelMap, no_ld);
//     cout << "Time new: " << t() << endl;

//     cout << "S1;" << endl;
//     std::vector<std::vector<std::vector<bool> > > copy;
//     for(int i = 0; i < _simulations.size(); i++){
//         //_simulations[i]->dump(ltss[i]->get_names());
//         copy.push_back(_simulations[i]->get_relation());
//     }

//     vector<SimulationRelation *> sim2;
//     cout << "Building LTSs and Simulation Relations" << endl;
//     vector<LabelledTransitionSystem *> ltss2;
//     for (auto a : _abstractions){
//         ltss2.push_back(a->get_lts(labelMap));
//         cout << "LTS built: " << ltss2.back()->size() << " states " << ltss2.back()->num_transitions() << " transitions " << _labels->get_size() << " num_labels"  << endl;
//         //Create initial goal-respecting relation
//         sim2.push_back(new SimulationRelationSimple(a));
//     }

//     LabelRelation TMPlabel2 (_labels);
//     if(!no_ld)
// 	TMPlabel2.init(ltss2, sim2, labelMap);

//     Timer t2;
//     compute_ld_simulation(_labels, ltss2, sim2, labelMap, no_ld);
//     cout << "Time old: " << t2() << endl;

//     /*cout << "S2;" << endl;
//    for(int i = 0; i < sim2.size(); i++){ 
//        sim2[i]->dump(ltss[i]->get_names()); 
//        }*/
//     cout << sim2.size() << endl;
//     for(int i = 0; i < sim2.size(); i++){
//         cout << "Check LTS " << i << endl;
//         for(int s = 0; s < ltss[i]->size(); s++){
//             for(int t = 0; t < ltss[i]->size(); t++){
//                 if(copy[i][s][t] != sim2[i]->simulates(s, t)){
//                     cout << "Error: different simulations    ";
//                     cout <<ltss[i]->get_names() [t] << " <= " <<ltss[i]->get_names() [s];
//                     cout << "   Old: " << sim2[i]->simulates(s, t) << " new: " << copy[i][s][t]  << endl;
//                 }
//             }
//         }
//     }
//     cout << "Checked!" << endl;


//     for(int i = 0; i < sim2.size(); i++){
//         for(int l = 0; l < TMPlabel.get_num_labels(); ++l){
//             for(int l2 = 0; l2 < TMPlabel.get_num_labels(); ++l2){
//                 if(TMPlabel.dominates(l, l2, i) != TMPlabel2.dominates(l, l2, i)){
//                     cout << "Different label relation for " << l << " and " << l2 << " in "<< i << ": " << TMPlabel.dominates(l, l2, i)  << "    " << TMPlabel2.dominates(l, l2, i)  << endl;

//                 }
//             }
//         }
//     }

//     cout << "Checked labels!" << endl;


//     // for (int i = 0; i < ltss.size(); i++){
//     // 	 SimulationRelationEfficient s1 (_abstractions[i]);
//     // 	 SimulationRelationSimple s2 (_abstractions[i]);
//     // 	 for(int s = 0; s < ltss[i]->size(); s++){
//     // 	     for(int t = 0; t < ltss[i]->size(); t++){
//     // 		 if(s1.simulates(s, t) != s2.simulates(s, t)){
//     // 		     cout << "Error: different simulations    ";
//     // 		     cout <<ltss[i]->get_names() [t] << " <= " <<ltss[i]->get_names() [s];
//     // 		     cout << "   Old: " << s1.simulates(s, t) << " new: " << s2.simulates(s, t) << endl;
//     // 		 }
//     // 	     }
//     // 	 }
//     //  }

//     exit(0);
// }


// If lts_id = -1 (default), then prunes in all ltss. If lts_id > 0,
// prunes transitions dominated in all in all LTS, but other
// transitions are only checked for lts_id
int LDSimulation::prune_subsumed_transitions(const LabelMap & labelMap,
					     LabelRelation & label_dominance, 
					     const vector<LabelledTransitionSystem *> & ltss, 
					     int lts_id){ 
    cout << "number of transitions before pruning:" << endl;
    for (auto abs : abstractions) {
        abs->statistics(false);
    }
    int num_pruned_transitions = 0;

    //a) prune transitions of labels that are completely dominated by
    //other in all LTS
    vector <int> labels_id;
    label_dominance.get_labels_dominated_in_all(labels_id);
    for (auto abs : abstractions){
	for (int l : labels_id){
	    num_pruned_transitions += abs->prune_transitions_dominated_label_all(labelMap.get_old_id(l));
	    label_dominance.kill_label(l);
	}
    }

    //b) prune transitions dominated by noop in a transition system
    for (int l = 0; l < label_dominance.get_num_labels(); l++){
        int lts = label_dominance.get_dominated_by_noop_in(l);
        if(lts >= 0 && (lts == lts_id || lts_id == -1)){
            // the index of the LTS and its corresponding abstraction should always be the same -- be careful about
            // this in the other code!
            num_pruned_transitions += abstractions[lts]->
		prune_transitions_dominated_label_noop(lts, ltss, 
						       simulations, label_dominance, 
						       labelMap, labelMap.get_old_id(l));
        }
    }

    //c) prune transitions dominated by other transitions
    for (int lts = 0; lts < abstractions.size(); lts++) {
	if(lts_id != -1 && lts != lts_id) continue; 
        Abstraction* abs = abstractions[lts];
        const auto & is_rel_label = abs->get_relevant_labels();
        //l : Iterate over relevant labels
        for (int l = 0; l < is_rel_label.size(); ++l){
            if(!is_rel_label[l]) continue;
            int label_l = labelMap.get_id(l);
            for (int l2 = l; l2 < is_rel_label.size(); ++l2){
                if(!is_rel_label[l2]) continue;
                int label_l2 = labelMap.get_id(l2);
                //l2 : Iterate over relevant labels
                /* PIET-edit: after talking to Alvaro, it seems that the only critical case, i.e., having two labels l and l'
                 * with l >= l' and l' >= l, and two transitions, s--l->t and s--l'->t' with t >= t' and t' >= t, cannot occur
                 * if we first run simulation-shrinking. So, if we make sure that it is run before irrelevance pruning, we
                 * should have no problems here.
                 */
//                if (l != l2 && label_dominance.dominates(label_l, label_l2, lts) && label_dominance.dominates(label_l2, label_l, lts)) {
//                    cerr << "Error: two labels dominating each other in all but one LTS; this CANNOT happen!" << endl;
//                    cerr << l << "; " << l2 << "; " << label_l << "; " << label_l2 << endl;
//                    exit(1);
//                }
                if (label_dominance.dominates(label_l2, label_l, lts)
                        && label_dominance.dominates(label_l, label_l2, lts)) {
                    num_pruned_transitions +=
			abs->prune_transitions_dominated_label_equiv(lts, ltss, simulations, label_dominance, labelMap, l, l2);
                } else if (label_dominance.dominates(label_l2, label_l, lts)) {
                    num_pruned_transitions +=
                            abs->prune_transitions_dominated_label(lts, ltss, simulations, label_dominance, labelMap, l, l2);
                } else if (label_dominance.dominates(label_l, label_l2, lts)) {
                    num_pruned_transitions +=
			abs->prune_transitions_dominated_label(lts, ltss, simulations, label_dominance, labelMap,l2, l);
                }
            }
        }
    }
    remove_dead_labels(abstractions);

    return num_pruned_transitions;
}

//void LDSimulation::compute_ld_simulation_after_merge(vector<Abstraction *> & all_abstractions,
//        vector<SimulationRelation *> & all_simulations,
//        const pair<int, int> & next_systems) {
//    Abstraction *new_abstraction = all_abstractions.back();
//    vector <Abstraction *> _abstractions; //vector with only the non-null pointers
//    vector <Abstraction *> _simulations;  //vector with only the non-null simulations
//    for (auto a : all_abstractions) if(a) _abstractions.push_back(a);
//    for (auto s : all_simulations) if(s) _simulations.push_back(s);
//
//    //b) Reset simulations
//    //for(auto & sim :_ simulations){
//    //sim->reset();
//    //}
//    //c) Initialize simulation from previous simulations
//    SimulationRelation * new_sim = new SimulationRelation(new_abstraction,
//            simulations[next_systems.first],
//            simulations[next_systems.second]);
//    all_simulations.push_back(new_sim);
//    _simulations.push_back(new_sim);
//}


void LDSimulation::dump_options() const {
    merge_strategy->dump_options();
    labels->dump_options();
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

void LDSimulation::initialize() {
    Timer timer;
    cout << "Initializing simulation heuristic..." << endl;
    dump_options();
    verify_no_axioms();

    if(limit_absstates_merge > 1){
        if(use_mas){
            build_abstraction();
        }else{
            build_factored_systems();
        }
    }else{
        Abstraction::build_atomic_abstractions(abstractions, labels.get());
    }


    cout << "Computing simulation..." << endl;
    if(skip_simulation){
        for (auto a : abstractions){
            simulations.push_back(new SimulationRelationIdentity(a));
        }
    }else{
        if (intermediate_simulations) {
            // remove the intermediately built simulations, and create the final ones from scratch
            vector<SimulationRelation *>().swap(simulations);
        }
        compute_ld_simulation();
    }



    cout << endl;
    cout << "Done initializing simulation heuristic [" << timer << "]"
            << endl;

    cout << "Final abstractions: " << abstractions.size() << endl;
    for (auto abs : abstractions) {
        abs->normalize();
        const vector<int> & varset = abs->get_varset();
        cout << varset.size() << " variables " << abs->size() << " states " << abs->total_transitions() << " transitions " << abs->get_num_nonreduced_labels() << " / " << abs->get_num_labels() << " labels still alive" << endl;
        cout << "used variables:";
        for (auto var : varset) {
            cout << " " << var;
        }
        cout << endl;
    }

    int num_equi = num_equivalences();
    int num_sims = num_simulations();
    cout << "Total Simulations: " << num_sims + num_equi*2  << endl;
    cout << "Similarity equivalences: " << num_equi  << endl;
    cout << "Only Simulations: " << num_sims << endl;
    cout << "Useless vars: " << useless_vars.size() << endl;
    /*for(int i = 0; i < simulations.size(); i++){
    cout << "States after simulation: " << simulations[i]->num_states() << " " 
	 << simulations[i]->num_different_states() << endl;
	 }*/

    

    if (prune_dead_operators) {
	vector<bool> dead_labels_ops (labels->get_size(), false);
	vector<bool> dead_operators (g_operators.size(), false);
	for (auto abs : abstractions) abs->check_dead_operators(dead_labels_ops, dead_operators);
        int num_dead = 0;
        for (int i = 0; i < dead_operators.size(); i++) {
            if (dead_operators[i])
                num_dead++;
        }

	printf("Dead operators due to dead labels: %d / %d (%.2lf%%)\n",
	       num_dead, g_operators.size(),
	       ((double) num_dead / g_operators.size()) * 100);

	if(!Abstraction::store_original_operators){
	    /*cout << "Dead Operators due to dead labels: " << num_dead << " / "
	      << g_operators.size() << " i.e., "
	      << ((double) num_dead / g_operators.size() * 100) << "%"
	      << endl;*/
	    for (int i = 0; i < g_operators.size(); i++){
		if (dead_operators[i]){
		    //	cout << g_operators[i].get_name() << " is dead." << endl;
		    g_operators[i].set_dead();
		}
	    }	
	}
    }

    if (Abstraction::store_original_operators) {
        boost::dynamic_bitset<> required_operators(g_operators.size());
        for (int i = 0; i < labels->get_size(); i++) {
            if (dead_labels[i] || labels->is_label_reduced(i)) {
                continue;
            }
            boost::dynamic_bitset<> required_operators_for_label;
            for (auto abs : abstractions) {
                if (!abs->get_relevant_labels()[i])
                    continue;
                const vector<AbstractTransition> & transitions = abs->get_transitions_for_label(i);
                //cout << transitions.size() << " " << abs->get_relevant_labels()[i] << endl;
                boost::dynamic_bitset<> required_operators_for_abstraction(g_operators.size());
                for (int j = 0; j < transitions.size(); j++) {
                    required_operators_for_abstraction |= transitions[j].based_on_operators;
                }
                if (required_operators_for_label.size() == 0) {
                    required_operators_for_label = required_operators_for_abstraction;
                } else {
                    required_operators_for_label &= required_operators_for_abstraction;
                }
            }
            //cout << endl;
            required_operators |= required_operators_for_label;
        }
        printf("Dead operators detected by storing original operators: %d / %d (%.2lf%%)\n",
                g_operators.size() - required_operators.count(),
                g_operators.size(),
                ((double) g_operators.size() - required_operators.count())
                        / g_operators.size() * 100);
        /*cout << "Dead Operators detected by storing original operators: "
                << (g_operators.size() - required_operators.count()) << " / "
                << g_operators.size() << " i.e., "
                << ((double) (g_operators.size() - required_operators.count())
                        / g_operators.size() * 100) << "%" << endl;*/

	    for (int i = 0; i < g_operators.size(); i++){
		if (!required_operators[i]){
		//	cout << g_operators[i].get_name() << " is dead." << endl;
		    g_operators[i].set_dead();
		}
	    }
    }


    if(compute_final_abstraction) {
	final_abstraction = unique_ptr<Abstraction> (complete_heuristic());

	cout << "Done initializing merge and shrink heuristic [" << timer << "]"
	     << endl;
    }



    /*for (int i = 0; i < g_operators.size(); i++) {
        if (required_operators[i])
            cout << g_operators[i].get_name() << endl;
    }*/
    /*for (auto abs : abstractions) {
        abs->statistics(true);
    }*/
}

int LDSimulation::num_equivalences() const {
    int res = 0;
    for(int i = 0; i < simulations.size(); i++){
        res += simulations[i]->num_equivalences();
    }
    return res;
}


int LDSimulation::num_simulations() const {
    int res = 0;
    for(int i = 0; i < simulations.size(); i++){
        res += simulations[i]->num_simulations(true);
    }
    return res;
}

void LDSimulation::precompute_dominated_bdds(SymVariables * vars){
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominated_bdds();
    }
}

void LDSimulation::precompute_dominating_bdds(SymVariables * vars){
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominating_bdds();
    }
}


bool LDSimulation::pruned_state(const State &state) const {
    for(auto sim : simulations) {
        if(sim->pruned(state)){
            return true;
        }
    }
    return false;
}

int LDSimulation::get_cost(const State &state) const {
    int cost = 0;
    if(final_abstraction) {
	cost = final_abstraction->get_cost(state);
    }else{
	for(auto sim : simulations) {
	    int new_cost = sim->get_cost(state);
	    if (new_cost == -1) return -1;
	    cost = max (cost, new_cost);
	}
    }
    return cost;
}


BDD LDSimulation::getSimulatedBDD(SymVariables * vars, const State &state ) const{
    BDD res = vars->oneBDD();
    try{
        for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
            res *= (*it)->getSimulatedBDD(state);
        }
    }catch(BDDError e){
        return vars->zeroBDD();
    }
    return res;
}

BDD LDSimulation::getSimulatingBDD(SymVariables * vars, const State &state ) const{
    BDD res = vars->oneBDD();
    try{
        for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
            res *= (*it)->getSimulatingBDD(state);
        }
    }catch(BDDError e){
        return vars->zeroBDD();
    }

    return res;
}


double LDSimulation::get_percentage_simulations(bool ignore_equivalences) const {
    double percentage = 1;
    for (auto sim : simulations){
        percentage *= sim->get_percentage_simulations(false);
    }
    if(ignore_equivalences){
        percentage -= get_percentage_equivalences();
    } else {
        percentage -= get_percentage_equal();
    }
    return percentage;
}


double LDSimulation::get_percentage_equal() const {
    double percentage = 1;
    for (auto sim : simulations){
        percentage *= 1/(sim->num_states()*sim->num_states());
    }
    return percentage;
}


double LDSimulation::get_percentage_equivalences() const {
    double percentage = 1;
    for (auto sim : simulations){
        percentage *= sim->get_percentage_equivalences();
    }
    return percentage;
}


void LDSimulation::add_options_to_parser(OptionParser &parser){
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


    parser.add_option<bool>("efficient_simulation",
            "Use the efficient method for simulation",
            "false");

    parser.add_option<bool>("efficient_lts",
            "Use the efficient method for LTS representation",
            "false");


    parser.add_option<bool>("skip_simulation",
            "Skip doing the simulation algorithm",
            "false");

    parser.add_option<bool>("nold_simulation",
            "Perform only simulation but without label dominance",
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

    parser.add_option<bool>("compute_final_abstraction",
            "Continue the mas process after having computing the simulation relations in order to compute a MaS heuristic",
            "false");


    parser.add_option<bool>("shrink_after_merge",
                            "If true, performs the shrinking after merge instead of before",
                            "false");

    parser.add_option<ShrinkStrategy *>("shrink_strategy",
					"shrink strategy; ", 
	"shrink_bisimulation");
}


//Returns a optimized variable ordering that reorders the variables
//according to the standard causal graph criterion
void LDSimulation::getVariableOrdering(vector <int> & var_order){
    vector<vector<int> > partitions;
    vector<int> partition_var (g_variable_domain.size(), 0);
    cout << "Init partitions"<< endl;
    vector<int> partition_order;

    for(auto a : abstractions){
        const vector<int> & vs = a->get_varset();
        for(int v : vs){
            partition_var[v] = partitions.size();
        }

        partition_order.push_back(partitions.size());
        partitions.push_back(vs);
    }

    cout << "Create IG partitions"<< endl;

    //First optimize the ordering of the abstractions
    InfluenceGraph ig_partitions (partitions.size());
    for(int v = 0; v < g_variable_domain.size(); v++){
        for (int v2 : g_causal_graph->get_successors(v)){
            if(partition_var[v] != partition_var[v2]){
                ig_partitions.set_influence(partition_var[v], partition_var[v2]);
            }
        }
    }
    cout << "Optimize partitions ordering " << endl;
    ig_partitions.get_ordering(partition_order);

    cout << "Partition ordering: ";
    for(int v : partition_order) cout << v << " ";
    cout  << endl;

    vector<int> partition_begin;
    vector<int> partition_size;

    for(int i : partition_order){
        partition_begin.push_back(var_order.size());
        partition_size.push_back(partitions[i].size());
        for(int v : partitions[i]) var_order.push_back(v);
    }

    InfluenceGraph ig_vars (g_variable_domain.size());
    for(int v = 0; v < g_variable_domain.size(); v++){
        for (int v2 : g_causal_graph->get_successors(v)){
            ig_vars.set_influence(v, v2);
        }
    }

    ig_vars.optimize_variable_ordering_gamer(var_order, partition_begin, partition_size);

    var_order.insert(end(var_order), begin(useless_vars), end(useless_vars));

    cout << "Var ordering: ";
    for(int v : var_order) cout << v << " ";
    cout  << endl;
}




