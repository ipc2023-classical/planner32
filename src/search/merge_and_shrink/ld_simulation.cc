#include "ld_simulation.h"

#include "abstraction.h"
#include "shrink_bisimulation.h"
#include "simulation_simple.h"
#include "simulation_identity.h"
#include "simulation_complex.h"
#include "simulation_complex_nold.h"
#include "merge_strategy.h"
#include "labelled_transition_system.h"
#include "opt_order.h"
#include "../globals.h"
#include "../causal_graph.h"
#include "label_reducer.h"
#include <boost/dynamic_bitset.hpp>
#include "label_relation.h"
#include "label_relation_identity.h"
#include "label_relation_noop.h"

using namespace std;

std::ostream & operator<<(std::ostream &os, const LabelDominanceType & type){
    switch(type){
    case LabelDominanceType::NONE: return os << "none";
    case LabelDominanceType::NOOP: return os << "noop";
    case LabelDominanceType::NORMAL: return os << "normal";
    default:
	std::cerr << "Name of LabelDominanceType not known";
	exit(-1);
    }
} 
const std::vector<std::string> LabelDominanceTypeValues {
    "NONE", "NOOP", "NORMAL"
	};

std::ostream & operator<<(std::ostream &os, const SimulationType & type){
    switch(type){
    case SimulationType::NONE: return os << "none";
    case SimulationType::SIMPLE: return os << "simple";
    case SimulationType::COMPLEX: return os << "complex";
    default:
	std::cerr << "Name of SimulationType not known";
	exit(-1);
    }
} 
const std::vector<std::string> SimulationTypeValues {
    "NONE", "SIMPLE", "COMPLEX"
	};

LDSimulation::LDSimulation(bool unit_cost, const Options &opts, OperatorCost cost_type) : 
    labels (new Labels(unit_cost, opts, cost_type)) {
}

LDSimulation::~LDSimulation(){
    //Abstractions should use shared_ptr to avoid a leak here
    for(auto abs : abstractions){
        delete abs;
    }
}

unique_ptr<DominanceRelation> LDSimulation::create_dominance_relation(SimulationType simulation_type, 
								      LabelDominanceType label_dominance_type, 
								      int switch_off_label_dominance) {
    switch(simulation_type){
    case SimulationType::NONE:
	return unique_ptr<DominanceRelation>(new DominanceRelationIdentity<LabelRelationIdentity>(labels.get()));
    case SimulationType::SIMPLE: 
	switch(label_dominance_type){
	case LabelDominanceType::NONE:
	    return unique_ptr<DominanceRelation>(new DominanceRelationSimple<LabelRelationIdentity>(labels.get())); 
	case LabelDominanceType::NOOP: 
	    return unique_ptr<DominanceRelation>(new DominanceRelationSimple<LabelRelationNoop>(labels.get())); 
	case LabelDominanceType::NORMAL: 
	    if (labels->get_size() > switch_off_label_dominance) {
		return unique_ptr<DominanceRelation>(new DominanceRelationSimple<LabelRelationNoop>(labels.get()));
	    }
	    return unique_ptr<DominanceRelation>(new DominanceRelationSimple<LabelRelation>(labels.get())); 
	}

    case SimulationType::COMPLEX: 
	switch(label_dominance_type){
	case LabelDominanceType::NONE:
	    return unique_ptr<DominanceRelation>(new DominanceRelationComplexNoLD<LabelRelationIdentity>(labels.get())); 
	case LabelDominanceType::NOOP: 
	    return unique_ptr<DominanceRelation>(new DominanceRelationComplex<LabelRelationNoop>(labels.get())); 
	case LabelDominanceType::NORMAL: 	   
	    if (labels->get_size() > switch_off_label_dominance) {
		return unique_ptr<DominanceRelation>(new DominanceRelationComplex<LabelRelationNoop>(labels.get()));
	    }
 
	    return unique_ptr<DominanceRelation>(new DominanceRelationComplex<LabelRelation>(labels.get())); 
	}
    }
    cerr << "Error: unkown type of simulation relation or label dominance" << endl;
    exit(-1);
}

void LDSimulation::init_atomic_abstractions() {
    Abstraction::build_atomic_abstractions(abstractions, labels.get());
    if(!useless_vars.empty()) remove_useless_atomic_abstractions(abstractions);
}

void LDSimulation::init_factored_systems(const std::vector<std::vector<int> > & partition_vars) {
    for (const auto & factor : partition_vars) {
        PDBAbstraction * abs_factor = new PDBAbstraction(labels.get(), factor);
        abstractions.push_back(abs_factor);
        abs_factor->normalize();
        abs_factor->compute_distances();
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

int LDSimulation::remove_useless_abstractions(vector<Abstraction *> & abstractions) {
    remove_dead_labels(abstractions);
    if(dominance_relation) dominance_relation->remove_useless();
    int removed_abstractions = 0;
    for(int i =0; i < abstractions.size(); ++i){
	if(abstractions[i] && abstractions[i]->is_useless()){
	    useless_vars.insert(end(useless_vars), begin(abstractions[i]->get_varset()), 
				end(abstractions[i]->get_varset()));
	    abstractions[i]->release_memory();
	    delete abstractions[i];
	    abstractions[i] = nullptr;
	    
	    removed_abstractions++;
	}
    }
    return removed_abstractions;			    
}


double LDSimulation::estimated_memory_MB (vector<Abstraction * > all_abstractions) const {
    double total_mem = 0;
    for (auto abs : all_abstractions) {
	if(abs) total_mem += abs->memory_estimate()/(1024.0*1024.0);
    }
    cout << "Total mem: " << total_mem << endl;
    return total_mem;
}

// Just main loop copied from merge_and_shrink heuristic 
void LDSimulation::complete_heuristic(MergeStrategy * merge_strategy, ShrinkStrategy * shrink_strategy, 
				      bool shrink_after_merge, int limit_seconds, int limit_memory_kb,
				      bool prune_dead_operators, bool use_expensive_statistics, 
				      std::vector<std::unique_ptr<Abstraction> > & res) const {
    Timer t_mas;
    cout << "Complete heuristic Initialized with " << abstractions.size() << " abstractions" << endl;
    //Insert atomic abstractions in the first g_variable_domain
    //variables, filling with nullptr. Then composite abstractions
    vector<Abstraction *> all_abstractions(g_variable_domain.size(), nullptr);
    int remaining_abstractions = 0;
    for (auto a : abstractions) {
	remaining_abstractions ++;
	if (a->get_varset().size() == 1) {
	    all_abstractions[*(a->get_varset().begin())] = a->clone();
	}else{
	    all_abstractions.push_back(a->clone());
	}
    }
    labels->reset_relevant_for(all_abstractions);

    vector<int> used_vars; 
    for(int i = 0; i < g_variable_domain.size(); i++) {
	if(!all_abstractions[i]) used_vars.push_back(i);
    }

    merge_strategy->init(all_abstractions);
    merge_strategy->remove_useless_vars (used_vars);
    
    if(abstractions.size() > 1){
	labels->reduce(make_pair(0, 1), all_abstractions); 
// With the reduction methods we use here, this should just apply label reduction on all abstractions
    }

    while (!merge_strategy->done() && remaining_abstractions > 1 && 
	   t_mas() < limit_seconds && get_peak_memory_in_kb() < limit_memory_kb ) {

	cout << endl << "Remaining: " << remaining_abstractions <<
	    " time: " << t_mas() << "/" << limit_seconds  << "s" <<  
	    " memory: " << get_peak_memory_in_kb()  << "/" << limit_memory_kb << " KB" << endl; 

	remaining_abstractions--;
        pair<int, int> next_systems = merge_strategy->get_next(all_abstractions);
        int system_one = next_systems.first;
        int system_two = next_systems.second;
	DEBUG_MAS(cout << " NEXT SYSTEMS: " << system_one <<  " " << system_two << endl;);
        assert(system_one != system_two);


        Abstraction *abstraction = all_abstractions[system_one];
        assert(abstraction);
        Abstraction *other_abstraction = all_abstractions[system_two];
        assert(other_abstraction);


	if(shrink_strategy && !shrink_after_merge){
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
	    if (!abstraction->is_solvable()) {
		exit_with(EXIT_UNSOLVABLE);
	    }
	    if (!other_abstraction->is_solvable()) {
		exit_with(EXIT_UNSOLVABLE);
	    }

	    shrink_strategy->shrink_before_merge(*abstraction, *other_abstraction);
	    // TODO: Make shrink_before_merge return a pair of bools
	    //       that tells us whether they have actually changed,
	    //       and use that to decide whether to dump statistics?
	    // (The old code would print statistics on abstraction iff it was
	    // shrunk. This is not so easy any more since this method is not
	    // in control, and the shrink strategy doesn't know whether we want
	    // expensive statistics. As a temporary aid, we just print the
	    // statistics always now, whether or not we shrunk.)
	    cout << "M1: "; abstraction->statistics(use_expensive_statistics);
	    cout << "M2: ";  other_abstraction->statistics(use_expensive_statistics);

	    if (!reduced_labels) {
		labels->reduce(make_pair(system_one, system_two), all_abstractions);
	    }
	    abstraction->normalize();
	    other_abstraction->normalize();

	    abstraction->compute_distances();
	    other_abstraction->compute_distances();

	    DEBUG_MAS(
		if (!reduced_labels) {
		    // only print statistics if we just possibly reduced labels
		    other_abstraction->statistics(use_expensive_statistics);
		    abstraction->statistics(use_expensive_statistics);
		});
	}else{
	    abstraction->normalize();
	    other_abstraction->normalize();
	}



        Abstraction *new_abstraction = new CompositeAbstraction(labels.get(),
                                                                abstraction,
                                                                other_abstraction);

        abstraction->release_memory();
        other_abstraction->release_memory();

        cout << "Merged: "; new_abstraction->statistics(use_expensive_statistics);
	
        all_abstractions[system_one] = 0;
        all_abstractions[system_two] = 0;
        all_abstractions.push_back(new_abstraction);

	/* this can help pruning unreachable/irrelevant states before starting on label reduction
	 * problem before: label reduction ran out of memory if unreachable/irrelevant states not
	 * pruned beforehand (at least in some instances)
	 * possible downside: Too many transitions here
	 */
	new_abstraction->compute_distances();
	if (!new_abstraction->is_solvable()) {
	    exit_with(EXIT_UNSOLVABLE);
	}
	
	if(shrink_strategy && shrink_after_merge){
            labels->reduce(make_pair(all_abstractions.size() - 1, 
				     all_abstractions.size() - 1), all_abstractions);
	    new_abstraction->normalize();
	    shrink_strategy->shrink(*new_abstraction, numeric_limits<int>::max(), true);
	    assert (new_abstraction->is_solvable());
	}
    }

    for (size_t i = 0; i < all_abstractions.size(); ++i) {
        if (all_abstractions[i]) {
	    all_abstractions[i]->compute_distances();

	    cout << "Final: "; all_abstractions[i]->statistics(use_expensive_statistics);
	    
	    if (!all_abstractions[i]->is_solvable()) 
		exit_with(EXIT_UNSOLVABLE);

            res.push_back(unique_ptr<Abstraction> (all_abstractions[i]));
        }
    }    

    if(prune_dead_operators) prune_dead_ops(all_abstractions);

    for (size_t i = 0; i < all_abstractions.size(); ++i) {
        if (all_abstractions[i]) {
	    all_abstractions[i]->release_memory();
	}
    }
}


int LDSimulation::remove_useless_atomic_abstractions(std::vector<Abstraction* > & abss) const {
    int total = 0;
    for (int i = 0; i < abss.size(); ++i) {
	if (abss[i]) {
	    const vector <int> & abs_varset = abss[i]->get_varset();
	    if(abs_varset.size() == 1) {
		if(std::find(begin(useless_vars), end(useless_vars), abs_varset[0]) != end(useless_vars)){
		    total ++;
		    abss[i] = nullptr;
		}
	    }
	}
    }
    return total;
}

void LDSimulation::build_abstraction(MergeStrategy * merge_strategy,  int limit_absstates_merge, 
				     int limit_transitions_merge, bool original_merge, 
				     ShrinkStrategy * shrink_strategy, bool forbid_lr, 
				     int limit_seconds, int limit_memory_kb,  
				     bool intermediate_simulations, bool incremental_simulations, 
				     SimulationType simulation_type, 
				     LabelDominanceType label_dominance_type, 
				     int switch_off_label_dominance, bool complex_lts, 
				     bool apply_subsumed_transitions_pruning, bool apply_label_dominance_reduction, 
				     bool apply_simulation_shrinking, bool preserve_all_optimal_plans, bool use_expensive_statistics ) {

    // TODO: We're leaking memory here in various ways. Fix this.
    //       Don't forget that build_atomic_abstractions also
    //       allocates memory.
    Timer t;

    int remaining_abstractions = 0;
    vector<Abstraction *> all_abstractions;
    if(abstractions.empty()) {
	// vector of all abstractions. entries with 0 have been merged.
	all_abstractions.reserve(g_variable_domain.size() * 2 - 1);
	Abstraction::build_atomic_abstractions(all_abstractions, labels.get());
	remaining_abstractions = all_abstractions.size();

	if(!useless_vars.empty()) 
	    remaining_abstractions -= remove_useless_atomic_abstractions(abstractions);
    } else {
	all_abstractions.resize(g_variable_domain.size(), nullptr);
	for (auto a : abstractions) {
	    remaining_abstractions++;
	    if (a->get_varset().size() == 1) {
		all_abstractions[*(a->get_varset().begin())] = a->clone();
	    }else{
		all_abstractions.push_back(a->clone());
	    }
	}
	labels->reset_relevant_for(all_abstractions);
	abstractions.clear();
    }

    vector<int> used_vars; 
    for(int i = 0; i < g_variable_domain.size(); i++) {
	if(!all_abstractions[i]) used_vars.push_back(i);
    }

    merge_strategy->init(all_abstractions);
    merge_strategy->remove_useless_vars (used_vars);


    // compute initial simulations, based on atomic abstractions
    if (intermediate_simulations) {
	if(!forbid_lr){
	    DEBUG_MAS(cout << "Reduce labels: " << labels->get_size() << " t: " << t() << endl;);
	    // With the reduction methods we use here, this should
	    // just apply label reduction on all abstractions
	    labels->reduce(make_pair(0, 1), all_abstractions);
	    DEBUG_MAS(cout << "Normalize: " << t() << endl;);
	    for (auto abs : all_abstractions) {
		if(abs){
		    abs->normalize();
		    DEBUG_MAS(abs->statistics(use_expensive_statistics););
		}
	    }
	}
        for (size_t i = 0; i < all_abstractions.size(); ++i) {
            if (all_abstractions[i]) {
                abstractions.push_back(all_abstractions[i]);
            }
        }

        // initialize simulations
        compute_ld_simulation(simulation_type, label_dominance_type, switch_off_label_dominance, 
			      complex_lts,
			      apply_subsumed_transitions_pruning, 
			      apply_label_dominance_reduction, 
			      apply_simulation_shrinking, preserve_all_optimal_plans);

    } else if (shrink_strategy) {
	if(!forbid_lr){
	    DEBUG_MAS(cout << "Reduce labels: " << labels->get_size() << " t: " << t() << endl;);
	    labels->reduce(make_pair(0, 1), all_abstractions); // With the reduction methods we use here, this should just apply label reduction on all abstractions
	    DEBUG_MAS(cout << "Normalize: " << t() << endl;);
	    for (auto abs : all_abstractions) {
		if(abs){
		    abs->normalize();
		    DEBUG_MAS(abs->statistics(use_expensive_statistics););
		}
	    }
	}
        // do not use bisimulation shrinking on atomic abstractions if simulations are used
        DEBUG_MAS(cout << "Bisimulation-shrinking atomic abstractions..." << endl;);
        for (size_t i = 0; i < all_abstractions.size(); ++i) {
	    if(all_abstractions[i]){
		all_abstractions[i]->compute_distances();
		if (!all_abstractions[i]->is_solvable()) {
		    exit_with(EXIT_UNSOLVABLE);
		}
		shrink_strategy->shrink_atomic(*all_abstractions[i]);
	    }
        }
    }


    remaining_abstractions -= remove_useless_abstractions(all_abstractions);

    DEBUG_MAS(cout << "Merging abstractions..." << endl;);

    merge_strategy->remove_useless_vars (useless_vars);
    while (!merge_strategy->done() && t() <= limit_seconds && 
	   get_peak_memory_in_kb() < limit_memory_kb && 
	   remaining_abstractions > 1) {
	
	cout << endl << "Remaining: " << remaining_abstractions <<
	    " time: " << t() << "/" << limit_seconds  << "s" <<  
	    " memory: " << get_peak_memory_in_kb()  << "/" << limit_memory_kb << " KB" << endl; 

        pair<int, int> next_systems = original_merge ? merge_strategy->get_next(all_abstractions) : 
	    merge_strategy->get_next(all_abstractions, limit_absstates_merge, limit_transitions_merge);
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

	if(original_merge){
	    if((limit_absstates_merge &&  
		abstraction->size() * other_abstraction->size() > limit_absstates_merge) || 
	       (limit_transitions_merge && abstraction->estimate_transitions(other_abstraction) > limit_transitions_merge)) {
		break;
	    }
	}
        DEBUG_MAS(cout << "Merge: " << t() << endl;);

        cout << "M1: "; abstraction->statistics(use_expensive_statistics);
        cout << "M2: "; other_abstraction->statistics(use_expensive_statistics);

        //TODO: Could be improved by updating simulation when doing incremental simulation
        CompositeAbstraction *new_abstraction = new CompositeAbstraction(labels.get(),
									 abstraction,
									 other_abstraction);

        abstraction->release_memory();
        other_abstraction->release_memory();

	remaining_abstractions --;
        cout << "Merged: "; new_abstraction->statistics(use_expensive_statistics);


        all_abstractions[system_one] = 0;
        all_abstractions[system_two] = 0;
        all_abstractions.push_back(new_abstraction);

        // Note: we do not reduce labels several times for the same abstraction
        bool reduced_labels = false;
        if (shrink_strategy && shrink_strategy->reduce_labels_before_shrinking()) {
	    remove_dead_labels(all_abstractions);
	    if(!forbid_lr){
		DEBUG_MAS(cout << "Reduce labels: " << labels->get_size() << " t: " << t() << endl;);
		//labels->reduce(make_pair(system_one, system_two), all_abstractions);
		if(remaining_abstractions == 1){
		    labels->reduce_to_cost();
		}else{
		    labels->reduce(make_pair(0, 1), all_abstractions); // With the reduction methods we use here, this should just apply label reduction on all abstractions
		}
		reduced_labels = true;
	    }
            DEBUG_MAS(cout << "Normalize: " << t() << endl;);
							     
            /*abstraction->normalize();
	      other_abstraction->normalize();
	      abstraction->statistics(use_expensive_statistics);
	      other_abstraction->statistics(use_expensive_statistics);*/
            new_abstraction->normalize();
            DEBUG_MAS(new_abstraction->statistics(use_expensive_statistics););
        }

        DEBUG_MAS(cout << "Compute distances: " << t() << endl;);
        // distances need to be computed before shrinking
        //abstraction->compute_distances();
        //other_abstraction->compute_distances();
        new_abstraction->compute_distances();
	if (!new_abstraction->is_solvable()){
	    exit_with(EXIT_UNSOLVABLE);
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
	    
	} else {
	    DEBUG_MAS(cout << "Normalize: " << t() << endl;);
	    //abstraction->normalize();
	    //other_abstraction->normalize();
	    new_abstraction->normalize();
	}

	if(shrink_strategy){
            /* PIET-edit: Well, we actually want to have bisim shrinking AFTER merge */
	    DEBUG_MAS(cout << "Shrink: " << t() << endl;);
            //shrink_strategy->shrink_before_merge(*abstraction, *other_abstraction);
            shrink_strategy->shrink(*new_abstraction, new_abstraction->size(), true); /* force bisimulation shrinking */

	    new_abstraction->normalize();
            //abstraction->statistics(use_expensive_statistics);
            //other_abstraction->statistics(use_expensive_statistics);
            DEBUG_MAS(new_abstraction->statistics(use_expensive_statistics););
        }


	new_abstraction->compute_distances();

        if (!reduced_labels) {
            // only print statistics if we just possibly reduced labels
            //other_abstraction->statistics(use_expensive_statistics);
            //abstraction->statistics(use_expensive_statistics);
            DEBUG_MAS(new_abstraction->statistics(use_expensive_statistics););
        }

        DEBUG_MAS(cout << "Next it: " << t() << endl;);
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
            
            if (incremental_simulations) {		
		dominance_relation->
		    init_incremental(new_abstraction, 
				     abstraction->get_simulation_relation(),
				     other_abstraction->get_simulation_relation());

	    }		
	    compute_ld_simulation(simulation_type, label_dominance_type, 
				  switch_off_label_dominance, complex_lts,
				  apply_subsumed_transitions_pruning, 
				  apply_label_dominance_reduction, 
				  apply_simulation_shrinking, preserve_all_optimal_plans,
				  incremental_simulations);
	    
        }

	remaining_abstractions -= remove_useless_abstractions(all_abstractions);
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
            DEBUG_MAS(all_abstractions[i]->statistics(use_expensive_statistics););
            abstractions.push_back(all_abstractions[i]);
	    // all_abstractions[i]->release_memory();
        }
    }

    DEBUG_MAS(cout << "Partition: " << endl;
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
	);
}


void LDSimulation::compute_ld_simulation(SimulationType simulation_type, 
					 LabelDominanceType label_dominance_type, 
					 int switch_off_label_dominance, bool complex_lts,
					 bool apply_subsumed_transitions_pruning, 
					 bool apply_label_dominance_reduction, 
					 bool apply_simulation_shrinking, bool preserve_all_optimal_plans, 
					 bool incremental_step ) {
    if(!dominance_relation) {
	dominance_relation = std::move(create_dominance_relation(simulation_type, 
								 label_dominance_type, 
								 switch_off_label_dominance));
    }
    
    LabelMap labelMap (labels.get());


    vector<LabelledTransitionSystem *> ltss_simple;
    vector<LTSComplex *> ltss_complex;
    // Generate LTSs and initialize simulation relations
    DEBUG_MSG(cout << "Building LTSs and Simulation Relations:";);
    for (auto a : abstractions){
        a->compute_distances();
	if (!a->is_solvable()){
	    exit_with(EXIT_UNSOLVABLE);
	}
        int lts_size, lts_trs;
        if(complex_lts){
            ltss_complex.push_back(a->get_lts_complex(labelMap));
            lts_size= ltss_complex.back()->size();
            lts_trs= ltss_complex.back()->num_transitions();
        }else{
            ltss_simple.push_back(a->get_lts(labelMap));
            lts_size= ltss_simple.back()->size();
            lts_trs= ltss_simple.back()->num_transitions();
        }
        DEBUG_MSG(cout << " " << lts_size << " (" << lts_trs << ")";);
    }
    DEBUG_MSG(cout << endl;);

    if (!incremental_step) {
	dominance_relation->init(abstractions);
    }

    // Algorithm to compute LD simulation 
    if(complex_lts){
        dominance_relation->compute_ld_simulation(ltss_complex, labelMap,
						  incremental_step);
    }else{
        dominance_relation->compute_ld_simulation(ltss_simple, labelMap,
						  incremental_step);
    }


    
    if (apply_subsumed_transitions_pruning) {
	Timer t;
	int lts_id = incremental_step ? dominance_relation->size() - 1 : -1;

	DEBUG_MAS(cout << "number of transitions before pruning:" << endl;
		  for (auto abs : abstractions) {
		      abs->statistics(false);
		  });
        int num_pruned_trs = dominance_relation->prune_subsumed_transitions(abstractions, labelMap, ltss_simple, lts_id/*TODO: Hack lts_complex will not work ever */, preserve_all_optimal_plans);

	remove_dead_labels(abstractions);


        if(num_pruned_trs){
	    std::cout << num_pruned_trs << " transitions pruned from LTS " <<  lts_id << ". ";
	}

        //_labels->prune_irrelevant_labels();
    }
   
    if (apply_label_dominance_reduction) {
	set<int> dangerous_LTSs;
	//labels->reduce(make_pair(0, 1), abstractions);

	labels->reduce(labelMap, *dominance_relation, dangerous_LTSs);
	DEBUG_MAS(cout << "Labels reduced. Dangerous for: " << dangerous_LTSs.size() << endl;);
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
		if(!dangerous_LTSs.count(dominance_relation->size() - 1)){
		    //TODO: HACK HACK we should refractor a bit this. 
		    dominance_relation->get_simulations().back()->shrink();
		}
	    } else {
		//cout << "Shrink all" << endl;
		for (int i = 0; i < dominance_relation->size(); i++) {
		    if(!dangerous_LTSs.count(i)){
			(*dominance_relation)[i].shrink();
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


    if(!abstractions.empty()) 
	cout << abstractions.back()->get_num_nonreduced_labels() << " / " << 
	    abstractions.back()->get_num_labels() << " labels still alive. " << endl; 
    DEBUG_MAS(cout << "Final LTSs: ";);
    for (auto abs : abstractions) {
        abs->normalize();
	DEBUG_MAS(cout << abs->size() << " (" << abs->total_transitions() << ") ";);
    }
    DEBUG_MAS(cout << endl<< endl;);
}

void LDSimulation::compute_final_simulation(SimulationType simulation_type, 
					    LabelDominanceType label_dominance_type, 
					    int switch_off_label_dominance, 
					    bool intermediate_simulations, bool complex_lts, 
					    bool apply_subsumed_transitions_pruning, 
					    bool apply_label_dominance_reduction, 
					    bool apply_simulation_shrinking, 
					    bool preserve_all_optimal_plans) {
    cout << "Computing simulation..." << endl;
    if (!dominance_relation){
	dominance_relation = std::move(create_dominance_relation(simulation_type, 
								 label_dominance_type, 
								 switch_off_label_dominance));
    } else if (intermediate_simulations) {
	// remove the intermediately built simulations, and create the final ones from scratch
	dominance_relation->clear();
    }
    compute_ld_simulation(simulation_type, label_dominance_type, switch_off_label_dominance, 
			  complex_lts,
			  apply_subsumed_transitions_pruning, 
			  apply_label_dominance_reduction, 
			  apply_simulation_shrinking, preserve_all_optimal_plans);    

    cout << endl;
    cout << "Done initializing simulation heuristic [" << g_timer << "]"
	 << endl;

    cout << "Final abstractions: " << abstractions.size() << endl;
    for (auto abs : abstractions) {
        abs->normalize();
        const vector<int> & varset = abs->get_varset();
        cout << "   " << varset.size() << " variables " << abs->size() << " states " << abs->total_transitions() << " transitions " << endl;
        DEBUG_MAS(cout << "used variables:";
        for (auto var : varset) {
            cout << " " << var;
        }
		  cout << endl;);
    }
    

    dominance_relation->dump_statistics(false);   
    if(!useless_vars.empty()) cout << "Useless vars: " << useless_vars.size() << endl;

}

void LDSimulation::prune_dead_ops (const vector<Abstraction*> & all_abstractions) const {
    vector<bool> dead_labels_ops (labels->get_size(), false);
    vector<bool> dead_operators (g_operators.size(), false);
    for (auto abs : all_abstractions) if(abs) abs->check_dead_operators(dead_labels_ops, dead_operators);
    int num_dead = 0;
    int were_dead = 0;
    for (int i = 0; i < dead_operators.size(); i++) {
	if(!g_operators[i].is_dead()) {
	    if (dead_operators[i])
		num_dead++;
	} else {
	    were_dead ++;
	}
    }

    printf("Dead operators due to dead labels: %d (new %d) / %lu (%.2lf%%)\n",
	   were_dead + num_dead, num_dead,  g_operators.size(),
	   ((double) (num_dead + were_dead) / g_operators.size()) * 100);

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
    } else {
        boost::dynamic_bitset<> required_operators(g_operators.size());
        for (int i = 0; i < labels->get_size(); i++) {
            if (dead_labels[i] || labels->is_label_reduced(i) || 
		(i < g_operators.size() && g_operators[i].is_dead())) {
                continue;
            }

            boost::dynamic_bitset<> required_operators_for_label;
	    bool irrelevant_for_all_abstractions = true;
            for (auto abs : all_abstractions) {
                if (!abs || !abs->get_relevant_labels()[i])
                    continue;
		irrelevant_for_all_abstractions = false;
                const vector<AbstractTransition> & transitions = abs->get_transitions_for_label(i);
                const auto & t_ops = abs->get_transition_ops_for_label(i);

                boost::dynamic_bitset<> required_operators_for_abstraction(g_operators.size());
                for (int j = 0; j < transitions.size(); j++) {
                    required_operators_for_abstraction |= t_ops[j];
                }
                if (required_operators_for_label.size() == 0) {

                    required_operators_for_label = required_operators_for_abstraction;
                } else {
                    required_operators_for_label &= required_operators_for_abstraction;
                }
            }
	    if(!irrelevant_for_all_abstractions)  required_operators |= required_operators_for_label;
        }
        printf("Dead operators detected by storing original operators: %lu / %lu (%.2lf%%)\n",
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
}


int LDSimulation::get_cost(const State &state) const {
    return dominance_relation->get_cost(state);
}


//Returns a optimized variable ordering that reorders the variables
//according to the standard causal graph criterion
void LDSimulation::getVariableOrdering(vector <int> & var_order){
    if(abstractions.empty()) return;
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

    // cout << "Var ordering: ";
    // for(int v : var_order) cout << v << " ";
    // cout  << endl;
}



void LDSimulation::release_memory() {
    for (size_t i = 0; i < abstractions.size(); ++i) {
        if (abstractions[i]) {
	    abstractions[i]->release_memory();
	}
    }
}



