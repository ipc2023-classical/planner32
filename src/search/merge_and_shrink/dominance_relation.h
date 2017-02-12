#ifndef MERGE_AND_SHRINK_DOMINANCE_RELATION_H
#define MERGE_AND_SHRINK_DOMINANCE_RELATION_H

#include <vector>
#include <memory>
#include "abstraction.h" 
#include "labels.h" 
#include "../sym/sym_variables.h"
#include "simulation_relation.h"

class LabelledTransitionSystem;
class LTSComplex;

/*
 * Class that represents the collection of simulation relations for a
 * factored LTS. Uses unique_ptr so that it owns the simulations and
 * it cannot be copied away.
 */
class DominanceRelation {
protected:
    std::vector<std::unique_ptr<SimulationRelation> > simulations;

    virtual std::unique_ptr<SimulationRelation> init_simulation (Abstraction * _abs) = 0;

    virtual std::unique_ptr<SimulationRelation> 
	init_simulation_incremental (CompositeAbstraction * _abs, 
				     const SimulationRelation & simrel_one, 
				     const SimulationRelation & simrel_two) = 0;


public: 
    //Methods to use the simulation 
    bool pruned_state(const State &state) const;
    int get_cost(const State &state) const;
    bool dominates(const State &t, const State & s) const;



    void init (const std::vector<Abstraction *> & abstractions);

    void init_incremental (CompositeAbstraction * _abs, 
			   const SimulationRelation & simrel_one, 
			   const SimulationRelation & simrel_two);
    
    virtual void compute_ld_simulation (std::vector<LabelledTransitionSystem *> & _ltss,
				       const LabelMap & labelMap, 
					bool incremental_step, bool dump) = 0;   

    virtual void compute_ld_simulation (std::vector<LTSComplex *> & _ltss,
				       const LabelMap & labelMap, 
				       bool incremental_step, bool dump) = 0;   


    virtual bool propagate_transition_pruning
	(int lts_id, const std::vector<LabelledTransitionSystem *> & ltss, 
	 int src, int label_id, int target) const = 0;
    
    virtual int prune_subsumed_transitions(std::vector<Abstraction *> & abstractions, 
					   const LabelMap & labelMap,
					   const std::vector<LabelledTransitionSystem *> & ltss, 
					   int lts_id, bool preserve_all_optimal_plans ) = 0;


    virtual EquivalenceRelation* get_equivalent_labels_relation
	(const LabelMap & labelMap, 
	 std::set<int> & dangerous_LTSs) const = 0;
   
    void remove_useless();


    //Statistics of the factored simulation 
    void dump_statistics (bool expensive_statistics) const;
    int num_equivalences() const;
    int num_simulations() const;   
    double num_st_pairs() const;
    double num_states_problem() const;
    
    //Computes the probability of selecting a random pair s, s' such
    //that s simulates s'.
    double get_percentage_simulations(bool ignore_equivalences) const;
    //Computes the probability of selecting a random pair s, s' such that
    //s is equivalent to s'.
    double get_percentage_equivalences() const;
    double get_percentage_equal() const;


    //Methods to obtain the BDD representation for pruning
    void precompute_dominated_bdds(SymVariables * vars);
    void precompute_dominating_bdds(SymVariables * vars);

    BDD getSimulatedBDD(SymVariables * vars, const State &state) const;
    BDD getSimulatingBDD(SymVariables * vars, const State &state) const;
    BDD getIrrelevantStates(SymVariables * vars) const;


    //Methods to access the underlying simulation relations
    const std::vector<std::unique_ptr<SimulationRelation> > & get_simulations () const{
	return simulations;
    }

    int size () const {
	return simulations.size();
    }
    
    SimulationRelation & operator[](int index) {
        return *(simulations[index]);
    }

    const SimulationRelation & operator[](int index) const {
        return *(simulations[index]);
    }

    void clear(){
	std::vector<std::unique_ptr<SimulationRelation>>().swap(simulations);
    }
};

/*
 * Class that represents the collection of simulation relations for a
 * factored LTS. Uses unique_ptr so that it owns the simulations and
 * it cannot be copied away.
 */
template <typename LR>
class DominanceRelationLR : public DominanceRelation {
    LR label_dominance;


    virtual void update(int lts_id, const LabelledTransitionSystem * lts, 
			const LR & label_dominance, 
			SimulationRelation & simrel) = 0;

    virtual void update(int lts_id, 
			const LTSComplex * lts, 
			const LR & label_dominance, 
			SimulationRelation & simrel) = 0;

    bool propagate_label_domination(int lts_id, 
				    const LabelledTransitionSystem * lts,
				    int l, int l2, SimulationRelation & simrel) const;

    template<typename LTS>
	void compute_ld_simulation_template(std::vector<LTS *> & _ltss,
					    const LabelMap & labelMap, 
					    bool incremental_step, 
					    bool dump) {
	assert(_ltss.size() == simulations.size());
	Timer t;

	/* if(!nold_simulation){ */
	/*     label_dominance.init(_ltss, *this, labelMap); */
	/* }else{ */
	/*   label_dominance.init_identity(_ltss.size(), labelMap); */
	/* } */

	int total_size = 0, max_size = 0, total_trsize = 0, max_trsize = 0;
	for (auto lts : _ltss) {
	    max_size = std::max(max_size, lts->size());
	    max_trsize = std::max(max_trsize, lts->num_transitions());
	    total_size +=  lts->size();
	    total_trsize += lts->num_transitions();
	}
	std::cout << "Compute LDSim on " << _ltss.size() << " LTSs."
		  << " Total size: " << total_size  
		  << " Total trsize: " << total_trsize  
		  << " Max size: " << max_size
		  << " Max trsize: " << max_trsize 
		  << std::endl;
	
        label_dominance.init(_ltss, *this, labelMap);
	
	std::cout << "  Init LDSim in " << t() << "s: " << std::flush;
	do{
	    //label_dominance.dump();
	    if (incremental_step) {
		// Should be enough to just update the last (i.e., the new) simulation here.
		update(simulations.size() - 1, _ltss.back(), 
		       label_dominance, *(simulations.back()));
		
	    } else {
		for (int i = 0; i < simulations.size(); i++){
		    update(i, _ltss[i], label_dominance, *(simulations[i]));
		    //_dominance_relation[i]->dump(_ltss[i]->get_names());
		}
	    }
	    std::cout << " " << t() << "s" << std::flush;
	    //return; //PIET-edit: remove this for actual runs; just here for debugging the complex stuff
	}while(label_dominance.update(_ltss, *this));
	std::cout << std::endl << "LDSim computed " << t() << std::endl;
	if(dump) {
	    for(int i = 0; i < _ltss.size(); i++){
		//_ltss[i]->dump();
		simulations[i]->dump(_ltss[i]->get_names());
	    }
	}
	//label_dominance.dump_equivalent();
	//label_dominance.dump_dominance();
	//exit(0);
	//}
    }

  public: 
DominanceRelationLR(Labels * labels) : label_dominance(labels) 
    {}

   
    virtual void compute_ld_simulation(std::vector<LabelledTransitionSystem *> & _ltss,
				       const LabelMap & labelMap, 
				       bool incremental_step, bool dump){
	compute_ld_simulation_template(_ltss, labelMap, incremental_step, dump);
    }  

    virtual void compute_ld_simulation(std::vector<LTSComplex *> & _ltss,
				       const LabelMap & labelMap, 
				       bool incremental_step, bool dump){
	compute_ld_simulation_template(_ltss, labelMap, incremental_step, dump);
    }  
 
    virtual bool propagate_transition_pruning(int lts_id,
					      const std::vector<LabelledTransitionSystem *> & ltss, 
					      int src, int label_id, int target) const{
	return label_dominance.propagate_transition_pruning(lts_id, ltss, *this, src, label_id, target);
    }


// If lts_id = -1 (default), then prunes in all ltss. If lts_id > 0,
// prunes transitions dominated in all in all LTS, but other
// transitions are only checked for lts_id
    virtual int prune_subsumed_transitions(std::vector<Abstraction *> & abstractions, 
					   const LabelMap & labelMap,
					   const std::vector<LabelledTransitionSystem *> & ltss, 
					   int lts_id, bool preserve_all_optimal_plans){
	int num_pruned_transitions = 0;
	    
	//a) prune transitions of labels that are completely dominated by
	//other in all LTS
	if(!preserve_all_optimal_plans) {
	    std::vector <int> labels_id = label_dominance.get_labels_dominated_in_all();
	    for (auto abs : abstractions){
		for (int l : labels_id){
		    num_pruned_transitions += abs->prune_transitions_dominated_label_all(labelMap.get_old_id(l));
		    label_dominance.kill_label(l);
		}
	    }
	}

	//b) prune transitions dominated by noop in a transition system
	for (int l = 0; l < label_dominance.get_num_labels(); l++){
	    int lts = label_dominance.get_dominated_by_noop_in(l);
	    if(lts >= 0 && (lts == lts_id || lts_id == -1)){
		// the index of the LTS and its corresponding abstraction should always be the same -- be careful about
		// this in the other code!
		//std::cout << "Abs pointer: " << l << " dominated by noop in " << lts << "   " << std::endl;
	    
		num_pruned_transitions += abstractions[lts]->
		    prune_transitions_dominated_label_noop(lts, ltss, 
							   *this, 
							   labelMap, labelMap.get_old_id(l));
	    }
	}

	if(!preserve_all_optimal_plans) {

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
				abs->prune_transitions_dominated_label_equiv(lts, ltss, *this, labelMap, l, l2);
			} else if (label_dominance.dominates(label_l2, label_l, lts)) {
			    num_pruned_transitions +=
				abs->prune_transitions_dominated_label(lts, ltss, *this, labelMap, l, l2);
			} else if (label_dominance.dominates(label_l, label_l2, lts)) {
			    num_pruned_transitions +=
				abs->prune_transitions_dominated_label(lts, ltss, *this, labelMap,l2, l);
			}
		    }
		}
	    }
	}
	return num_pruned_transitions;
    }
   
    virtual EquivalenceRelation* get_equivalent_labels_relation(const LabelMap & labelMap, 
								std::set<int> & dangerous_LTSs) const {
	return label_dominance.get_equivalent_labels_relation(labelMap, dangerous_LTSs);
    }    
};


#endif
