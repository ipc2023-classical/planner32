#ifndef MERGE_AND_SHRINK_NUMERIC_DOMINANCE_RELATION_H
#define MERGE_AND_SHRINK_NUMERIC_DOMINANCE_RELATION_H

#include <vector>
#include <memory>
#include "abstraction.h" 
#include "labels.h" 
#include "numeric_simulation_relation.h"
#include "numeric_label_relation.h"

class LabelledTransitionSystem;
class LTSComplex;

/*
 * Class that represents the collection of simulation relations for a
 * factored LTS. Uses unique_ptr so that it owns the simulations and
 * it cannot be copied away.
 */
class NumericDominanceRelation {

    NumericLabelRelation label_dominance;

protected:
    std::vector<std::unique_ptr<NumericSimulationRelation> > simulations;

    std::unique_ptr<NumericSimulationRelation> init_simulation (Abstraction * _abs);

    void update(int lts_id, const LabelledTransitionSystem * lts, 
		const NumericLabelRelation & label_dominance, NumericSimulationRelation & simrel);

    bool propagate_label_domination(int lts_id, 
				    const LabelledTransitionSystem * lts,
				    int l, int l2, NumericSimulationRelation & simrel) const;

    template<typename LTS>
	void compute_ld_simulation_template(std::vector<LTS *> & _ltss, const LabelMap & labelMap) {
	assert(_ltss.size() == simulations.size());
	Timer t;

	std::cout << "Compute numLDSim on " << _ltss.size() << " LTSs." << std::endl;
	
        label_dominance.init(_ltss, *this, labelMap);
	
	std::cout << "  Init numLDSim in " << t() << "s: " << std::flush;
	do{
	    //label_dominance.dump();
		for (int i = 0; i < simulations.size(); i++){
		    update(i, _ltss[i], label_dominance, *(simulations[i]));
		    //_dominance_relation[i]->dump(_ltss[i]->get_names());
		}
	    std::cout << " " << t() << "s" << std::flush;
	}while(label_dominance.update(_ltss, *this));
	std::cout << std::endl << "LDSim computed " << t() << std::endl;
	//for(int i = 0; i < _ltss.size(); i++){
	//_ltss[i]->dump();
	//	_dominance_relation[i]->dump(_ltss[i]->get_names());
	//}
	//label_dominance.dump_equivalent();
	//label_dominance.dump_dominance();
	//exit(0);
	//}
    }

public:

   NumericDominanceRelation(Labels * labels) : label_dominance(labels) 
    {}
 
    //Methods to use the dominance relation 
    bool pruned_state(const State &state) const;
    int get_cost(const State &state) const;
    bool dominates(const State &t, const State & s) const;

    void init (const std::vector<Abstraction *> & abstractions);
    
    void compute_ld_simulation (std::vector<LabelledTransitionSystem *> & _ltss, 
				const LabelMap & labelMap, 
				bool incremental_step);   


    //Methods to access the underlying simulation relations
    const std::vector<std::unique_ptr<NumericSimulationRelation> > & get_simulations () const{
	return simulations;
    }

    int size () const {
	return simulations.size();
    }
    
    NumericSimulationRelation & operator[](int index) {
        return *(simulations[index]);
    }

    const NumericSimulationRelation & operator[](int index) const {
        return *(simulations[index]);
    }    
};


#endif
