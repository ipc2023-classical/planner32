#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_NOLD_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_NOLD_H

#include <queue>
#include <iostream>
#include "simulation_complex.h"

template <typename LR> 
class ComplexNoLDDominanceRelation : public ComplexDominanceRelation<LR> {

    template<typename LTS> void init (int lts_id, const LTS * lts,
				      const LR & label_dominance, 
				      std::queue <Block *> & blocksToUpdate);

    /* template<typename LTS> void update_sim_nold (int lts_id, const LTS * lts, */
    /* 					    const LR & label_dominance,  */
    /* 					    SimulationRelation & simrel); */
    
    void update_sim_nold (int lts_id, const LTSComplex * lts,
			  const LR & label_dominance, 
			  SimulationRelation & simrel);
    
 public:
    ComplexNoLDDominanceRelation(Labels * labels) : ComplexDominanceRelation<LR>(labels) {}
  

    virtual void update(int /*lts_id*/, const LabelledTransitionSystem * /*lts*/,
			const LR & /*label_dominance*/, SimulationRelation & ){
	//update_sim(lts_id, lts, label_dominance);
    }
    virtual void update(int lts_id, const LTSComplex * lts,
			const LR & label_dominance,
			SimulationRelation & simrel){
	update_sim_nold(lts_id, lts, label_dominance, simrel);
    }

    virtual bool propagate_label_domination(int /*lts_id*/, 
					    const LabelledTransitionSystem * /*lts*/,
					    const LR & /*label_dominance*/, 
					    int /*l*/, int /*l2*/, SimulationRelation & /*simrel*/) const{
	std::cerr << "Error: ComputeSimulationRelationComplexNoLD::propagate_label_domination not implemented yet" << std::endl;
	std::exit(-1);
	return false;
    }
};

#endif
