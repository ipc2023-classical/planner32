#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_NOLD_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_NOLD_H

#include <queue>
#include <iostream>
#include "simulation_complex.h"

template <typename LR> 
class DominanceRelationComplexNoLD : public DominanceRelationComplex<LR> {

    template<typename LTS> void init (int lts_id, const LTS * lts,
				      const LR & label_dominance, 
				      std::queue <Block *> & blocksToUpdate);

    /* template<typename LTS> void update_sim_nold (int lts_id, const LTS * lts, */
    /* 					    const LR & label_dominance,  */
    /* 					    SimulationRelation & simrel); */
    
    /* void update_sim_nold (int lts_id, const LTSComplex * lts, */
    /*     		  const LR & label_dominance,  */
    /*     		  SimulationRelation & simrel); */
    
 public:
    DominanceRelationComplexNoLD(Labels * labels) : DominanceRelationComplex<LR>(labels) {}
  

    virtual void update(int /*lts_id*/, const LabelledTransitionSystem * /*lts*/,
			const LR & /*label_dominance*/, SimulationRelation & ){
	//update_sim(lts_id, lts, label_dominance);
    }
    /* virtual void update(int lts_id, const LTSComplex * lts, */
    /*     		const LR & label_dominance, */
    /*     		SimulationRelation & simrel){ */
    /*     update_sim_nold(lts_id, lts, label_dominance, simrel); */
    /* } */
};

#endif
