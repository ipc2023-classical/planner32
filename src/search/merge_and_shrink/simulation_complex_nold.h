#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_NOLD_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_NOLD_H

#include <queue>
#include "simulation_complex.h"

class SimulationRelationComplexNoLD : public SimulationRelationComplex{

    template<typename LTS> void init (int lts_id, const LTS * lts,
				      const LabelRelation & label_dominance, 
		   std::queue <Block *> & blocksToUpdate);
    
 public:
    SimulationRelationComplexNoLD(Abstraction * _abs);

    

    virtual void update(int /*lts_id*/, const LabelledTransitionSystem * /*lts*/,
			const LabelRelation & /*label_dominance*/){
	//update_sim(lts_id, lts, label_dominance);
    }
    virtual void update(int lts_id, const LTSComplex * lts,
			const LabelRelation & label_dominance){
	update_sim(lts_id, lts, label_dominance);
    }

    template<typename LTS> void update_sim (int lts_id, const LTS * lts,
				   const LabelRelation & label_dominance);
};

#endif
