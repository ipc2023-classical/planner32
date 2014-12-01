#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_EFFICIENT_NOLD_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_EFFICIENT_NOLD_H

#include <queue>
#include "simulation_efficient.h"

class SimulationRelationEfficientNoLD : public SimulationRelationEfficient{

    template<typename LTS> void init (int lts_id, const LTS * lts,
				      const LabelRelation & label_dominance, 
		   std::queue <Block *> & blocksToUpdate);
    
 public:
    SimulationRelationEfficientNoLD(const Abstraction * _abs);

    

    virtual void update(int /*lts_id*/, const LabelledTransitionSystem * /*lts*/,
			const LabelRelation & /*label_dominance*/){
	//update_sim(lts_id, lts, label_dominance);
    }
    virtual void update(int lts_id, const LTSEfficient * lts,
			const LabelRelation & label_dominance){
	update_sim(lts_id, lts, label_dominance);
    }

    template<typename LTS> void update_sim (int lts_id, const LTS * lts,
				   const LabelRelation & label_dominance);
};

#endif
