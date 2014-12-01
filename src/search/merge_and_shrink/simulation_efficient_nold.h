#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_EFFICIENT_NOLD_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_EFFICIENT_NOLD_H

#include <queue>
#include "simulation_efficient.h"

class SimulationRelationEfficientNoLD : public SimulationRelationEfficient{
 public:
    SimulationRelationEfficientNoLD(const Abstraction * _abs);

    virtual void init (int lts_id, const LTSEfficient * lts,
		   const LabelRelation & label_dominance, 
		   std::queue <Block *> & blocksToUpdate);

    virtual void update(int lts_id, const LTSEfficient * lts,
			const LabelRelation & label_dominance);
};

#endif
