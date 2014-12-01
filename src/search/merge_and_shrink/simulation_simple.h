#ifndef MERGE_AND_SHRINK_SIMULATION_SIMPLE_H
#define MERGE_AND_SHRINK_SIMULATION_SIMPLE_H

#include "simulation_relation.h"

class Labels;

// First implementation of a simulation relation. 
class SimulationRelationSimple : public SimulationRelation{

 public:
    SimulationRelationSimple(const Abstraction * _abs);

    /*
     * THIS IMPLEMENTATION IS VERY INNEFICIENT
     * ONLY TO BE USED AS A PROOF OF CONCEPT
     */
    virtual void update(int lts_id, const LabelledTransitionSystem * lts,
			const LabelRelation & label_dominance){
	update_sim(lts_id, lts, label_dominance);
    }
    virtual void update(int lts_id, const LTSEfficient * lts,
			const LabelRelation & label_dominance){
	update_sim(lts_id, lts, label_dominance);
    }

    template<typename LTS> void update_sim (int lts_id, const LTS * lts,
				   const LabelRelation & label_dominance);
};

#endif
