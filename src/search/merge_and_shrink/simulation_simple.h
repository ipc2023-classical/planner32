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
    void update(int lts_id, 
		const LabelledTransitionSystem * lts, 
		const LabelRelation & label_dominance);

    virtual void update(int , const LTSEfficient * ,
			const LabelRelation & ){}

};

#endif
