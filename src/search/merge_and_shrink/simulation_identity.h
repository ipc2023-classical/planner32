#ifndef MERGE_AND_SHRINK_SIMULATION_IDENTITY_H
#define MERGE_AND_SHRINK_SIMULATION_IDENTITY_H

#include "simulation_relation.h"


// First implementation of a simulation relation. 
class SimulationRelationIdentity : public SimulationRelation{

 public:
    SimulationRelationIdentity(const Abstraction * _abs);
    void update(int , 
		const LabelledTransitionSystem * , 
		const LabelRelation & ){}

    virtual void update(int , const LTSEfficient * ,
			const LabelRelation & ){}

};

#endif
