#ifndef MERGE_AND_SHRINK_SIMULATION_IDENTITY_H
#define MERGE_AND_SHRINK_SIMULATION_IDENTITY_H

#include "dominance_relation.h"
#include "simulation_relation.h"
#include "abstraction.h"


// First implementation of a simulation relation. 
template <typename LR> 
class DominanceRelationIdentity : public DominanceRelationLR<LR> {

 public:
    DominanceRelationIdentity(Labels * labels) : 
    DominanceRelationLR<LR>(labels) {}

    virtual std::unique_ptr<SimulationRelation> init_simulation (Abstraction * _abs) {
	std::unique_ptr<SimulationRelation> simrel = 
	    std::unique_ptr<SimulationRelation> (new SimulationRelation(_abs));
	simrel->init_identity ();
	return simrel;
    }

    virtual std::unique_ptr<SimulationRelation> 
	init_simulation_incremental (CompositeAbstraction * _abs, 
				     const SimulationRelation & /*simrel_one*/, 
				     const SimulationRelation & /*simrel_two*/) {
	return init_simulation (static_cast<Abstraction *> (_abs));
    }
    

    virtual void update(int , 
			const LabelledTransitionSystem * , 
			const LR & , SimulationRelation & ){}

    /* virtual void update(int , const LTSComplex * , */
    /*     		const LR &, */
    /*     		SimulationRelation & ){} */

    virtual bool propagate_label_domination(int /*lts_id*/, 
					    const LabelledTransitionSystem * /*lts*/,
					    const LR & /*label_dominance*/, 
					    int /*l*/, int /*l2*/, 
					    SimulationRelation & ) const{
	return true;
    }

};

#endif
