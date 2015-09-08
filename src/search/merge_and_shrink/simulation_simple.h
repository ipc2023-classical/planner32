#ifndef MERGE_AND_SHRINK_SIMULATION_SIMPLE_H
#define MERGE_AND_SHRINK_SIMULATION_SIMPLE_H

#include "compute_simulation.h"

class Labels;
 
/*
 * First implementation of a simulation relation.
 * THIS IMPLEMENTATION IS VERY INNEFICIENT
 * ONLY TO BE USED AS A PROOF OF CONCEPT
 */
class ComputeSimulationRelationSimple : public ComputeSimulationRelation{
    template<typename LTS> void update_sim (int lts_id, const LTS * lts,
					    const LabelRelation & label_dominance, 
					    SimulationRelation & simrel);
 public:
  ComputeSimulationRelationSimple() : ComputeSimulationRelation() {}

    virtual std::unique_ptr<SimulationRelation> init_simulation (Abstraction * _abs);

    virtual std::unique_ptr<SimulationRelation> 
	init_simulation_incremental (CompositeAbstraction * _abs, 
				     const SimulationRelation & simrel_one, 
				     const SimulationRelation & simrel_two);
 
   virtual void update(int lts_id, const LabelledTransitionSystem * lts,
		       const LabelRelation & label_dominance,
		       SimulationRelation & simrel){
       update_sim(lts_id, lts, label_dominance, simrel);
    }
    virtual void update(int lts_id, const LTSComplex * lts,
			const LabelRelation & label_dominance,
			SimulationRelation & simrel){
	update_sim(lts_id, lts, label_dominance, simrel);
    }

    bool propagate_label_domination(int lts_id, const LabelledTransitionSystem * lts,
				    const LabelRelation & label_dominance, 
				    int l, int l2, SimulationRelation & simrel) const;    
};

#endif
