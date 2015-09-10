#ifndef MERGE_AND_SHRINK_COMPUTE_SIMULATION_H
#define MERGE_AND_SHRINK_COMPUTE_SIMULATION_H

#include <vector>
#include <memory>

class Labels;
class SimulationRelation;
class LabelRelation;
class LTSComplex;
class LabelledTransitionSystem;
class Abstraction;
class CompositeAbstraction;
class FactoredSimulation; 

/* 
 * Abstract class for simulation relation algorithms
 */
class ComputeSimulationRelation {

    virtual std::unique_ptr<SimulationRelation> init_simulation (Abstraction * _abs) = 0;


    virtual std::unique_ptr<SimulationRelation> 
	init_simulation_incremental (CompositeAbstraction * _abs, 
				     const SimulationRelation & simrel_one, 
				     const SimulationRelation & simrel_two) = 0;

public:

    ComputeSimulationRelation() {}
    virtual ~ComputeSimulationRelation() {}

    void init (FactoredSimulation & simulations, 
	       const std::vector<Abstraction *> & abstractions);


    void init_simulation_incremental (FactoredSimulation & simulations, 
					  CompositeAbstraction * _abs, 
					  const SimulationRelation & simrel_one, 
					  const SimulationRelation & simrel_two);


    virtual void update(int lts_id, const LabelledTransitionSystem * lts, 
			const LabelRelation & label_dominance, 
			SimulationRelation & simrel) = 0;

    virtual void update(int lts_id, 
			const LTSComplex * lts, 
			const LabelRelation & label_dominance, 
			SimulationRelation & simrel) = 0;

    virtual bool propagate_label_domination(int lts_id, 
					    const LabelledTransitionSystem * lts,
					    const LabelRelation & label_dominance, 
					    int l, int l2, SimulationRelation & simrel) const = 0;
};

#endif
