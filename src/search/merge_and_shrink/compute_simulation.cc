#include "compute_simulation.h"

#include "dominance_relation.h"
#include "simulation_relation.h"

using namespace std;

void ComputeSimulationRelation::init (DominanceRelation & simulations, 
				      const std::vector<Abstraction *> & abstractions){ 
    simulations.clear();
    for (auto abs : abstractions){
	simulations.simulations.push_back(init_simulation(abs)); 
    }
}



void ComputeSimulationRelation::init_simulation_incremental (DominanceRelation & simulations, 
							     CompositeAbstraction * new_abs, 
							     const SimulationRelation & simrel_one, 
							     const SimulationRelation & simrel_two){
    simulations.simulations.push_back(init_simulation_incremental(new_abs, simrel_one, simrel_two));
    

    simulations.simulations.erase(std::remove_if(begin(simulations.simulations),
						 end(simulations.simulations),
						 [&](unique_ptr<SimulationRelation> & ptr){
						     return ptr.get() == (&simrel_one) || 
							 ptr.get() == (&simrel_two);
						 }), end(simulations.simulations));
}
