#include "dominance_relation.h"

#include "simulation_relation.h" 
#include "abstraction.h" 
#include "labels.h" 

using namespace std;


void DominanceRelation::init (const std::vector<Abstraction *> & abstractions){ 
    simulations.clear();
    for (auto abs : abstractions){
	simulations.push_back(init_simulation(abs)); 
    }
}



void DominanceRelation::init_incremental (CompositeAbstraction * new_abs, 
					  const SimulationRelation & simrel_one, 
					  const SimulationRelation & simrel_two){
    simulations.push_back(init_simulation_incremental(new_abs, simrel_one, simrel_two));
    

    simulations.erase(std::remove_if(begin(simulations),
				     end(simulations),
				     [&](unique_ptr<SimulationRelation> & ptr){
					 return ptr.get() == (&simrel_one) || 
					     ptr.get() == (&simrel_two);
				     }), end(simulations));
}


void DominanceRelation::remove_useless() {
    simulations.erase(std::remove_if(begin(simulations), end(simulations),
				     [&](unique_ptr<SimulationRelation> & sim){
					 return sim->get_abstraction()->is_useless();
				     }), end(simulations));
    
}

double DominanceRelation::get_percentage_simulations(bool ignore_equivalences) const {
    double percentage = 1;
    for (auto & sim : simulations){
        percentage *= sim->get_percentage_simulations(false);
    }
    if(ignore_equivalences){
        percentage -= get_percentage_equivalences();
    } else {
        percentage -= get_percentage_equal();
    }
    return percentage;
}


double DominanceRelation::get_percentage_equal() const {
    double percentage = 1;
    for (auto & sim : simulations){
        percentage *= 1/(sim->num_states()*sim->num_states());
    }
    return percentage;
}


double DominanceRelation::get_percentage_equivalences() const {
    double percentage = 1;
    for (auto & sim : simulations){
        percentage *= sim->get_percentage_equivalences();
    }
    return percentage;
}


BDD DominanceRelation::getSimulatedBDD(SymVariables * vars, const State &state ) const{
    BDD res = vars->oneBDD();
    try{
        for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
            res *= (*it)->getSimulatedBDD(state);
        }
    }catch(BDDError e){
        return vars->zeroBDD();
    }
    return res;
}

BDD DominanceRelation::getSimulatingBDD(SymVariables * vars, const State &state ) const{
    BDD res = vars->oneBDD();
    try{
        for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
            res *= (*it)->getSimulatingBDD(state);
        }
    }catch(BDDError e){
        return vars->zeroBDD();
    }

    return res;
}

BDD DominanceRelation::getIrrelevantStates(SymVariables * vars) const{
    BDD res = vars->zeroBDD();
    try{
        for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
            res += (*it)->getIrrelevantStates(vars);
        }
    }catch(BDDError e){
        return vars->zeroBDD();
    }
    return res;
}

void DominanceRelation::precompute_dominated_bdds(SymVariables * vars){
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominated_bdds();
    }
}

void DominanceRelation::precompute_dominating_bdds(SymVariables * vars){
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominating_bdds();
    }
}

int DominanceRelation::num_equivalences() const {
    int res = 0;
    for(int i = 0; i < simulations.size(); i++){
        res += simulations[i]->num_equivalences();
    }
    return res;
}

int DominanceRelation::num_simulations() const {
    int res = 0;
    for(int i = 0; i < simulations.size(); i++){
        res += simulations[i]->num_simulations(true);
    }
    return res;
}


void DominanceRelation::dump_statistics() const {
    int num_equi = num_equivalences();
    int num_sims = num_simulations();
    cout << "Total Simulations: " << num_sims + num_equi*2  << endl;
    cout << "Similarity equivalences: " << num_equi  << endl;
    cout << "Only Simulations: " << num_sims << endl;
    /*for(int i = 0; i < simulations.size(); i++){
      cout << "States after simulation: " << simulations[i]->num_states() << " " 
      << simulations[i]->num_different_states() << endl;
      }*/
}



bool DominanceRelation::pruned_state(const State &state) const {
    for(auto & sim : simulations) {
        if(sim->pruned(state)){
            return true;
        }
    }
    return false;
}


int DominanceRelation::get_cost(const State &state) const{
    int cost = 0;
    for(auto & sim : simulations) {
	int new_cost = sim->get_cost(state);
	if (new_cost == -1) return -1;
	cost = max (cost, new_cost);
    }
    return cost;
}
