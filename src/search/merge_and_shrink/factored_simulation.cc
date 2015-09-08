#include "factored_simulation.h"

#include "simulation_relation.h" 
#include "abstraction.h" 

using namespace std;


void FactoredSimulation::remove_useless() {
    simulations.erase(std::remove_if(begin(), end(),
				     [&](unique_ptr<SimulationRelation> & sim){
					 return sim->get_abstraction()->is_useless();
				     }), end());
    
}

double FactoredSimulation::get_percentage_simulations(bool ignore_equivalences) const {
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


double FactoredSimulation::get_percentage_equal() const {
    double percentage = 1;
    for (auto & sim : simulations){
        percentage *= 1/(sim->num_states()*sim->num_states());
    }
    return percentage;
}


double FactoredSimulation::get_percentage_equivalences() const {
    double percentage = 1;
    for (auto & sim : simulations){
        percentage *= sim->get_percentage_equivalences();
    }
    return percentage;
}





BDD FactoredSimulation::getSimulatedBDD(SymVariables * vars, const State &state ) const{
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

BDD FactoredSimulation::getSimulatingBDD(SymVariables * vars, const State &state ) const{
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

BDD FactoredSimulation::getIrrelevantStates(SymVariables * vars) const{
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

void FactoredSimulation::precompute_dominated_bdds(SymVariables * vars){
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominated_bdds();
    }
}

void FactoredSimulation::precompute_dominating_bdds(SymVariables * vars){
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominating_bdds();
    }
}

int FactoredSimulation::num_equivalences() const {
    int res = 0;
    for(int i = 0; i < simulations.size(); i++){
        res += simulations[i]->num_equivalences();
    }
    return res;
}

int FactoredSimulation::num_simulations() const {
    int res = 0;
    for(int i = 0; i < simulations.size(); i++){
        res += simulations[i]->num_simulations(true);
    }
    return res;
}


void FactoredSimulation::dump_statistics() const {
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


