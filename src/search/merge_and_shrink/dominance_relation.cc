#include "dominance_relation.h"

#include "simulation_relation.h" 
#include "abstraction.h" 
#include "labels.h" 
#include "labelled_transition_system.h" 

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
    Timer t;
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominated_bdds();
    }
    cout << "Precomputed dominated BDDs: " << t() << endl;
}

void DominanceRelation::precompute_dominating_bdds(SymVariables * vars){
    Timer t;
    for(auto & sim : simulations){
        sim->precompute_absstate_bdds(vars);
        sim->precompute_dominating_bdds();
    }
    cout << "Precomputed dominating BDDs: " << t() << endl;
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

double DominanceRelation::num_st_pairs() const {
    double res = 1;
    for(int i = 0; i < simulations.size(); i++){
        res *= simulations[i]->num_simulations(false);
    }
    return res;
}


double DominanceRelation::num_states_problem() const {
    double res = 1;
    for(int i = 0; i < simulations.size(); i++){
        res *= simulations[i]->num_states();
    }
    return res;
}


void DominanceRelation::dump_statistics(bool expensive) const {
    int num_equi = num_equivalences();
    int num_sims = num_simulations();

    int num_vars = 0;
    int num_vars_with_simulations = 0;
    for(int i = 0; i < simulations.size(); i++){
        if(simulations[i]->num_simulations(true) > 0){
            num_vars_with_simulations ++;
        }
        num_vars++;
    }

    cout << "Total Simulations: " << num_sims + num_equi*2  << endl;
    cout << "Similarity equivalences: " << num_equi  << endl;
    cout << "Only Simulations: " << num_sims << endl;
    cout << "Simulations Found in " << num_vars_with_simulations << " out of " << num_vars << " variables" << endl;
    
    if(expensive){
	double num_pairs = num_st_pairs();
	double problem_size = num_states_problem();
    
	cout << "Total st pairs: " << num_pairs  << endl;
	cout << "Percentage st pairs: " << num_pairs/(problem_size*problem_size)  << endl;
    }
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


bool DominanceRelation::dominates(const State &t, const State & s) const {
    for(auto & sim : simulations) {
	if (!sim->simulates(t, s)) {
	    return false;
	}
    }
    return true;
}


//l does not dominate l2 anymore, check if this changes the simulation relation
template <typename LR> 
bool DominanceRelationLR<LR>::propagate_label_domination(int lts_id, 
						       const LabelledTransitionSystem * lts,
						       int l, int l2, 
						       SimulationRelation & simrel) const {
    for (int s = 0; s < lts->size(); s++) {
	for (int t = 0; t < lts->size(); t++) { //for each pair of states t, s
	    if (s != t && simrel.simulates(t, s)) {
		//Check if really t simulates s //for each transition s--l->s':
		// a) with noop t >= s' and l dominated by noop?
		// b) exist t--l'-->t', t' >= s' and l dominated by l'?
		bool not_simulates_anymore = lts->applyPostSrc(s, [&](const LTSTransition & trs) {
			const vector<int> & labels_trs = lts->get_labels (trs.label_group);
			vector<bool> is_label_simulated (labels_trs.size(), false);
			int num_labels_trs_simulated = 0;

			if(lts->get_group_label(l2) == trs.label_group) {
			    if(labels_trs.size() == 1) {
				return false;
			    }
			    for(size_t i = 0; i <  labels_trs.size(); ++i) {
				if (labels_trs[i] == l2) {
				    is_label_simulated[i] = true;
				    num_labels_trs_simulated ++;				    
				}
			    }
			}
			    //cout << "Dominated by noop!" << endl;
			if(simrel.simulates (t, trs.target)) {
			    for(size_t i = 0; i <  labels_trs.size(); ++i) {
				if(label_dominance.dominated_by_noop(labels_trs[i], lts_id)) {
				    is_label_simulated [i] =  true;
				    num_labels_trs_simulated ++;
				}
			    }
			}

			if(num_labels_trs_simulated == labels_trs.size()) {
			    return false;
			}
			bool all_simulated = lts->applyPostSrc(t,[&](const LTSTransition & trt) {
 				if(simrel.simulates(trt.target, trs.target)) {
				    const vector<int> & labels_trt = lts->get_labels (trt.label_group);
				    for(size_t i = 0; i <  labels_trs.size(); ++i) {
					if(!is_label_simulated[i]) {				
					    for (int label_trt : labels_trt) {
						if(label_dominance.dominates(label_trt, labels_trs[i], lts_id)) {
						    is_label_simulated [i] =  true;
						    if (++num_labels_trs_simulated == labels_trs.size()) {
					return true;
				    }
						    break;
						}
					    }
					}
				    }					    
				}
				assert(num_labels_trs_simulated < labels_trs.size());
				    return false;
				});
			
			return !all_simulated;
		    });

		if(not_simulates_anymore) return false;
	    }
	}
    }
    return true;    
}
