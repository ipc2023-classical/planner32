#include "dominance_relation.h"

#include "simulation_relation.h" 
#include "abstraction.h" 

using namespace std;


bool DominanceRelation::propagate_transition_pruning(int lts_id, const vector<LabelledTransitionSystem *> & ltss, 
				  int src, int label_id, int target) const {
    return label_dominance.propagate_transition_pruning(lts_id, ltss, *this, src, label_id, target);
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


// If lts_id = -1 (default), then prunes in all ltss. If lts_id > 0,
// prunes transitions dominated in all in all LTS, but other
// transitions are only checked for lts_id
int DominanceRelation::prune_subsumed_transitions(vector<Abstraction *> & abstractions,
						   const LabelMap & labelMap,
						   const vector<LabelledTransitionSystem *> & ltss, 
						   int lts_id)  { 
    int num_pruned_transitions = 0;

    //a) prune transitions of labels that are completely dominated by
    //other in all LTS
    vector <int> labels_id = label_dominance.get_labels_dominated_in_all();
    for (auto abs : abstractions){
	for (int l : labels_id){
	    num_pruned_transitions += abs->prune_transitions_dominated_label_all(labelMap.get_old_id(l));
	    label_dominance.kill_label(l);
	}
    }

    //b) prune transitions dominated by noop in a transition system
    for (int l = 0; l < label_dominance.get_num_labels(); l++){
        int lts = label_dominance.get_dominated_by_noop_in(l);
        if(lts >= 0 && (lts == lts_id || lts_id == -1)){
            // the index of the LTS and its corresponding abstraction should always be the same -- be careful about
            // this in the other code!
	    cout << "Abs pointer: " << l << " dominated by noop in " << lts << "   " << endl;
	    
            num_pruned_transitions += abstractions[lts]->
		prune_transitions_dominated_label_noop(lts, ltss, 
						       *this, 
						       labelMap, labelMap.get_old_id(l));
        }
    }

    //c) prune transitions dominated by other transitions
    for (int lts = 0; lts < abstractions.size(); lts++) {
	if(lts_id != -1 && lts != lts_id) continue; 
        Abstraction* abs = abstractions[lts];
        const auto & is_rel_label = abs->get_relevant_labels();
        //l : Iterate over relevant labels
        for (int l = 0; l < is_rel_label.size(); ++l){
            if(!is_rel_label[l]) continue;
            int label_l = labelMap.get_id(l);
            for (int l2 = l; l2 < is_rel_label.size(); ++l2){
                if(!is_rel_label[l2]) continue;
                int label_l2 = labelMap.get_id(l2);
                //l2 : Iterate over relevant labels
                /* PIET-edit: after talking to Alvaro, it seems that the only critical case, i.e., having two labels l and l'
                 * with l >= l' and l' >= l, and two transitions, s--l->t and s--l'->t' with t >= t' and t' >= t, cannot occur
                 * if we first run simulation-shrinking. So, if we make sure that it is run before irrelevance pruning, we
                 * should have no problems here.
                 */
//                if (l != l2 && label_dominance.dominates(label_l, label_l2, lts) && label_dominance.dominates(label_l2, label_l, lts)) {
//                    cerr << "Error: two labels dominating each other in all but one LTS; this CANNOT happen!" << endl;
//                    cerr << l << "; " << l2 << "; " << label_l << "; " << label_l2 << endl;
//                    exit(1);
//                }
                if (label_dominance.dominates(label_l2, label_l, lts)
                        && label_dominance.dominates(label_l, label_l2, lts)) {
                    num_pruned_transitions +=
			abs->prune_transitions_dominated_label_equiv(lts, ltss, *this, labelMap, l, l2);
                } else if (label_dominance.dominates(label_l2, label_l, lts)) {
                    num_pruned_transitions +=
                            abs->prune_transitions_dominated_label(lts, ltss, *this, labelMap, l, l2);
                } else if (label_dominance.dominates(label_l, label_l2, lts)) {
                    num_pruned_transitions +=
			abs->prune_transitions_dominated_label(lts, ltss, *this, labelMap,l2, l);
                }
            }
        }
    }

    return num_pruned_transitions;
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
