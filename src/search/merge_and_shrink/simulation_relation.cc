#include "simulation_relation.h"

#include "label_relation.h"

#include <cstdlib>
#include <iostream>
#include <ext/slist>
#include "abstraction.h"

using namespace std;

SimulationRelation::SimulationRelation(Abstraction * _abs) : abs(_abs){
    int num_states = abs->size();
    const std::vector <bool> & goal_states = abs->get_goal_states();
    relation.resize(num_states);
    for(int i = 0; i < num_states; i++){
        relation[i].resize(num_states, true);
        if(!goal_states[i]){
            for (int j = 0; j < num_states; j++){
                if (goal_states[j]){
                    relation[i][j] = false;
                }
            }
        }
    }
    // set a pointer in the abstraction to this simulation relation
    abs->set_simulation_relation(this);
}

SimulationRelation::~SimulationRelation() {
    // make sure that the abstraction still works fine, even when this simulation relation is deleted
    abs->set_simulation_relation(0);
}

void SimulationRelation::apply_shrinking_to_table(const vector<int> & abstraction_mapping) {
    cout << "reducing simulation size from " << relation.size() << " to " << abs->size() << endl;
    int new_states = abs->size();
    vector<vector<bool> > new_relation(new_states);
    for (int i = 0; i < new_states; i++) {
        new_relation[i].resize(new_states, false);
    }
    int old_states = abstraction_mapping.size();
    for (int i = 0; i < old_states; i++) {
        int new_i = abstraction_mapping[i];
        if (new_i == Abstraction::PRUNED_STATE)
            continue;
        for (int j = 0; j < old_states; j++) {
            int new_j = abstraction_mapping[j];
            if (new_j == Abstraction::PRUNED_STATE)
                continue;
            new_relation[new_i][new_j] = relation[i][j];
        }
    }
    vector<vector<bool> >().swap(relation);
    relation.swap(new_relation);
}

void SimulationRelation::reset() {
    int num_states = abs->size();
    const std::vector <bool> & goal_states = abs->get_goal_states();
    fixed_relation.resize(num_states);
    for(int i = 0; i < num_states; i++){
        fixed_relation[i].resize(num_states, false);
        for (int j = 0; j < num_states; j++){
            if(relation[i][j]){
                fixed_relation[i][j] = true;
            }else if(goal_states[i] || !goal_states[j]){
                relation[i][j] = true;
            }
        }
    }
}

void SimulationRelation::dump(const vector<string> & names) const{ 
    cout << "SIMREL:" << endl;
    for(int j = 0; j < relation.size(); ++j){
        for(int i = 0; i < relation.size(); ++i){
            if(simulates(j, i) && i != j){
                if(simulates(i, j)){
                    if (j < i){
                        cout << names[i] << " <=> " << names[j] << endl;
                    }
                }else{
                    cout << names[i] << " <= " << names[j] << endl;
                }
            }
        }
    }
}

BDD SimulationRelation::getSimulatedBDD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return zeroBDD;
    else return dominated_bdds[absstate];
}

BDD SimulationRelation::getSimulatingBDD(const State & state) const{
    assert(!dominated_bdds.empty());
    int absstate = abs->get_abstract_state(state);
    if(absstate == -1) return zeroBDD;
    else return dominating_bdds[absstate];
}


void SimulationRelation::precompute_absstate_bdds(SymVariables * vars){
    abs->getAbsStateBDDs(vars, abs_bdds);
    zeroBDD = vars->zeroBDD();
}

void SimulationRelation::precompute_dominated_bdds(){

    for (int i = 0; i < abs->size(); i++){
        dominated_bdds.push_back(zeroBDD);
    }
    for(int i = 0; i < relation.size(); ++i){
        for(int j = 0; j < relation.size(); ++j){
            if (i == j && !simulates(i, j)){
                cerr << "Assertion error: simulation relation is not reflexive" << endl;
                exit(0);
            }
            if(simulates(i, j)){
                dominated_bdds[i] += abs_bdds[j];
            }
        }
    }
}

void SimulationRelation::precompute_dominating_bdds(){
    for (int i = 0; i < abs->size(); i++){
        dominating_bdds.push_back(zeroBDD);
    }
    for(int i = 0; i < relation.size(); ++i){
        for(int j = 0; j < relation.size(); ++j){
            if (i == j && !simulates(j, i)){
                cerr << "Assertion error: simulation relation is not reflexive" << endl;
                exit(0);
            }
            if(simulates(j, i)){
                dominating_bdds[i] += abs_bdds[j];
            }
        }
    }
}



int SimulationRelation::num_equivalences() const{
    int num = 0;
    std::vector<bool> counted (relation.size(), false);
    for(int i = 0; i < counted.size(); i++){
        if(!counted[i]){
            for(int j = i + 1; j < relation.size(); j++){
                if(similar(i, j)){
                    counted [j] = true;
                }
            }
        }else{
            num++;
        }
    }
    return num;
}

int SimulationRelation::num_simulations(bool ignore_equivalences) const{
    int res = 0;
    if(ignore_equivalences){
        std::vector<bool> counted (relation.size(), false);
        for(int i = 0; i < relation.size(); ++i){
            if(!counted[i]){
                for(int j = i+1; j < relation.size(); ++j){
                    if(similar(i, j)){
                        counted[j] = true;
                    }
                }
            }
        }
        for(int i = 0; i < relation.size(); ++i){
            if(!counted[i]){
                for(int j = i+1; j < relation.size(); ++j){
                    if(!counted[j]){
                        if(!similar(i, j) && (simulates(i, j) || simulates(j, i))){
                            res ++;
                        }
                    }
                }
            }
        }
    }else {
        for(int i = 0; i < relation.size(); ++i)
            for(int j = 0; j < relation.size(); ++j)
                if(simulates(i, j))
                    res++;
    }
    return res;
}

int SimulationRelation::num_different_states() const{
    int num = 0;
    std::vector<bool> counted (abs_bdds.size(), false);
    for(int i = 0; i < counted.size(); i++){
        if(!counted[i]){
            num++;
            for(int j = i + 1; j < relation.size(); j++){
                if(similar(i, j)){
                    counted [j] = true;
                }
            }
        }
    }
    return num;
}

//Computes the probability of selecting a random pair s, s' such
//that s simulates s'.
double SimulationRelation::get_percentage_simulations(bool ignore_equivalences) const {
    double num_sims = num_simulations (ignore_equivalences);
    double num_states = (ignore_equivalences ? num_different_states() : relation.size());
    return num_sims/(num_states*num_states);
}

//Computes the probability of selecting a random pair s, s' such that
//s is equivalent to s'.
double SimulationRelation::get_percentage_equivalences() const{
    double num_eq = 0;
    double num_states = relation.size();
    for(int i = 0; i < relation.size(); ++i)
        for(int j = 0; j < relation.size(); ++j)
            if(similar(i, j))
                num_eq++;
    return num_eq/(num_states*num_states);
}

void SimulationRelation::shrink() {
    std::vector<__gnu_cxx::slist<int> > equivRel;
    equivRel.reserve(relation.size());
    std::vector<bool> already_in (relation.size(), false);
    for (int i = 0; i < already_in.size(); i++) {
        // was already added due to being similar to some other state
        if (already_in[i])
            continue;
        already_in[i] = true;
        __gnu_cxx::slist<int> newEquiv;
        newEquiv.push_front(i);
        for (int j = i + 1; j < relation.size(); j++) {
            if (similar(i, j)) {
                already_in[j] = true;
                newEquiv.push_front(j);
            }
        }
        equivRel.push_back(newEquiv);
    }
    if (abs->size() != equivRel.size()) {
        cout << "Size for applying simulation shrinking: " << equivRel.size() << "; was: " << abs->size() << endl;
        //for (auto er : equivRel) {
        //    for (auto eq : er) {
        //        cout << eq << " ";
        //    }
        //    cout << endl;
        //}
        abs->apply_abstraction(equivRel);
        abs->normalize();

        // PIET-edit: This is now automatically handled by the actual shrinking
        /*vector<vector<bool> > newRelation(equivRel.size());
        for (int i = 0; i < newRelation.size(); i++) {
            newRelation[i].resize(equivRel.size());
            int old_i = equivRel[i].front();
            for (int j = 0; j < equivRel.size(); j++) {
                int old_j = equivRel[j].front();
                newRelation[i][j] = relation[old_i][old_j];
            }
        }
        vector<vector<bool> >().swap(relation);
        relation.swap(newRelation);*/
    } else {
        cout << "Simulation shrinking did not shrink anything" << endl;
    }
}
