#include "simulation_relation.h"

#include "label_relation.h"

#include <cstdlib>
#include <iostream>
#include "abstraction.h"

using namespace std;

SimulationRelation::SimulationRelation(const Abstraction * _abs, int num_states, const vector<bool> & goal_states) : abs(_abs){
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
}

void SimulationRelation::dump(const vector<string> & names) const{ 
  cout << "SIMREL:" << endl;
  for(int i = 0; i < relation.size(); ++i){
    for(int j = 0; j < relation.size(); ++j){    
      if(simulates(j, i) && i != j){
	cout << names[i] << " <= " << names[j] << endl;
      }
    }
  }
}

BDD SimulationRelation::getSimulatedBDD(const State & state) const{
  return dominated_by_bdds[abs->get_abstract_state(state)];
}

void SimulationRelation::precompute_dominated_bdds(SymVariables * vars){
  abs->getAbsStateBDDs(vars, abs_bdds);
  cout << "xx" << endl;
  for (int i = 0; i < abs->size(); i++){
    dominated_by_bdds.push_back(vars->zeroBDD());
  }
  for(int i = 0; i < relation.size(); ++i){
    for(int j = 0; j < relation.size(); ++j){    
      if(simulates(i, j)){
	dominated_by_bdds[i] += abs_bdds[j];
      }
    }
  }
}

