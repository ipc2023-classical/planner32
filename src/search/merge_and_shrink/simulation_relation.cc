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
  std::vector<bool> counted (abs_bdds.size(), false);
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

int SimulationRelation::num_simulations() const{
  int res = 0;
  std::vector<bool> counted (abs_bdds.size(), false);
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

