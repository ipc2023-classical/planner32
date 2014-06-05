#include "simulation_relation.h"

#include "label_relation.h"

#include <cstdlib>
#include <iostream>
using namespace std;

SimulationRelation::SimulationRelation(int num_states, const vector<bool> & goal_states){
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
