#include "simulation_identity.h"

#include "labelled_transition_system.h"

#include "label_relation.h"

#include <cstdlib>
#include <iostream>
#include "abstraction.h"

using namespace std;

SimulationRelationIdentity::SimulationRelationIdentity(const Abstraction * _abs) : SimulationRelation(_abs){
    for(int i = 0; i < relation.size(); i++){
	for(int j = 0; j < relation[i].size(); j++){
	    relation[i][j] = (i==j);
	}
    }
}
