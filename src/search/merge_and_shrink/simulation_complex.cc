#include "simulation_complex.h"

#include <queue> 
#include "../debug.h" 
#include "labelled_transition_system.h" 
#include "label_relation.h" 
#include "abstraction.h" 
#include "simulation_relation.h" 
#include "lts_complex.h" 
#include "label_relation_identity.h"
#include "label_relation_noop.h"

using namespace std;

template <typename LR> 
unique_ptr<SimulationRelation> ComplexDominanceRelation<LR>::init_simulation (Abstraction * abs) { 
    //cout << "Generate complex relation" << endl;
    int num_states = abs->size();
    //Generate a partition with two blocks, goal and non-goal (or N, based on computed_distances)
    Qp.resize(num_states, 0);
    Qp_block.resize(num_states, 0);
    Qp_pos.resize(num_states, 0);
    int goal_id =0;
    int nongoal_id = num_states -1;
    const vector <bool> & goal_states = abs->get_goal_states();
    for(int i = 0; i < num_states; i++){
	if(goal_states[i]){
	    Qp_pos[i] = goal_id;
	    Qp[goal_id++] = i;
	    Qp_block[i] = 0;
	}else{
	    Qp_pos[i] = nongoal_id;
	    Qp[nongoal_id--] = i;
	    Qp_block[i] = 1;
	}
    }
    Block * goalBlock = new Block(0, 0, goal_id - 1);
    partition.push_back(unique_ptr<Block> (goalBlock));	    
    if(nongoal_id < num_states - 1){
	Block * ngoalBlock = new Block(1, nongoal_id+1, num_states -1);
	partition.push_back(unique_ptr<Block> (ngoalBlock));
	ngoalBlock->add_rel(goalBlock);
    }

    std::unique_ptr<SimulationRelation> res (new SimulationRelation(abs));
    res->init_goal_respecting();
    return std::move (res);
}

//Splits this block between itself and a new block
unique_ptr<Block> Block::split(int index) {  
    unique_ptr<Block> newb {new Block(index, node.b, node.b + splitCount-1, rel, notRel, relCount)};
    node.b += splitCount;
    return newb;
}

void Block::dump  (const std::vector<int> & Qp){
    cout << index << " (";
    for(int i = node.b; i <= node.e; ++i) cout << " " << Qp[i]; 
    cout << ")";
}

// void Block::set_notin_pre_rel(int sl){
//     relCount[sl]--;
// }

template class ComplexDominanceRelation<LabelRelation>;
template class ComplexDominanceRelation<LabelRelationIdentity>;
template class ComplexDominanceRelation<LabelRelationNoop>;
