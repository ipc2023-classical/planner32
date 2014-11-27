#include "simulation_efficient.h"

#include <queue> 

using namespace std;



void Block::add_notRel(Block * block, const LTSEfficient * lts, const vector<int> &  Qp) {
    notRel.push_back(block->node);

    //Block does not simulate me anymore, so I should decrement
    //relCount of all the transitions to block!
    for(int i = block->node.b; i <= block->node.e; ++i){
        //cout << "Decrease transitions going to: " << Qp[i] << endl;
        lts->applyPreTarget(Qp[i], [&](const LTSTransitionEfficient & t){
            int sl = lts->get_pos_qa_post(t.label, t.src);
            //cout << "Decreasing " << sl << " because of s=" << t.src << " l=" << t.label << endl;
            relCount[sl]--;
            return false;
        });
    }	  

}



//Splits this block between itself and a new block
unique_ptr<Block> Block::split(int index) {  
    unique_ptr<Block> newb {new Block(index, node.b, node.b + splitCount-1, rel, notRel, relCount)};
    node.b += splitCount;
    return newb;
}

SimulationRelationEfficient::SimulationRelationEfficient(const Abstraction * _abs) 
: SimulationRelation (_abs){
    //cout << "Generate efficient relation" << endl;
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
    //cout << "Generated efficient relation" << endl;
}


void LabelData::add_to_remove_init(const Qa & qa,
        const LabelRelation &label_dominance,
        int lts_id){
    //cout << "Add to remove " << qa.label << " " << qa.state << endl;
    for(int l = 0; l < label_dominance.get_num_labels(); ++l){

        // We dont care about labels dominated by noop (will be skipped later)
        if(!label_dominance.dominated_by_noop(l, lts_id)){
            if(label_dominance.dominates(qa.label, l, lts_id)){
                //cout << "Adding to " << l << " because " << qa.label << " dominates " << l << " in " << lts_id <<  endl;
                remove_by_label[l].insert(qa.state);
            }
        }
    }
}



void SimulationRelationEfficient::init(int lts_id, const LTSEfficient * lts,
        const LabelRelation &label_dominance,
        queue <Block *> & blocksToUpdate) {

    // cout << "Init efficient relation" << endl;
    // cout << "Qp: "; for (auto q : Qp) cout << " " << q; cout << endl;
    // cout << "Qp_block: "; for (auto q : Qp) cout << " " << Qp_block[q]; cout << endl;
    // cout << "Qp_pos: "; for (auto q : Qp_pos) cout << " " << q; cout << endl;

    //Initialize relCounts. For every qa, set the number of transitions q-l>x with l >= a
    const vector<Qa> & qas = lts->get_qa_post();

    int num_labels = label_dominance.get_num_labels();
    LabelData label_data (num_labels);

    //Initialize relCounts and l.Remove
    //cout << "Init rel counts" << endl; 

    vector<int> relCount (qas.size(), 0); 
    for (auto & qa : qas){ 
        relCount[qa.index] += qa.size();

        //If a state does not have any outgoing transition labelled
        //with l >= a and a is not dominated by noop, it cannot
        //possibly simulate any other state with a
        label_data.add_to_remove_init(qa, label_dominance, lts_id);
    }    

    //cout << "Initialized relCount: ";     for(int i : relCount) cout << i << " "; cout << endl;

    //Set relCount to all blocks (initially all they have the same)
    for(auto & b : partition) b->initRelCounts(relCount); 


    //Perform split
    for(int l = 0; l < num_labels; l++) {
        //Skip labels dominated by noop (all states have them)
        if(label_dominance.dominated_by_noop(l, lts_id)){
            //cout << "Skip label " << l << " because is dominated by noop" << endl;
            continue;
        }
        //cout << "Split wrt label " << l << ": " << g_operators[l].get_name() << endl;

        vector<pair<Block *, Block *> > splitCouples;
        vector<Block *> blocksInRemove;
        split_block(label_data.get_remove(l), splitCouples, blocksInRemove);
        for(auto c : blocksInRemove){
            for(auto & d : partition){
                if(!d->markInRemove){
                    c->rel.erase(d->index);
                }
            }
        }
        for(auto c : blocksInRemove){ //Unset the mark
            c->markInRemove = false;
        }
    }

    for (auto & c : partition){
        /*cout << " Push c: " << c->index << endl;
	c->dump(Qp); cout << " relCount: "; 
	for(int i : c->relCount) cout << i;
	cout << endl;*/

        for (auto & aux : partition){
            if (!c->in_rel(aux.get())){
                //cout << "Adding block "; aux->dump(Qp);
                //cout << " as not rel of  "; c->dump(Qp); cout << endl;
                c->add_notRel(aux.get(), lts, Qp);
            }
        }
        if(!c->notRel.empty() && unmarked(c->markS))
            blocksToUpdate.push(c.get());

        /*c->dump(Qp); cout << " relCount: ";
	for(int i : c->relCount) cout << i;
	cout << endl;*/

    }



    /*cout << "Init efficient relation, done." << endl;
    lts->dump_names();
    cout << "Qp: "; for (auto q : Qp) cout << " " << q; cout << endl;
    cout << "Qp_block: "; for (auto q : Qp) cout << " " << Qp_block[q]; cout << endl;
    cout << "Qp_pos: "; for (auto q : Qp_pos) cout << " " << q; cout << endl;
    for(auto & b : partition){
    	b->dump(Qp); cout << " has as rel: " << endl;
    	for(auto idC : b->rel){ partition[idC]->dump(Qp); cout << endl;}
	}*/
}


void  SimulationRelationEfficient::update(int lts_id, const LTSEfficient * lts,
        const LabelRelation & label_dominance){

    Timer t;
    //cout << "Update efficient relation" << endl;
    //label_dominance.dump();
    LabelData label_data (label_dominance.get_num_labels());
    set<int> alph; //TODO: Not using the data structure suggested in the paper

    queue <Block *> blocksToUpdate; 

    //TODO: if info init_label(info) else
    init(lts_id, lts, label_dominance, blocksToUpdate);
    //if(true) return;

    while (!blocksToUpdate.empty()){
        //queue <Block *> copy(blocksToUpdate);
        // cout << " To update: ";
        // while(!copy.empty()){
        //     Block * b = copy.front();
        //     copy.pop();
        //     cout << "  ";
        //     b->dump(Qp);
        // }
        // cout << endl;

        Block * b = blocksToUpdate.front();
        blocksToUpdate.pop();
        b->markS = false;

        cout << "Selected block: "; b->dump(Qp); cout << endl;
        //Block b has notRel => c does not dominate b anymore. Check
        //transitions to c because they do not simulate transitions to
        //b anymore.

        /*cout << " relCount: ";
	for(int i : b->relCount) cout << i;
	cout << endl;*/
        // for(auto & x : partition){
        //     x->dump(Qp); cout << " has as rel: " << endl;
        //     for(auto idC : x->rel){ partition[idC]->dump(Qp); cout << endl;}
        // }

        assert(alph.empty()); // assert for all labels, PreB and Remove is empty

        // for (int l = 0; l < label_dominance.num_labels(); ++l){ //For each affected label
        //     if(!label_data.get_remove(l).empty()){
        // 	cout << "Remove is not empty: " << l << endl;
        // 	for (int a : label_data.get_remove(l) ) cout << a << " ";
        // 	cout << endl;
        // 	exit(0);
        //     }

        //     if(!label_data.get_preB(l).empty()){
        // 	cout << "PreB is not empty: " << l << endl;
        // 	for (int a : label_data.get_preB(l) ) cout << a << " ";
        // 	cout << endl;
        // 	exit(0);
        //     }
        // }

        //Set of pairs <src, label> to check we'll split states
        //depending on whether they have -- l' --> B or not
        set<pair<int, int> > slToCheck;

        //Compute label.remove with a single pass through the transitions
        for (auto & t : lts->get_transitions_pre()){
            //cout << "    Check T: " << t.src << " -- " << t.label << " --> " << t.target << endl;
            // cout << " Remove: ";
            // for (int a : label_data.get_remove(t.label)) cout << a << " ";
            // cout << endl;

            /*cout << "Is " << t.target << " in not rel of block " <<
	      b->index << "?   " << b->in_notRel(Qp_pos[t.target]) << endl;*/

            if(b->in_notRel(Qp_pos[t.target])){
                //cout << t << " does not support " << b->index << " anymore" << endl;

                // For every label l <= tl, s.t. exist a transition to block B
                for(int l = 0; l < label_dominance.get_num_labels(); ++l){
                    if(label_dominance.dominates(t.label, l, lts_id) &&
                            exists_transition_to(lts, l, b)){
                        slToCheck.insert(pair<int, int> (t.src, l));
                    }
                }
            }
        }

        //cout << "Check all those states that are not better than the block" << endl;
        //For every state in not rel and every label dominated by noop
        //with a transition to B
        for(int l = 0; l < label_dominance.get_num_labels(); ++l){
            if(label_dominance.dominated_by_noop(l, lts_id)){
                //cout << "Checking label " << l << " because is noop" << endl;

                if(exists_transition_to(lts, l, b)){
                    b->applyNotRel([&](int pos){ //For every concrete state in notRel(b) (s does not simulate b anymore)
                        int s = Qp[pos];
                        slToCheck.insert(pair<int, int> (s, l));
                    });
                }
            }
        }

        for(auto & p : slToCheck){
            int s = p.first; int l = p.second;
            //cout << "Check sl: " << s << " " << l << endl;
            if(!get_in_pre_rel(b, s, l, lts_id, lts, label_dominance)){
                alph.insert(l);
                label_data.add_to_remove(l, s);
                //cout << "Add " << s << " to remove of " << l << endl;
            }
        }

        b->clear_notRel(); //We have already refined with respect to (B, B.notRel)
        //cout << "Computed remove: " << endl;
        for (auto l : alph){ //Populate a.PreB
            //for each state s in B
            for(int i = b->node.b; i <= b->node.e; ++i){
                int s = Qp[i];
                //for each transition pointing s with l
                if(!lts->hasQaPre(l, s)) continue;
                const Qa & qa = lts->get_qa_pre(l, s);
                for(int j = qa.b; j <= qa.e; ++j){
                    const LTSTransitionEfficient & t = lts->get_transitions_pre()[j];
                    label_data.add_to_preB(t);  //add predecessor to PreB
                }
            }
        }

        //cout << " Computed preB" << endl;

        for (auto l : alph){ //For each affected label
            // cout << "Check label: " << l << " (" << g_operators[l].get_name() << ") and block: "; b->dump(Qp); cout <<" that has remove: ";
            // for (int a :  label_data.get_remove(l)) cout << a << " "; cout  << " and pre b: ";
            // for (int a :  label_data.get_preB(l)) cout << a << " "; cout  << endl;

            // for (int a :  label_data.get_remove(l)){
            // 	for (int b :  label_data.get_preB(l)){
            // 	    if(a == b){
            // 		cout << a << " is in remove and preB" << endl;
            // 		exit(0);
            // 	    }
            // 	}
            // }

            vector<pair<Block *, Block *> > splitCouples;
            vector<Block *> blocksInRemove;

            split_block(label_data.get_remove(l), splitCouples, blocksInRemove);

            for(const auto & p : splitCouples){
                cout <<"Blocks: " << p.first->index << " (";
                for(int i = p.first->node.b; i <= p.first->node.e; ++i) cout << " " << Qp[i];
                cout <<") and " << p.second->index << " (";
                for(int i = p.second->node.b; i <= p.second->node.e; ++i) cout << " " << Qp[i];
                cout << ") have been splitted" << endl;

                p.first->remove_rel(p.second);
                p.first->add_notRel(p.second, lts, Qp);

                /*for(auto & x : partition){
		     x->dump(Qp); cout << " has as rel: " << endl;
		     for(auto idC : x->rel){ partition[idC]->dump(Qp); cout << endl;}
		     }*/


                if(unmarked(p.first->markS)) blocksToUpdate.push(p.first);
            }

            //cout << "processed couples splitted." << endl;
            for(auto  d : blocksInRemove) {
                d->markInRemove = false;
                for (int s : label_data.get_preB(l)){
                    Block * c = partition[Qp_block[s]].get();
                    if (c->in_rel(d)){
                        d->dump(Qp);  cout << " does not simulate "; c->dump(Qp); cout << endl;
                        c->remove_rel(d);
                        c->add_notRel(d, lts, Qp);
                        if(unmarked(c->markS)) blocksToUpdate.push(c);
                    }
                }
            }

            //cout << "processed blocks in remove." << endl;
        }
        label_data.cleanup();
        alph.clear();
    }

    cout << "Done update efficient relation: " << t() << endl;
    for(int s = 0; s < relation.size(); s++){
        Block * bs = partition[Qp_block[s]].get();
        for(int t = 0; t < relation.size(); t++){
            if(s != t && !bs->in_rel(Qp_block[t])){//t does not simulate s
                relation[t][s] = false;
            }
        }
    }
    cout << "Updated relation" << endl;
}

void SimulationRelationEfficient::split_block(const set<int> & set_remove,
        vector<pair<Block *, Block *> > & splitCouples,
        vector<Block *> & blocksInRemove){
    // cout << "SplitBlock" << endl;
    // cout << "Remove set: "; for (int i : set_remove) cout  << i; cout << endl;
    // cout << "Qp: "; for (auto q : Qp) cout << " " << q; cout << endl;
    // cout << "Qp_block: "; for (auto q : Qp) cout << " " << Qp_block[q]; cout << endl;
    // cout << "Qp_pos: "; for (auto q : Qp_pos) cout << " " << q; cout << endl;

    set<Block *> touched; //TODO: we"re not using the same data structure than in the paper
    for (int r : set_remove){
        Block * c = partition[Qp_block[r]].get();
        touched.insert(c);

        int oldpos = Qp_pos[r];
        int newpos = c->get_newpos();

        int rp = Qp[newpos];
        Qp[newpos] = r; Qp_pos[r] = newpos;
        Qp[oldpos] = rp; Qp_pos[rp] = oldpos;
    }

    for (auto c : touched){
        if (c->splitCount == c->size()){
            if(unmarked(c->markInRemove))
                blocksInRemove.push_back(c);
        } else {
            int d_id = partition.size();
            partition.push_back(move(c->split(d_id)));

            Block * d = partition[d_id].get();
            if(unmarked(d->markInRemove))
                blocksInRemove.push_back(d);

            splitCouples.push_back(pair<Block *, Block *>(c, d));

            //Update block of each state
            for(int pos = d->node.b; pos <= d->node.e; ++pos){
                Qp_block[Qp[pos]] = d_id;
            }
        }
        c->splitCount = 0;
    }

    //Update the relation for other blocks
    for(const auto & p : splitCouples){
        for(auto & e : partition){
            if (e->in_rel(p.first)){
                e->add_rel(p.second);
            }
        }
    }
}


void Block::dump  (const std::vector<int> & Qp){
    cout << index << " (";
    for(int i = node.b; i <= node.e; ++i) cout << " " << Qp[i]; 
    cout << ")";
}



// void Block::set_notin_pre_rel(int sl){
//     relCount[sl]--;
// }


//Given a transition, s->tl->x, check whether s still has a transition
//s - l -> t' with l >= tl and t' in rel(B), i. e., if for any block C
//in rel(B), and label l, C.relCount(sl) > 0
bool SimulationRelationEfficient::get_in_pre_rel(Block * b, int src, int label, int lts_id,
        const LTSEfficient * lts,
        const LabelRelation & label_dominance) const {
    //cout << "Check whether " << src << " is in preRel" << b->index << " l" << label  << endl;
    //Check whether s noop is enough For each r --l-> B s. l is
    //dominated by noop and r does not simulate B anymore, check
    //whether r has some r -l'->B with l' >= l
    if(label_dominance.dominated_by_noop(label, lts_id) 
            && b->in_rel(Qp_block[src])){
        return true; //Accept this transition because is dominated by noop
    }

    // For every label l >= tl
    for(int l = 0; l < label_dominance.get_num_labels(); ++l){
        if(label_dominance.dominates(l, label, lts_id)){
            //cout << "Checking transitions with label " << l
            // 	 << " state " << src << " hasQa? " <<  lts->hasQaPost(l, src) << endl;
            if(lts->hasQaPost(l, src)){ //Get sl, (if exists)
                int sl = lts->get_qa_post(l,src).index;
                //cout << "sl index " << sl << endl;
                //For every block C in rel(B)
                for(int idC : b->rel){
                    //cout << "Have a transition labelled with " << l << " to ";
                    //partition[idC]->dump(Qp);
                    /*cout << " >= "; b->dump(Qp); cout << "?" << endl;
		     cout << " relCount: "; 
		     for(int i : partition[idC]->relCount) cout << i;
		     cout << endl;*/

                    if(partition[idC]->hasTransitionFrom(sl)){
                        //cout << "YES" << endl;
                        return true; //Found!
                    }
                }
            }
        }
    }
    //cout << " NOT!! " << src << endl;
    return false;
}
