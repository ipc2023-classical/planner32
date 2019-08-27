#include "simulation_complex_nold.h"

#include <queue> 
#include "../debug.h" 
#include "label_relation.h" 
#include "label_relation_identity.h"
#include "simulation_relation.h" 
#include "label_relation_noop.h"

using namespace std;

template <typename LR>  template <typename LTS> 
void DominanceRelationComplexNoLD<LR>::init(int /*lts_id*/, const LTS * lts,
        const LR &label_dominance,
        queue <Block *> & blocksToUpdate) {

    //Initialize relCounts. For every qa, set the number of transitions q-l>x with l >= a
    const vector<Qa> & qas = lts->get_qa_post();

    int num_labels = label_dominance.get_num_labels();
    LabelData<LR> label_data (num_labels);

    //Initialize relCounts and l.Remove
    vector<int> relCount (qas.size(), 0); 
    for (auto & qa : qas){ 
        relCount[qa.index] += qa.size();

        //If a state does not have any outgoing transition labelled
        //with l >= a and a is not dominated by noop, it cannot
        //possibly simulate any other state with a
        label_data.add_to_remove(qa.label, qa.state);
    }    

    //cout << "Initialized relCount: ";     for(int i : relCount) cout << i << " "; cout << endl;

    //Set relCount to all blocks (initially all they have the same)
    for(auto & b : DominanceRelationComplex<LR>::partition) b->initRelCounts(relCount); 

    //Perform split
    for(int l = 0; l < num_labels; l++) {
        vector<pair<Block *, Block *> > splitCouples;
        vector<Block *> blocksInRemove;
        DominanceRelationComplex<LR>::split_block(label_data.get_remove(l), splitCouples, blocksInRemove);
        for(auto c : blocksInRemove){
            for(auto & d : DominanceRelationComplex<LR>::partition){
                if(!d->markInRemove){
                    c->rel.erase(d->index);
                }
            }
        }
        for(auto c : blocksInRemove){ //Unset the mark
             c->markInRemove = false;
        }
    }

    for (auto & c : DominanceRelationComplex<LR>::partition){
        /*cout << " Push c: " << c->index << endl;
	c->dump(Qp); cout << " relCount: "; 
	for(int i : c->relCount) cout << i;
	cout << endl;*/

        for (auto & aux : DominanceRelationComplex<LR>::partition){
            if (!c->in_rel(aux.get())){
                //cout << "Adding block "; aux->dump(Qp);
                //cout << " as not rel of  "; c->dump(Qp); cout << endl;
                c->add_notRel(aux.get(), lts, DominanceRelationComplex<LR>::Qp);
            }
        }
        if(!c->notRel.empty() && DominanceRelationComplex<LR>::unmarked(c->markS))
            blocksToUpdate.push(c.get());
        /*c->dump(Qp); cout << " relCount: ";
	for(int i : c->relCount) cout << i;
	cout << endl;*/
    }

    /*cout << "Init complex relation, done." << endl;
    lts->dump_names();
    cout << "Qp: "; for (auto q : Qp) cout << " " << q; cout << endl;
    cout << "Qp_block: "; for (auto q : Qp) cout << " " << Qp_block[q]; cout << endl;
    cout << "Qp_pos: "; for (auto q : Qp_pos) cout << " " << q; cout << endl;
    for(auto & b : partition){
    	b->dump(Qp); cout << " has as rel: " << endl;
    	for(auto idC : b->rel){ partition[idC]->dump(Qp); cout << endl;}
	}*/
}

// template <typename LR> 
// void DominanceRelationComplexNoLD<LR>::update_sim_nold (int lts_id, const LTSComplex * lts,
// 							    const LR & label_dominance, 
// 							    SimulationRelation & simrel) {
//     Timer t;
//     //cout << "Update complex relation" << endl;
//     //label_dominance.dump();
//     LabelData<LR> label_data (label_dominance.get_num_labels());
//     set<int> alph; //TODO: Not using the data structure suggested in the paper

//     queue <Block *> blocksToUpdate; 

//     //TODO: if info init_label(info) else
//     init(lts_id, lts, label_dominance, blocksToUpdate);

//     while (!blocksToUpdate.empty()){
//         //queue <Block *> copy(blocksToUpdate);
//         // cout << " To update: ";
//         // while(!copy.empty()){
//         //     Block * b = copy.front();
//         //     copy.pop();
//         //     cout << "  ";
//         //     b->dump(Qp);
//         // }
//         // cout << endl;

//         Block * b = blocksToUpdate.front();
//         blocksToUpdate.pop();
//         b->markS = false;

//         DEBUG_MSG(cout << "Selected block: "; b->dump(DominanceRelationComplex<LR>::Qp); cout << endl;);
//         //Block b has notRel => c does not dominate b anymore. Check
//         //transitions to c because they do not simulate transitions to
//         //b anymore.

//         /*cout << " relCount: ";
// 	for(int i : b->relCount) cout << i;
// 	cout << endl;*/
//         // for(auto & x : partition){
//         //     x->dump(Qp); cout << " has as rel: " << endl;
//         //     for(auto idC : x->rel){ partition[idC]->dump(Qp); cout << endl;}
//         // }

//         assert(alph.empty()); // assert for all labels, PreB and Remove is empty

//         // for (int l = 0; l < label_dominance.num_labels(); ++l){ //For each affected label
//         //     if(!label_data.get_remove(l).empty()){
//         // 	cout << "Remove is not empty: " << l << endl;
//         // 	for (int a : label_data.get_remove(l) ) cout << a << " ";
//         // 	cout << endl;
//         // 	exit(0);
//         //     }

//         //     if(!label_data.get_preB(l).empty()){
//         // 	cout << "PreB is not empty: " << l << endl;
//         // 	for (int a : label_data.get_preB(l) ) cout << a << " ";
//         // 	cout << endl;
//         // 	exit(0);
//         //     }
//         // }

//         //Compute label.remove with a single pass through the transitions
//         for (auto & t : lts->get_transitions_pre()){
//             if(b->in_notRel(DominanceRelationComplex<LR>::Qp_pos[t.target]) &&
// 	       !b->hasTransitionFrom(lts->get_qa_post(t.label,t.src).index)) {
// 		alph.insert(t.label);
// 		label_data.add_to_remove(t.label, t.src);   
// 	    }
//         }

//         b->clear_notRel(); //We have already refined with respect to (B, B.notRel)
//         //cout << "Computed remove: " << endl;
//         for (auto l : alph){ //Populate a.PreB
//             //for each state s in B
//             for(int i = b->node.b; i <= b->node.e; ++i){
//                 int s = DominanceRelationComplex<LR>::Qp[i];
//                 //for each transition pointing s with l
//                 if(!lts->hasQaPre(l, s)) continue;
//                 const Qa & qa = lts->get_qa_pre(l, s);
//                 for(int j = qa.b; j <= qa.e; ++j){
//                     const LTSTransition & t = lts->get_transitions_pre()[j];
//                     label_data.add_to_preB(t);  //add predecessor to PreB
//                 }
//             }
//         }

//         //cout << " Computed preB" << endl;

//         for (auto l : alph){ //For each affected label
//             // cout << "Check label: " << l << " (" << g_operators[l].get_name() << ") and block: "; b->dump(Qp); cout <<" that has remove: ";
//             // for (int a :  label_data.get_remove(l)) cout << a << " "; cout  << " and pre b: ";
//             // for (int a :  label_data.get_preB(l)) cout << a << " "; cout  << endl;

//             // for (int a :  label_data.get_remove(l)){
//             // 	for (int b :  label_data.get_preB(l)){
//             // 	    if(a == b){
//             // 		cout << a << " is in remove and preB" << endl;
//             // 		exit(0);
//             // 	    }
//             // 	}
//             // }

//             vector<pair<Block *, Block *> > splitCouples;
//             vector<Block *> blocksInRemove;

//             DominanceRelationComplex<LR>::split_block(label_data.get_remove(l), splitCouples, blocksInRemove);

//             for(const auto & p : splitCouples){
//                 DEBUG_MSG(cout <<"Blocks: " << p.first->index << " (";
// 			  for(int i = p.first->node.b; i <= p.first->node.e; ++i) cout << " " << DominanceRelationComplex<LR>::Qp[i];
// 			  cout <<") and " << p.second->index << " (";
// 			  for(int i = p.second->node.b; i <= p.second->node.e; ++i) cout << " " << DominanceRelationComplex<LR>::Qp[i];
// 			  cout << ") have been splitted" << endl;
// 			  );

//                 p.first->remove_rel(p.second);
//                 p.first->add_notRel(p.second, lts, DominanceRelationComplex<LR>::Qp);

//                 /*for(auto & x : partition){
// 		     x->dump(Qp); cout << " has as rel: " << endl;
// 		     for(auto idC : x->rel){ partition[idC]->dump(Qp); cout << endl;}
// 		     }*/


//                 if(DominanceRelationComplex<LR>::unmarked(p.first->markS)) blocksToUpdate.push(p.first);
//                 else if (DominanceRelationComplex<LR>::unmarked(p.second->markS)) blocksToUpdate.push(p.second);
// 	    }

//             //cout << "processed couples splitted." << endl;
//             for(auto  d : blocksInRemove) {
//                 d->markInRemove = false;
//                 for (int s : label_data.get_preB(l)){
//                     Block * c = DominanceRelationComplex<LR>::partition[DominanceRelationComplex<LR>::Qp_block[s]].get();
//                     if (c->in_rel(d)){
//                         DEBUG_MSG(d->dump(DominanceRelationComplex<LR>::Qp);  cout << " does not simulate "; c->dump(DominanceRelationComplex<LR>::Qp); cout << endl;);
//                         c->remove_rel(d);
//                         c->add_notRel(d, lts, DominanceRelationComplex<LR>::Qp);
//                         if(DominanceRelationComplex<LR>::unmarked(c->markS)) blocksToUpdate.push(c);
//                     }
//                 }
//             }

//             //cout << "processed blocks in remove." << endl;
//         }
//         label_data.cleanup();
//         alph.clear();
//     }

//     //cout << "Done update complex relation: " << t() << endl;
//     for(int s = 0; s < simrel.num_states(); s++){
//         Block * bs = DominanceRelationComplex<LR>::partition[DominanceRelationComplex<LR>::Qp_block[s]].get();
//         for(int t = 0; t < simrel.num_states(); t++){
//             if(s != t && !bs->in_rel(DominanceRelationComplex<LR>::Qp_block[t])){//t does not simulate s
//                 simrel.remove(t, s);
//             }
//         }
//     }
//     //cout << "Updated relation" << endl;
// }

template class DominanceRelationComplexNoLD<LabelRelation>;
template class DominanceRelationComplexNoLD<LabelRelationIdentity>;
template class DominanceRelationComplexNoLD<LabelRelationNoop>;
