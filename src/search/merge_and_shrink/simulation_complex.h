#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_COMPLEX_H

/* Complex algorithm to compute a simulation relation.  It is based on
 * (Cece 2013): "Three Simulation Algorithms for Labelled Transition
 * Systems". NOT RECOMENDED: in practice is much slower than
 * simulation_simple.
 */

#include <set>
#include <queue>
#include <memory>
#include <iostream>
#include "dominance_relation.h"
#include "lts_complex.h" 
#include "../debug.h" 

class LTSComplex;

class Qa;

//Implementation of a list with constant time to check 
//class SetBits { 
//    std::vector<bool> bits;
//    int relCount;
//};
typedef std::set<int> SetBits;

class PartitionNode {
 public:
    int b, e; //List of states represented as a pair of positions

    PartitionNode (int b_, int e_) : b(b_), e(e_) {}
    PartitionNode (const PartitionNode & o) = default;
    
    bool contains (int pos) const {
	return b <= pos && pos <= e;
    }
};   

class Block {
 public:
    int index; //Position of the block in the partition

    //List of abstract states 
    PartitionNode node;

    SetBits rel;  // C s.t. B <= C
    std::vector<PartitionNode> notRel; 

    int splitCount; 

    bool markS; //Marks whether the block have been inserted in S (blocksToUpdate)
    bool markInRemove; //Marks whether it is in blocksInRemove after split

    //Counter for each sl
    std::vector<int> relCount;
 private:     
    //Constructor used for split
    Block (int index_, int b, int e, SetBits rel_, 
	   const std::vector<PartitionNode> & notRel_,
	   const std::vector<int> &relCount_) :
    index(index_), node(b, e), 
	 rel(rel_), notRel(notRel_), splitCount(0), 
	markS(false),  markInRemove(false),relCount(relCount_){
    }

 public:
 Block(int index_, int b, int e) : 
    index(index_), node(b, e),
	 splitCount(0), 
	 markS(false),  markInRemove(false) {
	rel.insert(index_); //Simulates itself
    }
    
    void initRelCounts(const std::vector<int> & qas){
	relCount = qas;
    }

    int size() const {
	return node.e - node.b + 1;
    }
    std::unique_ptr<Block> split(int index); //Splits this block between itself and a new block

    int get_newpos(){
	return node.b + splitCount++;
    }

    //void set_notin_pre_rel(int sl);

    bool in_notRel(int pos){
	for(const auto & n : notRel){
	    if(n.contains(pos)) 
		return true;
	}
	return false;
    }

    void applyNotRel(std::function<void(int pos_sNotRel)> && f) const {
	for(const auto & n : notRel){
	    for(int i = n.b; i <= n.e; ++i){
		f(i);
	    }	  
	}
    }   

    //Block does not simulate me anymore, so I should decrement
    //relCount of all the transitions to block!
    template<typename LTS> void add_notRel(Block * block, const LTS * lts, 
					   const std::vector<int> & Qp)  {
	notRel.push_back(block->node);

	//Block does not simulate me anymore, so I should decrement
	//relCount of all the transitions to block!
	for(int i = block->node.b; i <= block->node.e; ++i){
	    //cout << "Decrease transitions going to: " << Qp[i] << endl;
	    lts->applyPreTarget(Qp[i], [&](const LTSTransition & t){
		    int sl = lts->get_pos_qa_post(t.label, t.src);
		    //cout << "Decreasing " << sl << " because of s=" << t.src << " l=" << t.label << endl;
		    relCount[sl]--;
		    return false;
		});
	}
    }

    void clear_notRel(){
	notRel.clear();
    }

    bool in_rel(Block * block){
	return rel.count(block->index);
    }

    bool in_rel(int b_index){
	return rel.count(b_index);
    }

    void add_rel(Block * b){
	rel.insert(b->index);
    }

    void remove_rel(Block * b){
	rel.erase(b->index);
    }

    bool hasTransitionFrom (int sl) {
	return relCount[sl] > 0;
    }

    void dump(const std::vector<int> & qp);
};

//Class that helps to handle the functionality needed for Remove and
//preB
template <typename LR>
class LabelData {
    //Set of states that should be removed (because they have a
    //transitino that other states do not have)
    //SetBits hashsl; //Only one hash for remove and preB because they are disjoint

    //For each label, list of sl items
    std::vector<std::set<int>> remove_by_label;
    
    //All the states that have a transition to B with this label.
    std::vector<std::set<int> > preB_by_label; 

 public:
    LabelData(int num_labels){
	remove_by_label.resize(num_labels);
	preB_by_label.resize(num_labels);
    }

    void cleanup(){
	//std::cout << "CLEANING UP" << std::endl;
	//hashsl.clear();
	for (auto & l : remove_by_label) l.clear();
	for (auto & l : preB_by_label) l.clear();
    }

    void add_to_remove(int l, int s/*, int sl*/){
	/* if (!hashsl.count(sl)){ */
	/*     remove_by_label[l].push_back(s); */
	/*     hashsl.insert(sl); */
	/* } */
	remove_by_label[l].insert(s);
    }

    void add_to_remove_init(const Qa & qa,
			    const LR &label_dominance, int lts_id) {
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

    const std::set<int> & get_remove(int label) {
	return remove_by_label[label];	
    }

    void add_to_preB(const LTSTransition & t){
	/* if (!hashsl.count(t.sl)){ */
	/*     preB_by_label[t.label].push_back(t.src); */
	/*     hashsl.insert(t.sl); */
	/* } */
	preB_by_label[t.label].insert(t.src);
    }

    const std::set<int> & get_preB(int label) const {
	return preB_by_label[label];
    }
};

template <typename LR> 
class ComplexDominanceRelation : public DominanceRelationLR<LR> {
 protected:
    //By now we assume that the partition is unitary... we can improve
    //this later with EquivalenceRelation
    //Temporary structures to compute the simulation
    std::vector <std::unique_ptr<Block>> partition;
    std::vector<int> Qp;   //List of abstract states sorted according to the partition 
    std::vector<int> Qp_block; //Assigns each state to the corresponding block 
    std::vector<int> Qp_pos; // Dual of Qp (to reorder Qp in the split)

    void split_block(const std::set<int> & set_remove, 
		     std::vector<std::pair<Block *, Block *> > & splitCouples, 
		     std::vector<Block *> & blocksInRemove) {
	// cout << "SplitBlock" << endl;
	// cout << "Remove set: "; for (int i : set_remove) cout  << i; cout << endl;
	// cout << "Qp: "; for (auto q : Qp) cout << " " << q; cout << endl;
	// cout << "Qp_block: "; for (auto q : Qp) cout << " " << Qp_block[q]; cout << endl;
	// cout << "Qp_pos: "; for (auto q : Qp_pos) cout << " " << q; cout << endl;

	std::set<Block *> touched; //TODO: we"re not using the same data structure than in the paper
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

		splitCouples.push_back(std::pair<Block *, Block *>(c, d));

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
    

//Given a transition, s->tl->x, check whether s still has a transition
//s - l -> t' with l >= tl and t' in rel(B), i. e., if for any block C
//in rel(B), and label l, C.relCount(sl) > 0
    template<typename LTS>
	bool get_in_pre_rel(Block * b, int src, int label, int lts_id,
			    const LTS * lts,
			    const LR & label_dominance) const {
	//cout << "Check whether " << src << " is in preRel" << b->index << " l" << label  << std::endl;
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
		// 	 << " state " << src << " hasQa? " <<  lts->hasQaPost(l, src) << std::endl;
		if(lts->hasQaPost(l, src)){ //Get sl, (if exists)
		    int sl = lts->get_qa_post(l,src).index;
		    //cout << "sl index " << sl << std::endl;
		    //For every block C in rel(B)
		    if(b->hasTransitionFrom(sl)){
			//cout << "YES" << std::endl;
			return true; //Found!
		    }

		    // Bug solved: Checking b->hasTransitionFrom(sl)
		    // already checks whether sl has a transition to any
		    // block better than B. The loop causes a severe bug.
		    // for(int idC : b->rel){
		    //     //cout << "Have a transition labelled with " << l << " to ";
		    //     //partition[idC]->dump(Qp);
		    //     /*cout << " >= "; b->dump(Qp); std::cout << "?" << std::endl;
		    //      std::cout << " relCount: "; 
		    //      for(int i : partition[idC]->relCount) std::cout << i;
		    //      std::cout << std::endl;*/

		    //     if(partition[idC]->hasTransitionFrom(sl)){
		    //         //cout << "YES" << std::endl;
		    //         return true; //Found!
		    //     }
		    // }
		}
	    }
	}
	//cout << " NOT!! " << src << std::endl;
	return false;
    }

    template<typename LTS> void init (int lts_id, const LTS * lts,
				      const LR & label_dominance, 
				      std::queue <Block *> & blocksToUpdate) {

   

	// std::cout << "Init complex relation" << std::endl;
	// std::cout << "Qp: "; for (auto q : Qp) std::cout << " " << q; std::cout << std::endl;
	// std::cout << "Qp_block: "; for (auto q : Qp) std::cout << " " << Qp_block[q]; std::cout << std::endl;
	// std::cout << "Qp_pos: "; for (auto q : Qp_pos) std::cout << " " << q; std::cout << std::endl;

    //Initialize relCounts. For every qa, set the number of transitions q-l>x with l >= a
    const std::vector<Qa> & qas = lts->get_qa_post();

    int num_labels = label_dominance.get_num_labels();
    LabelData<LR> label_data (num_labels);

    //Initialize relCounts and l.Remove
    //cout << "Init rel counts" << std::endl; 

    std::vector<int> relCount (qas.size(), 0); 
    for (auto & qa : qas){ 
        relCount[qa.index] += qa.size();

        //If a state does not have any outgoing transition labelled
        //with l >= a and a is not dominated by noop, it cannot
        //possibly simulate any other state with a
        label_data.add_to_remove_init(qa, label_dominance, lts_id);
    }    

    //cout << "Initialized relCount: ";     for(int i : relCount) std::cout << i << " "; std::cout << std::endl;

    //Set relCount to all blocks (initially all they have the same)
    for(auto & b : partition) b->initRelCounts(relCount); 


    //Perform split
    for(int l = 0; l < num_labels; l++) {
        //Skip labels dominated by noop (all states have them)
        if(label_dominance.dominated_by_noop(l, lts_id)){
            //cout << "Skip label " << l << " because is dominated by noop" << std::endl;
            continue;
        }
        //cout << "Split wrt label " << l << ": " << g_operators[l].get_name() << std::endl;

        std::vector<std::pair<Block *, Block *> > splitCouples;
        std::vector<Block *> blocksInRemove;
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
        /*cout << " Push c: " << c->index << std::endl;
	c->dump(Qp); std::cout << " relCount: "; 
	for(int i : c->relCount) std::cout << i;
	cout << std::endl;*/

        for (auto & aux : partition){
            if (!c->in_rel(aux.get())){
                //cout << "Adding block "; aux->dump(Qp);
                //cout << " as not rel of  "; c->dump(Qp); std::cout << std::endl;
                c->add_notRel(aux.get(), lts, Qp);
            }
        }
        if(!c->notRel.empty() && unmarked(c->markS))
            blocksToUpdate.push(c.get());

        /*c->dump(Qp); std::cout << " relCount: ";
	for(int i : c->relCount) std::cout << i;
	cout << std::endl;*/

    }

    /*cout << "Init complex relation, done." << std::endl;
    lts->dump_names();
    std::cout << "Qp: "; for (auto q : Qp) std::cout << " " << q; std::cout << std::endl;
    std::cout << "Qp_block: "; for (auto q : Qp) std::cout << " " << Qp_block[q]; std::cout << std::endl;
    std::cout << "Qp_pos: "; for (auto q : Qp_pos) std::cout << " " << q; std::cout << std::endl;
    for(auto & b : partition){
    	b->dump(Qp); std::cout << " has as rel: " << std::endl;
    	for(auto idC : b->rel){ partition[idC]->dump(Qp); std::cout << std::endl;}
	}*/
}


    bool unmarked(bool & mark) const {
	if(!mark){
	    mark = true;
	    return true;
	}else{
	    return false;
	}
    }

    template <typename LTS>
    bool exists_transition_to (const LTS * lts,int label, Block * b){
	for(int sb = b->node.b; sb <= b->node.e; sb++){
	    if(lts->hasQaPre(label, Qp[sb])){
		//std::cout << " found transition to block!" << std::endl;
		return true;
	    }else{
		//std::cout << " there is no transition " << label << " to " << Qp[sb] << std::endl;
	    }
	}
	return false;
    }
    
 public:
    ComplexDominanceRelation (Labels * labels) : DominanceRelationLR<LR> (labels) {}

    virtual std::unique_ptr<SimulationRelation> init_simulation (Abstraction * _abs);

    virtual std::unique_ptr<SimulationRelation> 
	init_simulation_incremental (CompositeAbstraction * /*_abs*/, 
				     const SimulationRelation & /*simrel_one*/, 
				     const SimulationRelation & /*simrel_two*/) {
	std::cerr << "Error: ComputeSimulationRelationComplex::init_simulation_incremental not implemented yet" << std::endl;
	std::exit(-1);
	return false;
    }

    virtual void update(int /*lts_id*/, const LabelledTransitionSystem * /*lts*/,
			const LR & /*label_dominance*/, 
			SimulationRelation & /*simrel*/){
	//update_sim(lts_id, lts, label_dominance);
    }

    virtual void update(int lts_id, const LTSComplex * lts,
			const LR & label_dominance, 
			SimulationRelation & simrel){
	update_sim(lts_id, lts, label_dominance, simrel);
    }

    template<typename LTS> void update_sim (int lts_id, const LTS * lts,
					    const LR & label_dominance, 
					    SimulationRelation & simrel) {    
    Timer t;
    //cout << "Update complex relation" << std::endl;
    //label_dominance.dump();
    LabelData<LR> label_data (label_dominance.get_num_labels());
    std::set<int> alph; //TODO: Not using the data structure suggested in the paper

    std::queue <Block *> blocksToUpdate; 

    //TODO: if info init_label(info) else
    init(lts_id, lts, label_dominance, blocksToUpdate);
    //if(true) return;

    while (!blocksToUpdate.empty()){
        //queue <Block *> copy(blocksToUpdate);
        // std::cout << " To update: ";
        // while(!copy.empty()){
        //     Block * b = copy.front();
        //     copy.pop();
        //     std::cout << "  ";
        //     b->dump(Qp);
        // }
        // std::cout << std::endl;

        Block * b = blocksToUpdate.front();
        blocksToUpdate.pop();
        b->markS = false;

        DEBUG_MSG(std::cout << "Selected block: "; b->dump(Qp); std::cout << std::endl;);
        //Block b has notRel => c does not dominate b anymore. Check
        //transitions to c because they do not simulate transitions to
        //b anymore.

        /*cout << " relCount: ";
	for(int i : b->relCount) std::cout << i;
	cout << std::endl;*/
        // for(auto & x : partition){
        //     x->dump(Qp); std::cout << " has as rel: " << std::endl;
        //     for(auto idC : x->rel){ partition[idC]->dump(Qp); std::cout << std::endl;}
        // }

        assert(alph.empty()); // assert for all labels, PreB and Remove is empty

        // for (int l = 0; l < label_dominance.num_labels(); ++l){ //For each affected label
        //     if(!label_data.get_remove(l).empty()){
        // 	cout << "Remove is not empty: " << l << std::endl;
        // 	for (int a : label_data.get_remove(l) ) std::cout << a << " ";
        // 	cout << std::endl;
        // 	exit(0);
        //     }

        //     if(!label_data.get_preB(l).empty()){
        // 	cout << "PreB is not empty: " << l << std::endl;
        // 	for (int a : label_data.get_preB(l) ) std::cout << a << " ";
        // 	cout << std::endl;
        // 	exit(0);
        //     }
        // }

        //Set of pairs <src, label> to check we'll split states
        //depending on whether they have -- l' --> B or not
        std::set<std::pair<int, int> > slToCheck;

        //Compute label.remove with a single pass through the transitions
        for (auto & t : lts->get_transitions_pre()){
            //cout << "    Check T: " << t.src << " -- " << t.label << " --> " << t.target << std::endl;
            // std::cout << " Remove: ";
            // for (int a : label_data.get_remove(t.label)) std::cout << a << " ";
            // std::cout << std::endl;

            /*cout << "Is " << t.target << " in not rel of block " <<
	      b->index << "?   " << b->in_notRel(Qp_pos[t.target]) << std::endl;*/

            if(b->in_notRel(Qp_pos[t.target])){
                //cout << t << " does not support " << b->index << " anymore" << std::endl;

                // For every label l <= tl, s.t. exist a transition to block B
                for(int l = 0; l < label_dominance.get_num_labels(); ++l){
                    if(label_dominance.dominates(t.label, l, lts_id) &&
                            exists_transition_to(lts, l, b)){
                        slToCheck.insert(std::pair<int, int> (t.src, l));
                    }
                }
            }
        }

        //cout << "Check all those states that are not better than the block" << std::endl;
        //For every state in not rel and every label dominated by noop
        //with a transition to B
        for(int l = 0; l < label_dominance.get_num_labels(); ++l){
            if(label_dominance.dominated_by_noop(l, lts_id)){
                //cout << "Checking label " << l << " because is noop" << std::endl;

                if(exists_transition_to(lts, l, b)){
                    b->applyNotRel([&](int pos){ //For every concrete state in notRel(b) (s does not simulate b anymore)
                        int s = Qp[pos];
                        slToCheck.insert(std::pair<int, int> (s, l));
                    });
                }
            }
        }

        for(auto & p : slToCheck){
            int s = p.first; int l = p.second;
            //cout << "Check sl: " << s << " " << l << std::endl;
            if(!get_in_pre_rel(b, s, l, lts_id, lts, label_dominance)){
                alph.insert(l);
                label_data.add_to_remove(l, s);
                //cout << "Add " << s << " to remove of " << l << std::endl;
            }
        }

        b->clear_notRel(); //We have already refined with respect to (B, B.notRel)
        //cout << "Computed remove: " << std::endl;
        for (auto l : alph){ //Populate a.PreB
            //for each state s in B
            for(int i = b->node.b; i <= b->node.e; ++i){
                int s = Qp[i];
                //for each transition pointing s with l
                if(!lts->hasQaPre(l, s)) continue;
                const Qa & qa = lts->get_qa_pre(l, s);
                for(int j = qa.b; j <= qa.e; ++j){
                    const LTSTransition & t = lts->get_transitions_pre()[j];
                    label_data.add_to_preB(t);  //add predecessor to PreB
                }
            }
        }

        //cout << " Computed preB" << std::endl;

        for (auto l : alph){ //For each affected label
            // std::cout << "Check label: " << l << " (" << g_operators[l].get_name() << ") and block: "; b->dump(Qp); std::cout <<" that has remove: ";
            // for (int a :  label_data.get_remove(l)) std::cout << a << " "; std::cout  << " and pre b: ";
            // for (int a :  label_data.get_preB(l)) std::cout << a << " "; std::cout  << std::endl;

            // for (int a :  label_data.get_remove(l)){
            // 	for (int b :  label_data.get_preB(l)){
            // 	    if(a == b){
            // 		cout << a << " is in remove and preB" << std::endl;
            // 		exit(0);
            // 	    }
            // 	}
            // }

            std::vector<std::pair<Block *, Block *> > splitCouples;
            std::vector<Block *> blocksInRemove;

            split_block(label_data.get_remove(l), splitCouples, blocksInRemove);

            for(const auto & p : splitCouples){
                DEBUG_MSG(std::cout <<"Blocks: " << p.first->index << " (";
			  for(int i = p.first->node.b; i <= p.first->node.e; ++i) std::cout << " " << Qp[i];
			  std::cout <<") and " << p.second->index << " (";
			  for(int i = p.second->node.b; i <= p.second->node.e; ++i) std::cout << " " << Qp[i];
			  std::cout << ") have been splitted" << std::endl;
			  );

                p.first->remove_rel(p.second);
                p.first->add_notRel(p.second, lts, Qp);

                /*for(auto & x : partition){
		     x->dump(Qp); std::cout << " has as rel: " << std::endl;
		     for(auto idC : x->rel){ partition[idC]->dump(Qp); std::cout << std::endl;}
		     }*/


                if(unmarked(p.first->markS)) blocksToUpdate.push(p.first);
                else if (unmarked(p.second->markS)) blocksToUpdate.push(p.second);
	    }

            //cout << "processed couples splitted." << std::endl;
            for(auto  d : blocksInRemove) {
                d->markInRemove = false;
                for (int s : label_data.get_preB(l)){
                    Block * c = partition[Qp_block[s]].get();
                    if (c->in_rel(d)){
                        DEBUG_MSG(d->dump(Qp);  std::cout << " does not simulate "; c->dump(Qp); std::cout << std::endl;);
                        c->remove_rel(d);
                        c->add_notRel(d, lts, Qp);
                        if(unmarked(c->markS)) blocksToUpdate.push(c);
                    }
                }
            }

            //cout << "processed blocks in remove." << std::endl;
        }
        label_data.cleanup();
        alph.clear();
    }

    //cout << "Done update complex relation: " << t() << std::endl;
    for(int s = 0; s < simrel.num_states(); s++){
        Block * bs = partition[Qp_block[s]].get();
        for(int t = 0; t < simrel.num_states(); t++){
            if(s != t && !bs->in_rel(Qp_block[t])){ //t does not simulate s
                simrel.remove(t, s);
            }
        }
    }
    //cout << "Updated relation" << std::endl;
}


    virtual bool propagate_label_domination(int /*lts_id*/, 
					    const LabelledTransitionSystem * /*lts*/,
					    const LR & /*label_dominance*/, 
					    int /*l*/, int /*l2*/, SimulationRelation & /*simrel*/) const {
	std::cerr << "Error: ComputeSimulationRelationComplex::propagate_label_domination not implemented yet" << std::endl;
	std::exit(-1);
	return false;
    }
};

#endif
