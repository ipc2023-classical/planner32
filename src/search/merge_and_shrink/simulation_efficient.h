#ifndef MERGE_AND_SHRINK_SIMULATION_RELATION_EFFICIENT_H
#define MERGE_AND_SHRINK_SIMULATION_RELATION_EFFICIENT_H

#include <queue>
#include "simulation_relation.h"

#include "lts_efficient.h"

class Qa;
class LabelRelation;
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
					   const std::vector<int> & Qp);

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
			    const LabelRelation &label_dominance, int lts_id);


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


class SimulationRelationEfficient : public SimulationRelation {
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
		     std::vector<Block *> & blocksInRemove);

    template<typename LTS>
	bool get_in_pre_rel(Block * b, int src, int label, int lts_id,
			    const LTS * lts,
			    const LabelRelation & label_dominance) const;

    template<typename LTS> void init (int lts_id, const LTS * lts,
			     const LabelRelation & label_dominance, 
			     std::queue <Block *> & blocksToUpdate);
    
 public:
    SimulationRelationEfficient(const Abstraction * _abs);


    virtual void update(int /*lts_id*/, const LabelledTransitionSystem * /*lts*/,
			const LabelRelation & /*label_dominance*/){
	//update_sim(lts_id, lts, label_dominance);
    }
    virtual void update(int lts_id, const LTSEfficient * lts,
			const LabelRelation & label_dominance){
	update_sim(lts_id, lts, label_dominance);
    }

    template<typename LTS> void update_sim (int lts_id, const LTS * lts,
				   const LabelRelation & label_dominance);
  
   
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

};

#endif
