#ifndef MERGE_AND_SHRINK_VARIABLE_PARTITION_FINDER_H
#define MERGE_AND_SHRINK_VARIABLE_PARTITION_FINDER_H

#include <vector>
#include <set>
#include <map>

class VariablePartitionFinder {
 protected: 
    std::vector<std::vector<int> > partitions; //Variables contained in each cluster
    const int limit_size;
    virtual void find() = 0;

 public:
 VariablePartitionFinder(int limit) : limit_size(limit){}
~VariablePartitionFinder(){}
    void dump() const;

    const std::vector<std::vector<int> > & get_partition(){
	if(partitions.empty()){
	    find();
	}
	return partitions;
    }
};

class VariablePartitionGreedy : public VariablePartitionFinder {
std::vector<int> part_size; // size of each cluster
std::map<int, std::map<int, std::set<int>>> weights; //We store the index of the operators

    void init();
    void merge(int p1, int p2);

    //Selects the pair with maximum weight so that size is below the threshold
    std::pair<int, int> pick_parts() const;

 protected: 
    virtual void find();

 public: 
    VariablePartitionGreedy(int limit) : VariablePartitionFinder(limit) {}
};

/* class VariablePartitionGA : public VariablePartitionFinder { */
/*     const int max_iterations; */
    
/*     void find(int limit_size){ */
/* 	//Generate a random partition */
/* 	Partition best = gen_random(); */
/* 	double best_score = score(best); */
/* 	int num_it = 0; */
/* 	while(num_it < max_iterations){ */
/* 	    new_partition = modify_random(best); */
/* 	    double new_score = score(new_partition); */
/* 	    if(new_score < best_score){ */
/* 		best_score = new_score; */
/* 		best= new_partition; */
/* 	    } */
/* 	} */
/*     } */

/*     Partition gen_random(){ */
	
/*     } */

/*     Partition modify_random(Partition o){ */
	
/*     } */
/* }; */
#endif
