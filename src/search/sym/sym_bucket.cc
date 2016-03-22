#include "sym_bucket.h" 


void removeZero(Bucket & bucket) {
    bucket.erase(remove_if(begin(bucket), end(bucket),
			   [] (BDD & bdd){ return bdd.IsZero();}),
		 end(bucket));
}

void copyBucket(const Bucket & bucket, Bucket & res){
    if(!bucket.empty()){
	res.insert(end(res), begin(bucket), end(bucket));
    }
}

void moveBucket(Bucket & bucket, Bucket & res){
    copyBucket(bucket, res);
    Bucket().swap(bucket);
}

int nodeCount(const Bucket & bucket)  {
    int sum = 0;
    for(const BDD & bdd : bucket){
	sum += bdd.nodeCount();
    }
    return sum;
}


bool extract_states(Bucket & list, 
		    const Bucket & pruned, 
		    Bucket & res) {
    assert (!pruned.empty());

    bool somethingPruned = false;
    for(auto & bddList : list){
	BDD prun = pruned[0]*bddList;

	for(int i = 1; i < pruned.size(); ++i){
	    prun += pruned[i]*bddList;
	}

	if(!prun.IsZero()){
	    somethingPruned = true;
	    bddList -= prun;
	    res.push_back(prun);
	}
    }
    removeZero(list);
    return somethingPruned;
}

