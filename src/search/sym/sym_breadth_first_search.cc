#include "sym_breadth_first_search.h"

#include "sym_bucket.h"

using namespace std;

SymBreadthFirstSearch::SymBreadthFirstSearch(const SymParamsSearch & params): 
    SymExploration(params), estimation(params), parent(nullptr) {}


bool SymBreadthFirstSearch::init(SymManager * manager, bool forward){
    mgr = manager;
    fw = forward;
    //Ensure that the mgr of the original state space is initialized
    //(only to get the planner output cleaner)
    mgr->init();
  
    if(fw){
	open.push_back(mgr->getInitialState());
    }else{
	open.push_back(mgr->getGoal());
    }
    closedTotal = mgr->zeroBDD();
    
    return true;
}


bool SymBreadthFirstSearch::init(SymBreadthFirstSearch * other, SymManager * manager, 
				 int maxRelaxTime, int maxRelaxNodes){
    Timer t;
    bool success = true;
    
    mgr = manager;
    fw = other->isFW();
    parent = other;
    p.inheritParentParams(parent->p);
    closedTotal = mgr->zeroBDD();
    mgr->setTimeLimit(maxRelaxTime);

    try {
	for (BDD bdd : other->open) {
	    other->filterDuplicates(bdd);
	    bdd = mgr->shrinkExists (bdd, maxRelaxNodes);
	    open.push_back(bdd);
	}

	// if(other->closed.size() == 1) {
	//     BDD bdd = other->closed[0];
	//     bdd = mgr->shrinkForall(bdd, maxRelaxNodes);
	//     closed.push_back(bdd);
	// }
	closedTotal = mgr->shrinkForall(other->closedTotal, maxRelaxNodes);

	mgr->unsetTimeLimit();
    } catch(BDDError e) {
	mgr->unsetTimeLimit();
	for (int i = open.size(); i < other->open.size(); ++i) {
	    success = false;
	    BDD bdd = other->open[i];
	    other->filterDuplicates(bdd);
	    open.push_back(bdd);
	}
    }

    mgr->mergeBucket(open);
    
    return success;
}

BDD SymBreadthFirstSearch::getUnreachableStates () const {
    BDD res = closedTotal;
    const SymBreadthFirstSearch * aux = this;
    while (aux->parent){
	aux = aux->parent;
	res += aux->closedTotal; // TODO: Should we shrink here? 
    }
    return !res;
}




long SymBreadthFirstSearch::nextStepTime() const{
    return estimation.time();
}

long SymBreadthFirstSearch::nextStepNodes() const {    
    return estimation.nextNodes();
}

long SymBreadthFirstSearch::nextStepNodesResult() const {
    return std::max(0L, estimation.nodes());
}

BDD SymBreadthFirstSearch::pop (){
    mgr->mergeBucket(open);
    BDD res = open.front();
    open.erase(open.begin());
    filterDuplicates(res);
    // if (res.nodeCount() > 100000) {
    // 	BDD reduced = res.SubsetShortPaths(0, 100000, false);
    // 	cout << " Reduced from " << res.nodeCount() << " to " << reduced.nodeCount() << endl;
    // 	open.push_back(res*!reduced);
    // 	return reduced;
    // }
    return res;
}

bool SymBreadthFirstSearch::stepImage(int maxTime, int maxNodes){
    if(mgr->getAbstraction())
	cout << ">> Step: " << *(mgr->getAbstraction());
    else
	cout << ">> Step: original";
    cout << (fw ? " fw " : " bw ");
    cout << " frontierNodes: " << nodeCount(open) << " [" << open.size() << "]"  << " total time: " << g_timer 
	 << " total nodes: " << mgr->totalNodes() << " total memory: " << mgr->totalMemory()/1000000 << "M";
    
    mgr->init_transitions(); // Ensure that transitions have been initialized
    Timer step_time;

    BDD S = pop();
    cout << "  Expanding " << S.nodeCount() << endl;
    close(S);
    
    int nodesStep = S.nodeCount();
    //double statesStep = mgr->getVars()->numStates(S);
    map <int, vector<BDD> > resImage;
    mgr->setTimeLimit(maxTime);
    try{
	mgr->cost_image(fw, S, resImage, maxNodes);
	mgr->unsetTimeLimit();
    }catch(BDDError e){
	//Update estimation
	mgr->unsetTimeLimit();
	violated(TruncatedReason::IMAGE_COST, step_time(), maxTime, maxNodes);
	open.push_back(S);
	stats.add_image_time_failed(step_time());

	if(step_time()*1000.0 > p.maxStepTime){
	    double ratio = (double)p.maxStepTime/((double)step_time()*1000.0);
	    p.maxStepNodes = maxNodes*ratio;
	}else{
	    p.maxStepNodes = maxNodes*0.75; 
	}
	return false;
    }

    //Include new states in the open list 
    int stepNodes = nodesStep;
    
    for(auto & pairCostBDDs : resImage){ 
	// cout << "ResImage: "; 
	// for (BDD & bdd : pairCostBDDs.second) {
	//     cout << bdd.nodeCount() << " ";
	// }
	// cout << endl;
	mgr->mergeBucket(pairCostBDDs.second);
	

	//Check the cut (removing states classified, since they do not need to be included in open)
	// if (!isAbstracted()){
	//     checkCutOriginal(pairCostBDDs.second, cost); 
	//     exit_with(EXIT_PLAN_FOUND);
	// }

	mgr->filterMutex(pairCostBDDs.second, fw, false);
	filterDuplicates(pairCostBDDs.second);

	for(auto & bdd : pairCostBDDs.second){  
	    if(!bdd.IsZero()){
		stepNodes = max(stepNodes, bdd.nodeCount());
		open.push_back(bdd);
	    }
	}
    }
    
    estimation.stepTaken(1000*step_time(), stepNodes);
    estimation.nextStep(nodeCount(open));
    stats.add_image_time(step_time());
    
    assert(isAbstracted() || perfectHeuristic->hasEvalOrig());
    stats.step_time += step_time();
    return true;
}


void SymBreadthFirstSearch::violated(TruncatedReason reason, double ellapsed_seconds, int maxTime, int maxNodes){
    //DEBUG_MSG(
    cout << "Truncated in " << reason << ", took " << ellapsed_seconds << " s," << 
	" maxtime: " << maxTime << " maxNodes: " << maxNodes<< endl;
    //);
    int time = 1 + ellapsed_seconds*1000;
    estimation.violated(time, maxTime, maxNodes);
}

bool SymBreadthFirstSearch::isSearchableWithNodes(int maxNodes) const{
    //cout << "Is searchable? "<< estimation.nodes() << " " << maxNodes << endl;
    return /* ((fw && bdExp->getDir() != Dir::BW) || 
	      (!fw && bdExp->getDir() != Dir::FW)) &&*/
	estimation.nodes() <= maxNodes;
}

bool SymBreadthFirstSearch::isUseful(double /*ratio*/) const {
    return true;
}
