#include "sym_breadth_first_search.h"

#include "sym_bucket.h"

using namespace std;

SymBreadthFirstSearch::SymBreadthFirstSearch(const SymParamsSearch & params): 
    SymExploration(params), estimation(params) {}


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
    
    return true;
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
    BDD res = open.back();
    open.pop_back();
    return res;
}

bool SymBreadthFirstSearch::stepImage(int maxTime, int maxNodes){
    if(mgr->getAbstraction())
	cout << ">> Step: " << *(mgr->getAbstraction());
    else
	cout << ">> Step: original";
    cout << (fw ? " fw " : " bw ");
    cout << " frontierNodes: " << nodeCount(open) << " [" << open.size() << "]"  << " total time: " << g_timer 
	 << " total nodes: " << mgr->totalNodes() << " total memory: " << mgr->totalMemory() << endl;

    
    mgr->init_transitions(); // Ensure that transitions have been initialized
    Timer sTime;
    BDD S = pop();
    
    int nodesStep = S.nodeCount();
    //double statesStep = mgr->getVars()->numStates(S);
    Timer image_time;
    map <int, vector<BDD> > resImage;
    mgr->setTimeLimit(maxTime);
    try{
	mgr->cost_image(fw, S, resImage, maxNodes);
	mgr->unsetTimeLimit();
    }catch(BDDError e){
	//Update estimation
	mgr->unsetTimeLimit();
	violated(TruncatedReason::IMAGE_COST, image_time(), maxTime, maxNodes);
	open.push_back(S);
	stats.add_image_time_failed(image_time());

	if(sTime()*1000.0 > p.maxStepTime){
	    double ratio = (double)p.maxStepTime/((double)sTime()*1000.0);
	    p.maxStepNodes = maxNodes*ratio;
	}else{
	    p.maxStepNodes = maxNodes*0.75; 
	}
	return false;
    }

    //Include new states in the open list 
    int stepNodes = nodesStep;
    
    for(auto & pairCostBDDs : resImage){ 
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
    
    estimation.stepTaken(1000*image_time(), stepNodes);
    stats.add_image_time(image_time());
    
    

    assert(isAbstracted() || perfectHeuristic->hasEvalOrig());
    stats.step_time += sTime();
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
    return /* ((fw && bdExp->getDir() != Dir::BW) || 
	      (!fw && bdExp->getDir() != Dir::FW)) &&*/
	estimation.nodes() <= maxNodes;
}

bool SymBreadthFirstSearch::isUseful(double /*ratio*/) const {
    return true;
}
