#include "sym_uct_pdbs.h"

#include "sym_manager.h"
#include "symba_unsat.h"

#include "../globals.h"
#include "../rng.h"

using namespace std;


UCTNode::UCTNode(SymVariables * vars, const SymParamsMgr & mgrParams) :  pdb (nullptr),
									 mgr(new SymManager(vars, nullptr, mgrParams)), 
									 reward_fw(0), reward_bw(0), visits_fw(0), visits_bw(0) {
    for (int i = 0; i < g_variable_domain.size(); ++i) pattern.insert(i);
    children.resize(pattern.size());

}

UCTNode::UCTNode(SymVariables * vars, const SymParamsMgr & mgrParams, 
		 SymManager * omgr, AbsTRsStrategy absTRsStrategy,
		 const std::set<int> & pattern_) : pattern (pattern_), 
						   pdb(new SymPDB(vars, absTRsStrategy,  pattern_)), 
						   mgr(new SymManager(omgr, pdb.get(), mgrParams)), 
						   reward_fw(0), reward_bw(0), visits_fw(0), visits_bw(0) { 
    children.resize(pattern.size());
}


SymBreadthFirstSearch * UCTNode::initSearch(bool fw, const SymParamsSearch & searchParams) {
    if (fw) {
	assert(!fw_search);
	fw_search.reset(new SymBreadthFirstSearch (searchParams));
	fw_search->init(mgr.get(), fw);
	return fw_search.get();
    } else {
	assert(!bw_search);
	bw_search.reset(new SymBreadthFirstSearch (searchParams));
	bw_search->init(mgr.get(), fw);
	return bw_search.get();
    }
}



UCTNode * UCTNode::getChild (bool /*fw*/, SymBAUnsat * manager)  {
    int num_selected = g_rng.next(children.size());
    
    if (children[num_selected] == nullptr) {
	set<int> new_pattern (pattern);
	auto it = new_pattern.begin();
	for (int i = 0; i < num_selected; ++i) {
	    ++it;
	}
	new_pattern.erase(it);
	children [num_selected] = manager->getUCTNode(this, new_pattern);  
    }

    return children [num_selected]; 
} 


double UCTNode::uct_value(bool fw, int visits_parent, double UCT_C) const {
    if(fw) return reward_fw/visits_fw + UCT_C*sqrt(log(visits_parent)/visits_fw); 
    else   return reward_bw/visits_bw + UCT_C*sqrt(log(visits_parent)/visits_bw); 
}


SymBreadthFirstSearch * UCTNode::relax(SymBreadthFirstSearch * search, 
				       const SymParamsSearch & searchParams,
				       int maxRelaxTime, int maxRelaxNodes, 
				       double ratioRelaxTime, double ratioRelaxNodes, 
				       bool /*perimeterPDBs*/) {
    bool fw = search->isFW();
    assert(!getSearch(fw));
    
    unique_ptr<SymBreadthFirstSearch> newSearch (new SymBreadthFirstSearch(searchParams));
    newSearch->init(search, mgr.get(), maxRelaxTime, maxRelaxNodes);
    
    if (newSearch->nextStepNodes() < ratioRelaxNodes*search->nextStepNodes() && 
	newSearch->nextStepTime() < ratioRelaxTime*search->nextStepTime()) {
	cout << " Relaxed and searchable " << newSearch->nextStepNodes() << ", " << newSearch->nextStepTime() << 
	    " instead of " << search->nextStepNodes()  << " and"  << search->nextStepTime()  << endl;
	if (fw) fw_search = std::move(newSearch);
	else bw_search = std::move(newSearch); 
    } else {
	cout << " Relaxed but not searchable " << newSearch->nextStepNodes() << ", " << newSearch->nextStepTime() << 
	    " instead of " << search->nextStepNodes()  << " and"  << search->nextStepTime()  << endl;

    }

    return getSearch(fw);
}


