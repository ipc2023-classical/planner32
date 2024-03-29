#include "sym_uct_pdbs.h"

#include "sym_manager.h"
#include "symba_unsat.h"

#include "../globals.h"
#include "../rng.h"
#include "../debug.h"

#include "../legacy_causal_graph.h"

#include <algorithm>

using namespace std;


UCTNode::UCTNode(SymVariables * vars, 
		 const SymParamsMgr & mgrParams,
		 OperatorCost cost_type) : 
    pdb (nullptr), mgr(new SymManager(vars, nullptr, mgrParams, cost_type)),  
    reward_fw(0), reward_bw(0), visits_fw(0), visits_bw(0),
    redundant_fw(false), redundant_bw(false) {

    for (int i = 0; i < g_variable_domain.size(); ++i) pattern.insert(i);
}

UCTNode::UCTNode(const std::set<int> & pattern_) : 
    pattern (pattern_), pdb(nullptr), mgr(nullptr), 
    reward_fw(0), reward_bw(0), visits_fw(0), visits_bw(0), 
    redundant_fw(false), redundant_bw(false) { 

}

void UCTNode::init(SymVariables * vars, const SymParamsMgr & mgrParams, 
		   SymManager * omgr, AbsTRsStrategy absTRsStrategy) {
    pdb.reset(new SymPDB(vars, absTRsStrategy,  pattern));
    mgr.reset(new SymManager(omgr, pdb.get(), mgrParams)); 
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


void UCTNode::refine_pattern (const set<int> & pattern, vector<pair<int, set<int>>> & patterns, int relaxed_variable) const {    
    vector<set<int>> patterns_goal(pattern.size());
    vector<bool> pattern_goal_redundant (pattern.size(), false);
    vector<int> pattern_var (g_variable_domain.size(), -1);

    int g = -1;
    for (int v : pattern) {
	g++;
	vector<int> open;
	if (pattern_var[v] == -1) {
	    pattern_var[v] = g;
	    open.push_back(v);
	    patterns_goal[g].insert(v);
	} else {
	    int p = pattern_var[v];
	    pattern_goal_redundant[p] = true;
	    for (int vo : patterns_goal[p]) {
		pattern_var[vo] = g;
		patterns_goal[g].insert(vo);
	    }
	}

	for (int i = 0; i < open.size(); ++i) { 
	    for (int vpre : g_legacy_causal_graph->get_predecessors(open[i])) {
		if(!pattern.count(vpre)) continue;
		if (pattern_var[vpre] == -1) {
		    pattern_var[vpre] = g;
		    open.push_back(vpre);
		    patterns_goal[g].insert(vpre);
		} else if(pattern_var[vpre] != g) {
		    int p = pattern_var[vpre];
		    pattern_goal_redundant[p] = true;
		    for (int vo : patterns_goal[p]) {
			pattern_var[vo] = g;
			patterns_goal[g].insert(vo);
		    }
		}
		
	    }
	}
	if (patterns_goal[g].size() == pattern.size()) {
	    if (find_if(begin(patterns), end(patterns), [&](const pair<int, set<int>> & p1){
			return isSubset(patterns_goal[g], p1.second);}) == end(patterns)){
		patterns.erase(remove_if(begin(patterns), end(patterns), 
					 [&](const pair<int, set<int>>  & p1){return isSubset(p1.second, patterns_goal[g]);
					 }), end(patterns));
		patterns.push_back(pair<int, set<int>>(relaxed_variable, std::move(patterns_goal[g])));
	    }

	    return;
	}
    }

    for (int g = 0; g < pattern.size(); ++g) {
	if (!pattern_goal_redundant[g] && find_if(begin(patterns),
						  end(patterns),
						  [&](const pair<int, set<int>>  & p1){
						      return isSubset(patterns_goal[g], p1.second);}) == end(patterns)){
	    patterns.erase(remove_if(begin(patterns), end(patterns), 
				     [&](const pair<int, set<int>> & p1){return isSubset(p1.second, patterns_goal[g]);
				     }), end(patterns));

	    patterns.push_back(pair<int, set<int>>(relaxed_variable, std::move(patterns_goal[g])));
	}
    }
    
}


bool UCTNode::isSubset (const set<int> & p1, const set<int> & p2) {
    vector<int> res;
    std::set_difference(begin(p1), end(p1), begin(p2), end(p2), std::inserter(res, res.end()));
    return res.empty();
}
bool UCTNode::existsFinishedSuperset(const set<int> & pat, pair<bool, bool> & redundant, std::set<UCTNode *> & cache) {
    if(cache.count(this)) return false;
    cache.insert(this);

    if (isSubset(pat, pattern)) {
	if (fw_search && fw_search->finished()) redundant.first = true;
	if (bw_search && bw_search->finished()) redundant.second = true;
	if (redundant.first && redundant.second) return true;
	for (auto & c : children) {
	    if(c->existsFinishedSuperset(pat, redundant, cache)) return true;
	}
    } 
    
    return false;
}

void UCTNode::removeSubsets(const set<int> & pat, bool fw, set<UCTNode *> & cache) {
    if(cache.count(this)) return;
    cache.insert(this);
    if (isSubset(pattern, pat)) {
	if (fw && fw_search) {
	    redundant_fw = true;
	    //fw_search.reset(nullptr);
	}
	if (!fw && bw_search) {
	    redundant_bw = true;
	    //bw_search.reset(nullptr);
	}
    } 

    for (auto & c : children) {
	c->removeSubsets(pat, fw, cache);
	if (c->redundant_fw) children_fw.erase(find(begin(children_fw), end(children_fw), c), end(children_fw));
	if (c->redundant_bw) children_bw.erase(find(begin(children_bw), end(children_bw), c), end(children_bw));
    }
    
}


void UCTNode::initChildren (SymBAUnsat * manager)  {
    if (!children.empty()) return;
    
    vector<pair<int, set<int> > > patterns;
    for (int v : pattern) {
	set<int> new_pattern (pattern);
	new_pattern.erase(std::find(begin(new_pattern), end(new_pattern), v));

	refine_pattern (new_pattern, patterns, v);
    }

    for (const auto & pp : patterns) {
	int v = pp.first;
	const auto & p = pp.second;
	pair<bool, bool> redundant (false, false);
	if(manager->getRoot()->existsFinishedSuperset(p, redundant)) continue;
	    
	auto newChild = manager->getUCTNode(p);

	if(redundant.first) newChild->redundant_fw = true;
	if(redundant.second) newChild->redundant_bw = true;
	if(newChild) {
	    if (find(begin(children), 
		     end(children), newChild)  == end(children)) {
		children.push_back(newChild);
		//cout << "   " << *newChild << endl;
		children_fw.push_back(newChild);
		children_bw.push_back(newChild);
		children_var[newChild] = v;
	    }	       
	}
    }
}

std::ostream & operator<<(std::ostream &os, const UCTNode & node){
    os << "PDB (" << node.pattern.size() << "): ";
    for (int v : node.pattern){
	os << v << " ";
    }
    return os;
}

UCTNode * UCTNode::getRandomChild (bool fw)  {    
    const auto & children_list = (fw ? children_fw : children_bw);
    if(children_list.empty()) return nullptr;

    int num_selected = g_rng.next(children_list.size()); 
    return children_list[num_selected];	
}


double UCTNode::get_rave_value(int var, double UCT_C) const {
    return rave_reward[var] + UCT_C*sqrt(log(total_rave_visits)/(rave_visits[var] + 1)); 
}

UCTNode * UCTNode::getChild (bool fw, double UCT_C, double RAVE_K, UCTNode * root) const { 
    const auto & children_list = (fw ? children_fw : children_bw);
    if(children_list.empty()) return nullptr;
    vector<UCTNode *> best_children;
   
    if(!root->rave_reward.empty()) {
	double best_value = numeric_limits<double>::lowest();	 
	for(auto c : children_list) {
	    if(c->isExplored(fw)) continue;
	    double c_val = root->get_rave_value(children_var.at(c), UCT_C);
	    DEBUG_MSG(cout << "  Rave value " <<  children_var.at(c) << " is " << c_val ;);
	    if(best_value <= c_val) {
		if(best_value < c_val){
		    best_value = c_val;
		    best_children.clear();
		}

		best_children.push_back(c);
	    }
	}
	DEBUG_MSG(cout << endl;);
    } else {
	for(auto c : children_list) {
	    if(!c->isExplored(fw)) {
		best_children.push_back(c);
	    }
	}
    }

    
    if(!best_children.empty()) {
	int num_selected = g_rng.next(best_children.size()); 
	return best_children[num_selected];
    }	

    
    int total_visits = 0;
    for(auto c : children_list) {
	total_visits += (fw ? c->visits_fw : c->visits_bw);
    } 

    double best_value = numeric_limits<double>::lowest();
    for(auto c : children_list) {
	double c_val = 0;
	if(!rave_reward.empty()){
    
	    double RAVE_Q = rave_reward[children_var.at(c)];
	    double RAVE_B = sqrt(RAVE_K/(3*total_visits + RAVE_K));
    
	    c_val = c->uct_value(fw, total_visits, UCT_C, RAVE_B) + RAVE_B*RAVE_Q;
    
	}else{
	    c_val = c->uct_value(fw, total_visits, UCT_C, 0);
	}
	if(best_value <= c_val) {
	    if(best_value < c_val){
		best_value = c_val;
		best_children.clear();
	    }

	    best_children.push_back(c);
	}
    }
    assert(!best_children.empty());
	    
    int num_selected = g_rng.next(best_children.size()); 

    return best_children[num_selected];	
} 


    double UCTNode::uct_value(bool fw, int visits_parent, double UCT_C, double RAVE_B) const {
	// if(fw) cout << reward_fw << " " << visits_fw << endl;
	// else cout << reward_bw << " " << visits_bw << endl;
    if(fw) return (1.0-RAVE_B)*reward_fw + UCT_C*sqrt(log(visits_parent)/visits_fw); 
    else   return (1.0-RAVE_B)*reward_bw + UCT_C*sqrt(log(visits_parent)/visits_bw); 
}


SymBreadthFirstSearch * UCTNode::relax(std::shared_ptr<SymBreadthFirstSearch> search, 
				       const SymParamsSearch & searchParams,
				       int maxRelaxTime, int maxRelaxNodes, 
				       double ratioRelaxTime, double ratioRelaxNodes, 
				       bool /*perimeterPDBs*/) {
    bool fw = search->isFW();
    assert(!getSearch(fw));
    
    shared_ptr<SymBreadthFirstSearch> newSearch = make_shared<SymBreadthFirstSearch>(searchParams);
    newSearch->init(search, mgr.get(), maxRelaxTime, maxRelaxNodes);
    
    if (newSearch->nextStepNodes() < ratioRelaxNodes*search->nextStepNodes() && 
	newSearch->nextStepTime() < ratioRelaxTime*search->nextStepTime()) {
	// cout << " Relaxed and searchable " << newSearch->nextStepNodes() << ", " << newSearch->nextStepTime() << 
	//     " instead of " << search->nextStepNodes()  << " and"  << search->nextStepTime()  << endl;
	if (fw) fw_search = std::move(newSearch);
	else bw_search = std::move(newSearch); 
    } else {
	// cout << " Relaxed but not searchable " << newSearch->nextStepNodes() << ", " << newSearch->nextStepTime() << 
	//     " instead of " << search->nextStepNodes()  << " and"  << search->nextStepTime()  << endl;

    }

    return getSearch(fw);
}


void UCTNode::propagateNewDeadEnds(BDD bdd, bool isFW) {
    if(((isFW || pdb) && fw_search) || ((!isFW || pdb) && bw_search)){ 
	if (pdb){
	    bdd = mgr->shrinkForall(bdd);
	    if(!bdd.IsZero()) mgr->addDeadEndStates(!isFW, bdd);
	    
	}
	if(!bdd.IsZero()){
	    mgr->addDeadEndStates(!isFW, bdd);

	    if((isFW || pdb) && fw_search) {
		fw_search->notifyMutexes (bdd);
	    } 

	    if((!isFW || pdb) && bw_search) {
		bw_search->notifyMutexes (bdd);
	    } 
	}
    }
    
    for(auto c : children) {
	if (c) c->propagateNewDeadEnds(bdd, isFW);
    }
}

bool UCTNode::chooseDirection(double UCT_C) const {
    if(visits_fw == 0 && visits_bw == 0) return g_rng.next31()% 2;
    if(visits_bw == 0) return false;
    if(visits_fw == 0) return true;
    
    double rfw = reward_fw + UCT_C*sqrt(log((double)visits_fw + (double)visits_bw)/(double)visits_fw); 
    double rbw = reward_bw + UCT_C*sqrt(log((double)visits_fw + (double)visits_bw)/(double)visits_bw);

    DEBUG_MSG(cout << "Choosing " << visits_fw << " " << visits_bw << " " << rfw << " " <<  rbw << endl;);
    if (rfw == rbw) return g_rng.next31()% 2;
    return rfw > rbw;
}


void UCTNode::notifyReward (bool fw, double new_reward) {
    double & reward = (fw? reward_fw : reward_bw);
    int & n = (fw? visits_fw : visits_bw);
    //cout << "Reward " << reward  << " n: " << n << " new_reward: " << new_reward;
    reward = (reward *n + new_reward) / (n+1);
    //cout << " Result:  " << reward << endl;;
    n++;
}


void UCTNode::notifyRewardRAVE (bool fw, double new_reward, const std::set<int> & final_pattern) {
    //cout << "Notify Reward RAVE" << endl;
    notifyReward(fw, new_reward);

    if(rave_reward.empty()) {
	rave_reward.resize(g_variable_domain.size(), 0.0);
	rave_visits.resize(g_variable_domain.size(), 0);
	total_rave_visits = 0;
    }
    //Apply RAVE
    for (int v : pattern) {
	if(std::find (begin(final_pattern), end(final_pattern), v) != std::end(final_pattern))
	    continue;
	//vector<double> & rave_reward = (fw? rave_reward_fw : rave_reward_bw);
	//vector<int> & rave_visits = (fw? rave_visits_fw : rave_visits_bw);
    
	rave_reward[v] =  (rave_reward[v]*rave_visits[v] + new_reward) / (rave_visits[v]+1);
	rave_visits[v] ++;
	total_rave_visits++;
    }    
}
