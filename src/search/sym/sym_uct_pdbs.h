#ifndef SYM_UCT_PDBS_H
#define SYM_UCT_PDBS_H

#include "../option_parser.h"
#include "sym_abstraction.h"
#include "sym_enums.h" 
#include "sym_params.h"
#include "sym_pdb.h"
#include "sym_manager.h"

#include "sym_breadth_first_search.h"

class SymVariables;
class SymManager;
class SymBAUnsat;

class UCTNode {
private:
    std::set<int> pattern;
    
    std::unique_ptr<SymPDB> pdb;
    std::unique_ptr<SymManager> mgr;

    std::vector <UCTNode *> children, children_fw, children_bw; //Nodes more abstracted
    int num_children_explored_fw, num_children_explored_bw;

    //std::vector <UCTNode *> parents; //Nodes less abstracted

    std::unique_ptr<SymBreadthFirstSearch> fw_search, bw_search;
  
    double reward_fw, reward_bw;
    int visits_fw, visits_bw;
    bool redundant_fw, redundant_bw;

    bool isExplored (bool fw) const {
	return (fw && visits_fw > 0) || (!fw && visits_bw > 0);
    }

public:
    UCTNode(SymVariables * vars, const SymParamsMgr & mgrParams); // Constructor for the original state space
    UCTNode(const std::set<int> & pattern_); // Constructor for abstract state space

    UCTNode(const UCTNode & o) = delete;
    UCTNode(UCTNode &&) = default;
    UCTNode & operator=(const UCTNode& ) = delete;
    UCTNode & operator=(UCTNode &&) = default;
    ~UCTNode() = default; 

    void init(SymVariables * vars, const SymParamsMgr & mgrParams,
	      SymManager * omgr, AbsTRsStrategy absTRsStrategy);

    SymBreadthFirstSearch * initSearch(bool fw, 
				       const SymParamsSearch & searchParams);

    void initChildren(SymBAUnsat * manager);

    void refine_pattern (const std::set<int> & pattern, std::vector<std::set<int>> & patterns) const;


    SymBreadthFirstSearch * relax(SymBreadthFirstSearch * search, 
				  const SymParamsSearch & searchParams,
				  int maxRelaxTime, int maxRelaxNodes, 
				  double ratioRelaxTime, double ratioRelaxNodes, 
				  bool perimeterPDBs);


    SymBreadthFirstSearch * getSearch(bool fw){
	if (fw) return fw_search.get();
	else return bw_search.get(); 
    }


    UCTNode * getChild (bool fw, double UCT_C);

    double uct_value(bool fw, int visits_parent, double UCT_C) const;

    SymManager * getMgr() const {
	return mgr.get(); 
    }

    bool isAbstractable () const {
	return pattern.size() > 1;
    }

    const std::set<int> & getPattern() const {
	return pattern;
    }

    bool chooseDirection(double UCT_C) const;

    void propagateNewDeadEnds(BDD bdd, bool isFW);

    void notifyReward (bool fw, double numDeadEndsFound, const std::set<int> & pattern);

    bool existsFinishedSuperset(const std::set<int> & pat, 
				std::pair<bool, bool> & redundant, 
				std::set<UCTNode *> & cache);

    bool existsFinishedSuperset(const std::set<int> & pat, 
				std::pair<bool, bool> & redundant) {
	std::set<UCTNode *> cache;
	return existsFinishedSuperset(pat, redundant, cache);
    }


    void removeSubsets(const std::set<int> & pat, bool fw, 
		       std::set<UCTNode *> & cache);

    void removeSubsets(const std::set<int> & pat, bool fw) {
	std::set<UCTNode *> cache;
	removeSubsets(pat, fw, cache);
    }

    bool isSubset (const std::set<int> & p1, const std::set<int> & p2);


    friend std::ostream & operator<<(std::ostream &os, const UCTNode & node);
    

};

#endif
