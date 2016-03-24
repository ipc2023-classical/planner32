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

    std::vector <UCTNode *> children; //Nodes more abstracted
    std::vector <UCTNode *> parents; //Nodes less abstracted

    std::unique_ptr<SymBreadthFirstSearch> fw_search, bw_search;
  
    double reward_fw, reward_bw, visits_fw, visits_bw;
public:
    UCTNode(SymVariables * vars, const SymParamsMgr & mgrParams); // Constructor for the original state space
    UCTNode(SymVariables * vars, const SymParamsMgr & mgrParams,  
	    SymManager * omgr, AbsTRsStrategy absTRsStrategy, 
	    const std::set<int> & pattern_); // Constructor for abstract state space

    UCTNode(const UCTNode & o) = delete;
    UCTNode(UCTNode &&) = default;
    UCTNode & operator=(const UCTNode& ) = delete;
    UCTNode & operator=(UCTNode &&) = default;
    ~UCTNode() = default; 

    SymBreadthFirstSearch * initSearch(bool fw, 
				       const SymParamsSearch & searchParams);

    SymBreadthFirstSearch * relax(SymBreadthFirstSearch * search, 
				  const SymParamsSearch & searchParams,
				  int maxRelaxTime, int maxRelaxNodes, 
				  double ratioRelaxTime, double ratioRelaxNodes, 
				  bool perimeterPDBs);


    SymBreadthFirstSearch * getSearch(bool fw){
	if (fw) return fw_search.get();
	else return bw_search.get(); 
    }


    UCTNode * getChild (bool fw, SymBAUnsat * manager);

    double uct_value(bool fw, int visits_parent, double UCT_C) const;

    SymManager * getMgr() const {
	return mgr.get(); 
    }

    bool isAbstractable () const {
	return pattern.size() > 1;
    }
};

#endif
