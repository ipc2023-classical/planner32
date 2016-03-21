#include "sym_abs_minimization.h"

#include "sym_util.h"
#include "../globals.h"
#include "sym_transition.h"
#include "../debug.h"

using namespace std;

SymAbsMinimization::SymAbsMinimization(SymVariables * bdd_vars, AbsMinimizationType _type,
		     int _threshold, bool _safe, double _quality) : 
    SymAbstraction(bdd_vars, AbsTRsStrategy::SHRINK_AFTER_IMG), type(_type), 
    threshold(_threshold), safe(_safe), quality(_quality){ 
      for(int i = 0; i < g_variable_name.size(); i++){
	  absVars.insert(i);
      }
}

BDD SymAbsMinimization::shrinkExists(const BDD & bdd, int /*maxNodes*/) const{
    switch (type) {
    case AbsMinimizationType::APPROX: 
	return bdd.OverApprox(0, threshold, safe, quality);
    case AbsMinimizationType::REMAP_UNDER_APPROX:
	return bdd.RemapOverApprox(0, threshold, quality);
    case AbsMinimizationType::HEAVY_BRANCH:
	return bdd.SupersetHeavyBranch(0, threshold);
    case AbsMinimizationType::SHORT_PATHS:
	return bdd.SupersetShortPaths(0, threshold, safe);
    default:
	return bdd;
    }
}


BDD SymAbsMinimization::shrinkForall(const BDD & /*bdd*/, int /*maxNodes*/) const{
    return vars->zeroBDD();
}

BDD SymAbsMinimization::shrinkTBDD(const BDD & /*bdd*/, int /*maxNodes*/) const{
    cout << "Shrink TBD not implemented in AbsMinimization" << endl;
    exit(-1);
}

BDD SymAbsMinimization::getInitialState() const{
    return vars->getStateBDD(g_initial_state());
}

BDD SymAbsMinimization::getGoal() const{
    return vars->getPartialStateBDD(g_goal);
}

std::string SymAbsMinimization::tag() const{
    return "AbsMinimization";
}

void SymAbsMinimization::print(std::ostream &os, bool /*fullInfo*/) const {
    os << "AbsMinimization";
}

ADD SymAbsMinimization::getExplicitHeuristicADD(bool /*fw*/){
    return vars->getADD(0);
}
void SymAbsMinimization::getExplicitHeuristicBDD(bool /*fw*/, map<int, BDD> & res){
    res[0] = vars->zeroBDD();
}
