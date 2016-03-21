#ifndef SYM_ABS_MINIMIZATION_H
#define SYM_ABS_MINIMIZATION_H

#include "sym_abstraction.h"
#include "sym_variables.h"
#include <set>



/* 
 * Five different Minimization methods: 
 UnderApprox/OverApprox(int numVars, int threshold = 0, bool safe = false, double quality = 1.0) const;
 * RemapUnderApprox/RemapOverApprox(int numVars, int threshold = 0, double quality = 1.0)
 * BiasedUnderApprox/BiasedOverApprox(const BDD& bias, int numVars, int threshold = 0, double quality1 = 1.0, double quality0 = 1.0) const;
 * BDD SubsetHeavyBranch/SupersetHeavyBranch(int numVars, int threshold) const;
 * BDD SubsetShortPaths/SupersetShortPaths (int numVars, int threshold, bool hardlimit = false) const;
*/

class SymAbsMinimization : public SymAbstraction {
  AbsMinimizationType type;
  int threshold;
  bool safe;
  double quality;

 public:  
  SymAbsMinimization(SymVariables * bdd_vars, AbsMinimizationType _type,
		     int _threshold, bool _safe, double _quality);
  ~SymAbsMinimization(){}
  virtual BDD shrinkExists(const BDD & bdd, int maxNodes) const;
  virtual BDD shrinkForall(const BDD & bdd, int maxNodes) const;
  virtual BDD shrinkTBDD (const BDD & tBDD, int maxNodes) const;

  virtual ADD getExplicitHeuristicADD(bool fw);
  virtual void getExplicitHeuristicBDD(bool fw, std::map<int, BDD> & res);

  virtual BDD getInitialState() const;
  virtual BDD getGoal() const;
  virtual std::string tag() const;

  inline SymVariables * getVars() const {
    return vars;
  }

  virtual void print(std::ostream & os, bool fullInfo) const;
};



#endif
