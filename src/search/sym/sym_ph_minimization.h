#ifndef SYM_PH_MINIMIZATION_H
#define SYM_PH_MINIMIZATION_H

#include "sym_ph.h"
#include "sym_enums.h"

class SymHNode;

class SymPHMinimization : public SymPH {

  SymHNode * relaxation;

  AbsMinimizationType type;
  int threshold;
  bool safe;
  double quality;


/* 
 * Five different Minimization methods: 
 * UnderApprox/OverApprox(int numVars, int threshold = 0, bool safe = false, double quality = 1.0) const;
 * RemapUnderApprox/RemapOverApprox(int numVars, int threshold = 0, double quality = 1.0)
 * BiasedUnderApprox/BiasedOverApprox(const BDD& bias, int numVars, int threshold = 0, double quality1 = 1.0, double quality0 = 1.0) const;
 * BDD SubsetHeavyBranch/SupersetHeavyBranch(int numVars, int threshold) const;
 * BDD SubsetShortPaths/SupersetShortPaths (int numVars, int threshold, bool hardlimit = false) const;
 */

  
 public:
  SymPHMinimization(const Options & opts);
  virtual ~SymPHMinimization(){}

  virtual bool init();
  virtual SymBDExp * relax(SymBDExp * bdExp, SymHNode * iniHNode, Dir dir, int num_relaxations);

  virtual void dump_options() const;

  virtual bool relaxGetsHarder(){
    return false;
  }


};
#endif


