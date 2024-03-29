#ifndef SYM_PH_PDBS_H
#define SYM_PH_PDBS_H

#include "sym_ph.h"
#include "../merge_and_shrink/variable_order_finder.h"

typedef std::set<int> VarSet;

class SymHNode;

class SymPHPDBs : public SymPH {
 //We consider two relax strategies, one when relaxing the original
 //state space and other when relaxing abstract state spaces.
 const LinearPDBStrategy strategy, strategy_abstract;

 VariableOrderType var_strategy; //Strategy to select vars
 const bool randomize_strategy; //Whether to randomize the variable selection procedure

  std::map<VarSet, SymHNode *> generatedSets;
  
  //Statistics
  double time_relax;
 public:
  SymPHPDBs(const Options & opts);
  virtual ~SymPHPDBs(){}

  virtual bool init();
  virtual SymBDExp * relax(SymBDExp * bdExp, SymHNode * iniHNode, Dir dir, int num_relaxations);

  virtual void dump_options() const;
  virtual void statistics() const;

  virtual bool relaxGetsHarder(){
      return false;
  }
private:
  void getListAbstraction(SymBDExp * bdExp, SymHNode * hNode, std::vector<SymHNode *> & res);

  std::unique_ptr <SymBDExp> 
      select_binary_search(const std::vector <SymHNode *> & nodes, 
			   SymBDExp * bdExp, Dir dir, int num_relaxations);
  
  std::unique_ptr<SymBDExp> select_linear(const std::vector <SymHNode *> & nodes,
					  SymBDExp * bdExp, Dir dir, int num_relaxations);
  
  std::unique_ptr<SymBDExp> select_reverse(const std::vector <SymHNode *> & nodes,
					   SymBDExp * bdExp, Dir dir, int num_relaxations);
};
#endif
