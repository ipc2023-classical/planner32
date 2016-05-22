#include "sym_ph_minimization.h"

#include "sym_manager.h"
#include "sym_engine.h" 
#include "sym_hnode.h"
#include "sym_bdexp.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../globals.h"
#include "../rng.h"
#include "sym_pdb.h"
#include "../debug.h"
#include "sym_abs_minimization.h"

SymPHMinimization::SymPHMinimization(const Options &opts) : 
    SymPH(opts), relaxation(nullptr),
    type(AbsMinimizationType(opts.get_enum("type"))), 
    threshold(opts.get<int>("threshold")), 
    safe(opts.get<bool>("safe")), 
    quality(opts.get<double>("quality"))
{}

bool SymPHMinimization::init(){ 
  return true;
}


SymBDExp * SymPHMinimization::relax(SymBDExp * bdExp, SymHNode * iniHNode, 
				    Dir dir, int num_relaxations){ 

  if (!relaxation){
    unique_ptr <SymAbstraction> absMin = unique_ptr <SymAbstraction>(
	new SymAbsMinimization(vars, type, threshold, safe, quality)) ; 
    relaxation = iniHNode->getEngine()->createHNode(iniHNode, move(absMin), this);
  }else if (iniHNode == relaxation) {
    return nullptr; // No further relaxation possible 
  }
  unique_ptr<SymBDExp> newBDExp = createBDExp (dir, bdExp);
  if(!relax_in(bdExp, newBDExp, relaxation, num_relaxations)){
    return nullptr;
  }
  SymBDExp * res = addHeuristicExploration(bdExp, move(newBDExp));
  if(res){
    cout << ">> Abstracted exploration: " << *res << " total time: " << g_timer << endl;
    DEBUG_MSG(res->getHNode()->getAbstraction()->print(cout, true););
  }else{
    cout << ">> Abstracted not possible. total time: " << g_timer << endl;
  }
  return res;
}


static SymPH *_parse(OptionParser &parser) {
  SymPH::add_options_to_parser(parser, "SHRINK_AFTER_IMG", 1);

  parser.add_enum_option("type", AbsMinimizationTypeValues,"type of BDD minimization: ", "HEAVY_BRANCH");

  parser.add_option<int>("threshold", "value of threshold param for BDD minimization methods", "0");

  parser.add_option<bool>("safe", "value of safe param for BDD minimization methods", "true");

  parser.add_option<double>("quality", "value of quality param for BDD minimization methods", "1");

  
  Options opts = parser.parse();

  SymPH *policy = 0;
  if (!parser.dry_run()) {
    policy = new SymPHMinimization(opts);
  }

  return policy;
}

void SymPHMinimization::dump_options() const {
  SymPH::dump_options();
}

static Plugin<SymPH> _plugin("minimization", _parse);
