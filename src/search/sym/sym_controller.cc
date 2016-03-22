#include "sym_controller.h"

#include "../merge_and_shrink/opt_order.h"

#include "sym_manager.h"
#include "sym_bdexp.h"
#include "sym_hnode.h"
#include "../debug.h" 
#include "../globals.h"
#include "../option_parser.h"

using namespace std;

SymController::SymController(const Options &opts)
    : vars(new SymVariables()), mgrParams(opts), searchParams(opts), gamer_ordering(opts.get<bool>("gamer_ordering")){
  VariableOrderFinder vo (mgrParams.variable_ordering);
  vector <int> var_order; 
  if(gamer_ordering) {
      InfluenceGraph::compute_gamer_ordering(var_order);
  }else{
      while(!vo.done()){
	  var_order.push_back(vo.next());
      }
  }
  cout << "Sym variable order: ";
  for (int v : var_order) cout << v;
  cout << endl;
  
  vars->init(var_order, mgrParams);
  mgrParams.print_options();
  searchParams.print_options();
}

SymHNode * SymController::createHNode(SymHNode * node,
				      unique_ptr <SymAbstraction> && abs, 
				      SymPH * ph){
  SymHNode * newNode = new SymHNode(node, ph, move (abs));
  nodes.push_back(unique_ptr<SymHNode> (newNode));
  node->addChildren(newNode);
  newNode->addParent(node);
  return newNode;
}


void SymController::add_options_to_parser(OptionParser &parser, int maxStepTime, int maxStepNodes) {
  SymParamsMgr::add_options_to_parser(parser);
  SymParamsSearch::add_options_to_parser(parser,maxStepTime, maxStepNodes);
  parser.add_option<bool> ("gamer_ordering",
			  "Use Gamer ordering optimization", "true");
}
