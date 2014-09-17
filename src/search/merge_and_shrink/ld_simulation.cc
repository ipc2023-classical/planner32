#include "ld_simulation.h"

#include "abstraction.h"
#include "shrink_bisimulation.h"
#include "simulation_relation.h"
#include "merge_strategy.h"
#include "labelled_transition_system.h"

using namespace std;

LDSimulation::LDSimulation(bool unit_cost, const Options &opts, OperatorCost cost_type) : 
  use_expensive_statistics(opts.get<bool>("expensive_statistics")),
  limit_absstates_merge(opts.get<int>("limit_merge")),
  merge_strategy(opts.get<MergeStrategy *>("merge_strategy")),
  use_bisimulation(opts.get<bool>("use_bisimulation")), 
  intermediate_simulations(opts.get<bool>("intermediate_simulations")), 
  labels (new Labels(unit_cost, opts, cost_type)) //TODO: c++14::make_unique 
{ 
}

LDSimulation::LDSimulation(const Options &opts) : 
  use_expensive_statistics(opts.get<bool>("expensive_statistics")),
  limit_absstates_merge(opts.get<int>("limit_merge")),
  merge_strategy(opts.get<MergeStrategy *>("merge_strategy")),
  use_bisimulation(opts.get<bool>("use_bisimulation")), 
  intermediate_simulations(opts.get<bool>("intermediate_simulations")){
  //TODO: Copied from heuristic.cc. Move to the PlanningTask class
  //that someday will be added to FD.
  OperatorCost cost_type = OperatorCost::NORMAL;
  bool is_unit_cost = true;
  for (size_t i = 0; i < g_operators.size(); ++i) {
    if (get_adjusted_action_cost(g_operators[i], cost_type) != 1) {
      is_unit_cost = false;
      break;
    }
  }
  labels = unique_ptr<Labels> (new Labels(is_unit_cost, opts, cost_type)); //TODO: c++14::make_unique 
}

LDSimulation::~LDSimulation(){
  for(auto abs : abstractions){
    delete abs;
  }
  for(auto sim : simulations){
    delete sim;
  }
}


void LDSimulation::build_abstraction() {
  // TODO: We're leaking memory here in various ways. Fix this.
  //       Don't forget that build_atomic_abstractions also
  //       allocates memory.

  // vector of all abstractions. entries with 0 have been merged.
  vector<Abstraction *> all_abstractions;
  all_abstractions.reserve(g_variable_domain.size() * 2 - 1);
  Abstraction::build_atomic_abstractions(all_abstractions, labels.get());

  // vector of all simulations. only used when computing intermediate simulations
  // vector<SimulationRelation *> all_simulations;
  // if(intermediate_simulations){
  //   all_simulations.reserve(g_variable_domain.size() * 2 - 1);
  //   compute_ld_simulation(labels.get(), all_abstractions, all_simulations);
  // }

  unique_ptr<ShrinkStrategy> shrink_strategy;
  if(use_bisimulation){
    shrink_strategy = unique_ptr<ShrinkStrategy>(ShrinkBisimulation::create_default());
    cout << "Shrinking atomic abstractions..." << endl;
    for (size_t i = 0; i < all_abstractions.size(); ++i) {
      all_abstractions[i]->compute_distances();
      // if (!all_abstractions[i]->is_solvable())
      // 	return all_abstractions[i];
      shrink_strategy->shrink_atomic(*all_abstractions[i]);
    }
  }
    
  cout << "Merging abstractions..." << endl;

  while (!merge_strategy->done()) {
    pair<int, int> next_systems = merge_strategy->get_next(all_abstractions, 
							   limit_absstates_merge);
    int system_one = next_systems.first;
    if(system_one == -1){
      break; //No pairs to be merged under the limit
    }
    Abstraction *abstraction = all_abstractions[system_one];
    assert(abstraction);
    int system_two = next_systems.second;
    assert(system_one != system_two);
    Abstraction *other_abstraction = all_abstractions[system_two];
    assert(other_abstraction);

    // Note: we do not reduce labels several times for the same abstraction
    bool reduced_labels = false;
    if (shrink_strategy && shrink_strategy->reduce_labels_before_shrinking()) {
      labels->reduce(make_pair(system_one, system_two), all_abstractions);
      reduced_labels = true;
      abstraction->normalize();
      other_abstraction->normalize();
      abstraction->statistics(use_expensive_statistics);
      other_abstraction->statistics(use_expensive_statistics);
    }

    // distances need to be computed before shrinking
    abstraction->compute_distances();
    other_abstraction->compute_distances();
    // if (!abstraction->is_solvable())
    //     return abstraction;
    // if (!other_abstraction->is_solvable())
    //     return other_abstraction;

    if(shrink_strategy){
      shrink_strategy->shrink_before_merge(*abstraction, *other_abstraction);
      abstraction->statistics(use_expensive_statistics);
      other_abstraction->statistics(use_expensive_statistics);
    }

    if (shrink_strategy && !reduced_labels) {
      labels->reduce(make_pair(system_one, system_two), all_abstractions);
    }


    abstraction->normalize();
    other_abstraction->normalize();

    if (!reduced_labels) {
      // only print statistics if we just possibly reduced labels
      other_abstraction->statistics(use_expensive_statistics);
      abstraction->statistics(use_expensive_statistics);
    }

    //TODO: UPDATE SIMULATION WHEN DOING INCREMENTAL COMPUTATION 
    Abstraction *new_abstraction = new CompositeAbstraction(labels.get(),
							    abstraction,
							    other_abstraction);

    abstraction->release_memory();
    other_abstraction->release_memory();

    new_abstraction->statistics(use_expensive_statistics);

    all_abstractions[system_one] = 0;
    all_abstractions[system_two] = 0;
    all_abstractions.push_back(new_abstraction);

    // if(intermediate_simulations){
    // 	compute_ld_simulation(labels, all_abstractions, all_simulations);
    // }
  }

  for (size_t i = 0; i < all_abstractions.size(); ++i) {
    if (all_abstractions[i]) {
      abstractions.push_back(all_abstractions[i]);
      all_abstractions[i]->compute_distances();
      all_abstractions[i]->statistics(use_expensive_statistics);
      //all_abstractions[i]->release_memory();
    }
  }
}

void LDSimulation::compute_ld_simulation() {
    compute_ld_simulation(labels.get(), abstractions, simulations);
}

//Compute a simulation from 0
void LDSimulation::compute_ld_simulation(Labels * _labels, vector<Abstraction *> & _abstractions, 
					 vector<SimulationRelation *> & _simulations) {
    cout << "Building LTSs and Simulation Relations" << endl;
    vector<LabelledTransitionSystem *> ltss;
    for (auto a : _abstractions){
	ltss.push_back(a->get_lts());
	cout << "LTS built: " << ltss.back()->size() << " states " << ltss.back()->num_transitions() << " transitions " << _labels->get_size() << " num_labels"  << endl;
	//Create initial goal-respecting relation
	_simulations.push_back(new SimulationRelation(a));
    }

    compute_ld_simulation(_labels, ltss, _simulations);
}


void LDSimulation::compute_ld_simulation(Labels * _labels, vector<LabelledTransitionSystem *> & _ltss, 
					 vector<SimulationRelation *> & _simulations) {
    Timer t;
    LabelRelation label_dominance (_labels);
    label_dominance.init(_ltss, _simulations);
  do{
    cout << "Label dominance initialized: " << t() << endl;
    cout << "LDsimulation loop: ";
    //label_dominance.dump();
    for (int i = 0; i < _simulations.size(); i++){
	_simulations[i]->update(i, _ltss[i], label_dominance);
	//cout << "loooooooping" <<  t() << endl;

      //_simulations[i]->dump(_ltss[i]->get_names());
    }
    cout << " took " << t() << "s" << endl;
  }while(label_dominance.update(_ltss, _simulations));
}

// void LDSimulation::compute_ld_simulation_after_merge(vector<Abstraction *> & all_abstractions, 
// 						     vector<SimulationRelation *> & all_simulations, 
//  						     const pair<int, int> & next_systems) {
//     Abstraction *new_abstraction = all_abstractions.back();
//     vector <Abstraction *> _abstractions; //vector with only the non-null pointers
//     vector <Abstraction *> _simulations;  //vector with only the non-null simulations
//     for (auto a : all_abstractions) if(a) _abstractions.push_back(a);
//     for (auto s : all_simulations) if(s) _simulations.push_back(s);

//     //b) Reset simulations 
//     //for(auto & sim :_ simulations){
// 	//sim->reset();
// 	//}
//     //c) Initialize simulation from previous simulations
//     SimulationRelation * new_sim = new SimulationRelation(new_abstraction, 
// 							  simulations[next_systems.first],
// 							  simulations[next_systems.second]);
//     all_simulations.push_back(new_sim);
//     _simulations.push_back(new_sim);  
// }


void LDSimulation::dump_options() const {
  merge_strategy->dump_options();
  labels->dump_options();
  cout << "Expensive statistics: "
       << (use_expensive_statistics ? "enabled" : "disabled") << endl;

  if (use_expensive_statistics) {
    string dashes(79, '=');
    cerr << dashes << endl
	 << ("WARNING! You have enabled extra statistics for "
	     "merge-and-shrink heuristics.\n"
	     "These statistics require a lot of time and memory.\n"
	     "When last tested (around revision 3011), enabling the "
	     "extra statistics\nincreased heuristic generation time by "
	     "76%. This figure may be significantly\nworse with more "
	     "recent code or for particular domains and instances.\n"
	     "You have been warned. Don't use this for benchmarking!")
	 << endl << dashes << endl;
  }

}

void LDSimulation::initialize() {
  Timer timer;
  cout << "Initializing simulation heuristic..." << endl;
  dump_options();
  verify_no_axioms();
 
  if(limit_absstates_merge > 1){
    build_abstraction();
  }else{
    Abstraction::build_atomic_abstractions(abstractions, labels.get());
  }
   
  cout << "Computing simulation..." << endl;
  compute_ld_simulation();

  if(use_expensive_statistics){
    // for(int i = 0; i < simulations.size(); i++){ 
    //   simulations[i]->dump(ltss[i]->get_names()); 
    // } 
  }

  // for (auto l : ltss){
  //   delete l;
  // }
  //vector<LabelledTransitionSystem *>().swap(ltss);

  cout << "Done initializing simulation heuristic [" << timer << "]"
       << endl;
  int num_equi = num_equivalences();
  int num_sims = num_simulations();
  cout << "Total Simulations: " << num_sims + num_equi*2  << endl;
  cout << "Similarity equivalences: " << num_equi  << endl;
  cout << "Only Simulations: " << num_sims << endl;
  for(int i = 0; i < simulations.size(); i++){ 
    cout << "States after simulation: " << simulations[i]->num_states() << " " 
	 << simulations[i]->num_different_states() << endl;
  }
  //exit(0);
}

int LDSimulation::num_equivalences() const {
  int res = 0;
  for(int i = 0; i < simulations.size(); i++){ 
    res += simulations[i]->num_equivalences(); 
  }
  return res;  
}


int LDSimulation::num_simulations() const {
  int res = 0;
  for(int i = 0; i < simulations.size(); i++){ 
    res += simulations[i]->num_simulations(); 
  } 
  return res;  
}

void LDSimulation::precompute_dominated_bdds(SymVariables * vars){
  for(auto & sim : simulations){
    sim->precompute_absstate_bdds(vars);
    sim->precompute_dominated_bdds();
  }
}

void LDSimulation::precompute_dominating_bdds(SymVariables * vars){
  for(auto & sim : simulations){
    sim->precompute_absstate_bdds(vars);
    sim->precompute_dominating_bdds();
  }  
}


bool LDSimulation::pruned_state(const State &state) const {
  for(auto sim : simulations) {
    if(sim->pruned(state)){
      return true;
    }
  }
  return false;
}

int LDSimulation::get_cost(const State &state) const {
    int cost = 0;
    for(auto sim : simulations) {
	int new_cost = sim->get_cost(state);
	if (new_cost == -1) return -1;
	cost = max (cost, new_cost);
    }
    return cost;
}


BDD LDSimulation::getSimulatedBDD(SymVariables * vars, const State &state ) const{
  BDD res = vars->oneBDD();
  for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
    res *= (*it)->getSimulatedBDD(state); 
  }
  return res;
}

BDD LDSimulation::getSimulatingBDD(SymVariables * vars, const State &state ) const{
  BDD res = vars->oneBDD();
  for (auto it = simulations.rbegin(); it != simulations.rend(); it++){
    res *= (*it)->getSimulatingBDD(state); 
  }
  return res;
}

void LDSimulation::add_options_to_parser(OptionParser &parser){
  vector<string> label_reduction_method;
  label_reduction_method.push_back("NONE");
  label_reduction_method.push_back("OLD");
  label_reduction_method.push_back("TWO_ABSTRACTIONS");
  label_reduction_method.push_back("ALL_ABSTRACTIONS");
  label_reduction_method.push_back("ALL_ABSTRACTIONS_WITH_FIXPOINT");
  parser.add_enum_option("label_reduction_method", label_reduction_method,
			 "label reduction method: "
			 "none: no label reduction will be performed "
			 "old: emulate the label reduction as desribed in the "
			 "IJCAI 2011 paper by Nissim, Hoffmann and Helmert."
			 "two_abstractions: compute the 'combinable relation' "
			 "for labels only for the two abstractions that will "
			 "be merged next and reduce labels."
			 "all_abstractions: compute the 'combinable relation' "
			 "for labels once for every abstraction and reduce "
			 "labels."
			 "all_abstractions_with_fixpoaint: keep computing the "
			 "'combinable relation' for labels iteratively for all "
			 "abstractions until no more labels can be reduced.",
			 "ALL_ABSTRACTIONS_WITH_FIXPOINT");

  vector<string> label_reduction_system_order;
  label_reduction_system_order.push_back("REGULAR");
  label_reduction_system_order.push_back("REVERSE");
  label_reduction_system_order.push_back("RANDOM");
  parser.add_enum_option("label_reduction_system_order", label_reduction_system_order,
			 "order of transition systems for the label reduction methods "
			 "that iterate over the set of all abstractions. only useful "
			 "for the choices all_abstractions and all_abstractions_with_fixpoint "
			 "for the option label_reduction_method.", "RANDOM");
  parser.add_option<bool>("expensive_statistics",
			  "show statistics on \"unique unlabeled edges\" (WARNING: "
			  "these are *very* slow, i.e. too expensive to show by default "
			  "(in terms of time and memory). When this is used, the planner "
			  "prints a big warning on stderr with information on the performance impact. "
			  "Don't use when benchmarking!)",
			  "false");

  parser.add_option<int>("limit_merge",
			 "limit on the number of abstract states after the merge"
			 "By default: 1, does not perform any merge",
			 "1");


  parser.add_option<bool>("use_bisimulation",
			  "If activated, use bisimulation to shrink abstractions before computing the simulation",
			  "true");

  parser.add_option<bool>("intermediate_simulations",
			  "Compute intermediate simulations and use them for shrinking",
			  "false");

  parser.add_option<MergeStrategy *>(
				     "merge_strategy",
				     "merge strategy; choose between merge_linear and merge_dfp",
				     "merge_linear");

}
