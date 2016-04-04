#include "symba_unsat.h"

#include "sym_breadth_first_search.h"
#include "sym_uct_pdbs.h"

#include "../timer.h"
#include "../globals.h"
#include "../debug.h"
#include "../option_parser.h"
#include "../plugin.h"

SymBAUnsat::SymBAUnsat(const Options &opts) : 
    SearchEngine(opts), SymController(opts),
    searchDir(Dir(opts.get_enum("search_dir"))),
    abstractDir(Dir(opts.get_enum("abstract_dir"))),
    mgrParams(opts), searchParams(opts), 
    phTime(opts.get<double> ("ph_time")), phMemory(opts.get<double> ("ph_memory")), 
    maxRelaxTime(opts.get<int> ("max_relax_time")), 
    maxRelaxNodes(opts.get<int> ("max_relax_nodes")), 
    absTRsStrategy (AbsTRsStrategy(opts.get_enum("tr_st"))),
    perimeterPDBs (opts.get<bool>("perimeter_pdbs")), 
    ratioRelaxTime(opts.get<double> ("relax_ratio_time")), 
    ratioRelaxNodes(opts.get<double> ("relax_ratio_nodes")), 
    multiply_time_by (opts.get<double> ("multiply_time_by")), 
    num_fails_to_multiply_time(opts.get<int> ("num_fails_to_multiply_time")),
    time_fails_to_multiply_time(opts.get<int> ("time_fails_to_multiply_time")), 
    shouldAbstractRatio(opts.get<double> ("should_abstract_ratio")), 
    maxNumAbstractions(opts.get<int> ("max_abstractions")), 
    UCT_C(opts.get<double> ("uct_c")),
    rewardType (UCTRewardType(opts.get_enum("reward_type"))), 
    RAVE_K (opts.get<double> ("rave_k")),
    add_abstract_to_ongoing_searches_time (opts.get<int>( "add_abstract_to_ongoing_searches_time")),
    numAbstractions(0), num_iterations_without_reward (0),time_last_reward(g_timer()), 
    time_step_abstract(0), time_step_original(0), 
    time_select_exploration(0),  time_notify_mutex(0), 
    time_init(0) {  
}

UCTNode * SymBAUnsat::getUCTNode (const std::set<int> & pattern) {
    if (!nodesByPattern.count(pattern)) {
	UCTNode * newNode = new UCTNode(pattern);
	nodes.push_back(unique_ptr<UCTNode> (newNode));
	nodesByPattern[pattern] = newNode;

	return newNode;
    }

    return nodesByPattern[pattern];
    
}


void SymBAUnsat::initialize(){
    print_options();
   
    nodes.push_back(unique_ptr<UCTNode> (new UCTNode(vars.get(), mgrParams, cost_type)));
    
    if (searchDir != Dir::BW) {
	ongoing_searches.push_back(getRoot()->initSearch(true, searchParams));
    }

    if (searchDir != Dir::FW){
	ongoing_searches.push_back(getRoot()->initSearch(false, searchParams));
    }
    time_init = g_timer();
}

void SymBAUnsat::insertDeadEnds(BDD bdd, bool isFW) {
    if(isFW){
	dead_end_fw.push_back(bdd);
	getRoot()->getMgr()->mergeBucket(dead_end_fw);
    } else {
	dead_end_bw.push_back(bdd);
	getRoot()->getMgr()->mergeBucket(dead_end_bw);
    }

    //Propagate to symbolic managers
    getRoot()->propagateNewDeadEnds(bdd, isFW);   
}


UCTNode * SymBAUnsat::relax(UCTNode * node,
			    bool fw,  
			    std::vector<UCTNode *> & uct_trace, 
			    bool override_search) {
    Timer t_relax;
    shared_ptr<SymBreadthFirstSearch> searchToRelax = node->retrieveSearch(fw, override_search);
    SymManager * parentMgr = node->getMgr(); 
    
    while(node && node->isAbstractable()) {
	node->initChildren(this);
	
	if (rewardType == UCTRewardType::RAND) {
	    node = node->getRandomChild(fw);
	} else {
	    node = node->getChild(fw, UCT_C, RAVE_K, getRoot());
	}
	
	if(!node) break; 

	uct_trace.push_back(node);

	if (node->getSearch(fw)) {
	    searchToRelax = node->retrieveSearch(fw, override_search);
	    continue;
	}

	if(!node->getMgr()) {	    
	    node->init(vars.get(), mgrParams, parentMgr, absTRsStrategy);
	}
	parentMgr =node->getMgr(); 

	auto res = node->relax(searchToRelax, searchParams, 
			       maxRelaxTime, maxRelaxNodes,
			       ratioRelaxTime, ratioRelaxNodes, 
			       perimeterPDBs);


	if (res) {
	    node->getMgr()->addDeadEndStates(dead_end_fw, dead_end_bw);
	    time_select_exploration += t_relax();
	    return node;
	}
    }

    time_select_exploration += t_relax();
    return nullptr;
}
    

    
int SymBAUnsat::step(){
    
    SymBreadthFirstSearch * currentSearch = selectExploration();

    Timer timer; 
    if(currentSearch){
	currentSearch->step();
    
	if(!currentSearch->isAbstracted()){
	    time_step_original += timer();

	    if (currentSearch->finished()) {
		if (currentSearch->foundSolution()){
		    Plan empty_plan;
		    set_plan(empty_plan);
		    return SOLVED;
		} else {
		    return FAILED; 
		} 
	    }	
	}else{
	    if (currentSearch->finished()) {
		notifyFinishedAbstractSearch(currentSearch);
		auto it = std::find(begin(ongoing_searches), 
				    end(ongoing_searches), 
				    currentSearch);

		ongoing_searches.erase(it);
	    }

	    time_step_abstract += timer();	
	}
    }
    return IN_PROGRESS;
}

void SymBAUnsat::increase_num_iterations_without_reward(){
    num_iterations_without_reward ++;
    if(num_iterations_without_reward > num_fails_to_multiply_time && 
       g_timer() - time_last_reward > time_fails_to_multiply_time) {
	searchParams.maxStepNodesMin *= multiply_time_by;
	cout << "Increased node bound: " << searchParams.getMaxStepNodes() << endl;
	num_iterations_without_reward = 0;
	time_last_reward = g_timer();
    }
}

void SymBAUnsat::notifyFinishedAbstractSearch(SymBreadthFirstSearch * currentSearch, double time_spent, 
					      const vector<UCTNode *> & uct_trace){
    Timer t_notify;
    cout << "Finished " << ((currentSearch->isFW()? "fw" : "bw"))  << " in " << *(currentSearch->getAbstraction())<< ": " << flush;
    if (!currentSearch->foundSolution()){
	cout << "proved task unsolvable!" << endl;
	statistics(); 
	exit_with(EXIT_UNSOLVABLE); 
    } else {
	double multiplier = 1.0;
	BDD newDeadEnds = currentSearch->getUnreachableStates();
	//cout << " deadends " <<  vars->percentageNumStates(newDeadEnds) << flush; 
	try {
	    
	    newDeadEnds = getRoot()->getMgr()->filter_mutex(newDeadEnds, !currentSearch->isFW(), 1000000, true);
	    // cout << "  removed:  " << newDeadEnds.nodeCount() << endl;
	    
	} catch(BDDError e) {
	    cout << "  could not remove other dead ends. " << flush;
	    multiplier = 0.1; //Reduce reward by a factor of 10
	}

	newDeadEnds = newDeadEnds*vars->validStates();
	    
       	getRoot()->removeSubsets(currentSearch->getAbstraction()->getFullVars(), currentSearch->isFW());
	
	if (!newDeadEnds.IsZero()) {
	    cout << "  found dead ends: " << newDeadEnds.nodeCount() << endl; 

	    insertDeadEnds(newDeadEnds, currentSearch->isFW());
	    num_iterations_without_reward = 0;
	    time_last_reward = g_timer();
	} else {
	    cout <<  "  with no results "  << endl; 
	    increase_num_iterations_without_reward();
	}

	double reward = computeReward(newDeadEnds, time_spent)*multiplier;
	DEBUG_MSG(cout << "Reward: " << reward << endl;);

	for (UCTNode * node : uct_trace) {
	    if(RAVE_K) {
		node->notifyRewardRAVE(currentSearch->isFW(), reward, uct_trace.back()->getPattern());
	    }else{
		node->notifyReward(currentSearch->isFW(), reward);
	    }	    
	}
    }
    
    time_notify_mutex += t_notify();
} 


double SymBAUnsat::computeReward (const BDD & bdd, double time_spent) const {
    if(bdd.IsZero()) return 0;
    switch(rewardType) {
    case STATES_NODES:
	if (bdd.nodeCount() < 100000) {
	    return std::max<double>(1.0, 10*vars->percentageNumStates(bdd));
	} else {
	    return std::max <double> (0.25, (10000000 - bdd.nodeCount())/10000000.0);
	}
	
    case STATES: 
	return vars->percentageNumStates(bdd);
    case NODES:
	return  std::min(std::max <double> (0.5, bdd.nodeCount()/100000.0), 10.0);
    case STATES_TIME:
	return vars->percentageNumStates(bdd) * 1800.0/(1 + time_spent); 
    case NODES_TIME:
	return (bdd.nodeCount()/100000.0) * 1800.0/(1+ time_spent);
    case NONE: 
    case RAND:
	return 0;
    }
    return 0;
} 

SymBreadthFirstSearch * SymBAUnsat::selectExploration() {
    //DEBUG_PHPDBS(cout << "We have " << ongoing_searches.size() << " ongoing_searches" << endl;);
    //1) Look in already generated explorations => get the easiest one
    //(gives preference to shouldSearch abstractions)
    std::sort(begin(ongoing_searches), end(ongoing_searches),
    	      [this] (SymBreadthFirstSearch * e1, 
		      SymBreadthFirstSearch * e2){
    		  return e1->isBetter (*e2); 
    	      });

    for(auto exp : ongoing_searches){
	if(exp->isSearchable()){
	    return exp;
	}
    }
    
    //3) Ask hierarchy policies to generate new heuristics/explorations
    if (askHeuristic()) return nullptr;
    else return ongoing_searches.front(); 
}


bool SymBAUnsat::chooseDirection() const {
    if (abstractDir ==  Dir::FW) {
	return  true;
    }else if (abstractDir == Dir::BW) {
	return false;	
    }

    return getRoot()->chooseDirection(UCT_C);
}

bool SymBAUnsat::askHeuristic() {
    Timer t_gen_heuristic;
    //cout << "Ask heuristic" << endl;
    while(t_gen_heuristic() < phTime && 
	  vars->totalMemory() < phMemory &&
	  numAbstractions < maxNumAbstractions){	
	numAbstractions++;

	vector<UCTNode *> uct_trace;
	uct_trace.push_back(getRoot());

	bool fw = chooseDirection();
	    
	//1) Generate a new abstract exploration
	UCTNode * abstractNode = relax(getRoot(), fw, uct_trace, false);
	if(!abstractNode) {
	    increase_num_iterations_without_reward(); 
	    if(RAVE_K) for (auto node : uct_trace) node->notifyRewardRAVE(fw, 0, uct_trace.back()->getPattern());
	    else for (auto node : uct_trace) node->notifyReward(fw, 0);
	    continue;
	}

	SymBreadthFirstSearch * abstractExp = abstractNode->getSearch(fw);
	//2) Search the new exploration
	while(abstractExp &&
	      !abstractExp->finished() &&
	      vars->totalMemory() < phMemory && 
	      t_gen_heuristic() < phTime){
	    if (abstractExp->isSearchable()){
		Timer t_step;
		abstractExp->step();
		time_step_abstract += t_step();
		
		
		if(abstractExp->finished()) {
		    notifyFinishedAbstractSearch(abstractExp, t_gen_heuristic(), uct_trace); 
		    return true;
		}
	    } else {
		bool override_search =  g_timer() < add_abstract_to_ongoing_searches_time;
		
		if(!override_search) {
		    ongoing_searches.push_back(abstractExp); //Store ongoing searches
		} 
		//If we cannot continue the search, we relax it even more
		abstractNode = relax(abstractNode, fw, uct_trace, override_search);
		if(!abstractNode){
		    increase_num_iterations_without_reward();
		    if (RAVE_K) for (auto node : uct_trace) node->notifyRewardRAVE(fw, 0, uct_trace.back()->getPattern());
		    else for (auto node : uct_trace) node->notifyReward(fw, 0);
		    return true;
		}
		abstractExp = abstractNode->getSearch(fw);
	    }
	}
    }
    
    //I did not generate any heuristic
    return false;
}

void SymBAUnsat::print_options() const{
    cout << "SymBAUnsat* " << endl;
    cout << "   Search dir: " << searchDir <<   cout << endl;
    cout << "   Abstract dir: " << abstractDir <<   cout << endl;

    cout << "  Max num abstractions: " << maxNumAbstractions << endl;
    cout << "   Abs TRs Strategy: " << absTRsStrategy << endl;
    cout << "   PH time: " << phTime << ", memory: " << phMemory << endl;
    cout << "   Relax time: " << maxRelaxTime << ", nodes: " << maxRelaxNodes << endl;
    cout << "   Ratio relax time: " <<  ratioRelaxTime << ", nodes: " << ratioRelaxNodes << endl;
    cout << "   Perimeter Abstractions: " << (perimeterPDBs ? "yes" : "no") << endl;
    cout << "   ShouldAbstract ratio: " << shouldAbstractRatio << endl;

    mgrParams.print_options();
    searchParams.print_options();
}

void SymBAUnsat::statistics() const{
    cout << endl << "Total BDD Nodes: " << vars->totalNodes() << endl;
    cout << "Initialization time: " << time_init << endl;
    cout << "Total time spent in original search: " << time_step_original << "s" << endl;
    cout << "Total time spent in abstract searches: " << time_step_abstract<<"s" <<   endl;
    cout << "Total time spent relaxing: " << time_select_exploration << "s" << endl;
    cout << "Total time spent notifying mutexes: " << time_notify_mutex << "s" << endl;

    cout << "Total time: " << g_timer() << "s" << endl;

}

static SearchEngine *_parse(OptionParser &parser) {

    SearchEngine::add_options_to_parser(parser);
    SymController::add_options_to_parser(parser, 45e3, 1e7);

    parser.add_enum_option("search_dir", DirValues,
			   "search direction", "BIDIR");

    parser.add_enum_option("abstract_dir", DirValues,
			   "search direction in abstract searches", "BIDIR");

    SymParamsMgr::add_options_to_parser(parser);
    SymParamsSearch::add_options_to_parser(parser, 30e3, 1e7);

    parser.add_option<int>("max_abstractions",
			   "maximum number of calls to askHeuristic", "100000");

    parser.add_option<double>("ph_time", 
			      "allowed time to use the ph", "500");
    parser.add_option<double>("ph_memory",
			      "allowed memory to use the ph", to_string(3.0e9));

    parser.add_option<int>("max_relax_time",
			   "allowed time to relax the search", to_string(10e3));
    parser.add_option<int>("max_relax_nodes",
			   "allowed nodes to relax the search", to_string(10e7));
    parser.add_option<double>("relax_ratio_time",
			      "allowed time to accept the abstraction after relaxing the search.", 
			      "0.5");
    parser.add_option<double>("relax_ratio_nodes",
			      "allowed nodes to accept the abstraction after relaxing the search.", "0.5");


    parser.add_option<double>("multiply_time_by",
			      "", 
			      "1.5");

    parser.add_option<int>("num_fails_to_multiply_time",
			   "", "10");

    parser.add_option<int>("time_fails_to_multiply_time",
			   "", "100");

      
    parser.add_enum_option("tr_st", AbsTRsStrategyValues,
			   "abstraction TRs strategy", "IND_TR_SHRINK");
  
    parser.add_option<bool>("perimeter_pdbs",  "initializes explorations with the one being relaxed", "true");

    parser.add_option<double>("should_abstract_ratio",
			      "relax the search when has more than this estimated time/nodes· If it is zero, it abstract the current perimeter (when askHeuristic is called)", "0");

    parser.add_option<double>("ratio_increase", 
			      "maxStepTime is multiplied by ratio to the number of abstractions", "1.5");

    parser.add_option<double>("uct_c", 
			      "constant for uct formula", "1.0");

    parser.add_option<double>("rave_k", 
			      "constant for RAVE formula. Disabled by default", "0.0");

    parser.add_enum_option("reward_type", UCTRewardTypeValues,
			   "type of reward function", "RANDOM");

    parser.add_option<int>("add_abstract_to_ongoing_searches_time", 
			    "includes the abstract searches that were incomplete to the set of ongoing searches", 
			    "0");

    Options opts = parser.parse();

    SearchEngine *policy = 0;
    if (!parser.dry_run()) {
	policy = new SymBAUnsat(opts);
    }  
    return policy;
}

static Plugin<SearchEngine> _plugin("symba_unsat", _parse);
