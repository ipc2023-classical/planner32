#include "symba_unsat.h"

#include "sym_breadth_first_search.h"
#include "sym_uct_pdbs.h"

#include "../debug.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../globals.h"
#include "../rng.h"


SymBAUnsat::SymBAUnsat(const Options &opts) : 
    SearchEngine(opts), SymController(opts),
    searchDir(Dir(opts.get_enum("search_dir"))),
    mgrParams(opts), searchParams(opts), 
    phTime(opts.get<double> ("ph_time")), phMemory(opts.get<double> ("ph_memory")), 
    maxRelaxTime(opts.get<int> ("max_relax_time")), 
    maxRelaxNodes(opts.get<int> ("max_relax_nodes")), 
    absTRsStrategy (AbsTRsStrategy(opts.get_enum("tr_st"))),
    perimeterPDBs (opts.get<bool>("perimeter_pdbs")), 
    ratioRelaxTime(opts.get<double> ("relax_ratio_time")), 
    ratioRelaxNodes(opts.get<double> ("relax_ratio_nodes")),
    use_mutex_in_abstraction(opts.get<bool> ("use_mutex_in_abstraction")), 
    shouldAbstractRatio(opts.get<double> ("should_abstract_ratio")), 
    maxNumAbstractions(opts.get<int> ("max_abstractions")), 
    UCT_C(opts.get<double> ("uct_c")),
    numAbstractions(0), 
    time_step_abstract(0),
    time_step_original(0), 
    time_select_exploration(0) {  
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
   
    nodes.push_back(unique_ptr<UCTNode> (new UCTNode(vars.get(), mgrParams)));
    
    if (searchDir != Dir::BW) {
	ongoing_searches.push_back(getRoot()->initSearch(true, searchParams));
    }

    if (searchDir != Dir::FW){
	ongoing_searches.push_back(getRoot()->initSearch(false, searchParams));
    }
}

void SymBAUnsat::insertDeadEnds(BDD bdd, bool isFW) {
    if(isFW){
	dead_end_fw.push_back(bdd);
    } else {
	dead_end_bw.push_back(bdd);
    }

    //Propagate to symbolic managers
    getRoot()->propagateNewDeadEnds(bdd, isFW);   
}


std::pair<UCTNode *, bool> SymBAUnsat::relax() {
    bool fw = false;
    if(g_rng.next31()% 2) fw = true;

    return relax(getRoot(), fw);  
}

std::pair<UCTNode *, bool> SymBAUnsat::relax(UCTNode * node,
					     bool fw) {
    SymBreadthFirstSearch * searchToRelax = node->getSearch(fw);
    SymManager * parentMgr = node->getMgr(); 
    
    while(node && node->isAbstractable()) {
	node->initChildren(this);
	
	node = node->getChild(fw);
	
	if(!node) break; 

	if (node->getSearch(fw)) {
	    searchToRelax = node->getSearch(fw);
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
	    return std::pair<UCTNode *, bool>(node, fw);
	}
    }

    return std::pair<UCTNode *, bool>(nullptr, fw);
}
    

    
int SymBAUnsat::step(){
    Timer timer; 
    SymBreadthFirstSearch * currentSearch = selectExploration();

    time_select_exploration += timer.reset(); 
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


void SymBAUnsat::notifyFinishedAbstractSearch(SymBreadthFirstSearch * currentSearch){
    cout << "Finished abstract search: " << ((currentSearch->isFW()? "fw" : "bw"))  << " in " << *(currentSearch->getAbstraction())<< ": " ;
    if (!currentSearch->foundSolution()){
	cout << "proved task unsolvable!" << endl; 
	exit_with(EXIT_UNSOLVABLE); 
    } else {
	BDD newDeadEnds = currentSearch->getUnreachableStates();

	try {
	    newDeadEnds = getRoot()->getMgr()->filter_mutex(newDeadEnds, !currentSearch->isFW(), 100000, true);
	} catch(BDDError e) {
	    cout << "Could not remove other dead ends ";
	}

	getRoot()->removeSubsets(currentSearch->getAbstraction()->getFullVars(), currentSearch->isFW());

	if (!newDeadEnds.IsZero()) {
	    cout << "found " << vars->numStates(newDeadEnds) << " dead ends." << endl; 

	    insertDeadEnds(newDeadEnds, !currentSearch->isFW());
	
	//getRoot()->notifyReward (vars->numStates(newDeadEnds), 
	//			 currentSearch->getAbstraction()->getFullVars());

	} else {
	    cout <<  "  with no results "  << endl; 

	}
    }

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


bool SymBAUnsat::askHeuristic() {
    Timer t_gen_heuristic;
    //cout << "Ask heuristic" << endl;
    while(t_gen_heuristic() < phTime && 
	  vars->totalMemory() < phMemory &&
	  numAbstractions < maxNumAbstractions){	
	numAbstractions++;

	//1) Generate a new abstract exploration
	auto resRelax = relax();
	UCTNode * abstractNode = resRelax.first;
	if(!abstractNode) return false;
	SymBreadthFirstSearch * abstractExp = abstractNode->getSearch(resRelax.second);
	//2) Search the new exploration
	while(abstractExp &&
	      !abstractExp->finished() &&
	      vars->totalMemory() < phMemory && 
	      t_gen_heuristic() < phTime){
	    
	    if (abstractExp->isSearchable()){
		abstractExp->step();
		if(abstractExp->finished()) {
		    notifyFinishedAbstractSearch(abstractExp); 
		    return true;
		}
	    } else {
		ongoing_searches.push_back(abstractExp); //Store ongoing searches
		//If we cannot continue the search, we relax it even more
		resRelax = relax(resRelax.first, resRelax.second);
		abstractNode = resRelax.first;
		if(!abstractNode) return true;
		abstractExp = abstractNode->getSearch(resRelax.second);
	    }
	}
    }
    
    //I did not generate any heuristic
    return false;
}

void SymBAUnsat::print_options() const{
    cout << "SymBAUnsat* " << endl;
    cout << "   Search dir: " << searchDir <<   cout << endl;

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

    cout << "Total time spent in original search: " << time_step_original << endl;
    cout << "Total time spent in abstract searches: " << time_step_abstract<<  endl;
    cout << "Total time spent relaxing: " << time_select_exploration << endl;
}

static SearchEngine *_parse(OptionParser &parser) {

    SearchEngine::add_options_to_parser(parser);
    SymController::add_options_to_parser(parser, 45e3, 1e7);

    parser.add_enum_option("search_dir", DirValues,
			   "search direction", "BIDIR");

    SymParamsMgr::add_options_to_parser(parser);
    SymParamsSearch::add_options_to_parser(parser, 30e3, 1e7);

    parser.add_option<int>("max_abstractions",
			   "maximum number of calls to askHeuristic", "1000");

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
  
    parser.add_enum_option("tr_st", AbsTRsStrategyValues,
			   "abstraction TRs strategy", "IND_TR_SHRINK");
  
    parser.add_option<bool>("perimeter_pdbs",  "initializes explorations with the one being relaxed", "true");

    parser.add_option<bool>("use_mutex_in_abstraction", 
			    "uses mutex to prune abstract states in the abstraction procedure", "true");

    parser.add_option<double>("should_abstract_ratio", "relax the search when has more than this estimated time/nodes· If it is zero, it abstract the current perimeter (when askHeuristic is called)", "0");
    parser.add_option<double>("ratio_increase", 
			      "maxStepTime is multiplied by ratio to the number of abstractions", "2");

    parser.add_option<double>("uct_c", 
			      "constant for uct formula", "0.5");

    Options opts = parser.parse();


  
    SearchEngine *policy = 0;
    if (!parser.dry_run()) {
	policy = new SymBAUnsat(opts);
    }  
    return policy;
}

static Plugin<SearchEngine> _plugin("symba_unsat", _parse);
