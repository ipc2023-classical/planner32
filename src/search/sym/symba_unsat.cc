#include "symba_unsat.h"

#include "sym_ph.h" 
#include "sym_breadth_first_search.h"
#include "sym_hnode.h" 
#include "../debug.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../globals.h"
#include "../rng.h"

SymBAUnsat::SymBAUnsat(const Options &opts) : 
    SymEngine(opts),
    currentPH(0), time_step_abstract(0), time_step_original(0), time_select_exploration(0) {  
}

void SymBAUnsat::initialize(){
    print_options();
    SymEngine::initialize();

    if (searchDir != Dir::BW){
	unique_ptr<SymBreadthFirstSearch> fw_orig (new SymBreadthFirstSearch (getSearchParams()));
	fw_orig->init(originalStateSpace->getManager(), true);
	ongoing_searches.push_back(std::move(fw_orig));
    }

    if (searchDir != Dir::FW){
	unique_ptr<SymBreadthFirstSearch> bw_orig (new SymBreadthFirstSearch (getSearchParams()));
	bw_orig->init(originalStateSpace->getManager(), false);
	ongoing_searches.push_back(std::move(bw_orig));
    }
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
	    for(auto ph : phs){
		ph->operate(originalSearch);
	    }
	}else{
	    if (currentSearch->finished()) {
		if (!currentSearch->foundSolution()){
		    return FAILED; 
		} else {
		    BDD newDeadEnds = currentSearch->getUnreachableStates();
		    bool dirFw = currentSearch->isFW();

		    auto it = std::find_if(begin(ongoing_searches), end(ongoing_searches), [&](std::unique_ptr<SymBreadthFirstSearch> const& p) {
			    return p.get() == currentSearch;
			});
		    ongoing_searches.erase(it);
		    insertDeadEnds(newDeadEnds, !dirFw);
		}
	    }

	    time_step_abstract += timer();	
	}
    }
    return IN_PROGRESS;
}

SymBreadthFirstSearch * SymBAUnsat::selectExploration() {
    //DEBUG_PHPDBS(cout << "We have " << ongoing_searches.size() << " ongoing_searches" << endl;);
    //1) Look in already generated explorations => get the easiest one
    //(gives preference to shouldSearch abstractions)
    std::sort(begin(ongoing_searches), end(ongoing_searches),
    	      [this] (const unique_ptr<SymBreadthFirstSearch> & e1, const unique_ptr<SymBreadthFirstSearch> & e2){
    		  return e1->isBetter (*(e2.get())); 
    	      });

    for(auto & exp : ongoing_searches){
	if(exp->isSearchable()){
	    return exp.get();
	}
    }
    
    //Pick one exploration that seems "promising". For now: random
    /*int random = g_rng.next(ongoing_searches.size());
    SymBreadthFirstSearch * searchToAbstract = ongoing_searches[random].get();
    */
    //3) Ask hierarchy policies to generate new heuristics/explorations
    for(int i = 0; i < phs.size(); i++){ //Once per heuristic
	//bool didSomething = phs[currentPH]->askHeuristic(searchToAbstract);
	currentPH ++;
	if(currentPH >= phs.size()){
	    currentPH = 0;
	}
    }

    return ongoing_searches.front().get(); 
}

void SymBAUnsat::print_options() const{
    cout << "SymBAUnsat* " << endl;
    cout << "   Search dir: " << searchDir <<   cout << endl;
    for(auto ph : phs){
	ph->dump_options();
    }
}

void SymBAUnsat::statistics() const{
    SymEngine::statistics();
    
    cout << "Total time spent in original search: " << time_step_original << endl;
    cout << "Total time spent in abstract searches: " << time_step_abstract<<  endl;
    cout << "Total time spent relaxing: " << time_select_exploration << endl;
}

static SearchEngine *_parse(OptionParser &parser) {
    SymEngine::add_options_to_parser(parser);
    Options opts = parser.parse();
  
    SearchEngine *policy = 0;
    if (!parser.dry_run()) {
	policy = new SymBAUnsat(opts);
    }  
    return policy;
}

static Plugin<SearchEngine> _plugin("symba_unsat", _parse);
