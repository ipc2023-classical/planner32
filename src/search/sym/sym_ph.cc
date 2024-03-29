#include "sym_ph.h"

#include "sym_engine.h" 
#include "sym_bdexp.h"
#include "sym_manager.h"
#include "sym_hnode.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "sym_pdb.h"

SymPH::SymPH(const Options & opts) : 
    vars(nullptr), mgrParams(opts), searchParams(opts), 
    phTime(opts.get<double> ("ph_time")), phMemory(opts.get<double> ("ph_memory")), 
    maxRelaxTime(opts.get<int> ("max_relax_time")), 
    maxRelaxNodes(opts.get<int> ("max_relax_nodes")), 
    absTRsStrategy (AbsTRsStrategy(opts.get_enum("tr_st"))),
    perimeterPDBs (opts.get<bool>("perimeter_pdbs")), 
    relaxDir(RelaxDirStrategy(opts.get_enum("relax_dir"))),
    ratioRelaxTime(opts.get<double> ("relax_ratio_time")), 
    ratioRelaxNodes(opts.get<double> ("relax_ratio_nodes")),
    use_mutex_in_abstraction(opts.get<bool> ("use_mutex_in_abstraction")), 
    shouldAbstractRatio(opts.get<double> ("should_abstract_ratio")), 
    maxNumAbstractions(opts.get<int> ("max_abstractions")),
    numAbstractions(0), ignore_if_useful(false)
{
    dump_options();
}

bool SymPH::init(SymController * eng, SymVariables * v, SymManager * mgr){
    engine = eng;
    vars = v;
    if(use_mutex_in_abstraction){
	const auto & nmFw = mgr->getNotMutexBDDs(true);
	const auto & nmBw = mgr->getNotMutexBDDs(false);
	notMutexBDDs.insert(end(notMutexBDDs), begin(nmFw), end(nmFw));
	notMutexBDDs.insert(end(notMutexBDDs), begin(nmBw), end(nmBw));
    }
    return init();
}


bool SymPH::askHeuristic(SymBDExp * originalSearch, double allotedTime){
    DEBUG_PHPDBS(cout << "Ask heuristic" << endl;);
    Timer t_gen_heuristic;
    allotedTime = min<double> (allotedTime, phTime);

    SymBDExp * bdExp = originalSearch;
    assert (bdExp->getFw()->getClosed()->hasEvalOrig() || 
	    bdExp->getBw()->getClosed()->hasEvalOrig());
    
    if(!expPerimeters.empty()){  //Use the other perimeter instead
	DEBUG_PHPDBS(cout << ">> Reusing stored perimeter"<< endl;);
	bdExp = expPerimeters[0].get(); 
    }

  
    while(t_gen_heuristic() < allotedTime && 
	  vars->totalMemory() < phMemory && 
	  !originalSearch->isSearchable() && 
	  numAbstractions < maxNumAbstractions){	
	numAbstractions++;

	//1) Generate a new abstract exploration
	SymBDExp * abstractExp = relax(bdExp);

	int num_relaxations = 1;

	//2) Search the new exploration
	while(abstractExp &&
	      (!abstractExp->finishedMainDiagonal() || abstractExp->isUseful()) &&
	      vars->totalMemory() < phMemory && 
	      t_gen_heuristic() < allotedTime){

	    auto * selectedExp = abstractExp->selectBestDirection(false);
	    
	    if (selectedExp->isSearchable()){
		selectedExp->step();
	    } else {
		if (!abstractExp->isSearchable() && 
		    (abstractExp->getF() < originalSearch->getF() || 
		     !abstractExp->finishedMainDiagonal())){
		    // We skip this abstract state space if it is not
		    // searchable or we have not finished the main
		    // diagonal
		    abstractExp->desactivate();
		    numAbstractions --; //Does not count
		    //TODO: Here we should liberate all the memory
		    //used by this abstraction
		} 

		if(numAbstractions++ >= maxNumAbstractions ) {
		    break; 
		}
	
		//If we cannot continue the search, we relax it even more
		cout << "We cannot continue the search so we relax it even more" << endl;
		abstractExp = relax(bdExp, abstractExp->getHNode(), ++num_relaxations);			
	    }
	}
    
	string reason = !abstractExp ? ">> No abstract search remaining" : 
	    (abstractExp->finishedMainDiagonal() ? (abstractExp->finished() ? "totally explored" : 
						    (!abstractExp->isUseful() ? "is not useful" :
						     "f > concrete_f")) : 
	     "exhausted resources");

	cout << "Finished the exploration of the abstract state space (" <<  reason
	     << "): " << t_gen_heuristic() << "s spent of " << allotedTime << endl;
  
	//3) Add last heuristic 
	if(abstractExp){
	    DEBUG_PHPDBS(cout << "Set heuristic to abstract exp: " << *abstractExp << endl;);
	    bdExp->setHeuristic(*abstractExp);
      
	}else{
	    DEBUG_PHPDBS(cout << "Adding explicit heuristic to bw: "<< 
			 intermediate_heuristics_fw.size() << endl;);
	    //Add explicit heuristic
	    for(auto & inth : intermediate_heuristics_fw){
		originalSearch->getFw()->addHeuristic(make_shared<SymHeuristic>(*vars,inth));
	    }
    
	    DEBUG_PHPDBS(cout << "Adding explicit heuristic to fw: " <<
			 intermediate_heuristics_bw.size() << endl;);
	    for(auto & inth : intermediate_heuristics_bw){
		originalSearch->getBw()->addHeuristic(make_shared<SymHeuristic> (*vars,inth));
	    }
	}
    }

    DEBUG_PHPDBS(cout << *originalSearch << endl;);
    //I did not generate any heuristic
    return false;
}

void SymPH::operate(SymBDExp * originalSearch) {
    int nextStepNodes = max(originalSearch->getFw()->nextStepNodes(),
			    originalSearch->getBw()->nextStepNodes());
    if(expPerimeters.empty() && shouldAbstractRatio && 
       (nextStepNodes > searchParams.maxStepNodes*shouldAbstractRatio)){
	expPerimeters.push_back(createBDExp(getDir(originalSearch), originalSearch));
    }     
}

SymBDExp * SymPH::addHeuristicExploration(SymBDExp * oldExp,
					  unique_ptr<SymBDExp> && newExp){
    if(newExp){
	// Needed so that the abstract heuristic starts informing as
	// soon as possible (and to know whether it is useful)
	oldExp->setHeuristic(*newExp);  
	SymBDExp * ptr = newExp.get();
	SymHNode * hnode = newExp->getHNode();
	hnode->add_exploration(std::move(newExp));
	return ptr;
    }else{
	return nullptr;
    }
}

unique_ptr<SymBDExp> SymPH::createBDExp (Dir dir, SymBDExp * bdExp) const{
    return unique_ptr<SymBDExp> (new SymBDExp(bdExp, searchParams, dir));
}

bool SymPH::relax_in(SymBDExp * bdExp, unique_ptr<SymBDExp> & newExp, 
		     SymHNode * hNode, int num_relaxations) const{
  
    Timer relax_timer; //TODO: remove. Only used for debugging

    if(!hNode->empty()){ //Do not repeat the same hnode twice
	//hNode->notuseful_exploration(bdExp);
	return false;
    }

    if (bdExp->getFw() != newExp->getFw()->getParent() && 
	bdExp->getBw() != newExp->getBw()->getParent()){
	cerr << "Assertion error" << endl;
	exit(-1);
    }
  
    cout << ">> Abstract in hNode: " << *hNode << " total time: " << g_timer << endl;
    //I have received a hNode and does not have an exploration. Try.
    if(newExp->initFrontier(hNode, maxRelaxTime, maxRelaxNodes)){
	//Ok, I relaxed the frontier!
	//Check if it is useful
	DEBUG_PHPDBS(cout << "Frontier initialized. total time: " << g_timer << endl;);
	if(!ignore_if_useful && !newExp->isUsefulAfterRelax(searchParams.ratioUseful)){
	    DEBUG_PHPDBS(cout << " >> New exploration is not useful" << *hNode << " total time: " << g_timer << endl;);
	    DEBUG_PHPDBS(cout << "Time failed relaxation: " << relax_timer << endl;);
	    hNode->notuseful_exploration(bdExp);
	    return false;
	} else if(newExp->isSearchableAfterRelax(num_relaxations)){
	    DEBUG_PHPDBS(cout << "New exp is searchable. total time: " << g_timer << endl;);
	    if(!perimeterPDBs){
		newExp.reset(new SymBDExp(engine, searchParams, getDir(bdExp)));
		return newExp->initFrontier(hNode, maxRelaxTime, maxRelaxNodes) &&
		    newExp->initAll(maxRelaxTime, maxRelaxNodes);

	    }else if(newExp->initAll(maxRelaxTime, maxRelaxNodes)){
		DEBUG_PHPDBS(cout << "New exp initialized. total time: " << g_timer << endl;);
		DEBUG_PHPDBS(cout << "Time successful relaxation: " << relax_timer << endl;);

		return true; 
	    }else{
		DEBUG_PHPDBS(cout << " >> Could not initAll the new exploration. total time: " << g_timer << endl;);
	    }
	}else{
	    DEBUG_PHPDBS(cout << " >> Not searchable exploration: " <<
			 *(newExp.get())  << " total time: " << g_timer << endl;);
	    //return false; // If the exploration is not searchable, we do not say we have failed
	    // Changed so that we say that we have failed (simplify the abstraction strategy)
	}
    }else{
	DEBUG_PHPDBS(cout << " >> Could not initFrontier the new exploration. total time: " << g_timer << endl;);
    }

    DEBUG_PHPDBS(cout << "Time failed relaxation: " << relax_timer << endl;);
    hNode->failed_exploration(bdExp);
    return false;
}

void SymPH::add_options_to_parser(OptionParser & parser,
				  const string & default_tr_st, 
				  int abstraction_limit){
  SymParamsMgr::add_options_to_parser(parser);
  SymParamsSearch::add_options_to_parser(parser, 30e3, 1e7);

  parser.add_option<int>("max_abstractions",
			 "maximum number of calls to askHeuristic", to_string(abstraction_limit));

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
			    "0.75");
  parser.add_option<double>("relax_ratio_nodes",
			 "allowed nodes to accept the abstraction after relaxing the search.", "0.75");
  
  parser.add_enum_option("tr_st", AbsTRsStrategyValues,
			 "abstraction TRs strategy", default_tr_st);
  parser.add_enum_option("relax_dir", RelaxDirStrategyValues,
			 "direction allowed to relax",  "BIDIR");

			    parser.add_option<bool>("perimeter_pdbs",  "initializes explorations with the one being relaxed", "true");

  parser.add_option<bool>("use_mutex_in_abstraction", 
			  "uses mutex to prune abstract states in the abstraction procedure", "true");

			    parser.add_option<double>("should_abstract_ratio", "relax the search when has more than this estimated time/nodes· If it is zero, it abstract the current perimeter (when askHeuristic is called)", "0");
  parser.add_option<double>("ratio_increase", 
			    "maxStepTime is multiplied by ratio to the number of abstractions", "2");

}

void SymPH::dump_options() const {
  cout << "  Max num abstractions: " << maxNumAbstractions << endl;
  cout << "   Abs TRs Strategy: " << absTRsStrategy << endl;
  cout << "   PH time: " << phTime << ", memory: " << phMemory << endl;
  cout << "   Relax time: " << maxRelaxTime << ", nodes: " << maxRelaxNodes << endl;
  cout << "   Ratio relax time: " <<  ratioRelaxTime << ", nodes: " << ratioRelaxNodes << endl;
  cout << "   Perimeter Abstractions: " << (perimeterPDBs ? "yes" : "no") << endl;
  cout << "   Relax dir: " << relaxDir << endl;
  cout << "   ShouldAbstract ratio: " << shouldAbstractRatio << endl;

  mgrParams.print_options();
  searchParams.print_options();
}


//Select direction of the new BDExp based on relaxDir
Dir SymPH::getDir(SymBDExp * bdExp) const {
  switch (relaxDir){
  case RelaxDirStrategy::FW: 
    return Dir::FW;
  case RelaxDirStrategy::BW: 
    return Dir::BW;
  case RelaxDirStrategy::BIDIR: 
	if (!bdExp->getFw()->getClosed()->hasEvalOrig())
	    //FW search unfeasible. Use abstractions only in fw direction to inform bw search
	    return Dir::FW; 
	else if (!bdExp->getBw()->getClosed()->hasEvalOrig()) 
           //BW search unfeasible. Use abstractions only in bw direction to inform bw search
	    return Dir::BW;
	else return Dir::BIDIR;
  case RelaxDirStrategy::SWITCHBACK:
    if(bdExp->getDir() == Dir::FW){
      return Dir::BW;
    }else if(bdExp->getDir() == Dir::BW){
      return Dir::FW;	
    }else{
      cerr<< "Cannot use Switchback with bidirectional searches" << endl;
      exit(-1);
    }
  default: 
    cerr<< "Unkown RelaxDirStrategy" << endl;
    exit(-1);
  }
}


// double SymPH::getMaxStepTime() const{
//   return searchParams.maxStepTime * pow(ratioIncrease, ((numAbstractions+1)/5) - 1);
// }

// double SymPH::getMaxAfterRelaxTime() const{
//   return getMaxStepTime()*ratioRelaxTime;
// }

// double SymPH::getMaxStepNodes() const{
//   cout << "MAX AFTER RELAX: " <<  searchParams.maxStepNodes * (1 + pow(numAbstractions, ratioIncrease)) << endl;
//   return searchParams.maxStepNodes * (1 + pow(numAbstractions, ratioIncrease));
// }

// double SymPH::getMaxAfterRelaxNodes() const{
//   return getMaxStepNodes()*ratioRelaxNodes;
// }
