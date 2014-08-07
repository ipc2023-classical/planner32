#include "merge_criterion.h"

#include "../causal_graph.h"
#include "../globals.h"
#include "abstraction.h" 
#include "../option_parser.h"
#include "../option_parser_util.h"
#include "../plugin.h"
#include <cassert>
#include <cstdlib>
#include <vector>
#include <map>

using namespace std;

void MergeCriterionCG::init(){
    is_causal_predecessor.resize(g_variable_domain.size(), false);
}
void MergeCriterionCG::select_next(int var_no) {
    if(allow_incremental){
	const vector<int> &new_vars = g_causal_graph->get_eff_to_pre(var_no);
	for (int i = 0; i < new_vars.size(); ++i)
	    is_causal_predecessor[new_vars[i]] = true;
    }
}
void MergeCriterionCG::filter(const std::vector<Abstraction *> & /*all_abstractions*/,
			      std::vector <int> & vars, Abstraction * abstraction) {
    if(!abstraction) return;
    if (!allow_incremental){
	is_causal_predecessor.resize(g_variable_domain.size(), false);
	const vector<int> & varset = abstraction->get_varset();
	for(int v = 0; v < varset.size(); v++){
	    const vector<int> &new_vars = g_causal_graph->get_eff_to_pre(varset[v]);
	    for (int i = 0; i < new_vars.size(); ++i){
		is_causal_predecessor[new_vars[i]] = true;
	    }
	}
    }

    // cout << "variables Relevant: ";
    // for(int i = 0; i < is_causal_predecessor.size(); i++) if(is_causal_predecessor[i]) cout << i << " ";
    // cout << endl;
    MergeCriterion::filter(vars, is_causal_predecessor);
}


void MergeCriterionGoal::init(){
    is_goal_variable.resize(g_variable_domain.size(), false);
    for (int i = 0; i < g_goal.size(); ++i){
	is_goal_variable[g_goal[i].first] = true;
    }
}
void MergeCriterionGoal::select_next(int /*var_no*/) {}
void MergeCriterionGoal::filter(const std::vector<Abstraction *> &/*all_abstractions*/, 
				std::vector <int> & vars, Abstraction * /*abstraction*/) {
    MergeCriterion::filter(vars, is_goal_variable);
}

void MergeCriterionMinSCC::init(){
    is_causal_predecessor.resize(g_variable_domain.size(), false);
    const vector<vector<int> > & graph = 
	(complete_cg ? (reverse ? g_causal_graph->get_predecessors() : 
			          g_causal_graph->get_successors()) : 
	 (reverse ? g_causal_graph->get_eff_to_pre() : 
                    g_causal_graph->get_pre_to_eff()));

    scc = new SCC(graph);
    // for(int i = 0; i < scc->get_sccs().size(); i++){
    //   cout << "SCC ("<< scc->get_scc_layer()[i] << "): ";
    //   for(int j = 0; j < scc->get_sccs()[i].size(); j++){
    //     cout << scc->get_sccs()[i][j] << " ";
    //   }
    //   cout << endl;
    // }
}

void MergeCriterionMinSCC::select_next(int var_no) {
    if(allow_incremental){
	const vector<int> & new_vars = g_causal_graph->get_eff_to_pre(var_no);
	for (int i = 0; i < new_vars.size(); ++i)
	    is_causal_predecessor[new_vars[i]] = true;
    }
}

void MergeCriterionMinSCC::forbid_scc_descendants(int scc_index,
						  const vector<set<int> > & scc_graph, 
						  vector<bool> & forbidden_sccs) const{
    const set<int> & descendants = scc_graph[scc_index];
    for(set<int>::iterator it = descendants.begin(); 
	it != descendants.end(); ++it){
	if (!forbidden_sccs[*it]){
	    forbidden_sccs [*it] = true;
	    forbid_scc_descendants(*it, scc_graph, forbidden_sccs);
	}
    }
}

void MergeCriterionMinSCC::filter(const std::vector<Abstraction *> &/*all_abstractions*/, 
				  std::vector <int> & vars,
				  Abstraction * abstraction) {
    if(!abstraction) return;
    if (!allow_incremental){
	is_causal_predecessor.resize(g_variable_domain.size(), false);
        const vector<int> & varset = abstraction->get_varset();
	for(int v = 0; v < varset.size(); v++){
	    const vector<int> &new_vars = g_causal_graph->get_eff_to_pre(varset[v]);
	    for (int i = 0; i < new_vars.size(); ++i){
		is_causal_predecessor[new_vars[i]] = true;
	    }
	}
    }

    if(!MergeCriterion::filter(vars, is_causal_predecessor)){
	return; //No CG relevant vars => we do not prefer any variable over another
    }

    const vector<set<int> > & scc_graph = scc->get_scc_graph();
    const vector<int> & vars_scc = scc->get_vertex_scc();
    vector<bool> forbidden_sccs (scc_graph.size(), false);
    //In each SCC,select only the variable whose "level" is "closer to the root"
    //We consider a variable closer to the root if it has lower id 
    //If reverse is activated, we consider the opposite order
    map<int, int> vars_by_scc;
    //1) forbid all sccs pointed by scc_var
    for (int i = 0; i < vars.size(); i++){
	int var = vars[i];
	int scc_var = vars_scc [var];
	if(!forbidden_sccs[scc_var]){
	    forbid_scc_descendants(scc_var, scc_graph, forbidden_sccs);
	    if(!vars_by_scc.count(scc_var) || (!reverse && vars_by_scc[scc_var] > var) || (reverse && vars_by_scc[scc_var] < var)){
		vars_by_scc [scc_var] = var;
	    }
	}
    }

    //2) Filter all variables whose scc has been forbidden.
    vector<int> new_vars; 
    if(tie_by_level){ //For valid sccs, include the selected variable
	for(map<int, int>::iterator it = vars_by_scc.begin(); 
	    it != vars_by_scc.end(); ++it){
	    if(!forbidden_sccs[it->first]){
		new_vars.push_back(it->second);
	    }
	}
    }else{
	//For valid sccs, include all variables
	for (int i = 0; i < vars.size(); i++){ 
	    int var = vars[i];
	    int scc_var = vars_scc [vars[i]];
	    if(!forbidden_sccs[scc_var])
		new_vars.push_back(var);
	}
    }
    vars.swap(new_vars);
    /* cout << "Candidate variable after scc: ";
       for(int i = 0; i < vars.size(); i++)
       cout << vars[i] << " ";
       cout << endl;*/
}

MergeCriterionTRs::MergeCriterionTRs(const Options & opts) : 
    only_goals(opts.get<bool> ("goal")), 
    only_empty(opts.get<bool>("empty")){
}


void MergeCriterionTRs::init(){}
void MergeCriterionTRs::select_next(int /*var_no*/) {}
void MergeCriterionTRs::filter(const std::vector<Abstraction *> &all_abstractions, 
			       std::vector <int> & vars, Abstraction * abstraction) {
    if(abstraction){ //Check if abstraction exists, because the first variable is selected without abstraction    
	vector<int> score;
	abstraction->normalize();
	abstraction->count_transitions(all_abstractions, vars, 
				       only_empty, only_goals, score);
	//for(int i = 0; i < vars.size(); i++)
	//    cout << "ScoreTRs(" << only_empty << ", " << only_goals << ") " << vars[i] << ": " << score[vars[i]] << " " << endl;
	MergeCriterion::filter_best(vars, score, false);
    }
}

void MergeCriterionRelevant::init(){
    MergeCriterionCG::init();
    for (int i = 0; i < g_goal.size(); ++i){
	is_causal_predecessor[g_goal[i].first] = true;
    }
}

template <class T> 
static MergeCriterion *_parse(OptionParser & parser) {
    if (parser.dry_run()) {
	return 0;
    } else {
	return new T();
    }
}

//TODO: This should be removed, because it is equivalent to scc(level=false)
static MergeCriterion *_parse_scc_no_level(OptionParser & parser) {
    if (parser.dry_run()) {
	return 0;
    } else {
	return new MergeCriterionMinSCC(false, false);
    }
}


MergeCriterionMinSCC::MergeCriterionMinSCC(const Options & opts) : 
    reverse(opts.get<bool> ("reverse")), 
    tie_by_level(opts.get<bool>("level")), 
    complete_cg(opts.get<bool>("complete_cg")){
}
	    
static MergeCriterion *_parse_scc(OptionParser & parser) {
    parser.add_option<bool>("reverse", "reverse scc criterion", "false");
    parser.add_option<bool>("level", "use level or not in the scc criterion", "false");
    parser.add_option<bool>("complete_cg", "use the old or the new cg", "false");
    Options opts = parser.parse();

    if (parser.dry_run()) {
	return 0;
    } else {
	return new MergeCriterionMinSCC(opts);
    }
}

static MergeCriterion *_parse_tr(OptionParser & parser) {
    parser.add_option<bool>("goal", "only counts transitions leading to a goal state", "false");
    parser.add_option<bool>("empty", "only counts transitions that will become empty", "false");
    Options opts = parser.parse();

    if (parser.dry_run()) {
	return 0;
    } else {
	return new MergeCriterionTRs(opts);
    }
}

static Plugin<MergeCriterion> _plugin_cg("cg",   _parse<MergeCriterionCG>);
static Plugin<MergeCriterion> _plugin_goal("goal",   _parse<MergeCriterionGoal>);
static Plugin<MergeCriterion> _plugin_scc("scc",   _parse_scc);
static Plugin<MergeCriterion> _plugin_scc_no_level("scc_no_level",   _parse_scc_no_level);
static Plugin<MergeCriterion> _plugin_tr("tr", _parse_tr);
static Plugin<MergeCriterion> _plugin_rel("relevant",   _parse<MergeCriterionRelevant>);
