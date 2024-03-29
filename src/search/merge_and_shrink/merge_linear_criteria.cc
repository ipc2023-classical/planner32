#include "merge_linear_criteria.h"

#include "abstraction.h" 
#include "merge_criterion.h"

#include "../globals.h"
#include "../option_parser.h"
#include "../option_parser_util.h"
#include "../plugin.h"
#include <cassert>
#include <cstdlib>
#include <vector>
#include <map>

using namespace std;

MergeLinearCriteria::MergeLinearCriteria(const Options &opts) : 
    MergeStrategy(),
    criteria(opts.get_list<MergeCriterion *> ("criteria")), 
    order(MergeOrder(opts.get_enum("var_order"))) {

  int var_count = g_variable_domain.size();
  if (order == REVERSE_LEVEL) {
    for (int i = 0; i < var_count; ++i)
      remaining_vars.push_back(i);
  } else {
    for (int i = var_count - 1; i >= 0; i--)
      remaining_vars.push_back(i);
  }
  
  if (order == RANDOM){
    random_shuffle(remaining_vars.begin(),
		   remaining_vars.end());
  }

  for(int i = 0; i < criteria.size(); ++i){
    criteria[i]->init();
  }
}

void MergeLinearCriteria::init_strategy (const std::vector <Abstraction * > & abstractions){
    remaining_vars.clear();
    for (auto abs : abstractions) {
	if (abs) {
	    assert(abs->get_varset().size() == 1 || abs == abstractions.back()); 
	    if (abs->get_varset().size() == 1) remaining_vars.push_back(*(abs->get_varset().begin()));
	} 

    }

    if (order == REVERSE_LEVEL) {
	std::sort(begin(remaining_vars), end(remaining_vars));
    } else if (order == RANDOM){
	random_shuffle(remaining_vars.begin(),
		       remaining_vars.end());
    } else {
	std::sort(begin(remaining_vars), end(remaining_vars), std::greater<int>());
    }


    for(int i = 0; i < criteria.size(); ++i){
	criteria[i]->init();
    }
}


MergeLinearCriteria::~MergeLinearCriteria() {
}


pair<int, int> MergeLinearCriteria::get_next(const std::vector<Abstraction *> &all_abstractions,
					     int limit_abstract_states_merge, int min_limit_abstract_states_merge, 
					     int limit_transitions_merge ) {
    assert(!done());

    int first;
    if (all_abstractions.back()->is_atomic()) {
        first = next(all_abstractions);
	
        cout << "First variable: #" << first;
		for(int i = 0; i < g_fact_names[first].size(); ++i) 
	    	cout << " " << g_fact_names[first][i]; 
	assert(all_abstractions[first]);
	cout << endl;
    } else {
        // The most recent composite abstraction is appended at the end of
        // all_abstractions in merge_and_shrink.cc
        first = all_abstractions.size() - 1;
    }
    int second;
    do{
	second = next(all_abstractions, all_abstractions[first], 
		      limit_abstract_states_merge, min_limit_abstract_states_merge, limit_transitions_merge);
	if(second < 0){
	    if(remaining_vars.size() < 2){
		return pair<int,int> (-1, -1);
	    }
	    //Select another variable as first
	    first = next(all_abstractions);
	    cout << "First variable: #" << first;
	    for(int i = 0; i < g_fact_names[first].size(); ++i) 
	    	cout << " " << g_fact_names[first][i]; 
	    cout << endl;
	}
    }while(second < 0);

    cout << "Next variable: #" << second;
    for(int i = 0; i < g_fact_names[second].size(); ++i) 
		cout << " " << g_fact_names[second][i]; 
    cout << endl;

    assert(all_abstractions[first]);
    assert(all_abstractions[second]);
    --remaining_merges;
    return make_pair(first, second);
}

void MergeLinearCriteria::dump_strategy_specific_options() const {
    cout << "Linear merge criteria strategy: ";
    for(int i = 0; i < criteria.size(); i++){
      cout << criteria[i]->get_name() << "_";
    }
    switch(order){
    case LEVEL:
      cout << "LEVEL";
      break;
    case REVERSE_LEVEL:
      cout << "REVERSE_LEVEL";
      break;
    case RANDOM:
      cout << "RANDOM";
      break;
    default:
      cerr << "Unknown merge criterion.";
      abort();
    }

}

string MergeLinearCriteria::name() const {
    return "linear_criteria";
}

bool MergeLinearCriteria::is_linear() const {
    return true;
}

void MergeLinearCriteria::select_next(int var_no) {
  vector<int>::iterator position = find(remaining_vars.begin(), remaining_vars.end(), var_no);
  assert(position != remaining_vars.end());
  remaining_vars.erase(position);
  //selected_vars.push_back(var_no);
  for(int i = 0; i < criteria.size(); ++i){
    criteria[i]->select_next(var_no);
  }
}

int MergeLinearCriteria::next(const std::vector<Abstraction *> &all_abstractions,
			      Abstraction * abstraction, int limit_abstract_states_merge, int min_limit_abstract_states_merge, 
			      int limit_transitions_merge) {
  assert(!done());
  vector<int> candidate_vars (remaining_vars);

  //Remove candidate vars
  if(limit_abstract_states_merge > 0){
      int limit = limit_abstract_states_merge/abstraction->size();
      int min_limit = min_limit_abstract_states_merge/abstraction->size();
      candidate_vars.erase(remove_if(begin(candidate_vars), 
				     end(candidate_vars),
				     [&](int var){
					 return (!all_abstractions[var] ||
						 all_abstractions[var]->size() > limit  ||
						 (limit_transitions_merge && 
						  abstraction->estimate_transitions(all_abstractions[var]) > limit_transitions_merge && 
						  all_abstractions[var]->size() > min_limit));
				     }), end(candidate_vars));
  }

  if(candidate_vars.empty()) 
      return -1;

  //Apply the criteria in order, until its finished or there is only one remaining variable  
  for(int i = 0; candidate_vars.size() > 1 && 
	i < criteria.size(); ++i){
        criteria[i]->filter(all_abstractions, candidate_vars, abstraction);
  }
  assert(!candidate_vars.empty());
    
  cout << "Candidates: ";
  for(int i = 0; i < candidate_vars.size(); i++){
    cout << candidate_vars[i] << " ";
  }
  cout << endl;

  int var = candidate_vars[0];
  select_next(var);
  return var;
}

void MergeLinearCriteria::remove_useless_vars (const std::vector<int> & useless_vars) {
    for (int v : useless_vars) cout << "Remove var from merge consideration: " << v << endl;
    remaining_vars.erase(std::remove_if(begin(remaining_vars), end(remaining_vars), 
					[&](int i) {
					    return std::find (begin(useless_vars), 
							      end(useless_vars), 
							      i) != std::end(useless_vars);
					}), end(remaining_vars));
} 
// bool MergeLinearCriteria::reduce_labels_before_merge () const{
//     for (int i = 0; i < criteria.size(); ++i){
// 	if(criteria[i0.]->reduce_labels_before_merge()){
// 	    return true;
// 	}
//     }
//     return false;
// }

static MergeStrategy *_parse(OptionParser &parser) {

    vector<string> variable_orders;
    variable_orders.push_back("level");
    variable_orders.push_back("reverse_level");
    variable_orders.push_back("random");
    parser.add_enum_option("var_order", variable_orders,

			   "merge variable order for tie breaking", 
			    "RANDOM");

    parser.add_list_option<MergeCriterion *> 
	("criteria", "list of criteria for the merge linear strategy");

    Options opts = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return new MergeLinearCriteria(opts);
}

static Plugin<MergeStrategy> _plugin("merge_linear_criteria", _parse);
