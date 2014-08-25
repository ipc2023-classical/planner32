#ifndef MERGE_AND_SHRINK_MERGE_CRITERION_H
#define MERGE_AND_SHRINK_MERGE_CRITERION_H

// The new merge strategy is based on a list of criteria.  We start
// with the set of candidate variables and apply each criterion, that
// discards some variables, until only one variable is left. If more
// than one variable is left after applying all the criteria, the
// merge_order is used as final tie-breaking.
// 
// CG: prefer v to v' if v has causal influence over already
// merged variables and v' not.
//
// GOAL: prefer goal variables to non-goal variables
//
// Relevant: prefer relevant variables to non-relevant variables. A
// variable is relevant if a) it is a goal or b) it has a causal
// influence over already merged variables.
//
// MinSCC: Same as CG, but applying tie-breaking in case that more
// than one variable is causally relevant: prefer variable v to v' if
// SCC(v) has a path to SCC(v') if SCC(v) has a path to SCC(v').
// Optionally, only one variable is selected per SCC 
// (the one with smallest "level", i.e. closer to the SCC root)
//
// Empty: prefer variables that will make more labels of the current
// abstractionempty
// 
// Num: prefer variables that apper in more labels of the current
// abstraction.

class OptionParser;
class Options;

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include "../scc.h"

class Abstraction;

class MergeCriterion {
 protected:
    //Returns true in case that at least one variable fits the criterion
    bool filter(std::vector <int> & vars, 
		const std::vector<bool> & criterion) const {
	std::vector<int> aux; 
	for(int i =0; i < vars.size(); i++){
	    if(criterion[vars[i]]){
		aux.push_back(vars[i]);
	    }
	}
	if(!aux.empty()){
	    vars.swap(aux);
	    return true;
	}
	return false;
    }

    template<class T>
	void filter_best(std::vector <int> & vars, 
			 const std::vector<T> & criterion, 
			 bool minimize, double opt_margin=1, double opt_diff=0)  const {
	std::vector<int> aux; 

	if(opt_diff == 0 && opt_margin == 1){
	    T best = (minimize ? std::numeric_limits<T>::max() : 0.0); 
	    for(int i =0; i < vars.size(); i++){
		T score = criterion[vars[i]];
		if((minimize && score < best) || 
		   (!minimize && score > best)){
		    std::vector<int>().swap(aux);
		    best = score;
		}
		if(score == best){
		    aux.push_back(vars[i]);
		}
	    }
	}else{
	    T best = (minimize ? *(std::min_element(criterion.begin(), criterion.end())) :
		                  *(std::max_element(criterion.begin(), criterion.end())));
	    T bestcmp = (minimize ? std::max(best*opt_margin, best + opt_diff) :
			 std::min(best*opt_margin, best - opt_diff));
            std::cout << " (" << best << " " << bestcmp << ") ";
	    for(int i =0; i < vars.size(); i++){
		T score = criterion[vars[i]];
		if((minimize && score <= bestcmp) || 
		   (!minimize && score >= bestcmp)){
		    aux.push_back(vars[i]);
		}		    
	    }
	}
	if(!aux.empty()){
	    vars.swap(aux);
	}
    }

    bool allow_incremental;
 public:
 MergeCriterion() : allow_incremental(true){}
    virtual void init() = 0; 
    virtual void disable_incremental(){
	allow_incremental = false;
    }
    // Allows for incremental computation (currently used to compute
    // predecessor variables in the CG). However, it should only work if
    // we have not disabled this.
    virtual void select_next(int var_no) = 0; 
    virtual void filter(const std::vector<Abstraction *> &all_abstractions, 
			std::vector <int> & vars, Abstraction * abstraction)  = 0;
    virtual std::string get_name() const = 0;
    virtual bool reduce_labels_before_merge () const{
	return false;
    }
};

class MergeCriterionCG : public MergeCriterion {
 protected:
    std::vector<bool> is_causal_predecessor;
 public:
    virtual void init();
    virtual void select_next(int var_no);
    virtual void filter(const std::vector<Abstraction *> &all_abstractions, 
			std::vector <int> & vars, Abstraction * abstraction) ;
    virtual std::string get_name() const {
	return "CG";
    }
};

class MergeCriterionGoal : public MergeCriterion {
    std::vector<bool> is_goal_variable;
 public:
    virtual void init();
    virtual void select_next(int var_no);
    virtual void filter(const std::vector<Abstraction *> &all_abstractions, 
			std::vector <int> & vars, Abstraction * abstraction) ;
    virtual std::string get_name() const {
	return "GOAL";
    }
};

class MergeCriterionRelevant : public MergeCriterionCG {
    virtual void init();
    virtual std::string get_name() const {
	return "RELEVANT";
    }
};

class MergeCriterionMinSCC : public MergeCriterion {
    const bool reverse;
    const bool tie_by_level;
    // whether to use pre_to_eff or complete cg
    const bool complete_cg;  

    std::vector<bool> is_causal_predecessor;
    SCC * scc;

    void forbid_scc_descendants(int scc_index,
				const std::vector<std::set<int> > & scc_graph, 
				std::vector<bool> & forbidden_sccs) const;
 public:
 MergeCriterionMinSCC(bool reverse_ = false, 
		      bool tie_by_level_ = true, 
		      bool complete_cg_ = false) : 
    reverse (reverse_),
	tie_by_level(tie_by_level_), 
	complete_cg(complete_cg_),
	scc(NULL){
	}
    MergeCriterionMinSCC(const Options & opts);
    virtual ~MergeCriterionMinSCC(){
	delete scc; //TODO: Smart pointer
    }
    virtual void init();
    virtual void select_next(int var_no);
    virtual void filter(const std::vector<Abstraction *> &all_abstractions, 
			std::vector <int> & vars, Abstraction * abstraction) ;
    virtual std::string get_name() const {
	return "SCC";
    }
};

class MergeCriterionTRs : public MergeCriterion {
    const bool only_goals;
    const bool only_empty;

    // By default, this criterion discards all the variables that do
    // have the maximum number of shared TRs.  In practice, however,
    // if the maximum is 1000 and there is a variable with 9999, it
    // may make sense to let that as well.  We add two parameters to
    // control the amount of suboptimality allowed (in terms of a
    // factor or a difference)
    // allowed_value = min(best*opt_factor, best - opt_diff)
    const double opt_factor;
    const int opt_diff;
 public:
    MergeCriterionTRs(const Options & opts);
    virtual void init();
    virtual void select_next(int var_no);
    virtual void filter(const std::vector<Abstraction *> &all_abstractions, 
			std::vector <int> & vars, Abstraction * abstraction) ;
    virtual std::string get_name() const {
	std::stringstream ss;
	ss << "TRs(";
	if(only_goals) ss << "goals ";
	if(only_empty) ss << "empty";
        ss << ")";
        return ss.str();
    }
    virtual bool reduce_labels_before_merge () const{
	return true;
    }
};

#endif
