#ifndef MERGE_AND_SHRINK_ABSTRACTION_H
#define MERGE_AND_SHRINK_ABSTRACTION_H

#include "shrink_strategy.h"
#include "../utilities.h"
#include "../sym/sym_variables.h"
//#include <boost/any.hpp>
//#include <boost/dynamic_bitset.hpp>
#include <boost/dynamic_bitset.hpp>

#include <ext/slist>
#include <string>
#include <vector>

class EquivalenceRelation;
class Label;
class Labels;
class LabelMap;
class State;
class SimulationRelation;
class LabelledTransitionSystem;
class LTSEfficient;

typedef int AbstractStateRef;

struct AbstractTransition {
    AbstractStateRef src;
    AbstractStateRef target;
    boost::dynamic_bitset<> based_on_operators;

    AbstractTransition(AbstractStateRef src_, AbstractStateRef target_)
        : src(src_), target(target_) {
    }

    bool operator==(const AbstractTransition &other) const {
        return src == other.src && target == other.target;
    }

    bool operator!=(const AbstractTransition &other) const {
        return !(*this == other);
    }

    bool operator<(const AbstractTransition &other) const {
        return src < other.src || (src == other.src && target < other.target);
    }

    bool operator>=(const AbstractTransition &other) const {
        return !(*this < other);
    }
};

class Abstraction {
    friend class AtomicAbstraction;
    friend class CompositeAbstraction;
    friend class PDBAbstraction;

    friend class ShrinkStrategy; // for apply() -- TODO: refactor!
    friend class SimulationRelation; // for apply() -- TODO: refactor!
    friend class LDSimulation; // for setting store_original_operators -- TODO: refactor!

    static const int PRUNED_STATE;
    static const int DISTANCE_UNKNOWN;

    static bool store_original_operators;

    // There should only be one instance of Labels at runtime. It is created
    // and managed by MergeAndShrinkHeuristic. All abstraction instances have
    // a copy of this object to ease access to the set of labels.
    // Alvaro: Removed const to allow setting the abstraction as
    // irrelevant for some labels during normalize().
    Labels *labels;
    /* num_labels equals to the number of labels that this abstraction is
       "aware of", i.e. that have
       been incorporated into transitions_by_label. Whenever new labels are
       generated through label reduction, we do *not* update all abstractions
       immediately. This equals labels->size() after normalizing. */
    int num_labels;
    /* transitions_by_label and relevant_labels both have size of (2 * n) - 1
       if n is the number of operators, because when applying label reduction,
       at most n - 1 fresh labels can be generated in addition to the n
       original labels. */
    std::vector<std::vector<AbstractTransition> > transitions_by_label;
    std::vector<bool> relevant_labels;

    //TODO: Unify with transitions by label??
    std::unique_ptr<LabelledTransitionSystem> lts;
    std::unique_ptr<LTSEfficient> lts_efficient;

    // Alvaro: Information regarding the number of transitions by label
    // (needed by some merge_criterions, only computed on demand)
    // TODO: added as attribute in abstractions to avoid recomputation
    // when different criterions use this data. Move somewhere else?
    std::vector<int> num_transitions_by_label;
    std::vector<int> num_goal_transitions_by_label;

    int num_states;

    std::vector<int> init_distances;
    std::vector<int> goal_distances;
    std::vector<bool> goal_states;
    AbstractStateRef init_state;

    int max_f;
    int max_g;
    int max_h;

    bool transitions_sorted_unique;
    //Alvaro: substituted goal_relevant by number of goal variables =>
    //allow easy check of whether all the goal variables are relevant
    int goal_relevant_vars; 
    bool all_goals_relevant;

    mutable int peak_memory;

    void clear_distances();
    void compute_init_distances_unit_cost();
    void compute_goal_distances_unit_cost();
    void compute_init_distances_general_cost();
    void compute_goal_distances_general_cost();

    //Alvaro: Computes num_transitions_by_label and
    //num_goal_transitions_by_label
    void count_transitions_by_label();

    // are_transitions_sorted_unique() is used to determine whether the
    // transitions of an abstraction are sorted uniquely or not after
    // construction (composite abstraction) and shrinking (apply_abstraction).
    bool are_transitions_sorted_unique() const;

    void apply_abstraction(std::vector<__gnu_cxx::slist<AbstractStateRef> > &collapsed_groups);

    int unique_unlabeled_transitions() const;
protected:
    std::vector<int> varset;

    virtual void apply_abstraction_to_lookup_table(const std::vector<
                                                       AbstractStateRef> &abstraction_mapping) = 0;
    virtual int memory_estimate() const;
public:
    Abstraction(Labels *labels);
    virtual ~Abstraction();

    int total_transitions() const;
    int total_transition_operators() const;

    // Two methods to identify the abstraction in output.
    // tag is a convience method that upper-cases the first letter of
    // description and appends ": ";
    virtual std::string description() const = 0;
    std::string tag() const;

    //Returns a description of the abstract state 
    virtual std::string description(int s) const{
      return std::to_string(s);
    }

    static void build_atomic_abstractions(std::vector<Abstraction *> &result,
                                          Labels *labels);
    bool is_solvable() const;

    int get_cost(const State &state) const;
    int size() const;
    void statistics(bool include_expensive_statistics) const;
    const std::string & label_name (int l) const;

    int get_peak_memory_estimate() const;
    // NOTE: This will only return something useful if the memory estimates
    //       have been computed along the way by calls to statistics().
    // TODO: Find a better way of doing this that doesn't require
    //       a mutable attribute?

    bool are_distances_computed() const;
    void compute_distances();
    bool is_normalized() const;
    void normalize();
    void normalize2(); // Version of normalize handling the storing of original operators in the transitions
    EquivalenceRelation *compute_local_equivalence_relation() const;
    void release_memory();

    void dump_relevant_labels() const;
    void dump() const;

    // The following methods exist for the benefit of shrink strategies.
    int get_max_f() const;
    int get_max_g() const; // Not being used!
    int get_max_h() const;

    bool is_goal_state(int state) const {
        return goal_states[state];
    }

    int get_init_distance(int state) const {
        return init_distances[state];
    }

    int get_goal_distance(int state) const {
        return goal_distances[state];
    }

    // These methods should be private but is public for shrink_bisimulation
    int get_label_cost_by_index(int label_no) const;
    const std::vector<AbstractTransition> &get_transitions_for_label(int label_no) const;
    // This method is shrink_bisimulation-exclusive
    int get_num_labels() const;
    int get_num_nonreduced_labels() const;
    // These methods are used by non_linear_merge_strategy
    void compute_label_ranks(std::vector<int> &label_ranks);
    bool is_goal_relevant() const {
        return goal_relevant_vars > 0;
    }
    bool get_all_goal_vars_in() const{
	return all_goals_relevant;
    }
    // This is used by the "old label reduction" method
    const std::vector<int> &get_varset() const {
        return varset;
    }

    const std::vector <bool> & get_relevant_labels() const {
      return relevant_labels;
    }

    virtual AbstractStateRef get_abstract_state(const State &state) const = 0;
    virtual void getAbsStateBDDs(SymVariables * vars, std::vector<BDD> & abs_bdds) const = 0;

    LabelledTransitionSystem * get_lts(const LabelMap & labelMap);
    LTSEfficient * get_lts_efficient(const LabelMap & labelMap);

    //Alvaro: used by shrink_empty_labels
    const std::vector<bool> & get_goal_states() const {
	return goal_states;
    }

    //Alvaro: used by merge criterions. For each remaining
    //abstraction, counts the number of transitions in this
    //transitions labelled with a relevant label for that
    //abstraction. If only_empty is activated, only the transitions
    //that are relevant for 2 abstractions are counted (this
    //abstraction and the other). If only_goal is activated, only
    //transitions leading to a goal state are counted.
    void count_transitions(const std::vector<Abstraction *> &all_abstractions, 
			   const std::vector<int> & remaining, bool only_empty,
			   bool only_goal, std::vector<int> & result);
    
    bool is_own_label(int label_no);


    // Methods to prune irrelevant transitions Prune all the
    // transitions of a given label (it is completely dominated by
    // some other label or by noop)
    int prune_transitions_dominated_label_all(int label_no);

    // Prune all the transitions of label_no such that exist a better
    // transition for label_no_by
    int prune_transitions_dominated_label(int label_no, int label_no_by,
					  SimulationRelation & rel);
    // Prune all the transitions of label_no such that exist a better
    // transition for label_no_by
    int prune_transitions_dominated_label_equiv(int label_no, int label_no_by,
                      SimulationRelation & rel);
    
    //Prune all the transitions dominated by noop 
    int prune_transitions_dominated_label_noop(int label_no,
					       SimulationRelation & rel);

    int estimate_transitions(const Abstraction * other) const;

    bool check_dead_labels(std::vector<bool> & dead_labels, std::vector<bool> & dead_operators) const;
};

class AtomicAbstraction : public Abstraction {
    int variable;
    std::vector<AbstractStateRef> lookup_table;
protected:
    virtual std::string description() const;
    virtual std::string description(int s) const;
    virtual void apply_abstraction_to_lookup_table(
        const std::vector<AbstractStateRef> &abstraction_mapping);
    virtual int memory_estimate() const;
public:
    AtomicAbstraction(Labels *labels, int variable);
    virtual ~AtomicAbstraction();

    virtual AbstractStateRef get_abstract_state(const State &state) const;
    virtual void getAbsStateBDDs(SymVariables * vars, std::vector<BDD> & abs_bdds) const;
};

class CompositeAbstraction : public Abstraction {
    Abstraction *components[2];
    std::vector<std::vector<AbstractStateRef> > lookup_table;
protected:
    virtual std::string description() const;
    virtual std::string description(int s) const;
    virtual void apply_abstraction_to_lookup_table(
        const std::vector<AbstractStateRef> &abstraction_mapping);
    virtual int memory_estimate() const;
public:
    CompositeAbstraction(Labels *labels, Abstraction *abs1, Abstraction *abs2);
    virtual ~CompositeAbstraction();

    virtual AbstractStateRef get_abstract_state(const State &state) const;
    virtual void getAbsStateBDDs(SymVariables * vars, std::vector<BDD> & abs_bdds) const;
};

class PDBAbstraction : public Abstraction {
    // List of variables of the abstraction (ordered to perform ranking)
    std::vector<int> pattern; 
    std::vector<AbstractStateRef> lookup_table;
 private: 

    template <typename T>
	int rank(const T &state) const {
	int res = 0;
	for(int i = 0; i < pattern.size(); ++i){
	    int v = pattern[i];
	    if (res) res *= g_variable_domain[pattern[i]];
	    res += state[v];
	}
	return res;
    }

    BDD unrankBDD (SymVariables * vars, int id) const ;

    void insert_transitions(std::vector <int> & pre_vals, std::vector <int> & eff_vals,
			    int label_no, int pos);

    void insert_goals(std::vector <int> & goal_vals, int pos);


protected:
    virtual std::string description() const;
    virtual std::string description(int s) const;
    virtual void apply_abstraction_to_lookup_table(
        const std::vector<AbstractStateRef> &abstraction_mapping);
    virtual int memory_estimate() const;
public:
    PDBAbstraction(Labels *labels, const std::vector<int> & pattern_);
    virtual ~PDBAbstraction();

    virtual AbstractStateRef get_abstract_state(const State &state) const;
    virtual void getAbsStateBDDs(SymVariables * vars, std::vector<BDD> & abs_bdds) const;
};

#endif
