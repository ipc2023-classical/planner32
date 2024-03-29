#ifndef MERGE_AND_SHRINK_SHRINK_MERGE_STRATEGY_H
#define MERGE_AND_SHRINK_SHRINK_MERGE_STRATEGY_H

#include <string>
#include <vector>

class Abstraction;

class MergeStrategy {
protected:
    int total_merges; 
    int remaining_merges;


    virtual void dump_strategy_specific_options() const = 0;
public:
    MergeStrategy();
    virtual ~MergeStrategy() {}

    void dump_options() const;

    bool done() const {
        return remaining_merges == 0;
    }

    void init (const std::vector <Abstraction * > & abstractions) {
	total_merges = abstractions.size() - 1;
	remaining_merges = abstractions.size() - 1 ;
	init_strategy(abstractions);
    }

    virtual void init_strategy (const std::vector <Abstraction * > & ) = 0;

    virtual void remove_useless_vars (const std::vector<int> & ) {}



    // implementations of get_next should decrease remaining_merges by one
    // every time they return a pair of abstractions which are merged next.
    // Alvaro: Merge strategies have now a limit on the size of the
    // merge.  If specified (> 0), the pair returned should fit the
    // constraint: a1.size()*a2.size()<=limit
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions,
					 int limit_abstract_states_merge = 0, int min_limit_abstract_states_merge = 0, 
					 int limit_transitions_merge = 0) = 0;
    virtual std::string name() const = 0;

    //Alvaro: Added to know whether a merge strategy is linear or not
    // (more general than old way comparing the strategy name)
    virtual bool is_linear() const = 0;
};

#endif
