#ifndef MERGE_AND_SHRINK_SHRINK_MERGE_STRATEGY_H
#define MERGE_AND_SHRINK_SHRINK_MERGE_STRATEGY_H

#include <string>
#include <vector>

class Abstraction;

class MergeStrategy {
protected:
    int remaining_merges;
    virtual void dump_strategy_specific_options() const = 0;
public:
    MergeStrategy();
    virtual ~MergeStrategy() {}

    void dump_options() const;

    bool done() const {
        return remaining_merges == 0;
    }
    // implementations of get_next should decrease remaining_merges by one
    // every time they return a pair of abstractions which are merged next.
    // Alvaro: Merge strategies have now a limit on the size of the
    // merge.  If specified (> 0), the pair returned should fit the
    // constraint: a1.size()*a2.size()<=limit
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions, int limit_abstract_states_merge = 0) = 0;
    virtual std::string name() const = 0;

    //Alvaro: Added to know whether a merge strategy is linear or not
    // (more general than old way comparing the strategy name)
    virtual bool is_linear() const = 0;
};

#endif
