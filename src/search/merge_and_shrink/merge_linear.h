#ifndef MERGE_AND_SHRINK_SHRINK_MERGE_LINEAR_H
#define MERGE_AND_SHRINK_SHRINK_MERGE_LINEAR_H

#include "merge_strategy.h"

#include "variable_order_finder.h"

class Options;

class MergeLinear : public MergeStrategy {
    VariableOrderFinder order;
    bool need_first_index;
protected:
    virtual void dump_strategy_specific_options() const;
public:
    explicit MergeLinear(const Options &opts);
    virtual ~MergeLinear() {}

    // Alvaro: Merge strategies have now a limit on the size of the
    // merge.  If specified (> 0), the pair returned should fit the
    // constraint: a1.size()*a2.size()<=limit
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions, 
					 int limit_abstract_states_merge = 0);
    virtual std::string name() const;
};

#endif
