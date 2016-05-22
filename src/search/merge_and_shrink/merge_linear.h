#ifndef MERGE_AND_SHRINK_SHRINK_MERGE_LINEAR_H
#define MERGE_AND_SHRINK_SHRINK_MERGE_LINEAR_H

#include "merge_strategy.h"

#include "variable_order_finder.h"

class Options;
//Alvaro: Merge linear will behave as a non-linear merge in case that
// limit_abstract_states_merge is set
class MergeLinear : public MergeStrategy {
    VariableOrderFinder order;
    bool need_first_index;
protected:
    virtual void dump_strategy_specific_options() const;
public:
    explicit MergeLinear(const Options &opts);
    virtual ~MergeLinear() {}

    virtual void init_strategy (const std::vector <Abstraction * > & ) {}


    // Alvaro: Merge strategies have now a limit on the size of the
    // merge.  If specified (> 0), the pair returned should fit the
    // constraint: a1.size()*a2.size()<=limit
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions, 
					 int limit_abstract_states_merge = 0, 
					 int limit_transitions_merge = 0);
    virtual std::string name() const;
    virtual bool is_linear() const;
};

#endif
