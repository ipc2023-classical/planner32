#ifndef MERGE_AND_SHRINK_SHRINK_OWN_LABELS_H
#define MERGE_AND_SHRINK_SHRINK_OWN_LABELS_H

#include "shrink_strategy.h"
#include <list>
#include <vector>

class Options;

class ShrinkOwnLabels : public ShrinkStrategy {
    const bool perform_sg_shrinking;
public:
    ShrinkOwnLabels(const Options &opts);
    virtual ~ShrinkOwnLabels();

    virtual std::string name() const;
    virtual void dump_strategy_specific_options() const;

    virtual bool reduce_labels_before_shrinking() const;

    virtual void shrink(Abstraction &abs,
			int target, bool force = false);
    virtual void shrink_atomic(Abstraction &abs);
    virtual void shrink_before_merge(Abstraction &abs1, 
				     Abstraction &abs2);

    static ShrinkOwnLabels *create_default();
};

#endif
