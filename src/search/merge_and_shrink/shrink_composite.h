#ifndef MERGE_AND_SHRINK_SHRINK_COMPOSITE_H
#define MERGE_AND_SHRINK_SHRINK_COMPOSITE_H

#include "shrink_strategy.h"

class Options;
class Signature;

class ShrinkComposite : public ShrinkStrategy {
    const std::vector<ShrinkStrategy *> strategies;

public:
    ShrinkComposite(const Options &opts);
    virtual ~ShrinkComposite();

    virtual std::string name() const;
    virtual void dump_strategy_specific_options() const;

    virtual bool reduce_labels_before_shrinking() const;

    virtual void shrink(Abstraction &abs, int target, bool force = false);
    virtual void shrink_atomic(Abstraction &abs);
    virtual void shrink_before_merge(Abstraction &abs1, Abstraction &abs2);
};

#endif
