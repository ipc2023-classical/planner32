#ifndef MERGE_AND_SHRINK_SHRINK_BISIMULATION_H
#define MERGE_AND_SHRINK_SHRINK_BISIMULATION_H

#include <limits>
#include "shrink_strategy.h"

class Options;
class Signature;

class ShrinkBisimulation : public ShrinkStrategy {
    enum AtLimit {
        RETURN,
        USE_UP
    };

    /*
      threshold: Shrink the abstraction iff it is larger than this
      size. Note that this is set independently from max_states, which
      is the number of states to which the abstraction is shrunk.
    */

    const bool greedy;
    const int threshold;
    const bool group_by_h;
    const AtLimit at_limit;

    // Goal states in abstractions with all goal variables will always
    // remain goal -> we can ignore all outgoing transitions, as we
    // can never leave this state; should help aggregating more
    // states, as all goal states should become bisimilar (can also
    // increase the abstraction size by making more abstract states
    // reachable)
    const bool aggregate_goals;

    void compute_abstraction(
        Abstraction &abs,
        int target_size,
        EquivalenceRelation &equivalence_relation);

    int initialize_groups(const Abstraction &abs,
                          std::vector<int> &state_to_group);
    void compute_signatures(
        const Abstraction &abs,
        std::vector<Signature> &signatures,
        std::vector<int> &state_to_group);
public:
    ShrinkBisimulation(const Options &opts);
    virtual ~ShrinkBisimulation();

    virtual std::string name() const;
    virtual void dump_strategy_specific_options() const;

    virtual bool reduce_labels_before_shrinking() const;

    virtual void shrink(Abstraction &abs, int target, bool force = false);
    virtual void shrink_atomic(Abstraction &abs);
    virtual void shrink_before_merge(Abstraction &abs1, Abstraction &abs2);

    static ShrinkStrategy *create_default(bool aggregate_goals = false, int limit_states = std::numeric_limits<int>::max());
};

#endif
