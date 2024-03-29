#ifndef PRUNE_HEURISTIC_H
#define PRUNE_HEURISTIC_H

#include "heuristic.h"
#include "operator_cost.h"

class State;
class OptionParser;
class Options;
class BDD;
class SymTransition;
class SymManager;
class SearchProgress; 

class PruneHeuristic : public Heuristic {
protected:
    virtual int compute_heuristic(const State &state) = 0;

public:
    PruneHeuristic(const Options &options);
    virtual ~PruneHeuristic();
    virtual void initialize(bool force_initialization = false) = 0;

    //Methods for pruning explicit search
    virtual void prune_applicable_operators(const State &/*state*/, int /*g*/, std::vector<const Operator *> & /*operators*/, SearchProgress & /*search_progress*/) {}
    virtual bool prune_generation(const State &state, int g, const State &parent, int action_cost) = 0;
    virtual bool prune_expansion (const State &state, int g) = 0;
    virtual bool is_dead_end(const State &state) = 0;

    virtual void print_statistics() {}
    //Methods for pruning symbolic search. Return the BDD without
    //pruned states.
    /* virtual BDD prune_generation(const BDD &bdd, int g) = 0; */
    /* virtual BDD prune_expansion (const BDD &bdd, int g) = 0; */

    //Returns a TR that can be used to generate dominated/dominating
    //state sets
    //virtual SymTransition * getTR(SymManager * mgr) = 0;

    static void add_options_to_parser(OptionParser &parser);
    static Options default_options();
};

#endif
