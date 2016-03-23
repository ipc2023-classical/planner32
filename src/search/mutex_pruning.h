#ifndef MUTEX_PRUNING_H
#define MUTEX_PRUNING_H

#include <memory>
#include "prune_heuristic.h"
#include "sym/sym_params.h"


class SymVariables;

class MutexPruning : public PruneHeuristic {
protected:
      const SymParamsMgr mgrParams; //Parameters for SymManager configuration.

      std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  
//std::unique_ptr<SymManager> mgr;    //The symbolic manager to handle mutex BDDs

      std::vector<BDD> mutex_bdds;

public:
    MutexPruning(const Options &options);
    virtual ~MutexPruning(){}
    virtual void initialize();

    virtual int compute_heuristic(const State &/*state*/){return 0;}
    //Methods for pruning explicit search
    virtual bool prune_generation(const State &/*state*/, int /*g*/){return false;}
    virtual bool prune_expansion (const State &/*state*/, int /*g*/){return false;}
    virtual bool is_dead_end(const State &state);

    virtual void print_statistics();
};

#endif
