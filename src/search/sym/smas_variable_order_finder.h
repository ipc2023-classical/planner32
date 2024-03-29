//WARNING: This is outdated and does not even use the right causal graph!!!
// Definetively not recommended to use it
#ifndef SMAS_VARIABLE_ORDER_FINDER_H
#define SMAS_VARIABLE_ORDER_FINDER_H

#include <vector>
enum VarOrder {
    MERGE_LINEAR_CG_GOAL_LEVEL,
    MERGE_LINEAR_CG_GOAL_RANDOM,
    MERGE_LINEAR_GOAL_CG_LEVEL,
    MERGE_LINEAR_RANDOM,
    MERGE_DFP,
    MERGE_LINEAR_LEVEL,
    MERGE_LINEAR_REVERSE_LEVEL,
    MAX_MERGE_STRATEGY
};

class SMASVariableOrderFinder {
    const VarOrder merge_strategy;
    std::vector<int> selected_vars;
    std::vector<int> remaining_vars;
    std::vector<bool> is_goal_variable;
    std::vector<bool> is_causal_predecessor;

    void select_next(int position, int var_no);
public:
    SMASVariableOrderFinder(VarOrder merge_strategy, bool is_first = true);
    ~SMASVariableOrderFinder();
    bool done() const;
    int next();
};

#endif
