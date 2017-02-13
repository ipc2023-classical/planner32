#ifndef MERGE_AND_SHRINK_MERGE_LINEAR_CRITERIA_H
#define MERGE_AND_SHRINK_MERGE_LINEAR_CRITERIA_H

#include "merge_strategy.h"
#include <vector>

class Options;
class MergeCriterion;

enum MergeOrder {
  LEVEL,
  REVERSE_LEVEL, 
  RANDOM
};

class MergeLinearCriteria : public MergeStrategy {
  const std::vector <MergeCriterion *> criteria;
  const MergeOrder order;

  //  std::vector<int> selected_vars;
  std::vector<int> remaining_vars;

  void select_next(int var_no);

  //Selects an abstraction based on the criterions
  virtual int next(const std::vector<Abstraction *> &all_abstractions, 
		   Abstraction * abstraction = 0, int limit_abstract_states_merge = 0, int min_limit_abstract_states_merge = 0, 
		   int limit_transitions_merge = 0);
protected:
    virtual void dump_strategy_specific_options() const;

    virtual void init_strategy (const std::vector <Abstraction * > & abstractions);
 public:
  explicit MergeLinearCriteria(const Options &opts);
  virtual ~MergeLinearCriteria();

  virtual void remove_useless_vars (const std::vector<int> & useless_vars);

  virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions, 
				       int limit_abstract_states_merge, int min_limit_abstract_states_merge, 
				       int limit_transitions_merge = 0);
  virtual std::string name() const;
  virtual bool is_linear() const;  
  // virtual bool reduce_labels_before_merge () const;
};

#endif
