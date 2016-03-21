#ifndef SYM_HEURISTIC_GENERATOR_H
#define SYM_HEURISTIC_GENERATOR_H

class SymHeuristic;

class SymHeuristicGenerator {
 public:
    virtual void getSymHeuristic(SymVariables * svars, 
				 std::vector<std::shared_ptr<SymHeuristic> > & heurs)  = 0;
};

#endif
