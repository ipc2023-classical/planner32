#ifndef MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H
#define MERGE_AND_SHRINK_SIMULATION_HEURISTIC_H

#include "../heuristic.h"

class Abstraction;
class Labels;
class SimulationRelation;

class SimulationHeuristic : public Heuristic {
    Labels *labels;

    std::vector<Abstraction *> abstractions;
    std::vector<SimulationRelation *> simulations;    

    void dump_options() const;
protected:
    virtual void initialize();
    virtual int compute_heuristic(const State &state);
public:
    SimulationHeuristic(const Options &opts);
    ~SimulationHeuristic();
};

#endif
