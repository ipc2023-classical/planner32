#include "sym_exploration.h"

using namespace std;

SymExploration::SymExploration(const SymParamsSearch & params):
    mgr(nullptr), p(params), fw(true) {}



void SymExploration::statistics() const {
    cout << "Exp " << (fw ? "fw" : "bw") << " time: " << stats.step_time << "s (img:" << 
	stats.image_time <<	"s, heur: " << stats.time_heuristic_evaluation << 
	"s) in " << stats.num_steps_succeeded  << " steps ";

    cout << "  ";
} 
