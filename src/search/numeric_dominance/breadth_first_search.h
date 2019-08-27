#ifndef NUMERIC_DOMINANCE_BREADTH_FIRST_SEARCH_H
#define NUMERIC_DOMINANCE_BREADTH_FIRST_SEARCH_H

#include <vector>
#include <limits>

template <typename T> 
void breadth_first_search_reachability_distances_one(const std::vector<std::vector<int> > &graph,
						     int initial_state,
						     std::vector<T> & distances,
						     std::vector<int> & reachable_states) {
    int num_states = graph.size();
    distances.resize(num_states);
    std::fill(distances.begin(), distances.end(), std::numeric_limits<int>::max());

    int increase_distance_in = 0;
    int current_distance = 0;

    distances[initial_state] = 0;
    reachable_states.push_back(initial_state);
    for(size_t i = 0; i < reachable_states.size(); ++i) {
	if (increase_distance_in == i) {
	    current_distance ++;
	    increase_distance_in = reachable_states.size();
	}
	for (int successor : graph[reachable_states[i]]) {
	    if (distances[successor] > current_distance) {
		distances[successor] = current_distance;
		reachable_states.push_back(successor);
	    }
	}
    }
}



#endif






