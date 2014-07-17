#include "scc.h"

#include "equivalence_relation.h"

using namespace std;

void SCC::dfs_equivalence(const vector<vector<int> > &graph, 
	     vector<bool> &is_goal, 
	     int vertex, int & current_dfs_number,
	     vector<int> & stack,
	     vector<int> & stack_indices,
	     vector<int> & dfs_numbers,
	     vector<int> & dfs_minima,
             std::vector<__gnu_cxx::slist<int> >  & final_sccs) {
    int vertex_dfs_number = current_dfs_number++;
    dfs_numbers[vertex] = dfs_minima[vertex] = vertex_dfs_number;
    stack_indices[vertex] = stack.size();
    stack.push_back(vertex);

    const vector<int> &successors = graph[vertex];
    for (int i = 0; i < successors.size(); i++) {
	int succ = successors[i];
	int succ_dfs_number = dfs_numbers[succ];
	if (succ_dfs_number == -1) {
	    dfs_equivalence(graph, is_goal, succ, current_dfs_number, 
			    stack, stack_indices, dfs_numbers, dfs_minima, final_sccs);
	    dfs_minima[vertex] = min(dfs_minima[vertex], dfs_minima[succ]);
	} else if (succ_dfs_number < vertex_dfs_number && stack_indices[succ] != -1) {
	    dfs_minima[vertex] = min(dfs_minima[vertex], succ_dfs_number);
	}
	if(is_goal[succ]){
	    is_goal[vertex] = true;
	}
    }

    if (dfs_minima[vertex] == vertex_dfs_number) {
	int stack_index = stack_indices[vertex];
	final_sccs.push_back(__gnu_cxx::slist<int>());
	__gnu_cxx::slist<int> & scc = final_sccs.back();
	for (int i = stack_index; i < stack.size(); i++) {
	    scc.push_front(stack[i]);
	    stack_indices[stack[i]] = -1;
	}
	stack.erase(stack.begin() + stack_index, stack.end());
    }
}


void SCC::compute_scc_equivalence(const vector<vector<int> > &graph, 
				  vector<bool> &is_goal, 
				  std::vector<__gnu_cxx::slist<int> >  & result){
    int node_count = graph.size();

    // The following three are indexed by vertex number.
    vector<int> dfs_numbers(node_count, -1);
    vector<int> dfs_minima (node_count, -1);
    vector<int> stack_indices (node_count, -1);
    vector<int> stack; // This is indexed by the level of recursion.
    stack.reserve(node_count);
    int current_dfs_number = 0;

    for (int i = 0; i < node_count; i++)
	if (dfs_numbers[i] == -1)
	    dfs_equivalence(graph, is_goal, i, current_dfs_number, 
			    stack, stack_indices, dfs_numbers, dfs_minima, result);
    //reverse(result.begin(), result.end());
    //reverse(scc_goal.begin(), scc_goal.end());
}
