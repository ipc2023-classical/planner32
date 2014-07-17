#ifndef SCC_H
#define SCC_H

#include <ext/slist>

/*
  This implements Tarjan's linear-time algorithm for finding the maximal
  strongly connected components. It takes time proportional to the sum
  of the number of vertices and arcs.

  Call SCC with a graph represented as a vector of vectors, where
  graph[i] is the vector of successors of vertex i.

  Method compute_scc_equivalence() returns a vector of strongly
  connected components, each of which is a list of vertices (ints).
  This is a partitioning of all vertices where each SCC is a maximal
  subset such that each node in an SCC is reachable from all other
  nodes in the SCC.  Note that the derived graph where each SCC is a
  single "supernode" is necessarily acyclic. The SCCs returned by
  get_result() are in a topological sort order with regard to this
  derived DAG.
*/

#include <vector>
#include <algorithm>
class EquivalenceRelation;

//TODO: It could be convinient to use template functions here to allow
//returning the SCC result in different ways (vector<vector<int> >,
//EquivalenceRelation, or even vector <int> state_to_group. However,
//it is not so easy to do so since we need to change the whole
//insertion of sccs into the result
class SCC {
    static void dfs_equivalence(
    const std::vector<std::vector<int> > &graph, 
	std::vector<bool> &is_goal, 
	int vertex, int & current_dfs_number,
	std::vector<int> & stack,
	std::vector<int> & stack_indices,
	std::vector<int> & dfs_numbers,
	std::vector<int> & dfs_minima,
	std::vector<__gnu_cxx::slist<int> > & final_sccs);
public:
static void compute_scc_equivalence(const std::vector<std::vector<int> > &graph, 
				    std::vector<bool> &is_goal, 
					std::vector<__gnu_cxx::slist<int> > & result);
};

#endif
