#include "opt_order.h"

#include <ostream>
#include "../globals.h"

using namespace std;

void InfluenceGraph::get_ordering(vector <int> & ordering) const{
    long value_optimization_function = optimize_variable_ordering_gamer(ordering, 50000);
    cout << "Value: " <<  value_optimization_function << endl;

    for (int counter = 0; counter < 20; counter++) {
	vector <int> new_order; 
	randomize(ordering, new_order); //Copy the order randomly
	long new_value = optimize_variable_ordering_gamer(new_order, 50000); 

	if (new_value < value_optimization_function){
	    value_optimization_function = new_value;
	    ordering.swap(new_order);
	    cout << "New value: " <<  value_optimization_function << endl;
	}
    }
}


void InfluenceGraph::randomize(vector <int> & ordering, vector<int> & new_order) {
    for(int i = 0; i < ordering.size(); i++){
	int rnd_pos = g_rng.next(ordering.size() - i);
	int pos = -1;
	do{       
	    pos++;
	    bool found;
	    do {
		found = false;
		for(int j = 0; j < new_order.size(); j++){
		    if(new_order[j] == ordering[pos]){
			found = true;
			break;
		    }
		} 
		if(found) pos++;
	    }while(found);
	}while(rnd_pos-- > 0);
	new_order.push_back(ordering[pos]);
    }
}



long InfluenceGraph::optimize_variable_ordering_gamer(vector <int> & order,
						      int iterations) const{
  long totalDistance = compute_function(order);
   
  long oldTotalDistance = totalDistance;
  //Repeat iterations times
  for (int counter = 0; counter < iterations; counter++) {
    //Swap variable
    int swapIndex1 = g_rng.next(order.size());
    int swapIndex2 = g_rng.next(order.size());
    if(swapIndex1 == swapIndex2) continue;

    //Compute the new value of the optimization function
    for(int i = 0; i < order.size(); i++){
      if(i == swapIndex1 || i == swapIndex2) continue;

      if(influence(order[i], order[swapIndex1]))
	totalDistance += (- (i - swapIndex1)*(i - swapIndex1)
			  + (i - swapIndex2)*(i - swapIndex2));
	
      if(influence(order[i], order[swapIndex2]))
	totalDistance += (- (i - swapIndex2)*(i - swapIndex2)
			  + (i - swapIndex1)*(i - swapIndex1));	
    }

    //Apply the swap if it is worthy
    if (totalDistance < oldTotalDistance){
      int tmp = order[swapIndex1];
      order[swapIndex1] = order[swapIndex2];
      order[swapIndex2] = tmp;
      oldTotalDistance = totalDistance;
      
      /*if(totalDistance != compute_function(order)){
	cerr << "Error computing total distance: " << totalDistance << " " << compute_function(order) << endl;
	exit(-1);
      }else{
	cout << "Bien: " << totalDistance << endl;
      }*/
    }else{
      totalDistance = oldTotalDistance;
    }
  }
//  cout << "Total distance: " << totalDistance << endl;
  return totalDistance;
}



long InfluenceGraph::compute_function(const std::vector <int> & order) const {
    long totalDistance = 0;
    for (int i = 0; i < order.size() - 1; i++) {
	for (int j = i+1; j < order.size(); j++) {
	    if(influence(order[i], order[j])){
		totalDistance += (i-j)*(i -j);
	    }
	}
    }
    return totalDistance;
}


InfluenceGraph::InfluenceGraph(int num)  {
    influence_graph.resize(num);
    for (auto & i : influence_graph){
	i.resize(num, 0);
    }
}



void InfluenceGraph::optimize_variable_ordering_gamer(vector <int> & order,
						      vector <int> & partition_begin,
						      vector <int> & partition_sizes,
						      int iterations) const{
  long totalDistance = compute_function(order);
   
  long oldTotalDistance = totalDistance;
  //Repeat iterations times
  for (int counter = 0; counter < iterations; counter++) {
    //Swap variable
    int partition = g_rng.next(partition_begin.size()); 
    if(partition_sizes[partition] <= 1) continue;
    int swapIndex1 = partition_begin[partition] + g_rng.next(partition_sizes[partition]);
    int swapIndex2 = partition_begin[partition] + g_rng.next(partition_sizes[partition]);
    if(swapIndex1 == swapIndex2) continue;

    //Compute the new value of the optimization function
    for(int i = 0; i < order.size(); i++){
      if(i == swapIndex1 || i == swapIndex2) continue;

      if(influence(order[i], order[swapIndex1]))
	totalDistance += (- (i - swapIndex1)*(i - swapIndex1)
			  + (i - swapIndex2)*(i - swapIndex2));
	
      if(influence(order[i], order[swapIndex2]))
	totalDistance += (- (i - swapIndex2)*(i - swapIndex2)
			  + (i - swapIndex1)*(i - swapIndex1));	
    }

    //Apply the swap if it is worthy
    if (totalDistance < oldTotalDistance){
      int tmp = order[swapIndex1];
      order[swapIndex1] = order[swapIndex2];
      order[swapIndex2] = tmp;
      oldTotalDistance = totalDistance;
      
      /*if(totalDistance != compute_function(order)){
	cerr << "Error computing total distance: " << totalDistance << " " << compute_function(order) << endl;
	exit(-1);
      }else{
	cout << "Bien: " << totalDistance << endl;
      }*/
    }else{
      totalDistance = oldTotalDistance;
    }
  }
//  cout << "Total distance: " << totalDistance << endl;

}