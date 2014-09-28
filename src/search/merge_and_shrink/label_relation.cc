#include "label_relation.h"

#include "labels.h"
#include "simulation_relation.h"

using namespace std;

LabelRelation::LabelRelation(Labels * _labels) : labels (_labels){}

void LabelRelation::dump() const {
  for (int l = 0; l < dominates_in.size(); ++l){
    if (l < 10){
      cout << "l" << l << ": ";  dump(l);
    }else{
      cout << "l" << l << ":";  dump(l);
    }

  }
}

void LabelRelation::dump(int label) const {
  if(dominated_by_noop_in[label] >= 0 && dominated_by_noop_in[label] <= 9){
    cout << " Dominated by noop: " << dominated_by_noop_in[label] << ", labels: ";
  }else{
    cout << " Dominated by noop:" << dominated_by_noop_in[label] << ", labels: ";
  } 

  for (int l2 = 0; l2 < dominates_in[label].size(); ++l2){
    if(dominates_in[l2][label] >= 0 && dominates_in[l2][label] <= 9) cout <<" ";
    cout << dominates_in[l2][label] << " ";
  }
  cout  << endl;
}

void LabelRelation::prune_operators(){
    //cout << "We have " << labels->get_size() << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
    for (int l = 0; l < dominates_in.size(); ++l){
	//labels->get_label_by_index(l)->dump();	
	if (dominated_by_noop_in[l]== DOMINATES_IN_ALL){
	    cout << g_operators[l].get_name() << " is dominated by noop " << endl;
	}

	for (int l2 = 0; l2 < dominates_in.size(); ++l2){
	    if (l2 != l && dominates_in[l2][l] == DOMINATES_IN_ALL){
		cout << g_operators[l].get_name() << " is dominated by " << g_operators[l2].get_name() << endl;
	    }
	}
    }
}

void LabelRelation::get_labels_dominated_in_all(std::vector<int> & labels_dominated_in_all){
    //cout << "We have " << labels->get_size() << " labels "<< dominates_in.size() << " " << dominates_in[0].size()  << endl;
    for (int l = 0; l < dominates_in.size(); ++l){
	//labels->get_label_by_index(l)->dump();	
	if (dominated_by_noop_in[l]== DOMINATES_IN_ALL){
	    labels_dominated_in_all.push_back(l);
	    continue;
	}

	for (int l2 = 0; l2 < dominates_in.size(); ++l2){
	    if (l2 != l && dominates_in[l2][l] == DOMINATES_IN_ALL && 
		!dominates_in[l][l2] == DOMINATES_IN_ALL){
		labels_dominated_in_all.push_back(l);
		break;
	    }
	}
    }
}
