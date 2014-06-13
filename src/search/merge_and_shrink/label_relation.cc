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


