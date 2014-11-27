#include "variable_partition_finder.h"

#include "../globals.h"
#include "../operator.h"

using namespace std;

void VariablePartitionFinder:: dump()  const{
    cout << "Partition: " << endl;
    for(auto & p  : partitions){
        int size = 1;
        for(int v : p){
            cout << " " << v << "(" << g_fact_names[v][0] << ")";
            size *= g_variable_domain[v];
        }
        cout << " (" << size << ")" << endl;
    }
    cout << endl;
}
void VariablePartitionGreedy::init(){
    cout << "Init" << endl;
    partitions.resize(g_variable_domain.size());
    for (int i = 0; i < g_variable_domain.size(); i++){
        partitions[i].push_back(i);
    }
    cout << "Copy g_var domain" << endl;

    part_size = g_variable_domain;

    //Set weights according to the causal graph
    for(int op = 0; op < g_operators.size(); op++){
        std::set<int> pre_vars, eff_vars;
        g_operators[op].get_vars(pre_vars, eff_vars);
        for(int v : eff_vars){
            for(int v2 : pre_vars){
                weights[v][v2].insert(op);
                weights[v2][v].insert(op);
            }
            for(int v2 : eff_vars){
                weights[v][v2].insert(op);
                weights[v2][v].insert(op);
            }
        }
    }
}

void VariablePartitionGreedy::merge(int p1, int p2){
    part_size.push_back(part_size[p1]*part_size[p2]);
    part_size[p1] = 0; part_size[p2] = 0;
    int newp = partitions.size();
    partitions.push_back(partitions[p1]);
    partitions.back().insert(partitions.back().end(),
            begin(partitions[p2]),
            end(partitions[p2]));
    for(int i = 0; i < partitions.size(); i++){
        std::set_union(begin(weights[p1][i]), end(weights[p1][i]),
                begin(weights[p2][i]), end(weights[p2][i]),
                inserter(weights[newp][i], end(weights[newp][i])));


        std::set_union(begin(weights[i][p1]), end(weights[i][p1]),
                begin(weights[i][p2]), end(weights[i][p2]),
                inserter(weights[i][newp], end(weights[i][newp])));
    }
}


void VariablePartitionGreedy::find(){
    cout << "Find greedy" << endl;
    init();
    while(true) {
        //Greedily pick two partitions
        pair<int, int> parts = pick_parts();
        cout << "Picked: " << parts.first << " " << parts.second << endl;
        if(parts.first == -1) break;
        merge(parts.first, parts.second); //Merge both partitions
    }

    cout << "Find greedy loop done" << endl;
    vector<vector<int> > aux;
    for (int i = 0; i < partitions.size(); ++i){
        if(part_size[i])
            aux.push_back(partitions[i]);
    }
    partitions.swap(aux);
    cout << "Find greedy done: " << partitions.size() << " " << aux.size() << endl;
}


pair<int, int> VariablePartitionGreedy::pick_parts() const{
    auto best = pair<int, int> {-1, -1};
    int best_weight = -1;
    for(int i = 0; i < partitions.size(); ++i){
        if(!part_size[i]) continue;
        for(int j = i+1; j < partitions.size(); ++j){
            if(!part_size[j]) continue;
            int w = 0;
            if (weights.count(i) && weights.at(i).count(j)){
                weights.at(i).at(j).size();
            }
            //cout << w << " " << i << " " << j << endl;
            if(part_size[i]*part_size[j] < limit_size &&
                    w > best_weight){
                best_weight = w;
                best = pair<int, int> {i, j};
            }
        }
    }
    return best;
}
