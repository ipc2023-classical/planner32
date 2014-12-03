#ifndef MERGE_AND_SHRINK_LABELS_H
#define MERGE_AND_SHRINK_LABELS_H

#include "../operator_cost.h"

#include <vector>
#include <set>
#include <list>

class Abstraction;
class Label;
class LabelMap;
class LabelReducer;
class LabelRelation;
class Options;


/*
 The Labels class is basically a container class for the set of all
 labels used by merge-and-shrink abstractions.
 */
class Labels {
    const bool unit_cost;
    const LabelReducer *label_reducer;

    std::vector<Label *> labels;
public:
    Labels(bool unit_cost, const Options &options, OperatorCost cost_type);
    ~Labels();
    void reduce(std::pair<int, int> next_merge,
            const std::vector<Abstraction *> &all_abstractions);
    void reduce(const LabelMap & labelMap, const LabelRelation & label_dominance);
    // TODO: consider removing get_label_by_index and forwarding all required
    // methods of Label and giving access to them by label number.
    const Label *get_label_by_index(int index) const;
    bool is_label_reduced(int label_no) const;
    void dump() const;
    void dump_options() const;

    int get_size() const {
        return labels.size();
    }
    bool is_unit_cost() const {
        return unit_cost;
    }

    void set_relevant_for (int label_no, Abstraction * abstraction);
    void set_irrelevant_for (int label_no, Abstraction * abstraction);
    void set_irrelevant_for_all_labels (Abstraction * abstraction);

    const std::set<Abstraction *> & get_relevant_for (int label_no) const;

    void prune_irrelevant_labels();

    bool applies_perfect_label_reduction() const;
};


class LabelMap{
    //mapping from labels to labels for LTSs (hack to get rid of not useful labels)
    int num_valid_labels;
    std::vector<int> label_id;
    std::vector<int> old_label_id;
public:
    LabelMap(Labels * labels){
        num_valid_labels = 0;
        label_id.reserve(labels->get_size());
        old_label_id.reserve(labels->get_size());
        for(int i = 0; i < labels->get_size(); i++){
            if(labels->is_label_reduced(i)){
                label_id.push_back (-1);
            }else{
                old_label_id.push_back(i);
                label_id.push_back(num_valid_labels++);
            }
        }
    }

    int get_id(int i) const{
        return label_id[i];
    }
    int get_old_id(int i) const{
        return old_label_id[i];
    }

    int get_num_labels() const{
        return num_valid_labels;
    }
};

#endif
