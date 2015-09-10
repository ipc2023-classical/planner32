#include "labels.h"

#include "label.h"
#include "label_reducer.h"
#include "factored_simulation.h"

#include "../globals.h"
#include "../utilities.h"
#include "../debug.h"
#include "../equivalence_relation.h"

#include <algorithm>
#include <cassert>
#include <set>


using namespace std;

Labels::Labels(bool unit_cost_, const Options &options, OperatorCost cost_type)
    : unit_cost(unit_cost_) {
    label_reducer = new LabelReducer(options);
    if (!g_operators.empty())
        labels.reserve(g_operators.size() * 2 - 1);
    for (size_t i = 0; i < g_operators.size(); ++i) {
      DEBUG_MSG(cout << "OPERATOR " << i << ": " << g_operators[i].get_name() << endl;);
        labels.push_back(new OperatorLabel(i, get_adjusted_action_cost(g_operators[i], cost_type),
                                           g_operators[i].get_prevail(), g_operators[i].get_pre_post()));
    }
}

Labels::~Labels() {
    delete label_reducer;
}

void Labels::reduce(pair<int, int> next_merge,
                    const std::vector<Abstraction *> &all_abstractions) {
    label_reducer->reduce_labels(next_merge, all_abstractions, labels);
}

void Labels::reduce(const LabelMap & labelMap, const FactoredSimulation & dominance_relation, 
		    std::set<int> & dangerous_LTSs) {
    EquivalenceRelation* equiv_rel = dominance_relation.get_equivalent_labels_relation(labelMap, 
										    dangerous_LTSs);
    label_reducer->reduce_exactly(equiv_rel, labels);
    delete equiv_rel;
}

void Labels::reduce_to_cost(){
    label_reducer->reduce_labels_to_cost(labels);
}

const Label *Labels::get_label_by_index(int index) const {
    assert(index >= 0 && index < labels.size());
    return labels[index];
}

bool Labels::is_label_reduced(int label_no) const {
    return get_label_by_index(label_no)->is_reduced();
}

void Labels::dump() const {
    cout << "no of labels: " << labels.size() << endl;
    for (size_t i = 0; i < labels.size(); ++i) {
        const Label *label = labels[i];
        label->dump();
    }
}

void Labels::dump_options() const {
    label_reducer->dump_options();
}

void Labels::set_irrelevant_for(int label_no, Abstraction * abstraction){
    labels[label_no]->set_irrelevant_for(abstraction);
}

void Labels::set_irrelevant_for_all_labels (Abstraction * abstraction){    
    for (size_t i = 0; i < labels.size(); ++i) {
	labels[i]->set_irrelevant_for(abstraction);
    }
}

void Labels::set_relevant_for(int label_no, Abstraction * abstraction){
    labels[label_no]->set_relevant_for(abstraction);
}

const std::set<Abstraction *> & Labels::get_relevant_for (int label_no) const{
    return labels[label_no]->get_relevant_for();
}

void Labels::prune_irrelevant_labels(){
    set<int> ops;
    for (auto label : labels) {
        if (label->is_irrelevant()) {
            label->get_operators(ops);
        }
    }
    cout << ops.size() << " irrelevant operators." << endl;
    for (int op : ops) {
        cout << "Irrelevant operator: " << g_operators[op].get_name() << endl;
    }
}

bool Labels::applies_perfect_label_reduction() const {
    return label_reducer->applies_perfect_label_reduction();
}
