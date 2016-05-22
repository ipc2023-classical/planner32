#ifndef MERGE_AND_SHRINK_LABEL_H
#define MERGE_AND_SHRINK_LABEL_H

#include "../operator.h"

#include <vector>
#include <set>

class CompositeLabel;
class Abstraction;

/* This class implements labels as used by merge-and-shrink abstractions.
   It abstracts from the underlying regular operators and allows to store
   additional information associated with labels.
   NOTE: operators that are axioms are currently not supported! */

class Label {
    friend class CompositeLabel; // required for update_root
    int id;
    int cost;
    // prevail and pre_posts are references to those of one "canonical"
    // operator, which is the operator an OperatorLabel was built from or
    // the "first" label of all parent labels when constructing a CompositeLabel.
    const std::vector<Prevail> &prevail;
    const std::vector<PrePost> &pre_post;
    
    //Alvaro: Sets of abstraction the label is relevant for. This is
    //needed to compute the own-labels for an abstraction (those that
    //are only relevant for the abstraction).
    std::set<Abstraction *> relevant_for; 
protected:
    // root is a pointer to a composite label that this label has been reduced
    // to, if such a label exists, or to itself, if the label has not been
    // reduced yet.
    Label *root;

    Label(int id, int cost, const std::vector<Prevail> &prevail,
          const std::vector<PrePost> &pre_post);
    virtual ~Label() {}
    virtual void update_root(CompositeLabel *new_root) = 0;
public:
    const std::vector<Prevail> &get_prevail() const {return prevail; }
    const std::vector<PrePost> &get_pre_post() const {return pre_post; }
    int get_id() const {
        return id;
    }
    int get_cost() const {
        return cost;
    }
    bool is_reduced() const;
    virtual const std::vector<Label *> &get_parents() const = 0;
    void dump() const;

    //Alvaro: Methods to access relevant_for.
    void set_relevant_for(Abstraction * abstraction);
    void set_irrelevant_for(Abstraction * abstraction);
    bool is_relevant_for(Abstraction * abstraction) const;


    void reset_relevant_for (const std::vector<Abstraction *> &  abstractions); 


    const std::set<Abstraction *> & get_relevant_for () const{
	return relevant_for;
    }

    bool is_irrelevant() const;
    virtual void get_operators(std::set<int> & ops_id) const = 0;
};

class OperatorLabel : public Label {
  //std::string name;
    void update_root(CompositeLabel *new_root);
public:
    OperatorLabel(int id, int cost, const std::vector<Prevail> &prevail,
                  const std::vector<PrePost> &pre_post);

    const std::vector<Label *> &get_parents() const;
    virtual void get_operators(std::set<int> & ops_id) const{
	ops_id.insert(get_id());
    }
};

class CompositeLabel : public Label {
    std::vector<Label *> parents;

    void update_root(CompositeLabel *new_root);
public:
    CompositeLabel(int id, const std::vector<Label *> &parents);
    const std::vector<Label *> &get_parents() const {
        return parents;
    }

    virtual void get_operators(std::set<int> & ops_id) const{
	for (auto p : parents){
	    p->get_operators(ops_id);
	}
    }

};

#endif
