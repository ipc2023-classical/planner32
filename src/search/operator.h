#ifndef OPERATOR_H
#define OPERATOR_H

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "globals.h"
#include "state.h"

struct Prevail {
    int var;
    int prev;
    Prevail(std::istream &in);
    Prevail(int v, int p) : var(v), prev(p) {}

    bool is_applicable(const State &state) const {
        assert(var >= 0 && var < g_variable_name.size());
        assert(prev >= 0 && prev < g_variable_domain[var]);
        return state[var] == prev;
    }

    bool operator==(const Prevail &other) const {
        return var == other.var && prev == other.prev;
    }

    bool operator!=(const Prevail &other) const {
        return !(*this == other);
    }

    bool operator<(const Prevail &other) const {
        return var < other.var || (var == other.var && prev < other.prev);
    }

    void dump() const;
};

struct PrePost {
    int var;
    int pre, post;
    std::vector<Prevail> cond;
    PrePost() {} // Needed for axiom file-reading constructor, unfortunately.
    PrePost(std::istream &in);
    PrePost(int v, int pr, int po, const std::vector<Prevail> &co)
        : var(v), pre(pr), post(po), cond(co) {}

    bool is_applicable(const State &state) const {
        assert(var >= 0 && var < g_variable_name.size());
        assert(pre == -1 || (pre >= 0 && pre < g_variable_domain[var]));
        return pre == -1 || state[var] == pre;
    }

    bool does_fire(const State &state) const {
        for (int i = 0; i < cond.size(); i++)
            if (!cond[i].is_applicable(state))
                return false;
        return true;
    }

    bool operator<(const PrePost &other) const {
        //TODO: Conditions are ignored
        return var < other.var || (var == other.var && (pre < other.pre || (pre == other.pre && post < other.post)));
    }


    void dump() const;
};

class Operator {
    bool is_an_axiom;
    std::vector<Prevail> prevail;      // var, val
    std::vector<PrePost> pre_post;     // var, old-val, new-val, effect conditions
    std::string name;
    int cost;

    mutable bool marked = false; // Used for short-term marking of preferred operators
    bool dead = false;
public:
    Operator(std::istream &in, bool is_axiom);

    Operator (bool _is_an_axiom,
                 std::vector<Prevail> && _prevail,
                 std::vector<PrePost> && _pre_post,
                 std::string _name,
                 int _cost) : is_an_axiom(_is_an_axiom),
                              prevail(std::move(_prevail)),
                              pre_post (std::move(_pre_post)),
                              name(_name), cost(_cost) {
    }

    Operator(const Operator & op) = default;

    void dump() const;
    std::string get_name() const {return name; }

    bool is_axiom() const {return is_an_axiom; }

    void set_dead() {
	dead = true;
    }

    bool is_dead() const{
	return dead;
    }

    const std::vector<Prevail> &get_prevail() const {return prevail; }
    const std::vector<PrePost> &get_pre_post() const {return pre_post; }


    void get_vars(std::set<int> & pre_vars, std::set<int> & eff_vars) const {
	for (auto & p : prevail){
	    pre_vars.insert(p.var);
	}
	for (auto & p : pre_post){
	    eff_vars.insert(p.var);
	    if(p.pre != -1){
		pre_vars.insert(p.var);
	    }
	}
    }

    bool is_applicable(const State &state) const {
        for (int i = 0; i < prevail.size(); i++)
            if (!prevail[i].is_applicable(state))
                return false;
        for (int i = 0; i < pre_post.size(); i++)
            if (!pre_post[i].is_applicable(state))
                return false;
        return true;
    }

    bool is_marked() const {
        return marked;
    }
    void mark() const {
        marked = true;
    }
    void unmark() const {
        marked = false;
    }

    mutable bool marker1, marker2; // HACK! HACK!

    int get_cost() const {return cost; }
};

#endif
