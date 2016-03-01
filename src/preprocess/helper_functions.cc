#include <cstdlib>
#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <fstream>


#include "helper_functions.h"
#include "state.h"
#include "mutex_group.h"
#include "operator.h"
#include "axiom.h"
#include "variable.h"
#include "successor_generator.h"
#include "domain_transition_graph.h"

using namespace std;


static const int SAS_FILE_VERSION = 3;
static const int PRE_FILE_VERSION = SAS_FILE_VERSION;


void check_magic(istream &in, string magic) {
    string word;
    in >> word;
    if (word != magic) {
        cerr << "Failed to match magic word '" << magic << "'." << endl;
        cerr << "Got '" << word << "'." << endl;
        if (magic == "begin_version") {
            cerr << "Possible cause: you are running the preprocessor "
                 << "on a translator file from an " << endl
                 << "older version." << endl;
        }
        exit(1);
    }
}

void read_and_verify_version(istream &in) {
    int version;
    check_magic(in, "begin_version");
    in >> version;
    check_magic(in, "end_version");
    if (version != SAS_FILE_VERSION) {
        cerr << "Expected translator file version " << SAS_FILE_VERSION
             << ", got " << version << "." << endl;
        cerr << "Exiting." << endl;
        exit(1);
    }
}

void read_metric(istream &in, bool &metric) {
    check_magic(in, "begin_metric");
    in >> metric;
    check_magic(in, "end_metric");
}

void read_variables(istream &in, vector<Variable> &internal_variables,
                    vector<Variable *> &variables) {
    int count;
    in >> count;
    internal_variables.reserve(count);
    // Important so that the iterators stored in variables are valid.
    for (int i = 0; i < count; i++) {
        internal_variables.push_back(Variable(in));
        variables.push_back(&internal_variables.back());
    }
}

void read_mutexes(istream &in, vector<MutexGroup> &mutexes,
                  const vector<Variable *> &variables) {
    size_t count;
    in >> count;
    for (size_t i = 0; i < count; ++i)
        mutexes.push_back(MutexGroup(in, variables));
}

void read_goal(istream &in, const vector<Variable *> &variables,
               vector<pair<Variable *, int> > &goals) {
    check_magic(in, "begin_goal");
    int count;
    in >> count;
    for (int i = 0; i < count; i++) {
        int varNo, val;
        in >> varNo >> val;
        goals.push_back(make_pair(variables[varNo], val));
    }
    check_magic(in, "end_goal");
}

void dump_goal(const vector<pair<Variable *, int> > &goals) {
    cout << "Goal Conditions:" << endl;
    for (int i = 0; i < goals.size(); i++)
        cout << "  " << goals[i].first->get_name() << ": "
             << goals[i].second << endl;
}

void read_operators(istream &in, const vector<Variable *> &variables,
                    vector<Operator> &operators) {
    int count;
    in >> count;
    for (int i = 0; i < count; i++)
        operators.push_back(Operator(in, variables));
}

void read_axioms(istream &in, const vector<Variable *> &variables,
                 vector<Axiom> &axioms) {
    int count;
    in >> count;
    for (int i = 0; i < count; i++)
        axioms.push_back(Axiom(in, variables));
}

void read_preprocessed_problem_description(istream &in,
                                           bool &metric,
                                           vector<Variable> &internal_variables,
                                           vector<Variable *> &variables,
                                           vector<MutexGroup> &mutexes,
                                           State &initial_state,
                                           vector<pair<Variable *, int> > &goals,
                                           vector<Operator> &operators,
                                           vector<Axiom> &axioms) {
    read_and_verify_version(in);
    read_metric(in, metric);
    read_variables(in, internal_variables, variables);
    read_mutexes(in, mutexes, variables);
    initial_state = State(in, variables);
    read_goal(in, variables, goals);
    read_operators(in, variables, operators);
    read_axioms(in, variables, axioms);
}

void dump_preprocessed_problem_description(const vector<Variable *> &variables,
                                           const State &initial_state,
                                           const vector<pair<Variable *, int> > &goals,
                                           const vector<Operator> &operators,
                                           const vector<Axiom> &axioms) {
    cout << "Variables (" << variables.size() << "):" << endl;
    for (int i = 0; i < variables.size(); i++)
        variables[i]->dump();

    cout << "Initial State:" << endl;
    initial_state.dump();
    dump_goal(goals);

    for (int i = 0; i < operators.size(); i++)
        operators[i].dump();
    for (int i = 0; i < axioms.size(); i++)
        axioms[i].dump();
}

void dump_DTGs(const vector<Variable *> &ordering,
               vector<DomainTransitionGraph> &transition_graphs) {
    for (int i = 0; i < transition_graphs.size(); i++) {
        cout << "Domain transition graph for " << ordering[i]->get_name() << ":" << endl;
        transition_graphs[i].dump();
    }
}

void generate_cpp_input(bool /*solvable_in_poly_time*/,
                        const vector<Variable *> &ordered_vars,
                        const bool &metric,
                        const vector<MutexGroup> &mutexes,
                        const State &initial_state,
                        const vector<pair<Variable *, int> > &goals,
                        const vector<Operator> &operators,
                        const vector<Axiom> &axioms,
                        const SuccessorGenerator &sg,
                        const vector<DomainTransitionGraph> transition_graphs,
                        const CausalGraph &cg) {
    /* NOTE: solvable_in_poly_time flag is no longer included in output,
       since the planner doesn't handle it specially any more anyway. */

    ofstream outfile;
    outfile.open("output", ios::out);

    outfile << "begin_version" << endl;
    outfile << PRE_FILE_VERSION << endl;
    outfile << "end_version" << endl;

    outfile << "begin_metric" << endl;
    outfile << metric << endl;
    outfile << "end_metric" << endl;

    outfile << ordered_vars.size() << endl;
    for (int i = 0; i < ordered_vars.size(); i++)
        ordered_vars[i]->generate_cpp_input(outfile);

    outfile << mutexes.size() << endl;
    for (int i = 0; i < mutexes.size(); i++)
        mutexes[i].generate_cpp_input(outfile);

    int var_count = ordered_vars.size();
    outfile << "begin_state" << endl;
    for (int i = 0; i < var_count; i++)
        outfile << initial_state[ordered_vars[i]] << endl;  // for axioms default value
    outfile << "end_state" << endl;

    vector<int> ordered_goal_values;
    ordered_goal_values.resize(var_count, -1);
    for (int i = 0; i < goals.size(); i++) {
        int var_index = goals[i].first->get_level();
        ordered_goal_values[var_index] = goals[i].second;
    }
    outfile << "begin_goal" << endl;
    outfile << goals.size() << endl;
    for (int i = 0; i < var_count; i++)
        if (ordered_goal_values[i] != -1)
            outfile << i << " " << ordered_goal_values[i] << endl;
    outfile << "end_goal" << endl;

    outfile << operators.size() << endl;
    for (int i = 0; i < operators.size(); i++)
        operators[i].generate_cpp_input(outfile);

    outfile << axioms.size() << endl;
    for (int i = 0; i < axioms.size(); i++)
        axioms[i].generate_cpp_input(outfile);

    outfile << "begin_SG" << endl;
    sg.generate_cpp_input(outfile);
    outfile << "end_SG" << endl;

    for (int i = 0; i < var_count; i++) {
        outfile << "begin_DTG" << endl;
        transition_graphs[i].generate_cpp_input(outfile);
        outfile << "end_DTG" << endl;
    }

    outfile << "begin_CG" << endl;
    cg.generate_cpp_input(outfile, ordered_vars);
    outfile << "end_CG" << endl;

    outfile.close();
}


void generate_unsolvable_cpp_input() {
    ofstream outfile;
    outfile.open("output", ios::out);
    outfile << "begin_version" << endl;
    outfile << PRE_FILE_VERSION << endl;
    outfile << "end_version" << endl;

    outfile << "begin_metric" << endl << "0" << "end_metric" << endl;

    //variables
    outfile << "1" << endl << "begin_variable" << endl
            << "var1" << endl
            << "0" << endl
            << "2" << endl
	    << "val0" << endl
	    << "val1" << endl
	    << "end_variable" << endl;

    //Mutexes 
    outfile << "0" << endl;

    //Initial state and goal
    outfile << "begin_state" << endl << "0" << endl << "end_state" << endl;
    outfile << "begin_goal" << endl << "1" << endl << "0 1" << endl << "end_goal" << endl;

    //Operators
    outfile << "0" << endl;

    //Axioms
    outfile << "0" << endl;

    outfile << "begin_SG" << endl;
    outfile << "check 0" << endl;
    outfile << "end_SG" << endl;

    outfile << "begin_DTG" << endl << "0" << endl  << "0" << endl
     << "end_DTG" << endl;

    outfile << "begin_CG" << endl << "0" << endl << "end_CG" << endl;

    outfile.close();
}



void write_variables(vector<Variable *> & variables, string filename){
  ofstream vars_file; 
  vars_file.open(filename.c_str());
  write_variables(variables, vars_file);
  vars_file.close();
}



void write_variables(vector<Variable *> & variables, ofstream & file){
  for(int i = 0; i < variables.size(); ++i){
    Variable * var = variables[i];
    if(var->is_necessary()){
      for(int j = 0; j < var->get_range(); ++j){
	if(var->is_reachable(j)){
	  file << var->get_fact_name(j) << endl;
	}
      }
      file << endl;
    }
  }
  
}

void write_operators(vector<Operator> & ops, string filename){
  ofstream file; 
  file.open(filename.c_str());
  write_operators(ops, file);
  file.close();
}

void write_operators(vector<Operator> & ops, ofstream & file){
  for(int i = 0; i < ops.size(); ++i){
    if (!ops[i].is_redundant())
      file << ops[i].get_name() << endl;
  }
}


void write_mutexes(vector<Variable*> & variables, vector<MutexGroup> & m, string filename){
  ofstream mutexes_file; 
  mutexes_file.open(filename.c_str());
  write_mutexes(variables, m, mutexes_file);
  mutexes_file.close();
}

void write_mutexes(vector<Variable*> & variables, vector<MutexGroup> & m, ofstream & file){
  for(int i = 0; i < variables.size(); ++i){
    if(variables[i]->is_necessary()){
      MutexGroup g_v (variables[i]);
      g_v.generate_gamer_input(file);
    }
  }
  for(int i = 0; i < m.size(); ++i){
    const MutexGroup & mutex = m[i];
    if(!mutex.is_redundant()){
      mutex.generate_gamer_input(file);
    }
  }
}


