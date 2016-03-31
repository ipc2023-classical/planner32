#include "state.h"
#include "helper_functions.h"

class Variable;

State::State(istream &in, const vector<Variable *> &variables) {
    check_magic(in, "begin_state");
    for (Variable *var : variables) {
        int value;
        in >> value; //for axioms, this is default value
        values[var] = value;
    }
    check_magic(in, "end_state");
}

int State::operator[](Variable *var) const {
    return values.find(var)->second;
}

void State::dump() const {
    for (const auto &value : values)
        cout << "  " << value.first->get_name() << ": " << value.second << endl;
}

void State::remove_unreachable_facts(){
    map<Variable *, int> newvalues;
    for (auto it = values.begin();
       it != values.end(); ++it){
    if(it->first->is_necessary())
      newvalues[it->first] = it->first->get_new_id(it->second);
  }
  newvalues.swap(values);
}
