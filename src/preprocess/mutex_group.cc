#include "mutex_group.h"

#include "helper_functions.h"
#include "variable.h"

MutexGroup::MutexGroup(istream &in, const vector<Variable *> &variables) : dir(FW) {
  //Mutex groups detected in the translator are "fw" mutexes
    int size;
    check_magic(in, "begin_mutex_group");
    in >> size;
    for (size_t i = 0; i < size; ++i) {
        int var_no, value;
        in >> var_no >> value;
        facts.push_back(make_pair(variables[var_no], value));
    }
    check_magic(in, "end_mutex_group");
}

MutexGroup::MutexGroup(const vector<pair<int, int> > & f, 
		       const vector<Variable *> &variables, 
		       bool regression) {
  if(regression){
    dir = BW;
  }else{
    dir = FW;
  }
  for (size_t i = 0; i < f.size(); ++i) {
    int var_no = f[i].first;
    int value = f[i].second;
    facts.push_back(make_pair(variables[var_no], value));
  }
}

MutexGroup::MutexGroup(const Variable * var) : dir(FW){
  for(int i = 0; i < var->get_range(); ++i){
    facts.push_back(make_pair(var, i));
  }
}


int MutexGroup::get_encoding_size() const {
    return facts.size();
}

void MutexGroup::dump() const {
    cout << "mutex group of size " << facts.size() << ":" << endl;
    for (size_t i = 0; i < facts.size(); ++i) {
        const Variable *var = facts[i].first;
        int value = facts[i].second;
        cout << "   " << var->get_name() << " = " << value
             << " (" << var->get_fact_name(value) << ")" << endl;
    }
}

void MutexGroup::generate_cpp_input(ofstream &outfile) const {
  string groupname = "mutex" ;
  string dirname = dir == FW ? "fw" : "bw";  
  outfile << "begin_mutex_group" << endl <<  groupname << endl << dirname << endl
	  << facts.size() << endl;
  for (size_t i = 0; i < facts.size(); ++i) {
    outfile << facts[i].first->get_level()
	    << " " << facts[i].second << endl;
  }
  outfile << "end_mutex_group" << endl;
}

void MutexGroup::strip_unimportant_facts() {
    int new_index = 0;
    for (int i = 0; i < facts.size(); i++) {
      if (facts[i].first->get_level() != -1 && facts[i].first->is_necessary())
	facts[new_index++] = facts[i];
    }
    facts.erase(facts.begin() + new_index, facts.end());
}

bool MutexGroup::is_redundant() const {
    // Only mutex groups that talk about two or more different
    // finite-domain variables are interesting.
    for (int i = 1; i < facts.size(); ++i)
        if (facts[i].first != facts[i - 1].first)
            return false;
    return true;
}

void strip_mutexes(vector<MutexGroup> &mutexes) {
    int old_count = mutexes.size();
    int new_index = 0;
    for (int i = 0; i < mutexes.size(); i++) {
        mutexes[i].strip_unimportant_facts();
        if (!mutexes[i].is_redundant())
            mutexes[new_index++] = mutexes[i];
    }
    mutexes.erase(mutexes.begin() + new_index, mutexes.end());
    cout << mutexes.size() << " of " << old_count
         << " mutex groups necessary." << endl;
}


void MutexGroup::remove_unreachable_facts(){
  vector<pair<const Variable *, int> > newfacts;
  for (int i = 0; i < facts.size(); ++i){
    if(facts[i].first->is_necessary() && facts[i].first->is_reachable(facts[i].second)){
      newfacts.push_back(make_pair(facts[i].first, facts[i].first->get_new_id(facts[i].second)));
    }
  }
  newfacts.swap(facts);
}

void MutexGroup::get_mutex_group(vector<pair<int, int>> &invariant_group) const {
    invariant_group.reserve(facts.size());
    for (size_t j = 0; j < facts.size(); ++j) {
        int var = facts[j].first->get_level();
        int val = facts[j].second;
        invariant_group.push_back(make_pair(var, val));
    }
}

void MutexGroup::generate_gamer_input(ofstream &outfile) const {
  string groupname = "mutex_group" ;
  string dirname = dir == FW ? "fw" : "bw";
  outfile << "begin_" << groupname << "_" << dirname << endl;
  for (size_t i = 0; i < facts.size(); ++i) {
    string fname = facts[i].first->get_fact_name(facts[i].second);
      outfile << fname;
    if(fname == "<none of those>"){
      bool first = true;
      for(int j = 0; j < facts[i].first->get_range(); ++j){
	if(facts[i].second != j){
	  outfile << (first ? ":" : ";") << facts[i].first->get_fact_name(j);
	  first = false;
	}
      }
    }
    outfile << endl;
  }
  outfile << "end_" << groupname << "_" << dirname << endl << endl;
}

bool MutexGroup::hasPair(int var, int val) const{
  for(int i = 0; i < facts.size(); ++i){
    if(facts[i].first->get_level () == var && facts[i].second == val){
      return true;
    }
  }
  return false;
}
