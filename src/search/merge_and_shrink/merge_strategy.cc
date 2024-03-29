#include "merge_strategy.h"

#include "../globals.h"
#include "../plugin.h"

#include <iostream>

using namespace std;

MergeStrategy::MergeStrategy()
    : total_merges(g_variable_domain.size() - 1), 
      remaining_merges(g_variable_domain.size() - 1) {
    // There are number of variables many atomic abstractions and we have
    // to perform one less merges than this number until we have merged
    // all abstractions into one composite abstraction.
}

void MergeStrategy::dump_options() const {
    cout << "Merge strategy: " << name() << endl;
    dump_strategy_specific_options();
}



static MergeStrategy *_parse(OptionParser &) {
    return nullptr;
}

static Plugin<MergeStrategy> _plugin("none", _parse);
