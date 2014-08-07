#include "shrink_composite.h"

#include "abstraction.h"

#include "../option_parser.h"
#include "../plugin.h"

using namespace std;
ShrinkComposite::ShrinkComposite(const Options &opts)
    : ShrinkStrategy(opts), strategies(opts.get_list<ShrinkStrategy *>("strategies")) {    
}

ShrinkComposite::~ShrinkComposite() {
}

string ShrinkComposite::name() const {
    return "composite";
}

void ShrinkComposite::dump_strategy_specific_options() const {
    for(int i = 0; i < strategies.size(); ++i) {
	strategies[i]->dump_options();
    }
}

bool ShrinkComposite::reduce_labels_before_shrinking() const {
    for(int i = 0; i < strategies.size(); ++i) {
	if(strategies[i]->reduce_labels_before_shrinking()){
	    return true;
	}
    }
    return false;
}

void ShrinkComposite::shrink(Abstraction &abs, int target, bool force) {
    //TODO: This method does not make much sense in here. Should it be
    //removed from shrink_strategy?
    for(int i = 0; i < strategies.size(); ++i) {
	if(i!=0){
	    abs.normalize();
	    abs.compute_distances();
	}
	strategies[i]->shrink(abs, target, force);
    }   
}

void ShrinkComposite::shrink_atomic(Abstraction &abs) {
    for(int i = 0; i < strategies.size(); ++i) {
	if(i!=0){
	    abs.normalize();
	    abs.compute_distances();
	}
	strategies[i]->shrink_atomic(abs);
    }
}

void ShrinkComposite::shrink_before_merge(
    Abstraction &abs1, Abstraction &abs2) {
    for(int i = 0; i < strategies.size(); ++i) {
	if(i != 0) {
	    abs1.normalize();
	    abs2.normalize();
	    abs1.compute_distances();
	    abs2.compute_distances();
	}
	strategies[i]->shrink_before_merge(abs1, abs2);
    }
}

static ShrinkStrategy *_parse(OptionParser &parser) {
    ShrinkStrategy::add_options_to_parser(parser);
    parser.add_list_option<ShrinkStrategy *>("strategies");

    Options opts = parser.parse();

    if (parser.help_mode())
        return 0;

    ShrinkStrategy::handle_option_defaults(opts);

    opts.verify_list_non_empty<ShrinkStrategy *>("strategies");



    if (!parser.dry_run())
        return new ShrinkComposite(opts);
    else
        return 0;
}

static Plugin<ShrinkStrategy> _plugin("shrink_composite", _parse);
