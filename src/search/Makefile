# By default, build in 32-bit mode. Use BITWIDTH=64 (on command
# line or in environment) to build in 64-bit mode and BITWIDTH=native
# to use the native bitness of the OS.
#BITWIDTH ?= 32
BITWIDTH ?= 64

## Set USE_LP to 1 to enable linear programming stuff.
USE_LP=0

## Set LINK_RELEASE_STATICALLY to 0 or 1 (default) to disable/enable
## static linking of the executable in release mode.
## On OS X, this is unsupported and will be silently disabled.
LINK_RELEASE_STATICALLY=1

## On a supported operating system, there should be no need to override
## the OS setting. If the provided code does not work even though your
## operating system is a supported one, please report this as a bug.
OS=auto


HEADERS = \
          axioms.h \
          causal_graph.h \
          combining_evaluator.h \
          domain_transition_graph.h \
          eager_search.h \
          enforced_hill_climbing_search.h \
          equivalence_relation.h \
          exact_timer.h \
          g_evaluator.h \
          globals.h \
          heuristic.h \
          prune_heuristic.h \
          mutex_pruning.h \
          int_packer.h \
          ipc_max_heuristic.h \
          iterated_search.h \
          lazy_search.h \
          legacy_causal_graph.h \
          mutex_group.h \
          max_evaluator.h \
          operator.h \
          operator_cost.h \
          option_parser.h \
          option_parser_util.h \
          segmented_vector.h \
          per_state_information.h \
          pref_evaluator.h \
          relaxation_heuristic.h \
          rng.h \
          scc.h \
          search_engine.h \
          search_node_info.h \
          search_progress.h \
          search_space.h \
          state.h \
          state_id.h \
          state_registry.h \
          successor_generator.h \
          sum_evaluator.h \
          timer.h \
          utilities.h \
          weighted_evaluator.h \
          \
          open_lists/alternation_open_list.h \
          open_lists/open_list_buckets.h \
          open_lists/pareto_open_list.h \
          open_lists/standard_scalar_open_list.h \
          open_lists/tiebreaking_open_list.h \

## Each of the following "HEADERS += ..." constructs defines a
## "plugin" feature that can be enabled or disabled by simply
## commenting out the respective lines in the Makefile.

## As with all changes to the Makefile, you will have to force a
## rebuild after such a change. Deleting the executable and running
## "make" to relink is enough; no need to do a complete rebuild.

HEADERS += additive_heuristic.h
HEADERS += blind_search_heuristic.h
HEADERS += cea_heuristic.h
HEADERS += cg_heuristic.h cg_cache.h
HEADERS += ff_heuristic.h
HEADERS += goal_count_heuristic.h
HEADERS += hm_heuristic.h
HEADERS += lm_cut_heuristic.h
HEADERS += max_heuristic.h

HEADERS += merge_and_shrink/abstraction.h \
           merge_and_shrink/label.h \
           merge_and_shrink/labels.h \
           merge_and_shrink/label_reducer.h \
           merge_and_shrink/merge_and_shrink_heuristic.h \
           merge_and_shrink/merge_criterion.h \
           merge_and_shrink/merge_dfp.h \
           merge_and_shrink/merge_linear.h \
           merge_and_shrink/merge_linear_criteria.h \
           merge_and_shrink/merge_strategy.h \
           merge_and_shrink/shrink_bisimulation.h \
           merge_and_shrink/shrink_bucket_based.h \
           merge_and_shrink/shrink_composite.h \
           merge_and_shrink/shrink_fh.h \
           merge_and_shrink/shrink_own_labels.h \
           merge_and_shrink/shrink_random.h \
           merge_and_shrink/shrink_strategy.h \
           merge_and_shrink/variable_order_finder.h \
           merge_and_shrink/simulation_relation.h \
           merge_and_shrink/dominance_relation.h \
           merge_and_shrink/dominance_pruning_simulation.h \
           merge_and_shrink/alternative_label_relation.h \
           merge_and_shrink/label_relation.h \
           merge_and_shrink/label_relation_identity.h \
           merge_and_shrink/label_relation_noop.h \
           merge_and_shrink/labelled_transition_system.h \
           merge_and_shrink/simulation_simple.h \
           merge_and_shrink/simulation_identity.h \
           merge_and_shrink/ld_simulation.h \
           merge_and_shrink/variable_partition_finder.h \
           merge_and_shrink/opt_order.h \
	   merge_and_shrink/abstraction_builder.h \
#          merge_and_shrink/variable_partition_dominance.h
	 # merge_and_shrink/simulation_complex.h \
#           merge_and_shrink/simulation_complex_nold.h \

HEADERS += numeric_dominance/numeric_simulation_relation.h \
           numeric_dominance/numeric_dominance_relation.h \
           numeric_dominance/numeric_label_relation.h \
           numeric_dominance/numeric_dominance_pruning.h \
           numeric_dominance/int_epsilon.h \
           numeric_dominance/dijkstra_search_epsilon.h \
	   numeric_dominance/tau_labels.h \
	   numeric_dominance/breadth_first_search.h

HEADERS += landmarks/exploration.h \
           landmarks/h_m_landmarks.h \
           landmarks/lama_ff_synergy.h \
           landmarks/landmark_cost_assignment.h \
           landmarks/landmark_count_heuristic.h \
           landmarks/landmark_status_manager.h \
           landmarks/landmark_graph_merged.h \
           landmarks/landmark_graph.h \
           landmarks/landmark_factory.h \
           landmarks/landmark_factory_rpg_exhaust.h \
           landmarks/landmark_factory_rpg_sasp.h \
           landmarks/landmark_factory_zhu_givan.h \
           landmarks/util.h \

# HEADERS += learning/AODE.h \
#            learning/classifier.h \
#            learning/composite_feature_extractor.h \
#            learning/feature_extractor.h \
#            learning/maximum_heuristic.h \
#            learning/naive_bayes_classifier.h \
#            learning/PDB_state_space_sample.h \
#            learning/probe_state_space_sample.h \
#            learning/selective_max_heuristic.h \
#            learning/state_space_sample.h \
#            learning/state_vars_feature_extractor.h \

HEADERS += pdbs/canonical_pdbs_heuristic.h \
           pdbs/dominance_pruner.h \
           pdbs/match_tree.h \
           pdbs/max_cliques.h \
           pdbs/pattern_generation_edelkamp.h \
           pdbs/pattern_generation_haslum.h \
           pdbs/pdb_heuristic.h \
           pdbs/util.h \
           pdbs/zero_one_pdbs_heuristic.h \

HEADERS += sym/sym_controller.h \
	   sym/sym_engine.h \
           sym/symba.h \
           sym/symba_unsat.h \
           sym/sym_variables.h \
           sym/sym_manager.h \
           sym/sym_util.h \
           sym/sym_enums.h \
           sym/sym_abstraction.h \
           sym/sym_pdb.h \
           sym/sym_ph_minimization.h \
           sym/sym_abs_minimization.h \
           sym/sym_transition.h \
           sym/sym_params.h \
           sym/sym_heuristic.h \
           sym/sym_bucket.h \
           sym/sym_exploration.h \
	   sym/sym_astar.h \
           sym/sym_astar_closed.h \
           sym/sym_astar_open.h \
	   sym/sym_breadth_first_search.h \
           sym/sym_hnode.h \
           sym/sym_bdexp.h \
           sym/sym_ph.h \
           sym/sym_ph_pdbs.h \
           sym/sym_estimate.h \
           sym/sym_estimate_linear.h \
           sym/sym_solution.h \
           sym/spmas_heuristic.h \
	   sym/test/sym_test.h \
	   sym/test/sym_test_mutex.h \
	   sym/sym_prune_heuristic.h \
	   sym/sym_uct_pdbs.h \
           #sym/sym_ph_smas.h \
           sym/smas_shrink_bisimulation.h \
           sym/smas_shrink_bucket_based.h \
           sym/smas_shrink_fh.h \
           sym/smas_shrink_random.h \
           sym/smas_shrink_all.h \
           sym/smas_shrink_strategy.h \
           sym/smas_abs_state.h \
           sym/smas_shrink_state.h \
           sym/sym_smas.h \

SHELL = /bin/bash

SOURCES = planner.cc $(HEADERS:%.h=%.cc)
TARGET  = downward

ARGS_PROFILE = --search 'astar(lmcut())' < profile-input.pre


ifeq ($(OS), auto)
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
OS=osx
else ifeq ($(UNAME), Linux)
OS=linux
else
UNAME_O := $(shell uname -o)
endif

ifeq ($(UNAME_O), Cygwin)
## For now, we treat Linux and Cygwin equally in this Makefile, so we
## just call Cygwin "linux".
OS=linux
else ifeq ($(OS), auto)
$(warning OS detection failed -- setting to Linux and hoping for the best!)
OS=linux
endif

endif # OS == auto

ifeq ($(OS), osx)
## Disable static linking on OS X.
LINK_RELEASE_STATICALLY=0
endif


OBJECT_SUFFIX_RELEASE = .release
TARGET_SUFFIX_RELEASE = -release
OBJECT_SUFFIX_DEBUG   = .debug
TARGET_SUFFIX_DEBUG   = -debug
OBJECT_SUFFIX_PROFILE = .profile
TARGET_SUFFIX_PROFILE = -profile

OBJECTS_RELEASE = $(SOURCES:%.cc=.obj/%$(OBJECT_SUFFIX_RELEASE).o)
TARGET_RELEASE  = $(TARGET)$(TARGET_SUFFIX_RELEASE)

OBJECTS_DEBUG   = $(SOURCES:%.cc=.obj/%$(OBJECT_SUFFIX_DEBUG).o)
TARGET_DEBUG    = $(TARGET)$(TARGET_SUFFIX_DEBUG)

OBJECTS_PROFILE = $(SOURCES:%.cc=.obj/%$(OBJECT_SUFFIX_PROFILE).o)
TARGET_PROFILE  = $(TARGET)$(TARGET_SUFFIX_PROFILE)

CXX     = g++
DEPEND = g++ -MM

## CXXFLAGS, LDFLAGS, POSTLINKOPT are options for compiler and linker
## that are used for all three targets (release, debug, and profile).
## (POSTLINKOPT are options that appear *after* all object files.)

ifeq ($(BITWIDTH), 32)
BITWIDTHOPT = -m32
else ifeq ($(BITWIDTH), 64)
BITWIDTHOPT = -m64
else ifneq ($(BITWIDTH), native)
$(error Bad value for BITWIDTH)
endif

CXXFLAGS = -Iext
# if you have Mac OS, you may be getting tree.hh from /opt/local/include
ifeq ($(OS), osx)
  CXXFLAGS += -I/opt/local/include/
endif

CXXFLAGS += -g
CXXFLAGS += $(BITWIDTHOPT)
CXXFLAGS += -Wall -W -Wno-sign-compare -Wno-deprecated -pedantic -Werror

## The following lines contain workarounds for bugs when
## cross-compiling to 64 bit on 32-bit systems using gcc 4.4 or gcc
## 4.5 in some Ubuntu releases. (We don't usually cross-compile to
## 64-bit, but in some cases we do; e.g. we did for the IPC.) See
## http://stackoverflow.com/questions/4643197/missing-include-bits-cconfig-h-when-cross-compiling-64-bit-program-on-32-bit.
ifeq ($(OS), linux) # workarounds only valid for Linux...
HAVE_GCC_4_4 := $(shell expr "$$(gcc -dumpversion)" : \\\(4\.4\.\\\))
HAVE_GCC_4_5 := $(shell expr "$$(gcc -dumpversion)" : \\\(4\.5\.\\\))

ifdef HAVE_GCC_4_4
CXXFLAGS += -I/usr/include/c++/4.4/i686-linux-gnu
endif

ifdef HAVE_GCC_4_5
CXXFLAGS += -I/usr/include/c++/4.5/i686-linux-gnu
endif
endif

LDFLAGS += $(BITWIDTHOPT)

POSTLINKOPT =

## Additional specialized options for the various targets follow.
## In release mode, we link statically since this makes it more likely
## that local compiles will work on the various grids (gkigrid, Black
## Forest Grid, maia).
##
## NOTE: This precludes some uses of exceptions.
##        For details, see man gcc on -static-libgcc.

CXXFLAGS_RELEASE  = -O3 -DNDEBUG -fomit-frame-pointer
CXXFLAGS_DEBUG    = -O3
CXXFLAGS_PROFILE  = -pg -O2

LDFLAGS_RELEASE =
LDFLAGS_DEBUG    =
LDFLAGS_PROFILE  = -pg

POSTLINKOPT_RELEASE =
POSTLINKOPT_DEBUG   =
POSTLINKOPT_PROFILE =

ifeq ($(LINK_RELEASE_STATICALLY), 1)
LDFLAGS_RELEASE += -static -static-libgcc
endif

ifeq ($(OS), linux)
ifeq ($(LINK_RELEASE_STATICALLY), 0)
POSTLINKOPT_RELEASE += -lrt
else
POSTLINKOPT_RELEASE += -Wl,-Bstatic -lrt
endif
POSTLINKOPT_DEBUG  += -lrt
POSTLINKOPT_PROFILE += -lrt
# POSTLINKOPT_PROFILE += -lc_p
endif


## Define the default target up here so that the LP stuff below
## doesn't define a default target.

default: release

ifeq ($(USE_LP),1)

COIN_ROOT = lp/coin

$(COIN_ROOT):
	cd lp && ./setup

$(TARGET_RELEASE) $(OBJECTS_RELEASE): $(COIN_ROOT)
$(TARGET_DEBUG) $(OBJECTS_DEBUG): $(COIN_ROOT)
$(TARGET_PROFILE) $(OBJECTS_PROFILE): $(COIN_ROOT)

## We want to link the Linear Programming libraries statically since
## they are unlikely to be preinstalled on the grids we use for
## evaluation. Static linking is a bit tricky: we need to specify the
## libraries *after* the source files and in such an order that if A
## depends on B, A is listed before B. (In case of dependency cycles,
## we can and must list them multiple times.) The following set of
## libraries and their ordering have been determined experimentally
## and hence might break if we use more functions from the LP
## libraries. See
## http://ask.metafilter.com/117792/How-to-fix-C-static-linking-problems

COIN_LIBS = OsiClp Clp CoinUtils Osi

## We want to always link the COIN libraries statically, even if static
## linking is otherwise disabled. We accomplish this by
## using -Wl,-Bstatic before the COIN libs and -Wl,-Bdynamic
## afterwards (unless in release mode with static linking enabled). See
## http://ubuntuforums.org/showthread.php?t=491455

COIN_CXXFLAGS = -I$(COIN_ROOT)/include/coin -D USE_LP -D COIN_USE_CLP
ifeq ($(OS), osx)
COIN_LDFLAGS = -L$(COIN_ROOT)/lib $(COIN_LIBS:%=-l %)
else
COIN_LDFLAGS = -L$(COIN_ROOT)/lib -Wl,-Bstatic $(COIN_LIBS:%=-l %)
endif

CXXFLAGS += $(COIN_CXXFLAGS)

POSTLINKOPT_RELEASE += $(COIN_LDFLAGS)
POSTLINKOPT_DEBUG += $(COIN_LDFLAGS)
POSTLINKOPT_PROFILE += $(COIN_LDFLAGS)

ifeq ($(OS), linux)
ifeq ($(LINK_RELEASE_STATICALLY), 0)
POSTLINKOPT_RELEASE += -Wl,-Bdynamic
endif
POSTLINKOPT_DEBUG   += -Wl,-Bdynamic
POSTLINKOPT_PROFILE += -Wl,-Bdynamic
endif

endif # LP code

CUDD_ROOT = cudd-3.0.0/cudd/.libs/libcudd.a

$(CUDD_ROOT):
	cd cudd-3.0.0 &&  aclocal && autoheader && automake --add-missing && autoconf && ./configure --enable-obj --build=i686-pc-linux-gnu "CFLAGS=-m64 -Wall -Wextra -g -O3" "CXXFLAGS=-m64 -Wall -Wextra -std=c++0x -g -O3" "LDFLAGS=-m64" && make && make cplusplus



$(TARGET_RELEASE) $(OBJECTS_RELEASE): $(CUDD_ROOT)
$(TARGET_DEBUG) $(OBJECTS_DEBUG): $(CUDD_ROOT)
$(TARGET_PROFILE) $(OBJECTS_PROFILE): $(CUDD_ROOT)

CUDD_CCOPT = -Icudd-3.0.0/cudd -Icudd-3.0.0/cplusplus

CUDD_LINKOPT = cudd-3.0.0/cudd/.libs/libcudd.a \
   cudd-3.0.0/dddmp/.libs/libdddmp.a

#CUDD Includes
CXXFLAGS += $(CUDD_CCOPT)

POSTLINKOPT_RELEASE += $(CUDD_LINKOPT)
POSTLINKOPT_DEBUG += $(CUDD_LINKOPT)
POSTLINKOPT_PROFILE += $(CUDD_LINKOPT)


#end CUDD code

all: release debug profile

## Build rules for the release target follow.

release: $(TARGET_RELEASE)

$(TARGET_RELEASE): $(OBJECTS_RELEASE)
	$(CXX) $(LDFLAGS) $(LDFLAGS_RELEASE) $(OBJECTS_RELEASE) $(POSTLINKOPT) $(POSTLINKOPT_RELEASE) -o $(TARGET_RELEASE)

$(OBJECTS_RELEASE): .obj/%$(OBJECT_SUFFIX_RELEASE).o: %.cc
	@mkdir -p $$(dirname $@)
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_RELEASE) -c $< -o $@

## Build rules for the debug target follow.

debug: $(TARGET_DEBUG)

$(TARGET_DEBUG): $(OBJECTS_DEBUG)
	$(CXX) $(LDFLAGS) $(LDFLAGS_DEBUG) $(OBJECTS_DEBUG) $(POSTLINKOPT) $(POSTLINKOPT_DEBUG) -o $(TARGET_DEBUG)

$(OBJECTS_DEBUG): .obj/%$(OBJECT_SUFFIX_DEBUG).o: %.cc
	@mkdir -p $$(dirname $@)
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_DEBUG) -c $< -o $@

## Build rules for the profile target follow.

profile: $(TARGET_PROFILE)

$(TARGET_PROFILE): $(OBJECTS_PROFILE)
	$(CXX) $(LDFLAGS) $(LDFLAGS_PROFILE) $(OBJECTS_PROFILE) $(POSTLINKOPT) $(POSTLINKOPT_PROFILE) -o $(TARGET_PROFILE)

$(OBJECTS_PROFILE): .obj/%$(OBJECT_SUFFIX_PROFILE).o: %.cc
	@mkdir -p $$(dirname $@)
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_PROFILE) -c $< -o $@

## Additional targets follow.

PROFILE: $(TARGET_PROFILE)
	./$(TARGET_PROFILE) $(ARGS_PROFILE)
	gprof $(TARGET_PROFILE) | (cleanup-profile 2> /dev/null || cat) > PROFILE

clean:
	rm -rf .obj
	rm -f *~ *.pyc
	rm -f Makefile.depend gmon.out PROFILE core
	rm -f sas_plan

distclean: clean
	rm -f $(TARGET_RELEASE) $(TARGET_DEBUG) $(TARGET_PROFILE)

## NOTE: If we just call gcc -MM on a source file that lives within a
## subdirectory, it will strip the directory part in the output. Hence
## the for loop with the sed call.

Makefile.depend: $(CUDD_ROOT) $(SOURCES) $(HEADERS)
	rm -f Makefile.temp
	for source in $(SOURCES) ; do \
	    $(DEPEND) $(CXXFLAGS) $$source > Makefile.temp0; \
	    objfile=$${source%%.cc}.o; \
	    sed -i -e "s@^[^:]*:@$$objfile:@" Makefile.temp0; \
	    cat Makefile.temp0 >> Makefile.temp; \
	done
	rm -f Makefile.temp0 Makefile.depend
	sed -e "s@\(.*\)\.o:\(.*\)@.obj/\1$(OBJECT_SUFFIX_RELEASE).o:\2@" Makefile.temp >> Makefile.depend
	sed -e "s@\(.*\)\.o:\(.*\)@.obj/\1$(OBJECT_SUFFIX_DEBUG).o:\2@" Makefile.temp >> Makefile.depend
	sed -e "s@\(.*\)\.o:\(.*\)@.obj/\1$(OBJECT_SUFFIX_PROFILE).o:\2@" Makefile.temp >> Makefile.depend
	rm -f Makefile.temp

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include Makefile.depend
endif
endif

.PHONY: default all release debug profile clean distclean
