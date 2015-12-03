#ifndef MERGE_AND_SHRINK_DOMINANCE_PRUNING_SIMULATION_H
#define MERGE_AND_SHRINK_DOMINANCE_PRUNING_SIMULATION_H

#include "../prune_heuristic.h"
#include "../sym/sym_variables.h"
#include "../sym/sym_params.h"

class LDSimulation;
class SymVariables;
class SymManager;

enum class PruningDD {BDD_MAP, ADD, BDD, BDD_MAP_DISJ, SKYLINE_BDD_MAP, SKYLINE_BDD};
std::ostream & operator<<(std::ostream &os, const PruningDD & m);
extern const std::vector<std::string> PruningDDValues;

enum class PruningType {Expansion, Generation, None};
std::ostream & operator<<(std::ostream &os, const PruningType & m);
extern const std::vector<std::string> PruningTypeValues;

class DominancePruningSimulation : public PruneHeuristic {  
 protected:
  //Parameters to control the pruning
  const SymParamsMgr mgrParams; //Parameters for SymManager configuration.

  bool initialized;
  const bool remove_spurious_dominated_states;
  const bool insert_dominated;
  const PruningType pruning_type;

  /*
   * Three parameters help to decide whether to apply dominance
   * pruning or not. Dominance pruning is used until
   * min_insertions_desactivation are performed. At that moment, if
   * the ratio pruned/checked is lower than min_desactivation_ratio
   * the pruning is desactivated. If not, the pruning remains
   * activated until the planner finishes.
   */
  const int min_insertions_desactivation;
  const double min_desactivation_ratio;
  
  std::unique_ptr<SymVariables> vars; //The symbolic variables are declared here  
  std::unique_ptr<SymManager> mgr;    //The symbolic manager to handle mutex BDDs

  std::unique_ptr<LDSimulation> ldSimulation;

  bool all_desactivated;
  bool activation_checked;

  int states_inserted; //Count the number of states inserted
  int states_checked; //Count the number of states inserted
  int states_pruned; //Count the number of states pruned
  int deadends_pruned; //Count the number of dead ends detected

  void dump_options() const;

  /* Methods to help concrete classes */
  BDD getBDDToInsert(const State &state);

  //Methods to keep dominated states in explicit search
  //Check: returns true if a better or equal state is known
  virtual bool check (const State & state, int g) = 0;
  virtual void insert (const State & state, int g) = 0;

  /* //Methods to keep dominated states in symbolic search */
  /* virtual BDD check (const BDD & bdd, int g) = 0; */
  /* virtual void insert (const BDD & bdd, int g) = 0; */

  //void build_abstraction();

  inline bool is_activated() {
      if(!activation_checked && states_inserted > min_insertions_desactivation){
	  activation_checked = true;
	  all_desactivated = states_pruned == 0 || 
	      states_pruned < states_checked*min_desactivation_ratio;
	  std::cout << "Simulation pruning " << (all_desactivated ? "desactivated: " : "activated: ")
		    << states_pruned << " pruned " << states_checked << " checked " << 
	      states_inserted << " inserted " << deadends_pruned << " deadends " << std::endl;
      }
      return !all_desactivated;
  }


/*   inline bool insert_is_activated() const { */
/*       if (states_inserted > min_insertions && !activation_checked){ */
/* 	  all_activated = states_pruned >= states_inserted*min_insert_ratio; */
/*       } */
/*       return !all_activated; */
/*       //    return all_desactivated states_inserted < min_insertions ||  */
/*       //  states_pruned >= states_inserted*min_insert_ratio; */
/*   } */
/*   inline bool prune_is_activated() const { */
/*       return !all_activated; */
/*       //return states_inserted < min_insertions ||  */
/* //	  states_pruned >= states_inserted*min_pruning_ratio; */
/*   } */
/*   inline bool deadend_is_activated() const { */
/*       return !all_activated; */
/*       //return prune_is_activated(); */
/*       //|| states_inserted < min_deadends  */
/*       //  || deadends_pruned >= states_inserted*min_deadend_ratio; */
/*   } */

 public:
  virtual void initialize();

  //Methods for pruning explicit search
  virtual bool prune_generation(const State &state, int g);
  virtual bool prune_expansion (const State &state, int g);

  virtual bool is_dead_end(const State &state);

  virtual int compute_heuristic(const State &state);
  DominancePruningSimulation(const Options &opts);
  virtual ~DominancePruningSimulation();
};

class DominancePruningSimulationBDDMap : public DominancePruningSimulation {
  std::map<int, BDD> closed;
 public:
  DominancePruningSimulationBDDMap (const Options &opts) : 
  DominancePruningSimulation(opts)
  {}
  virtual ~DominancePruningSimulationBDDMap (){}

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g);
  virtual void insert (const State & state, int g);
};

class DominancePruningSimulationBDDMapDisj : public DominancePruningSimulation {
    std::map<int, std::vector<BDD> > closed;
 public:
  DominancePruningSimulationBDDMapDisj (const Options &opts) : 
  DominancePruningSimulation(opts)
  {}
  virtual ~DominancePruningSimulationBDDMapDisj (){}

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g);
  virtual void insert (const State & state, int g);
};


class DominancePruningSimulationBDD : public DominancePruningSimulation {
  BDD closed, closed_inserted;
  bool initialized;
 public:
  DominancePruningSimulationBDD (const Options &opts) : 
  DominancePruningSimulation(opts), initialized(false)
  {}
  virtual ~DominancePruningSimulationBDD (){}

  //Methods to keep dominated states in explicit search
  virtual bool check (const State & state, int g);
  virtual void insert (const State & state, int g);
};

class DominancePruningSimulationSkylineBDDMap : public DominancePruningSimulation {
    // We have the set of states inserted with each g-value that could
    // dominate each fluent
    std::vector<std::vector<std::map<int, BDD> > > closed;
    std::set<int> g_values;
public:
    DominancePruningSimulationSkylineBDDMap (const Options &opts) : 
    DominancePruningSimulation(opts)
    {}
    virtual ~DominancePruningSimulationSkylineBDDMap (){}
    
    //Methods to keep dominated states in explicit search
    virtual bool check (const State & state, int g);
    virtual void insert (const State & state, int g);
};

class DominancePruningSimulationSkylineBDD : public DominancePruningSimulation {
    std::vector<std::vector<BDD> > closed;
public:
    DominancePruningSimulationSkylineBDD (const Options &opts) : 
    DominancePruningSimulation(opts)
    {}
    virtual ~DominancePruningSimulationSkylineBDD (){}
    
    //Methods to keep dominated states in explicit search
    virtual bool check (const State & state, int g);
    virtual void insert (const State & state, int g);
};


/* class DominancePruningSimulationADD : public DominancePruningSimulation { */
/*   ADD closed; */
/*  public: */
/*   DominancePruningSimulationADD (const Options &opts,  */
/* 			  bool _insert_dominated) :  */
/*   DominancePruningSimulation(opts, _insert_dominated) */
/*   {} */
/*   virtual ~DominancePruningSimulationADD (){} */

/*     virtual bool check (const State & state, int g); */
/*     virtual void insert (const State & state, int g); */
/* }; */

#endif
