#ifndef SYM_CLOSED_H
#define SYM_CLOSED_H

#include "../debug.h"
#include "sym_manager.h"
#include "sym_abstraction.h"
#include <vector>
#include <set>
#include <map>


class SymSolution;
class Timer;
class SymClosed;

//Auxiliar class to denote an heuristic evaluation
class Evaluation {
public:
    SymExploration * exp;
    std::vector<BDD> bucket; // States to be evaluated
    int f, h; // h and f values needed to prune states

    Evaluation (const Evaluation & o) = default;
Evaluation(SymExploration * exploration, int fval, int hval) :
    exp(exploration), f(fval), h(hval){}


    // Update eval removing all states that cannot be
    // pruned with this search as heuristic
    void updateClosed (const BDD & S, int hClosed);
    
    void evaluate (SymClosed * c) const ;
};


class SymClosed /*: public SymHeuristic */ {  
    friend class Evaluation; //For using evaluate_abs_orig
private:
    SymManager * mgr; //Symbolic manager to perform bdd operations
    SymExploration * exploration;

    std::map<int, BDD> closed;   // Mapping from cost to set of states

    // Auxiliar BDDs for the number of 0-cost action steps
    // ALERT: The information here might be wrong
    // It is just used to extract path more quickly, but the information
    // here is an (admissible) estimation and this should be taken into account
    std::map<int, vector<BDD>> zeroCostClosed; 
    BDD closedTotal;             // All closed states.

    int hNotClosed, fNotClosed; // Bounds on h and g for those states not in closed
    std::map<int, BDD> closedUpTo;  // Disjunction of BDDs in closed  (auxiliar useful to take the maximum between several BDDs)
    std::set<int> h_values; //Set of h_values of the heuristic

    SymClosed * parent;
    // children: Related heuristics derived from relaxations of this
    // search. We need this because we will only use my children to
    // evaluate a state if I do not have a value (otherwise could be
    // non-admissible due to the initialization of the closed list).
    std::vector<SymClosed *> children; 


    //std::map<SymExploration *, Evaluation> evals; 
    // For now, searches only inform the abstract search, so it is only
    // needed to store a single Evaluation
    std::unique_ptr <Evaluation> evalOrig; 
    std::unique_ptr <Evaluation> evalOpposite; 

    // Statistics of heuristic evaluation
    //int num_calls_eval;
    //double time_eval_states, time_closed_states, time_pruned_states, time_prune_some, time_prune_all, time_prune_some_children;

    BDD evaluate_orig_orig(const BDD & bdd, int fVal, int hVal, SymExploration * expAsking, bool store_eval);
    BDD evaluate_abs_abs(const BDD & bdd, int fVal, int hVal, SymExploration * expAsking, bool store_eval);
    BDD evaluate_abs_orig(const BDD & bdd, int fVal, int hVal, SymExploration * expAsking, bool store_eval);
    BDD evaluate(const BDD & bdd, int hVal);

    void cleanEvalOrig();

    void checkEval(Evaluation * eval);

public:
    SymClosed();
    void init(SymExploration * exp, SymManager * manager); 
    void init(SymExploration * exp, SymManager * manager, const SymClosed & other);

    void insert (int h, const BDD & S); 
    void setHNotClosed(int h);
    void setFNotClosed(int f);
    void newHValue(int h_value);


    bool can_prune(int fVal, int hVal);
    BDD evaluate(const BDD & bdd, int fVal, int hVal, SymExploration * expAsking, bool store_eval);
 

    bool accept(int f, int h) const;
    void getNextF(int f, pair<int, int> & upper_bound) const;
    const std::set<int> & getHValues(); 
    
    void cleanEvals(bool orig) {
	if (orig) cleanEvalOrig();
	evalOpposite.reset();
    }


    void setOrigF(int f, int h);
    void setOppositeF(int f, int h);

    void getUsefulExps(std::vector<SymExploration *> & useful_exps) const;
    bool isUsefulAfterRelax(double ratio, 
			    const std::vector<BDD> & newFrontier) const; 
    
    bool isUseful(double ratio); 

    /* bool isUseful(const SymExploration & frontierOriginal,  */
    /* 		const std::vector<BDD> & frontierAbstract, */
    /* 		double ratio) { */
    /*     double rUseful = frontierOriginal.ratioUseful(frontierAbstract); */
    /*     return rUseful > 0 && rUseful >= ratio  ; */
    /* } */
  

    //Check if any of the states is closed.
    //In case positive, return a solution pair <f_value, S>
    SymSolution checkCut(const BDD & states, int g, bool fw) const;

    void extract_path(const BDD & cut, int h, bool fw,
		      std::vector <const Operator *> & path) const;

    inline BDD getClosed() const{ 
	return closedTotal;
    }

    inline BDD notClosed() const{ 
	return !closedTotal;
    }

    inline std::map<int, BDD> getClosedList() const{
	return closed;
    }


    inline int getHNotClosed() const{
	return hNotClosed;
    }

    inline int getFNotClosed() const{
	return fNotClosed;
    }

    inline SymExploration * getExploration(){
	return exploration;
    }

    void addChild(SymClosed * c); 

    void desactivate(){
	if(parent)
	    parent->children.erase( std::remove( std::begin(parent->children),
					     std::end(parent->children), this),
				    std::end(parent->children) );
	evalOrig = nullptr;
    }

     bool hasEvalOrig() const {
	 return evalOrig.get() != nullptr;
    }

    bool hasChildren() const{
	return !children.empty();
    }

    void write(const std::string & fname, std::ofstream & file) const;
    void init(SymExploration * exp, SymManager * manager, const string & fname, ifstream & file); 

    ADD getHeuristic(int previousMaxH = -1) const;
    void getHeuristic(std::vector<ADD> & heuristics,
		      std::vector <int> & maxHeuristicValues) const;

    void statistics () const; 

    friend std::ostream & operator<<(std::ostream &os, const SymClosed & c);
};

#endif // SYM_CLOSED


//  void extract_path(const State &state, std::vector <Operator *> & path);
/*  
    void storeHeuristic(const string & filename) const ;*/

//Inherited from Heuristic
//int getH (const BDD & bdd) const; 
//std::set<int> getHValues() const; 
//BDD getBDD(const BDD & bdd, int fVal, int hVal) const;
//BDD getBDDUpTo(const BDD & bdd, int fVal, int hVal);
//int getH(const State & state);

//virtual void notifyToExploration(SymExploration * exploration);

