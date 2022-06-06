#ifndef __MLBPPO_FORMULATION_H__
#define __MLBPPO_FORMULATION_H__


#include "problems.h"
#include "mipsolver.h"

template<typename> struct Instance;
template<typename> struct Solution;

class MLBPPOFormulation : public MIPFormulation<MLBPPO>
{
public:
	virtual void createDecisionVariables(IloEnv env, const Instance<MLBPPO>& inst);
	virtual void addConstraints(IloEnv env, IloModel model, const Instance<MLBPPO>& inst);
	virtual void addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBPPO>& inst);
	virtual void extractSolution(IloCplex cplex, const Instance<MLBPPO>& inst, Solution<MLBPPO>& sol);
private:
	// binary decision variables x_{ijk}: item j at level i is inserted into bin k at level i + 1
	// DOES NOT include top level bins, so x.size = amount of levels
	IloArray<IloArray<IloNumVarArray>> x;

	// binary decision variables y_{ij}: bin j at level i is used (=1) or not (=0)
	// INCLUDES items at level 0, so y[0] is not used. y.size = amount of levels + 1
	IloArray<IloNumVarArray> y;

	// binary decision variables p_{ijk}: if item i goes to bin k at level j + 1, then it is 1
	// DOES NOT include top level bins, so p.size = amount of levels
	IloArray<IloArray<IloNumVarArray>> p;
};


#endif
