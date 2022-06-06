#ifndef __MLBP_FORMULATION_H__
#define __MLBP_FORMULATION_H__


#include "problems.h"
#include "mipsolver.h"

template<typename> struct Instance;
template<typename> struct Solution;

class MLBPFormulation : public MIPFormulation<MLBP>
{
public:
	virtual void createDecisionVariables(IloEnv env, const Instance<MLBP>& inst);
	virtual void addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst);
	virtual void addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst);
	virtual void extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol);
private:
	// binary decision variables x_{ijk}: item j at level i is inserted into bin k at level i + 1
	// DOES NOT include top level bins, so x.size = amount of levels
	IloArray<IloArray<IloNumVarArray>> x;

	// binary decision variables y_{ij}: bin j at level i is used (=1) or not (=0)
	// INCLUDES items at level 0, so y[0] is not used. y.size = amount of levels + 1
	IloArray<IloNumVarArray> y;

	// virtual void addBinCost(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol, int level, int bin);
};


#endif // __MLBP_FORMULATION_H__
