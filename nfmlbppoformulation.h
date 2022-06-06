#ifndef __NFMLBPPO_FORMULATION_H__
#define __NFMLBPPO_FORMULATION_H__


#include "problems.h"
#include "mipsolver.h"

template<typename> struct Instance;
template<typename> struct Solution;

class NFMLBPPOFormulation : public MIPFormulation<MLBPPO>
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

	// flow types f_{ijkh}: item/bin at level i provides types of flow that is equal to the amount of items
	// to a bin k at level i + 1, if the flow type is h, then h = 1, 0 otherwise
	// Includes top level bins, which send flow to the sink, so f.size = amount of levels + 1
	IloArray<IloArray<IloArray<IloNumVarArray>>> f;
};


#endif
