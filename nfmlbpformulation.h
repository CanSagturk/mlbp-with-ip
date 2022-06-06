#ifndef __NF_MLBP_FORMULATION_H__
#define __NF_MLBP_FORMULATION_H__

#include <vector>

#include "problems.h"
#include "mipsolver.h"

template<typename> struct Instance;
template<typename> struct Solution;

class NFMLBPFormulation : public MIPFormulation<MLBP>
{
public:
	virtual void createDecisionVariables(IloEnv env, const Instance<MLBP>& inst);
	virtual void addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst);
	virtual void addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst);
	virtual void extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol);
private:
	virtual int indexFinder(int node, const Instance<MLBP>& inst);

	// binary decision variables x_{i], where edge i is 1 if used, 0 otherwise
	IloNumVarArray x;

	// integer array f_{i}, represents flow going through edge i
	IloNumVarArray f;

	// Starting node of all edges
	std::vector<int> origin;

	// Ending node of all edges
	std::vector<int> destination;

	// Size added by all edges
	std::vector<int> size;

	// Cost of all edges
	std::vector<int> cost;


	// Demand for all nodes. For the source, this is the amount of items. For items, -1. All else is 0
	std::vector<int> demand;

	// Capacity of all nodes
	std::vector<int> capacity;

	int nodeCount;
	int edgeCount;
};


#endif // __NF_MLBP_FORMULATION_H__
