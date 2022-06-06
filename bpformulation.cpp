#include "bpformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"


void BPFormulation::createDecisionVariables(IloEnv env, const Instance<BP>& inst)
{
	// decision variables x_{ij}
	x = IloArray<IloNumVarArray>(env, inst.n);  // create array of decision variables (inst.n elements)
	for (int i : inst.I)
		x[i] = IloNumVarArray(env, inst.n, 0, 1, ILOBOOL);
	//                               ^     ^  ^     ^
	//                               |     |  |     +------ variable type: ILOINT for integers
	//                               |     |  +------------ maximum value of variables
	//                               |     +--------------- minium value of variables
	//                               +--------------------- number of variables
	MIP_OUT(TRACE) << "created " << inst.n * inst.n << " x_{ij} variables" << std::endl;

	// decision variables y_j
	y = IloNumVarArray(env, inst.n, 0, 1, ILOBOOL);
	MIP_OUT(TRACE) << "created " << inst.n << " y_{i} variables" << std::endl;
}

void BPFormulation::addConstraints(IloEnv env, IloModel model, const Instance<BP>& inst)
{
	// each item must be inserted into exactly one bin
	for (int i : inst.I) {
		IloExpr sum(env);  // represents a linear expression of deicison variables and constants
		for (int j = 0; j < inst.n; j++)
			sum += x[i][j];   // cplex overloads +,-,... operators
		model.add(sum == 1);  // add constraint to model
		sum.end();  // IloExpr must always call end() to free memory!
	}
	MIP_OUT(TRACE) << "added " << inst.n << " constraints to enforce the packing of each item" << std::endl;

	// the size of the content of a bin must not exceed the bin's capacity
	for (int j = 0; j < inst.n; j++) {
		IloExpr sum(env);
		for (int i : inst.I)
			sum += x[i][j] * inst.s[i];
		model.add(sum <= y[j] * inst.smax);
		sum.end();
	}
	MIP_OUT(TRACE) << "added " << inst.n << " capacity constraints" << std::endl;

	// symmetry breaking constraints -> make sure that the j-th bin is used before the j+1-th bin is used
	for (int j = 0; j < inst.n-1; j++)
		model.add(y[j] >= y[j+1]);
	MIP_OUT(TRACE) << "added " << inst.n-1 << " symmetry breaking constraints" << std::endl;
}

void BPFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<BP>& inst)
{
	IloExpr sum(env);
	for (int j = 0; j < inst.n; j++)
		sum += y[j];
	model.add(IloMinimize(env, sum));
	sum.end();
}

void BPFormulation::extractSolution(IloCplex cplex, const Instance<BP>& inst, Solution<BP>& sol)
{
	// cplex.getValue(x) returns the assigned value of decision variable x
	sol.item_to_bins.assign(inst.n, -1);
	sol.total_bins = 0;
	for (int j = 0; j < inst.n; j++) {
		if (cplex.getValue(y[j]) < 0.5)
			break;  // valid due to symmetry breaking constraints

		sol.total_bins++;
		for (int i : inst.I)
			if (cplex.getValue(x[i][j]) > 0.5)
				sol.item_to_bins[i] = j;
	}
}

