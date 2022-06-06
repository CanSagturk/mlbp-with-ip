#include "mlbpformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"


void MLBPFormulation::createDecisionVariables(IloEnv env, const Instance<MLBP>& inst)
{
	// decision variables x_{ijk}
	x = IloArray<IloArray<IloNumVarArray>>(env, inst.m);

	// decision variables y_{ij}
	y = IloArray<IloNumVarArray>(env, inst.m + 1);

	// x: For each item/bin with index j at each level with index i, create an IloNumVarArray with size of bins at level i + 1
	// y: For each item/bin with index i, create an IloNumVarArray with size of bins at level i
	int createdVariablesX = 0;
	int createdVariablesY = 0;
	for (int i = 0; i <= inst.m; i++) {
		if (i < inst.m) {
			x[i] = IloArray<IloNumVarArray>(env, inst.n[i]);
			for (int j : inst.B[i]) {
				x[i][j] = IloNumVarArray(env, inst.n[i + 1], 0, 1, ILOBOOL);
				createdVariablesX += inst.n[i + 1];
			}
		}
		y[i] = IloNumVarArray(env, inst.n[i], 0, 1, ILOBOOL);
		createdVariablesY += inst.n[i];
	}
	CS_OUT(TRACE) << "created " << createdVariablesX << " x_{ijk} variables" << std::endl;
	CS_OUT(TRACE) << "created " << createdVariablesY << " y_{ij} variables" << std::endl;
}

void MLBPFormulation::addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	// each item must be inserted into exactly one bin
	for (int j : inst.B[0]) {
		IloExpr sum(env);  // represents a linear expression of deicison variables and constants
		for (int k = 0; k < inst.n[1]; k++)
			sum += x[0][j][k];   // cplex overloads +,-,... operators
		model.add(sum == 1);  // add constraint to model
		sum.end();  // IloExpr must always call end() to free memory!
	}
	CS_OUT(TRACE) << "added " << inst.n[0] << " constraints to enforce the packing of each item" << std::endl;

	// each bin x must be inserted into exactly one bin y if x is being used
	int createdBinUsageConstraints = 0;
	for (int i = 1; i < inst.m; i++) {
		for (int j : inst.B[i]) {
			IloExpr sum(env);
			for (int k : inst.B[i + 1]) {
				sum += x[i][j][k];
			}
			model.add((y[i][j] == 1 && sum == 1) || (y[i][j] == 0 && sum == 0)); // y[i][j] == sum
			sum.end();
		}
		createdBinUsageConstraints += inst.n[i];
	}
	CS_OUT(TRACE) << "added " << createdBinUsageConstraints << " constraints to enforce the packing of each item" << std::endl;

	// the size of the content of a bin must not exceed the bin's capacity
	// This (might be) is inefficient. (This is not the case. Only assignning is inefficient, not the model) 
	// A better way to do it would be to create a data structure x for each level
	// Go through all the levels, sum the sizes of the items/bins put into a bin k at the higher level 
	// Insert the sum into index k in x.
	// While on next level, add the constraints using the sums from x
	int createdCapacityConstraints = 0;
	for (int i : inst.M) {
		for (int k : inst.B[i]) { // index of the bin the item/bin was put into
			IloExpr sum(env);
			for (int j : inst.B[i - 1]) { // index of the item/bin
				sum += x[i - 1][j][k] * inst.s[i - 1][j];
			}
			model.add(sum <= y[i][k] * inst.w[i][k]);
			sum.end();
		}
		createdCapacityConstraints += inst.n[i];
	}
	CS_OUT(TRACE) << "added " << createdCapacityConstraints << " capacity constraints" << std::endl;

//	for (int i = 0; i < inst.m; i++) {
//		for (int j : inst.B[i]) {
//			for (int k : inst.B[i + 1]) if (inst.w[i + 1][k] < inst.s[i][j]) {
//				model.add(x[i][j][k] == 0);
//			}
//		}
//	}

	// TODO: symmetry breaking constraints -> make sure that the j-th bin is used before the j+1-th bin is used
	//for (int i : inst.M) {
	//	for (int j : inst.B[i]) {
	//		for (int k : inst.)
	//	}
	//}

	// Symmetry Breaking: If an item/bin has the same weight with another one,
	// Then don't allocate the later without allocating the first
	//for (int i = 1; i < inst.m; i++) {
	//	for (int first : inst.B[i]) {
	//		for (int second = first + 1; second < inst.n[i]; second++) if (inst.s[i][first] == inst.s[i][second]) {
	//			IloExpr sumFirst(env);
	//			IloExpr sumSecond(env);
	//			for (int k : inst.B[i + 1]) {
	//				sumFirst += x[i][first][k];
	//				sumSecond += x[i][second][k];
	//			}
	//			model.add(sumFirst >= y[i][first] * sumSecond);
	//		}
	//	}
	//}
}

void MLBPFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	IloExpr sum(env);
	for (int i : inst.M) {
		for (int j : inst.B[i]) {
			sum += y[i][j] * inst.c[i][j];
		}
	}
	model.add(IloMinimize(env, sum));
	sum.end();
}

void MLBPFormulation::extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol)
{
	sol.cost = 0;
	for (int i = 0; i <= inst.m; i++) {
		if (i == inst.m) {
			for (int j : inst.B[i]) {
				if (cplex.getValue(y[i][j]) > 0.5) {
					sol.cost += inst.c[i][j];
				}
			}
		}
		else {
			sol.item_bin_indexes[i].assign(inst.n[i], -1);
			for (int j : inst.B[i]) {
				for (int k : inst.B[i + 1]) {
					if (cplex.getValue(x[i][j][k]) > 0.5) {
						sol.item_bin_indexes[i][j] = k;
						continue;
					}
				}
				if (i != 0 && cplex.getValue(y[i][j]) > 0.5) {
					sol.cost += inst.c[i][j];
				}
			}
		}
	}
}
