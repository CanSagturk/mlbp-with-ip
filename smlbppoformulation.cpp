#include "smlbppoformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"


void SMLBPPOFormulation::createDecisionVariables(IloEnv env, const Instance<MLBPPO>& inst)
{
	// decision variables x_{ijk}
	x = IloArray<IloArray<IloNumVarArray>>(env, inst.m);

	// decision variables y_{ij}
	y = IloArray<IloNumVarArray>(env, inst.m + 1);

	p = IloArray<IloNumVarArray>(env, inst.n[0]);

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

	for (int item : inst.B[0]) {
		int maxBins = 0;
		for (int levelSize : inst.n) if (levelSize > maxBins)
			maxBins = levelSize;
		p[item] = IloNumVarArray(env, inst.m, 0, maxBins, ILOINT);
	}
}

void SMLBPPOFormulation::addConstraints(IloEnv env, IloModel model, const Instance<MLBPPO>& inst)
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

	// Tracing the last level item
	for (int item : inst.B[0]) {
		for (int bin : inst.B[1]) {
			model.add(IloIfThen(env, x[0][item][bin] == 1, p[item][0] == bin));
		}
	}

	// Tracing the item
	for (int item : inst.B[0]) {
		for (int i = 1; i < inst.m; i++) {
			for (int j : inst.B[i]) {
				for (int k : inst.B[i + 1]) {
					model.add(IloIfThen(env, p[item][i - 1] == j && x[i][j][k] == 1, p[item][i] == k));
				}
			}
		}
	}

	//for (int item : inst.B[0]) {
	//	for (int level = 0; level < inst.m; level++) {
	//		IloExpr sum(env);
	//		for (int bin : inst.B[level + 1]) {
	//			sum += p[item][level][bin];
	//		}
	//		model.add(sum == 1);
	//		sum.end();
	//	}
	//}

	// Precedence constraints
	for (std::pair<int, int> pair : inst.pos) {
		for (int topLevelIndex : inst.B[inst.m]) {
			model.add(p[pair.first][inst.m - 1] <= p[pair.second][inst.m - 1]);
		}
	}
}

void SMLBPPOFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBPPO>& inst)
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

void SMLBPPOFormulation::extractSolution(IloCplex cplex, const Instance<MLBPPO>& inst, Solution<MLBPPO>& sol)
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

	for (int item : inst.B[0]) {
		for (int level = 0; level < inst.m; level++) {
			for (int bin : inst.B[level + 1]) {
				sol.pActual[item][level] = cplex.getValue(p[item][level]);
			}
		}
	}
}
