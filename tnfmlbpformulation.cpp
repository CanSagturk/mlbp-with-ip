#include "tnfmlbpformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"


void TNFMLBPFormulation::createDecisionVariables(IloEnv env, const Instance<MLBP>& inst)
{
	// decision variables x_{ijk}
	x = IloArray<IloArray<IloNumVarArray>>(env, inst.m);

	// flow variables f_{ijk}
	f = IloArray<IloArray<IloNumVarArray>>(env, inst.m);

	// decision variables y_{ij}
	y = IloArray<IloNumVarArray>(env, inst.m + 1);

	// x: For each item/bin with index j at each level with index i, create an IloNumVarArray with size of bins at level i + 1
	// y: For each item/bin with index i, create an IloNumVarArray with size of bins at level i
	int createdVariablesX = 0;
	int createdVariablesY = 0;
	for (int i = 0; i <= inst.m; i++) {
		if (i < inst.m) {
			x[i] = IloArray<IloNumVarArray>(env, inst.n[i]);
			f[i] = IloArray<IloNumVarArray>(env, inst.n[i]);
			for (int j : inst.B[i]) {
				x[i][j] = IloNumVarArray(env, inst.n[i + 1], 0, 1, ILOBOOL);
				//if (i == 0)
				//	f[i][j] = IloNumVarArray(env, inst.n[i + 1], 0, 1, ILOINT);
				//else
				f[i][j] = IloNumVarArray(env, inst.n[i + 1], 0, inst.n[0], ILOINT);
				createdVariablesX += inst.n[i + 1];
			}
		}
		y[i] = IloNumVarArray(env, inst.n[i], 0, 1, ILOBOOL);
		createdVariablesY += inst.n[i];
	}
	CS_OUT(TRACE) << "created " << createdVariablesX << " x_{ijk} and f_{ijk} variables" << std::endl;
	CS_OUT(TRACE) << "created " << createdVariablesY << " y_{ij} variables" << std::endl;
}

void TNFMLBPFormulation::addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	// All items provide one flow to the network
	// FROM SNF
	for (int j : inst.B[0]) {
		IloExpr sum(env);
		for (int k : inst.B[1]) {
			sum += f[0][j][k];
		}
		model.add(sum == 1);
		sum.end();
	}

	// All bins except the top level consume no flow
	// FROM SNF
	for (int i = 1; i < inst.m; i++) {
		for (int j : inst.B[i]) {
			IloExpr sum(env);
			for (int k : inst.B[i - 1]) {
				sum -= f[i - 1][k][j];
			}
			for (int k : inst.B[i + 1]) {
				sum += f[i][j][k];
			}
			model.add(sum == 0);
			sum.end();
		}
	}

	// Top level bins get flow equal to the amount of items
	// FROM SNF
	IloExpr topLevelSum(env);
	for (int j : inst.B[inst.m]) {
		for (int k : inst.B[inst.m - 1]) {
			topLevelSum += f[inst.m - 1][k][j];
		}
	}
	model.add(topLevelSum == inst.n[0]);
	topLevelSum.end();

	// TODO: Might want to combine this with the first constraint
	// If a bin is not used, then there is no flow
	//for (int i : inst.M) {
	//	for (int k : inst.B[i]) {
	//		IloExpr sum(env);
	//		for (int j : inst.B[i - 1]) {
	//			sum += f[i - 1][j][k];
	//			// model.add(IloIfThen(env, f[i - 1][j][k] >= 1, x[i - 1][j][k] == 1));
	//		}
	//		model.add(IloIfThen(env, y[i][k] == 0, sum == 0)); // try: ((y[i][k] == 0 && sum == 0) || y[i][k] == 1)
	//		sum.end();
	//	}
	//}

	// each item must be inserted into exactly one bin
	// FROM STANDARD MLBP
	for (int j : inst.B[0]) {
		IloExpr sum(env);  // represents a linear expression of deicison variables and constants
		for (int k = 0; k < inst.n[1]; k++)
			sum += x[0][j][k];   // cplex overloads +,-,... operators
		model.add(sum == 1);  // add constraint to model
		sum.end();  // IloExpr must always call end() to free memory!
	}

	// each bin x must be inserted into exactly one bin y if x is being used
	// FROM STANDARD MLBP
	for (int i = 1; i < inst.m; i++) {
		for (int j : inst.B[i]) {
			IloExpr sum(env);
			for (int k : inst.B[i + 1]) {
				sum += x[i][j][k];
			}
			model.add((y[i][j] == 1 && sum == 1) || (y[i][j] == 0 && sum == 0)); // y[i][j] == sum
			sum.end();
		}
	}

	// the size of the content of a bin must not exceed the bin's capacity
	// FROM STANDARD MLBP
	for (int i : inst.M) {
		for (int k : inst.B[i]) { // index of the bin the item/bin was put into
			IloExpr sum(env);
			for (int j : inst.B[i - 1]) { // index of the item/bin
				sum += x[i - 1][j][k] * inst.s[i - 1][j];
			}
			model.add(sum <= y[i][k] * inst.w[i][k]);
			sum.end();
		}
	}

	// NEW
	// Makes it slower
	//for (int i = 0; i < inst.m; i++) {
	//	for (int j : inst.B[i]) {
	//		for (int k : inst.B[i + 1]) {
	//			model.add(IloIfThen(env, x[i][j][k] == 0, f[i][j][k] == 0));
	//			model.add(IloIfThen(env, f[i][j][k] == 0, x[i][j][k] == 0));

	//			model.add(f[i][j][k] <= x[i][j][k] * inst.n[0]);
	//		}
	//	}
	//}
}

void TNFMLBPFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst)
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

void TNFMLBPFormulation::extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol)
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
					if (cplex.getValue(x[i][j][k]) > 0.5) { // CHANGE x TO f TO SEE THAT NF IS NOT USED
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
