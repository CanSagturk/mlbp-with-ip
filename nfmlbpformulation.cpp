#include "nfmlbpformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"


void NFMLBPFormulation::createDecisionVariables(IloEnv env, const Instance<MLBP>& inst)
{
	for (int i = 0; i <= inst.m; i++) {
		if (i == inst.m) {
			edgeCount += inst.n[i];
		}
		else {
			edgeCount += inst.n[i] * inst.n[i + 1];
		}
		nodeCount += inst.n[i];
	}
	nodeCount++;

	x = IloNumVarArray(env, edgeCount, 0, 1, ILOBOOL);
	f = IloNumVarArray(env, edgeCount, 0, inst.n[0], ILOINT);
	origin.assign(edgeCount, -1);
	destination.assign(edgeCount, -1);

	size.assign(edgeCount, 0);
	capacity.assign(nodeCount, 0);
	cost.assign(edgeCount, 0);

	demand.assign(nodeCount, 0);

	int nodeIndex = 0;
	int edgeIndex = 0;

	for (int i = inst.m + 1; i > 0; i--) {
		if (i == inst.m + 1) {
			for (int k : inst.B[i - 1]) {
				origin[edgeIndex] = 0;
				destination[edgeIndex] = k + 1;
				cost[edgeIndex] = inst.c[inst.m][k];
				size[edgeIndex] = 0;
				edgeIndex++;
			}
			nodeIndex++;
			demand[0] = inst.n[0];
			capacity[0] = 0;
		}
		else {
			int nextNodeIndex = nodeIndex + inst.n[i];
			for (int j : inst.B[i]) {
				int currentNodeIndex = nodeIndex + j;

				demand[currentNodeIndex] = 0;
				capacity[currentNodeIndex] = inst.w[i][j];

				for (int k : inst.B[i - 1]) {
					origin[edgeIndex] = currentNodeIndex;
					destination[edgeIndex] = nextNodeIndex + k;
					size[edgeIndex] = inst.s[i - 1][k];
					if (i == 1)
						cost[edgeIndex] = 1;
					else
						cost[edgeIndex] = inst.c[i - 1][k];
					edgeIndex++;
				}
			}
			nodeIndex = nextNodeIndex;
		}
	}

	for (int i : inst.B[0]) {
		demand[nodeIndex + i] = -1;
		capacity[nodeIndex + i] = 0;
	}

	// int asd = inst.c[14323][45];
	CS_OUT(TRACE) << "created " << edgeIndex << " edges" << std::endl;
}

void NFMLBPFormulation::addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	for (IloInt i = 0; i < nodeCount; ++i) {
		IloNumExpr sum(env);
		for (IloInt j = 0; j < edgeCount; ++j) {
			if (origin[j] == i)
				sum += f[j];
			if (destination[j] == i)
				sum -= f[j];
		}
		model.add(sum >= demand[i]); // maybe == ?
		sum.end();
	}

	for (IloInt i = 0; i < edgeCount; ++i) {
		model.add(IloIfThen(env, x[i] == 0, f[i] == 0));
	}

	// Capacity Contraint
	for (int i = 0; i < nodeCount; ++i) {
		IloNumExpr sum(env);
		for (int j = 0; j < edgeCount; ++j) {
			if (origin[j] == i)
				sum += x[j] * size[j];
		}
		model.add(sum <= capacity[i]);
		sum.end();
	}
}

void NFMLBPFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	IloExpr sum(env);
	for (int i = 0; i < edgeCount; ++i) {
		sum += x[i] * cost[i];
	}
	model.add(IloMinimize(env, sum));
	sum.end();
}

void NFMLBPFormulation::extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol)
{
	sol.cost = 0;
	// sol.item_bin_indexes = { {1, 2}, {2} };

	for (IloInt e = 0; e < edgeCount; ++e) {
		if (cplex.getValue(x[e]) > 0.5) {
			sol.cost += cost[e];
		}
	}
	sol.cost -= inst.n[0];

	int nodeIndex = nodeCount - 1;

	for (int i = 0; i < inst.m; i++) {
		sol.item_bin_indexes[i].assign(inst.n[i], -1);
		for (int j : inst.B[i]) {
			bool edgeFound = false;
			for (int e = 0; e < edgeCount; e++) {
				if (destination[e] == nodeIndex) {
					CS_OUT(TRACE) << "Origin: " << origin[e] << ", Destination: " << destination[e] << ", Flow: " << cplex.getValue(f[e]) 
						<< ", Is Used: " << cplex.getValue(x[e]) << std::endl;
					if (cplex.getValue(x[e]) > 0.5) {
						sol.item_bin_indexes[i][j] = indexFinder(origin[e], inst); // indexFinder(origin[e], inst)
					}
					


					//if (edgeFound)
					//	sol.item_bin_indexes[i][j] = -5;
					//else {
					//	sol.item_bin_indexes[i][j] = origin[e];
					//	edgeFound = true;
					//}
				}
			}
			nodeIndex--;
		}
	}
}

int NFMLBPFormulation::indexFinder(int node, const Instance<MLBP>& inst) {
	int levelIterator = 1;

	for (int i = inst.m; levelIterator + inst.n[i] <= node; i--) {
		if (i < 0) {
			return INT16_MIN;
		}
		levelIterator += inst.n[i];
	}

	return (node - levelIterator);
}
