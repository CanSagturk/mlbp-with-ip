#ifndef __BP_FORMULATION_H__
#define __BP_FORMULATION_H__

#include "problems.h"
#include "mipsolver.h"

template<typename> struct Instance;
template<typename> struct Solution;

/**
 * MIP Formulation for the Bin Packing Problem.
 * Uses two kinds of binary decision variables:
 *
 *  *) x_{ij} -> item i is inserted into bin j
 *  *) y_j    -> bin j is used
 *
 *  Since we need for each possible used bin a decision variable,
 *  this formulation needs an upper bound on the number of bins
 *  in order to create enough decision variables.
 *
 *  We use a naive upper bound equal to the number of items.
 */
class BPFormulation : public MIPFormulation<BP>
{
public:
	// create all required decision variables
	virtual void createDecisionVariables(IloEnv env, const Instance<BP>& inst) override;

	// add all constraints to the model
	virtual void addConstraints(IloEnv env, IloModel model, const Instance<BP>& inst) override;

	// add objective function to the model
	virtual void addObjectiveFunction(IloEnv env, IloModel model, const Instance<BP>& inst) override;

	// derive solution from the cplex object
	virtual void extractSolution(IloCplex cplex, const Instance<BP>& inst, Solution<BP>& sol) override;

private:
	// binary decision variables x_{ij}: item i is inserted into bin j (=1) or not (=0)
	IloArray<IloNumVarArray> x;

	// binary decision variables y_j: bin j is used (=1) or not (=0)
	IloNumVarArray y;
};


#endif // __BP_FORMULATION_H__
