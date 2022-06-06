#ifndef __SOLUTION_VERIFIER_H__
#define __SOLUTION_VERIFIER_H__


#include <vector>
#include <string>

#include "problems.h"
#include "solution.h"


template<typename>
struct Instance;

/*
 * Generic Solution Verifier.
 * To write a solution verifier class for a specific problem,
 * speialize this class for the corresponding problem type.
 */
template<typename ProbT>
struct SolutionVerifier
{
	static bool verify(const Instance<ProbT>& inst, const Solution<ProbT>& sol, std::vector<std::string>* error_msg = nullptr);
};


/*****************************************************************************************/
/** Bin Packing Problem ******************************************************************/
/*****************************************************************************************/
template<>
struct SolutionVerifier<BP>
{
	static bool verify(const Instance<BP>& inst, const Solution<BP>& sol, std::vector<std::string>* error_msg = nullptr);
};


/*****************************************************************************************/
/** Multi-Level Bin Packing Problem ******************************************************/
/*****************************************************************************************/
template<>
struct SolutionVerifier<MLBP>
{
	static bool verify(const Instance<MLBP>& inst, const Solution<MLBP>& sol, std::vector<std::string>* error_msg = nullptr);
	static void streamErrorMessage(std::stringstream& messageStream, std::vector<std::string>* errorStream);
};


/*****************************************************************************************/
/** Class Constrained Multi-Level Bin Packing Problem ************************************/
/*****************************************************************************************/
template<>
struct SolutionVerifier<CCMLBP>
{
	static bool verify(const Instance<CCMLBP>& inst, const Solution<CCMLBP>& sol, std::vector<std::string>* error_msg = nullptr);
};


/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Conflict Constraints ****************************/
/*****************************************************************************************/
template<>
struct SolutionVerifier<MLBPCC>
{
	static bool verify(const Instance<MLBPCC>& inst, const Solution<MLBPCC>& sol, std::vector<std::string>* error_msg = nullptr);
};


/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Partial Orders **********************************/
/*****************************************************************************************/
template<>
struct SolutionVerifier<MLBPPO>
{
	static bool verify(const Instance<MLBPPO>& inst, const Solution<MLBPPO>& sol, std::vector<std::string>* error_msg = nullptr);
	static int topLevelFinder(int itemIndex, std::vector<std::vector<int>> item_bin_indexes);
};


/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Time Windows ************************************/
/*****************************************************************************************/
template<>
struct SolutionVerifier<MLBPTW>
{
	static bool verify(const Instance<MLBPTW>& inst, const Solution<MLBPTW>& sol, std::vector<std::string>* error_msg = nullptr);
};


/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Fragmentation Constraints ***********************/
/*****************************************************************************************/
template<>
struct SolutionVerifier<MLBPFC>
{
	static bool verify(const Instance<MLBPFC>& inst, const Solution<MLBPFC>& sol, std::vector<std::string>* error_msg = nullptr);
};

#endif // __SOLUTION_VERIFIER_H__
