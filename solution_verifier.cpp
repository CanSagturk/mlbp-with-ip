#include "solution_verifier.h"
#include "instance.h"
#include "solution.h"

#include <sstream>
#include <algorithm>
#include <set>



/*************************************************************************************************/
/* BP ********************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<BP>::verify(const Instance<BP>& inst, const Solution<BP>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	std::vector<int> levels(sol.total_bins, 0);  // size of each bin
	for (int i : inst.I) {
		int bin = sol.item_to_bins[i];
		if (bin >= sol.total_bins) {
			if (error_msg) {
				std::stringstream ss;
				ss << "Detected bin with id " << bin << " but solutions claims that there are only " << sol.total_bins << ".";
				error_msg->push_back(ss.str());
			}
			ret = false;
		} else if (bin < 0) {
			if (error_msg) {
				std::stringstream ss;
				ss << "Item " << i << " is not assigned to any bin.";
				error_msg->push_back(ss.str());
			}
			ret = false;
		} else {
			levels[bin] += inst.s[i];
		}
	}
	for (int bin = 0; bin < sol.total_bins; bin++) {
		if (levels[bin] > inst.smax) {
			if (error_msg) {
				std::stringstream ss;
				ss << "Bin " << bin << " of size " << levels[bin] << " exceeds maximum capaicty (" << inst.smax << ").";
				error_msg->push_back(ss.str());
			}
			ret = false;
		}
	}
	return ret;
}


/*************************************************************************************************/
/* MLBP ******************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBP>::verify(const Instance<MLBP>& inst, const Solution<MLBP>& sol, std::vector<std::string>* error_msg)
{
	std::vector<std::vector<int>> levels(inst.m + 1); // +1 since levels vector contains the unused items level as well

	for (int i : inst.M) { // Level 0 is not used since it represents items and not bins
		levels[i].assign(inst.n[i], 0);
	} 

	for (int i = 0; i < inst.m; i++) {
		for (int j : inst.B[i]) {
			if (sol.item_bin_indexes[i][j] > inst.n[i + 1]) {
				std::stringstream ss;
				ss << "Detected bin with id " << sol.item_bin_indexes[i][j] << " but solutions claims that there are only " << inst.n[i + 1] << ".";
				streamErrorMessage(ss, error_msg);
				return false;
			}
			// Check if all items are placed for first level only
			if (i == 0) {
				if (sol.item_bin_indexes[i][j] < 0) {
					std::stringstream ss;
					ss << "Item " << j << " is not assigned to any bin.";
					streamErrorMessage(ss, error_msg);
					return false;
				}
			}
			else if (levels[i][j] > 0 && sol.item_bin_indexes[i][j] < 0) {
				std::stringstream ss;
				ss << "Bin " << j << " at level " << i << " is not assigned to any bin despite being used.";
				streamErrorMessage(ss, error_msg);
				return false;
			}
			int itemBinIndex = sol.item_bin_indexes[i][j];
			if (itemBinIndex >= 0) {
				levels[i + 1][itemBinIndex] += inst.s[i][j];
			}
		}
	}

	for (int i : inst.M) {
		for (int j : inst.B[i]) {
			if (levels[i][j] > inst.w[i][j]) {
				std::stringstream ss;
				ss << "Bin " << j << " at level " << i << " is over its capacity limit." << std::endl;
				ss << "Capacity limit: " << inst.w[i][j] << ", Loaded " << ((i == 0) ? "item" : "bin") << " weight : " << levels[i][j];
				streamErrorMessage(ss, error_msg);
				return false;
			}
		}
	}

	return true;
}

void SolutionVerifier<MLBP>::streamErrorMessage(std::stringstream& messageStream, std::vector<std::string>* errorStream) {
	if (errorStream) {
		errorStream->push_back(messageStream.str());
	}
}





/*************************************************************************************************/
/* CCMLBP ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<CCMLBP>::verify(const Instance<CCMLBP>& inst, const Solution<CCMLBP>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	// TODO
	return ret;
}




/*************************************************************************************************/
/* MLBPCC ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPCC>::verify(const Instance<MLBPCC>& inst, const Solution<MLBPCC>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	// TODO
	return ret;
}




/*************************************************************************************************/
/* MLBPPO ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPPO>::verify(const Instance<MLBPPO>& inst, const Solution<MLBPPO>& sol, std::vector<std::string>* error_msg)
{
	std::vector<std::vector<int>> levels(inst.m + 1); // +1 since levels vector contains the unused items level as well

	for (int i : inst.M) { // Level 0 is not used since it represents items and not bins
		levels[i].assign(inst.n[i], 0);
	}

	for (int i = 0; i < inst.m; i++) {
		for (int j : inst.B[i]) {
			if (sol.item_bin_indexes[i][j] > inst.n[i + 1]) {
				std::stringstream ss;
				ss << "Detected bin with id " << sol.item_bin_indexes[i][j] << " but solutions claims that there are only " << inst.n[i + 1] << ".";
				SolutionVerifier<MLBP>::streamErrorMessage(ss, error_msg);
				return false;
			}
			// Check if all items are placed for first level only
			if (i == 0) {
				if (sol.item_bin_indexes[i][j] < 0) {
					std::stringstream ss;
					ss << "Item " << j << " is not assigned to any bin.";
					SolutionVerifier<MLBP>::streamErrorMessage(ss, error_msg);
					return false;
				}
			}
			else if (levels[i][j] > 0 && sol.item_bin_indexes[i][j] < 0) {
				std::stringstream ss;
				ss << "Bin " << j << " at level " << i << " is not assigned to any bin despite being used.";
				SolutionVerifier<MLBP>::streamErrorMessage(ss, error_msg);
				return false;
			}
			int itemBinIndex = sol.item_bin_indexes[i][j];
			if (itemBinIndex >= 0) {
				levels[i + 1][itemBinIndex] += inst.s[i][j];
			}
		}
	}

	for (int i : inst.M) {
		for (int j : inst.B[i]) {
			if (levels[i][j] > inst.w[i][j]) {
				std::stringstream ss;
				ss << "Bin " << j << " at level " << i << " is over its capacity limit." << std::endl;
				ss << "Capacity limit: " << inst.w[i][j] << ", Loaded " << ((i == 0) ? "item" : "bin") << " weight : " << levels[i][j];
				SolutionVerifier<MLBP>::streamErrorMessage(ss, error_msg);
				return false;
			}
		}
	}

	for (std::pair<int, int> pair : inst.pos) {
		int firstTop = topLevelFinder(pair.first, sol.item_bin_indexes);
		int secondTop = topLevelFinder(pair.second, sol.item_bin_indexes);
		if (firstTop > secondTop) {
			std::stringstream ss;
			ss << "Item " << pair.second << " cannot be before " << pair.first << ". Item " << pair.first << " bin: " << firstTop;
			ss << ", Item " << pair.second << " bin: " << secondTop << std::endl;
			SolutionVerifier<MLBP>::streamErrorMessage(ss, error_msg);
			return false;
		}
	}

	return true;
}
int SolutionVerifier<MLBPPO>::topLevelFinder(int itemIndex, std::vector<std::vector<int>> item_bin_indexes) 
{
	int index = itemIndex;
	for (std::vector<int> level : item_bin_indexes) {
		index = level[index];
	}
	return index;
}




/*************************************************************************************************/
/* MLBPTW ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPTW>::verify(const Instance<MLBPTW>& inst, const Solution<MLBPTW>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	// TODO
	return ret;
}




/*************************************************************************************************/
/* MLBPFC ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPFC>::verify(const Instance<MLBPFC>& inst, const Solution<MLBPFC>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	// TODO
	return ret;
}
