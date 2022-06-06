#include "solution.h"

#include "lib/util.h"
#include "users.h"

#include <algorithm>
#include <numeric>
#include <cassert>



/*****************************************************************************************/
/** Bin Packing Problem ******************************************************************/
/*****************************************************************************************/
Solution<BP>::Solution(const Instance<BP>& inst) : item_to_bins(inst.s.size())
{
	std::iota(item_to_bins.begin(), item_to_bins.end(), 0);
	total_bins = (int)inst.s.size();
}

std::ostream& operator<<(std::ostream& os, const Solution<BP>& sol)
{
	os << sol.item_to_bins;
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem ******************************************************/
/*****************************************************************************************/
Solution<MLBP>::Solution(const Instance<MLBP>& inst) : db(-1) // item_bin_indexes(inst.m), 
{
	for (int i = 0; i < inst.m; i++) {
		std::vector<int> next(inst.s[i]);
		std::iota(next.begin(), next.end(), 0);
		item_bin_indexes.push_back(next);
	}
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBP>& sol) {
	for (std::vector<int> arr : sol.item_bin_indexes) {
		os << arr;
	}
	return os;
}




/*****************************************************************************************/
/** Class Constrained Multi-Level Bin Packing Problem ************************************/
/*****************************************************************************************/
Solution<CCMLBP>::Solution(const Instance<CCMLBP>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<CCMLBP>& sol) {
	// TODO
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Conflict Constraints ****************************/
/*****************************************************************************************/
Solution<MLBPCC>::Solution(const Instance<MLBPCC>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPCC>& sol) {
	// TODO
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Partial Orders **********************************/
/*****************************************************************************************/
Solution<MLBPPO>::Solution(const Instance<MLBPPO>& inst) : Solution<MLBP>(inst)
{
	for (int i : inst.B[0]) {
		std::vector<int> levels(inst.m);
		pActual.push_back(levels);
	}
	for (int i = 0; i <= inst.m; i++) {
		std::vector<std::vector<int>> binHolder(inst.n[i]);
		for (int j : inst.B[i]) {
			std::vector<int> bin(0);
			binHolder.push_back(bin);
		}
		flow.push_back(binHolder);
	}
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPPO>& sol) {
	os << "Indexes: ";
	for (std::vector<int> arr : sol.item_bin_indexes) {
		os << arr;
	}
	os << ", P Variables: ";
	for (std::vector<int> arr : sol.pActual) {
		os << arr;
	}
	os << ", Flow leaving each bin: ";
	int levelCount = 0;
	int binCount = 0;
	for (std::vector<std::vector<int>> levels : sol.flow) {
		for (std::vector<int> bins : levels) {
			if (!bins.empty())
				os << "Level: " << levelCount << ", " << (levelCount == 0 ? "Item" : "Bin") << ": " << binCount << " => " << bins << std::endl;
			binCount++;
		}
		binCount = 0;
		levelCount++;
	}
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Time Windows ************************************/
/*****************************************************************************************/
Solution<MLBPTW>::Solution(const Instance<MLBPTW>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPTW>& sol) {
	// TODO
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Fragmentation Constraints ***********************/
/*****************************************************************************************/
Solution<MLBPFC>::Solution(const Instance<MLBPFC>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPFC>& sol) {
	// TODO
	return os;
}


