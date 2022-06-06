#ifndef __PROBLEMS_H__
#define __PROBLEMS_H__

#include <string_view>

/**
 * Defines for each problem an own problem type.
 * This problem type is used in generic algorithms.
 */

/* Bin Packing Problem */
struct BP
{
	constexpr static const std::string_view name = "Bin Packing Problem";
};

/* Multi-Level Bin Packing Problem */
struct MLBP
{
	constexpr static const std::string_view name = "Multi-Level Bin Packing Problem";
};

/* Class Constrained Multi-Level Bin Packing Problem */
struct CCMLBP
{
	constexpr static const std::string_view name = "Class Constrained Multi-Level Bin Packing Problem";
};

/* Multi-Level Bin Packing Problem with Conflict Constraints */
struct MLBPCC
{
	constexpr static const std::string_view name = "Multi-Level Bin Packing Problem with Conflict Constraints";
};

/* Multi-Level Bin Packing Problem with Partial Orders */
struct MLBPPO
{
	constexpr static const std::string_view name = "Multi-Level Bin Packing Problem with Partial Orders";
};

/* Multi-Level Bin Packing Problem with Time Windows */
struct MLBPTW
{
	constexpr static const std::string_view name = "Multi-Level Bin Packing Problem with Time Windows";
};

/* Multi-Level Bin Packing Problem with Fragmentation Constraints */
struct MLBPFC
{
	constexpr static const std::string_view name = "Multi-Level Bin Packing Problem with Fragmentation Constraints";
};


#endif // __PROBLEMS_H__
