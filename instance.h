#ifndef __INSTANCE_H__
#define __INSTANCE_H__


#include <string>
#include <vector>
#include <iostream>

#include "problems.h"


template<typename>
struct Solution;


/*
 * Generic problem instance class.
 * To write an instance class for a specific problem,
 * speialize this class for the corresponding problem type.
 */
template<typename ProbT>
struct Instance { };




/*****************************************************************************************/
/** Bin Packing Problem ******************************************************************/
/*****************************************************************************************/
template<>
struct Instance<BP>
{
	friend std::ostream& operator<<(std::ostream& os, const Instance<BP>& inst);

	Instance(const std::string& input_file);

	int objective(const Solution<BP>& sol) const;

	std::string filename;

	int n;              // number of items
	std::vector<int> I; // index set of items [0,...,n-1]
	std::vector<int> s; // size of itmes
	int smax;           // maximum bin size
};

std::ostream& operator<<(std::ostream& os, const Instance<BP>& inst);




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem ******************************************************/
/*****************************************************************************************/
template<>
struct Instance<MLBP>
{
	friend std::ostream& operator<<(std::ostream& os, const Instance<MLBP>& inst);

	Instance(const std::string& input_file);

	int objective(const Solution<MLBP>& sol) const;

	std::string filename;

	int m;               // number of levels
	std::vector<int> n;  // number of items/bins at each level; level 0: number of items, level > 0: number of bins

	std::vector<std::vector<int> > s;  // size of each item/bin at each level
	std::vector<std::vector<int> > w;  // capacity of each bin; index 0 is reserved for items and therefore always empty
	std::vector<std::vector<int> > c;  // cost of each bin; index 0 is reserved for items and therefore always empty

	std::vector<int> M;  // index set of levels 1...m
	std::vector<std::vector<int> > B;  // index set of items/bins for each level

protected:
	Instance() { }
	void read_instance(std::ifstream& file);
};

std::ostream& operator<<(std::ostream& os, const Instance<MLBP>& inst);




/*****************************************************************************************/
/** Class Constrained Multi-Level Bin Packing Problem ************************************/
/*****************************************************************************************/
template<>
struct Instance<CCMLBP> : public Instance<MLBP>
{
	friend std::ostream& operator<<(std::ostream& os, const Instance<CCMLBP>& inst);

	Instance(const std::string& input_file);

	int q;                   // number of classes
	std::vector<int> Q;      // bound of different classes within a bin for each level, index 0 is empty
	std::vector<int> kappa;  // class of each item
};

std::ostream& operator<<(std::ostream& os, const Instance<CCMLBP>& inst);




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Conflict Constraints ****************************/
/*****************************************************************************************/
template<>
struct Instance<MLBPCC> : public Instance<MLBP>
{
	friend std::ostream& operator<<(std::ostream& os, const Instance<MLBPCC>& inst);

	Instance(const std::string& input_file);

	std::vector<std::vector<int> > conflict;  // conflict matrix
};

std::ostream& operator<<(std::ostream& os, const Instance<MLBPCC>& inst);




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Partial Orders **********************************/
/*****************************************************************************************/
template<>
struct Instance<MLBPPO> : public Instance<MLBP>
{
	friend std::ostream& operator<<(std::ostream& os, const Instance<MLBPPO>& inst);

	Instance(const std::string& input_file);

	std::vector<std::pair<int, int> > pos; // partial orders between items
};




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Time Windows ************************************/
/*****************************************************************************************/
template<>
struct Instance<MLBPTW> : public Instance<MLBP>
{
	friend std::ostream& operator<<(std::ostream& os, const Instance<MLBPTW>& inst);

	Instance(const std::string& input_file);

	int objective(const Solution<MLBPTW>& sol) const;

	int p; // penalty factor
	std::vector<int> e;  // earliest starting time of each item
	std::vector<int> l;  // latest starting time of each item
};




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Fragmentation Constraints ***********************/
/*****************************************************************************************/
template<>
struct Instance<MLBPFC> : public Instance<MLBP>
{
	friend std::ostream& operator<<(std::ostream& os, const Instance<MLBPFC>& inst);

	Instance(const std::string& input_file);

	int objective(const Solution<MLBPFC>& sol) const;

	int q; // number of groups
	int p; // penalty factor
	std::vector<int> g; // group of each item
};

std::ostream& operator<<(std::ostream& os, const Instance<MLBPFC>& inst);


#endif // __INSTANCE_H__
