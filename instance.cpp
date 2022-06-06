#include "instance.h"
#include "solution.h"
#include "lib/util.h"
#include "users.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <cassert>



/*************************************************************************************************/
/* BP ********************************************************************************************/
/*************************************************************************************************/
Instance<BP>::Instance(const std::string& input_file) : filename(input_file)
{
	std::ifstream file(input_file);
	if (!file.is_open())
		throw std::runtime_error("Cannot open file!");

	file >> n >> smax;
	
	I.reserve(n);
	s.assign(n, -1);
	for (int i = 0; i < n; i++) {
		I.push_back(i);
		file >> s[i];
	}
}

int Instance<BP>::objective(const Solution<BP>& sol) const
{
	return sol.total_bins;
}


std::ostream& operator<<(std::ostream& os, const Instance<BP> &inst)
{
	os << "n=" << inst.s.size() << ", smax=" << inst.smax << ", s=" << inst.s;
	return os;
}




/*************************************************************************************************/
/* MLBP ******************************************************************************************/
/*************************************************************************************************/
Instance<MLBP>::Instance(const std::string& input_file) : filename(input_file)
{
	std::ifstream file(input_file);
	if (!file.is_open())
		throw std::runtime_error("Cannot open file!");

	read_instance(file);
}

void Instance<MLBP>::read_instance(std::ifstream& file)
{
	// reading data from the file
	file >> m;

	// read number of tools and bins for each level
	n.assign(m+1, 0);
	M.reserve(m);
	for (int i = 0; i < m+1; i++) {
		if (i != 0) M.push_back(i);
		file >> n[i];
	}

	// read item sizes and bin sizes for each level
	s.assign(m+1, std::vector<int>());
	B.assign(m+1, std::vector<int>());
	for (int i = 0; i < m+1; i++) {
		s[i].assign(n[i], 0);
		for (int j = 0; j < n[i]; j++) {
			B[i].push_back(j);
			file >> s[i][j];
		}
	}

	// read capacities of bins for each level
	w.assign(m+1, std::vector<int>());
	for (int i = 1; i < m+1; i++) {
		w[i].assign(n[i], 0);
		for (int j = 0; j < n[i]; j++)
			file >> w[i][j];
	}

	// read cost of each bin for each level
	c.assign(m+1, std::vector<int>());
	for (int i = 1; i < m+1; i++) {
		c[i].assign(n[i], 0);
		for (int j = 0; j < n[i]; j++)
			file >> c[i][j];
	}
}

int Instance<MLBP>::objective(const Solution<MLBP>& sol) const
{
	return sol.cost;
}

std::ostream& operator<<(std::ostream& os, const Instance<MLBP> &inst)
{
	os << "m=" << inst.m << ", n=" << inst.n;
	return os;
}




/*************************************************************************************************/
/* CCMLBP ****************************************************************************************/
/*************************************************************************************************/
Instance<CCMLBP>::Instance(const std::string& input_file) : Instance<MLBP>()
{
	filename = input_file;

	std::ifstream file(input_file);
	if (!file.is_open())
		throw std::runtime_error("Cannot open file!");

	read_instance(file);

	file >> q;
	Q.assign(m+1, 0);
	for (int i = 1; i < m+1; i++)
		file >> Q[i];
	kappa.assign(n[0], 0);
	for (int j = 0; j < n[0]; j++)
		file >> kappa[j];
}

std::ostream& operator<<(std::ostream& os, const Instance<CCMLBP> &inst)
{
	os << "m=" << inst.m << ", n=" << inst.n << ", q=" << inst.q << ", Q=" << inst.Q;
	return os;
}




/*************************************************************************************************/
/* MLBPCC ****************************************************************************************/
/*************************************************************************************************/
Instance<MLBPCC>::Instance(const std::string& input_file) : Instance<MLBP>()
{
	filename = input_file;

	std::ifstream file(input_file);
	if (!file.is_open())
		throw std::runtime_error("Cannot open file!");

	read_instance(file);

	conflict.assign(n[0], std::vector<int>(n[0], 0));
	std::string line;
	int j, q;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		iss >> j;
		while(iss >> q) {
			conflict[j][q] = 1;
			conflict[q][j] = 1;
		}
	}
}

std::ostream& operator<<(std::ostream& os, const Instance<MLBPCC> &inst)
{
	os << "m=" << inst.m << ", n=" << inst.n << "\n";
	return os;
}




/*************************************************************************************************/
/* MLBPPO ****************************************************************************************/
/*************************************************************************************************/
Instance<MLBPPO>::Instance(const std::string& input_file) : Instance<MLBP>()
{
	filename = input_file;

	std::ifstream file(input_file);
	if (!file.is_open())
		throw std::runtime_error("Cannot open file!");

	read_instance(file);

	int p = 0;
	file >> p;

	pos.reserve(p);
	int v1, v2;
	for (int i = 0; i < p; i++) {
		file >> v1 >> v2;
		pos.emplace_back(v1, v2);
	}
}

std::ostream& operator<<(std::ostream& os, const Instance<MLBPPO> &inst)
{
	os << "m=" << inst.m << ", n=" << inst.n << ", #pos=" << inst.pos.size();
	return os;
}




/*************************************************************************************************/
/* MLBPTW ****************************************************************************************/
/*************************************************************************************************/
Instance<MLBPTW>::Instance(const std::string& input_file) : Instance<MLBP>()
{
	filename = input_file;

	std::ifstream file(input_file);
	if (!file.is_open())
		throw std::runtime_error("Cannot open file!");

	read_instance(file);

	file >> p;

	e.assign(n[0], -1);
	l.assign(n[0], -1);
	for (int i = 0; i < n[0]; i++)
		file >> e[i];
	for (int i = 0; i < n[0]; i++)
		file >> l[i];
}

int Instance<MLBPTW>::objective(const Solution<MLBPTW>& sol) const
{
	int obj = 0;
	// TODO
	return obj;
}

std::ostream& operator<<(std::ostream& os, const Instance<MLBPTW> &inst)
{
	os << "m=" << inst.m << ", n=" << inst.n << ", pl=["
	   << *std::min_element(inst.e.begin(), inst.e.end()) << "," << *std::max_element(inst.l.begin(), inst.l.end()) << "], p=" << inst.p;
	return os;
}




/*************************************************************************************************/
/* MLBPFC ****************************************************************************************/
/*************************************************************************************************/
Instance<MLBPFC>::Instance(const std::string& input_file) : Instance<MLBP>()
{
	filename = input_file;

	std::ifstream file(input_file);
	if (!file.is_open())
		throw std::runtime_error("Cannot open file!");

	read_instance(file);

	file >> q >> p;
	g.reserve(n[0]);
	int gj;
	for (int i = 0; i < n[0]; i++) {
		file >> gj;
		g.push_back(gj-1);
	}
}

int Instance<MLBPFC>::objective(const Solution<MLBPFC>& sol) const
{
	int obj = 0;
	// TODO
	return obj;
}

std::ostream& operator<<(std::ostream& os, const Instance<MLBPFC> &inst)
{
	os << "m=" << inst.m << ", n=" << inst.n << ", q=" << inst.q << ", p=" << inst.p << "\n";
	os << "g=" << inst.g;
	return os;
}

