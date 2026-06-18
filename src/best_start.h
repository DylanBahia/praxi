#include <vector>
#include <list>
#include <tuple>
#include <limits>
#include "cm.h"
#include "n_future_normals.h"

cm best_start(const std::vector<double>& z,
			const int& n,
			const double& b,
			const std::list<cm>& cms,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars);