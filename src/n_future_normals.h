#include <vector>
#include <list>
#include <tuple>
#include <limits>
#include "cm.h"
#include "anom_cost.h"

cm n_future_normals(const std::vector<double>& z,
			const int& n,
			const double& b,
			const cm& ck,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars);