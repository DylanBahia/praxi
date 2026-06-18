#include <vector>
#include <list>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include "cm.h"
#include "best_start.h"

cm cm_normal(const double& z_t,
		const cm& cm_t,
		const int& t);

void cm_update(const std::vector<double>& z,
			const double& b,
			std::list<cm>& cms,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars);

cm dm_update(const std::vector<double>& z,
			const double& b,
			std::list<cm>& cms,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars);