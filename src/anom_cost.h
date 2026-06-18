#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <cmath>
#include "cm.h"


double gam_sum(const std::vector<double>& gam,
		const double& i,
		const double& p);

std::pair<double,double> anom_cost(std::tuple<double,double,double> coeffs,
		     	const std::vector<double>& z,
		     	const int& n,
			const int& coll_len,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars);