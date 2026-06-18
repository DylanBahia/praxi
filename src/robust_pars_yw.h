#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "q_squared.h"

double gamma_hat(const std::vector<double>& y,
			const int& n,
			const int& p);

std::vector<double> yule_walker(const std::vector<double>& covs,
				const int& p);

std::vector<double> robust_pars_yw(const std::vector<double>& y,
					const double& marg_var,
					const int& p,
					const int& n);