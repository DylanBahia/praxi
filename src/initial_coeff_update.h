#include <vector>
#include <list>
#include <utility>
#include <tuple>
#include <numeric>
#include <limits>
#include "cm.h"
#include "quad_sum.h"	


void initial_coeff_update(std::list<cm>& cms,
				const double& z_t,
				const int& t,
				const std::tuple<int, std::vector<double>, double, double>& ar_pars);
