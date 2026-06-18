#include <vector>
#include <list>
#include <tuple>
#include <utility>
#include <cmath> 
#include <iterator>
#include <numeric>
#include <limits>
#include "cm.h"
#include "quad_sum.h"

void coeff_update(std::list<cm>& cms,
			const std::vector<double>& z,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars);