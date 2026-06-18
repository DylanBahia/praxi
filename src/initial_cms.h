#include <vector>
#include <list>
#include <utility>
#include <tuple>
#include <cmath>
#include <limits>
#include "cm.h"
#include "initial_coeff_update.h"

std::list<cm> initial_cms(const std::vector<double>& z,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars);