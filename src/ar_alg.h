#include <vector>
#include <algorithm>
#include "robust_mu.h"
#include "robust_marginal_var.h"
#include "robust_pars_yw.h"
#include "recursion.h"



std::list<std::vector<int>> ar_alg(std::vector<double> y,
					const int& p,
					const double& b);
