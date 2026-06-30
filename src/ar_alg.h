#include <vector>
#include <algorithm>
#include "robust_mu.h"
#include "robust_marginal_var.h"
#include "robust_pars_yw.h"
#include "recursion.h"
#include "cm.h"


cm ar_alg(std::vector<double> y,
					const int& p,
					const double& b);
