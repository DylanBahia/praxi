#include <vector>
#include <list>
#include <utility>
#include <tuple>
#include <iostream>
#include <limits>
#include "cm.h"
#include "initial_cms.h"
#include "coeff_update.h"
#include "pruning.h"
#include "cm_update.h"






cm recursion(const std::vector<double>& z,
			const double& b,
			const int& n,
			const int& p,
			const std::vector<double>& gam,
			const double& mu_0,
			const double& sig_0
);
