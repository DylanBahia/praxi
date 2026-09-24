#include <vector>
#include "ar_alg.h"
#include "cm.h"
#include "Rcpp.h"

List ar_alg_call(std::vector<double> y,
					const int& p,
					const double& b);
