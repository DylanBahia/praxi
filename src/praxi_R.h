#include <vector>
#include "ar_alg.h"
#include "cm.h"
#include "Rcpp.h"

Rcpp::NumericMatrix ar_alg_call(std::vector<double> y,
					const int& p,
					const double& b);
