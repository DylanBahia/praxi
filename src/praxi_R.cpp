#include <vector>
#include "ar_alg.h

#include "Rcpp.h"	
using namespace Rcpp;


// [[Rcpp::export]]

std::list<std::vector<int>> praxi(std::vector<double> y,
					const int& p,
					const double& b){
					
	return ar_alg(y,p,b);
}
