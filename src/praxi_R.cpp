#include <vector>
#include "ar_alg.h"
#include "cm.h"

#include "Rcpp.h"	
using namespace Rcpp;


// [[Rcpp::export]]

Rcpp::NumericMatrix ar_alg_call(std::vector<double> y,
				       const int& p,
					const double& b){
					
	cm output = ar_alg(y,p,b);
	 
	Rcpp::NumericMatrix mat(output.pars.size(), 3);

	for (size_t i = 0; i < output.pars.size(); ++i) {
	    mat(i, 0) = output[i].anoms.first;
	    mat(i, 1) = output[i].anoms.second;
	    mat(i, 2) = output[i].pars;
	}

	return out;
}
