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
	 
	int n = output.pars.size();
	 
	Rcpp::NumericMatrix mat(n, 3);

	for (std::size_t i = 0; i < n; ++i) {
        	mat(i, 0) = std::get<0>(output.anoms[i]);
        	mat(i, 1) = std::get<1>(output.anoms[i]);
        	mat(i, 2) = output.pars[i];
    	}
    	
    	Rcpp::colnames(mat) = Rcpp::CharacterVector::create("Start", "End", "Mean");
    	
	return mat;
}
