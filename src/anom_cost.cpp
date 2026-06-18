#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <cmath>
#include "cm.h"


//computes the sum of the AR coefficients from the
//i^th coefficient to the last one.
double gam_sum(const std::vector<double>& gam,
		const double& i,
		const double& p){
	
	double val = 0;
	for(int j=i; j<p; j++){
		val = val + gam[j];
	}
	return val;
}

//computes the cost of the anomalous segment with n future normals.
//z is a vector of the last n observations.
//mu is the mean of the anomalous segment.
//gam is a vector of the AR parameters.
//coeffs is a vector of the coefficients of the
//quadratic beginning at time k (given in R code).
std::pair<double,double> anom_cost(std::tuple<double,double,double> coeffs,
		     	const std::vector<double>& z,
		     	const int& n,
			const int& coll_len,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars){

	int p = std::get<0>(ar_pars);
	std::vector<double> gam = std::get<1>(ar_pars);
	double mu_0 = std::get<2>(ar_pars);
	double sig_0 = std::get<3>(ar_pars);
	

	for(int i=0; i<p-n; i++){
		double a = 1-gam_sum(gam,0,std::min(coll_len+i,p));
		std::get<0>(coeffs) = std::get<0>(coeffs)+std::pow(a,2);
		std::get<1>(coeffs) = std::get<1>(coeffs)-2*z[i]*a;
	}


	for(int i=0; i<n; i++){
		double a = gam_sum(gam,i,std::min(coll_len+i,p));
		//double b = 1-gam_sum(gam,0,std::min(coll_len+i,p));
		//coeffs[0] = coeffs[0]+std::pow(a,2)-std::pow(b,2);
		//coeffs[1] = coeffs[1]+2*z[p-n+i]*(a+b);
		std::get<0>(coeffs) = std::get<0>(coeffs)+std::pow(a,2);
		std::get<1>(coeffs) = std::get<1>(coeffs)+2*z[p-n+i]*a;
	}
	
	double c = std::get<0>(coeffs);
	double d = std::get<1>(coeffs);
	double e = std::get<2>(coeffs);

	double cost = e - (std::pow(d,2))/(4*c);
	
	double mu = -(d/(2*c))*sig_0+mu_0;

	std::pair<double,double> cost_and_mu(cost,mu);

	return cost_and_mu;

}