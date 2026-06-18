#include <vector>
#include <list>
#include <tuple>
#include <utility>
#include <cmath>
#include <iterator>
#include <numeric>
#include <limits>
#include "cm.h"
#include "quad_sum.h"

void coeff_update(std::list<cm>& cms,
			const std::vector<double>& z,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars){

	int p = std::get<0>(ar_pars);
	std::vector<double> gam = std::get<1>(ar_pars);

	auto cond = std::prev(cms.end(),p);
	for(auto ck = cms.begin(); ck!= cond; ++ck){
		int dist = t-p-ck->k;
		double mu_coeff = 1 - std::accumulate(gam.begin(),gam.begin()+std::min(dist-1,p),0.0);
		std::pair<double,double> new_quad(mu_coeff,-z[t-p-1]);
		ck->coeffs = quad_sum(ck->coeffs,new_quad);
		std::get<2>(ck->coeffs) += std::pow(z[t-1],2);	

	}

	for(auto ck = cond; ck != cms.end(); ++ck){
		std::get<2>(ck->coeffs) += std::pow(z[t-1],2);
	}

}