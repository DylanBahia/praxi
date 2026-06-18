#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "q_squared.h"

double gamma_hat(const std::vector<double>& y,
			const int& n,
			const int& p){
	std::vector<double> u(y.begin(),y.end()-p);
	std::vector<double> v(y.begin()+p,y.end());
	std::vector<double> u_plus_v(u.size());
	std::vector<double> u_minus_v(u.size());
	std::transform(u.begin(),u.end(),v.begin(),u_plus_v.begin(),[](double a, double b){return a+b;});
	std::transform(u.begin(),u.end(),v.begin(),u_minus_v.begin(),[](double a, double b){return a-b;});
	return (q_squared(u_plus_v)-q_squared(u_minus_v))/4;
}

std::vector<double> yule_walker(const std::vector<double>& covs,
				const int& p){
	double vari = covs[0];
	std::vector<double> h(covs.begin()+1,covs.end());

	std::vector<double> v(p);
	v[0] = h[0]/vari;

	if(p>1){
		for(int i = 1; i<p; i++){ 
			double bet = (h[i]-std::inner_product(h.begin(),h.begin()+i,v.rbegin()+p-i,0.0))/(vari-std::inner_product(h.begin(),h.begin()+i,v.begin(),0.0));
			std::vector<double> z(i);
			for(int j = 0; j<i; j++){
				z[j] = v[j] - bet*v[i-j-1];
			}
			std::copy(z.begin(),z.end(),v.begin());
			v[i] = bet;
		}
	}
	double sig = std::sqrt(vari-std::inner_product(h.begin(),h.end(),v.begin(),0.0));

	std::vector<double> pars(p+1);
	pars[0] = sig;
	std::copy(v.begin(),v.end(),pars.begin()+1);

	return pars;
}

std::vector<double> robust_pars_yw(const std::vector<double>& y,
					const double& marg_var,
					const int& p,
					const int& n){
	std::vector<double> covs(p+1);
	covs[0] = marg_var;
	for(int i = 1; i<=p; i++){
		covs[i] = gamma_hat(y,n,i);
	}

	return yule_walker(covs,p);
}