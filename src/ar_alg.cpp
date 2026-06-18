#include <vector>
#include <algorithm>
#include "robust_mu.h"
#include "robust_marginal_var.h"
#include "robust_pars_yw.h"
#include "recursion.h"

std::list<std::vector<int>> ar_alg(std::vector<double> y,
					const int& p,
					const double& b){

	int n = y.size();
	
	//centering observations
	double mu_0 = robust_mu(y,n);
	transform(y.begin(),y.end(),y.begin(),[mu_0](double a){return a - mu_0;});
	
	//estimating AR model
	double marg_var = robust_marginal_var(y);
	std::vector<double> yw_pars = robust_pars_yw(y,marg_var,p,n);
	double sig_0 = yw_pars[0];
	std::vector<double> gam(yw_pars.begin()+1,yw_pars.end());
	
	//taking residuals
	std::vector<double> z = y;
	for(int i = p; i<n; i++){
		for(int j = 0; j<p; j++){
			z[i] = z[i] - gam[j]*y[i-j-1];
		}
	}
	
	//standardising residuals
	std::transform(z.begin(),z.end(),z.begin(),[sig_0](double a){return a/sig_0;});	

	return recursion(z,b,n,p,gam,mu_0,sig_0);	
}


