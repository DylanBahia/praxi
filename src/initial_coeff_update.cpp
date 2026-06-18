#include <vector>
#include <list>
#include <utility>
#include <tuple>
#include <numeric>
#include <limits>
#include "quad_sum.h"
#include "cm.h"

void initial_coeff_update(std::list<cm>& cms,
				const double& z_t,
				const int& t,
				const std::tuple<int, std::vector<double>, double, double>& ar_pars){
int p = std::get<0>(ar_pars);
std::vector<double> gam = std::get<1>(ar_pars);

for(cm& ck:cms){
	int dist = t-ck.k;
	double mu_coeff = 1 - std::accumulate(gam.begin(),gam.begin()+std::min(dist-1,p),0.0);
	std::pair<double,double> new_quad(mu_coeff,-z_t);
	ck.coeffs = quad_sum(ck.coeffs,new_quad);
}

}