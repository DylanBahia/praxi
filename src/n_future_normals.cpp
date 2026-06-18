#include <vector>
#include <list>
#include <tuple>
#include <limits>
#include "cm.h"
#include "anom_cost.h"

cm n_future_normals(const std::vector<double>& z,
			const int& n,
			const double& b,
			const cm& ck,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars){

	int k = ck.k;
	
	int coll_len = t-n-k;

	int p = std::get<0>(ar_pars);
	std::vector<double> sub_vec(z.begin()+t-p,z.begin()+t);
	std::pair<double,double> cost_and_mu = anom_cost(ck.coeffs,sub_vec,n,coll_len,ar_pars);
	double total_cost = b+cost_and_mu.first;

	cm candidate;
	candidate.cost = total_cost;
	candidate.k = t;
	std::vector<std::pair<int,int>> new_anoms = ck.anoms;
	new_anoms.push_back(std::make_pair(k+1,t-n));
	candidate.anoms = new_anoms;
	std::vector<double> new_pars = ck.pars;
	new_pars.push_back(cost_and_mu.second);
	candidate.pars = new_pars;
	candidate.coeffs = std::make_tuple(0,0,total_cost);
	candidate.mu_range = std::make_pair(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());

	return candidate;
}