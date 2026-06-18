#include <vector>
#include <list>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#include "cm.h"
#include "best_start.h"

cm cm_normal(const double& z_t,
		const cm& cm_t,
		const int& t){

	double c1_cost = cm_t.cost + std::pow(z_t,2);
	cm c1 = {c1_cost,t,cm_t.anoms,cm_t.pars,std::make_tuple(0,0,c1_cost),std::make_pair(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity())};

	return c1;
} 

void cm_update(const std::vector<double>& z,
			const double& b,
			std::list<cm>& cms,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars){

	cm ct = cm_normal(z[t-1],cms.back(),t);
	
	int p = std::get<0>(ar_pars);

	cm best_start_point = best_start(z,p,b,cms,t,ar_pars);

	if(ct.cost<=best_start_point.cost){
		cms.push_back(ct);
	}else{
		cms.push_back(best_start_point);
	} 	

}

cm dm_update(const std::vector<double>& z,
			const double& b,
			std::list<cm>& cms,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars){

	cm dm = cms.back();
	
	int p = std::get<0>(ar_pars);
	
	std::list<cm> best_starts;

	for(int i=0; i<p; i++){
		best_starts.push_back(best_start(z,i,b,cms,t,ar_pars));
	}
	
	cm best_start = *std::min_element(best_starts.begin(),best_starts.end(),[](const cm& a, const cm& b){return a.cost<b.cost;});

	if(dm.cost<=best_start.cost){
		return dm;
	}else{
		return best_start;
	} 	

}

			
			
