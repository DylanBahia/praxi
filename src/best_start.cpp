#include <vector>
#include <list>
#include <tuple>
#include <limits>
#include "cm.h"
#include "n_future_normals.h"

cm best_start(const std::vector<double>& z,
			const int& n,
			const double& b,
			const std::list<cm>& cms,
			const int& t,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars){

	cm best_start = {std::numeric_limits<double>::infinity(),0,{},{},{},{}};
	for(const cm& ck:cms){
		if(ck.k<=(t-1-n)){
			cm candidate = n_future_normals(z,n,b,ck,t,ar_pars);
			if(candidate.cost<=best_start.cost){
				best_start = candidate;
			}
		}else{
			return best_start;
		}
	}
	return best_start;	
}