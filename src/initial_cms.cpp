#include <vector>
#include <list>
#include <utility>
#include <tuple>
#include <cmath>
#include <limits>
#include "cm.h"
#include "initial_coeff_update.h"


std::list<cm> initial_cms(const std::vector<double>& z,
			const std::tuple<int, std::vector<double>, double, double>& ar_pars){

int p = std::get<0>(ar_pars);
std::vector<double> gam = std::get<1>(ar_pars);
std::list<cm> cms;
double log_like_val = 0;

cms.emplace_back(cm{log_like_val,p,{},{},std::make_tuple(0,0,log_like_val),std::make_pair(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity())});

for(int t = p+1; t<=2*p; t++){
	log_like_val = log_like_val+std::pow(z[t-1],2);

	initial_coeff_update(cms,z[t-1],t,ar_pars);

	cms.emplace_back(cm{log_like_val,p,{},{},std::make_tuple(0,0,log_like_val),std::make_pair(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity())});

	if(p>1){
		for(int t = 2*p+1; t<=3*p-1; t++){
			log_like_val = log_like_val+std::pow(z[t-1],2);
                        cms.emplace_back(cm{log_like_val,p,{},{},std::make_tuple(0,0,log_like_val),std::make_pair(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity())});
			
		}
	}
}

return(cms);

}