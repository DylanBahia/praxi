#include <vector>
#include <list>
#include <utility>
#include <tuple>
#include <iostream>
#include <limits>
#include "cm.h"
#include "initial_cms.h"
#include "coeff_update.h"
#include "pruning.h"
#include "cm_update.h"







cm recursion(const std::vector<double>& z,
			const double& b,
			const int& n,
			const int& p,
			const std::vector<double>& gam,
			const double& mu_0,
			const double& sig_0
){


std::tuple<int, std::vector<double>, double, double> ar_pars(p,gam,mu_0,sig_0);

std::list<cm> cms = initial_cms(z,ar_pars);

for(int t = 3*p; t<=n; t++){
	coeff_update(cms,z,t,ar_pars);
	pruning(cms,p);
	cm_update(z,b,cms,t,ar_pars);
}

cm final = dm_update(z,b,cms,n,ar_pars);


//for(auto i:final.anoms){
//	std::cout << i.first << std::endl;
//	std::cout << i.second << std::endl;
//}

//for(auto i:final.pars){
//	std::cout << i << std::endl;
//}

//std::list<std::vector<int>> output;

//for(auto i:final.anoms){
//	output.push_back({i.first,i.second});
//}

return final;

}


