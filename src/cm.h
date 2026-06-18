#ifndef cm_H
#define cm_H

struct cm {
	double cost;
	int k;
	std::vector<std::pair<int,int>> anoms;
	std::vector<double> pars;
	std::tuple<double,double,double> coeffs;
	std::pair<double,double> mu_range;
};

#endif