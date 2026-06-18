#include <tuple>
#include <utility>
#include <cmath>

std::pair<double,double> quad_compare(const std::tuple<double,double,double>& coeffs_one,
	 	  std::tuple<double,double,double> coeffs_two){

	std::get<0>(coeffs_two) -= std::get<0>(coeffs_one);
	std::get<1>(coeffs_two) -= std::get<1>(coeffs_one);
	std::get<2>(coeffs_two) -= std::get<2>(coeffs_one);	
	
	double a = std::get<0>(coeffs_two);
	double b = std::get<1>(coeffs_two);
	double c = std::get<2>(coeffs_two);

	double term = std::pow(b/(2*a),2) - c/a;	

	std::pair<double,double> bounds = {0,0};

	if(term>0){
		bounds.first = -term - b/(2*a);
		bounds.second = term - b/(2*a);
	}

	return bounds;
		

}

