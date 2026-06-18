#include <utility>
#include <tuple>
#include <cmath>

std::tuple<double,double,double> quad_sum(std::tuple<double,double,double> old_coeffs,
			     const std::pair<double,double>& new_coeffs){
	
	double d = new_coeffs.first;
	double e = new_coeffs.second;

	std::pair<double,double> operand(std::pow(d,2),2*d*e);
	
	std::get<0>(old_coeffs) += operand.first;
	std::get<1>(old_coeffs) += operand.second;
	
	return old_coeffs;
}
