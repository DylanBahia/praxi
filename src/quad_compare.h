#include <tuple>
#include <utility>
#include <cmath>


//given two quadratics, returns true if the second quadratic
//is always larger than the first.
//both quadratics of the form ax^2+bx+c.
std::pair<double,double> quad_compare(const std::tuple<double,double,double>& coeffs_one,
	 	    std::tuple<double,double,double> coeffs_two);