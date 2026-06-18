#include <utility>
#include <tuple>
#include <cmath>


//returns the coefficients of the sum of the quadratics.
//old quadratic is of the form ax^2 + bx + c; three coefficients.
//new quadratic is of the form (dx+e)^2; two coefficients.
std::tuple<double,double,double> quad_sum(std::tuple<double,double,double> old_coeffs,
			     const std::pair<double,double>& new_coeffs);