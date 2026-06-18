#include <vector>
#include <algorithm>

double robust_mu(std::vector<double> y,
			const int& n){
	sort(y.begin(),y.end());
	double median;
	
	if(n%2 == 0){
		median = (y[(n/2)-1]+y[n/2])/2.0;
	}else{
		median = y[(n-1)/2];
	}

	return median;
}