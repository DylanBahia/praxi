#include <vector>
#include <list>
#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>


double weighted_median(std::vector<double> A,
                       std::vector<long long> weight){
  long long n = A.size();
  long long left = 0;
  long long right = n-1;
  long long ind;
  
  while(true){
    ind = std::floor((left+right)/2);
    double sum_one = std::accumulate(weight.begin(),weight.begin()+ind,0.0);
    double sum_two = std::accumulate(weight.begin()+ind+1,weight.end(),0.0);
    
    if((sum_one+weight[ind])>=sum_two){
      if(sum_one<=(weight[ind]+sum_two)){
        break;
      }else{
        right = ind-1;
      }
    }else{
      left = ind+1;
    }
    
  }
  return(A[ind]);
}

double q_squared(const std::vector<double>& z){

	
	long long n = z.size();
	long long nchoose = (n*(n-1))/2;
	long long K = std::floor((nchoose+2)/4) + 1 + n + nchoose;

	std::vector<double> X(z.begin(),z.end());
	std::sort(X.begin(),X.end());

	std::vector<double> Y = X; 	
	std::reverse(Y.begin(),Y.end());
	std::transform(Y.begin(),Y.end(),Y.begin(),std::negate<double>());

	long long L = 0;
	long long R = std::pow(n,2);

	std::vector<long long> Lb(n);
	std::vector<long long> Rb(n,n);

	for(long long i=1; i<=n; i++){
		Lb[i-1] = n-i+2;
	}

	while(R-L>n){
		std::vector<double> A;
		std::vector<long long> weight;
		
		for(long long i=1; i<=n; i++){
			if(Lb[i-1]<=Rb[i-1]){
				A.push_back(X[i-1]+Y[std::floor((Lb[i-1]+Rb[i-1])/2)-1]);
				weight.push_back(Rb[i-1]-Lb[i-1]+1);
			}
		}
		double a_m = weighted_median(A,weight);	
		std::vector<long long> P(n,0);
		std::vector<long long> Q(n,0);
		//compute P
		long long j = 0;
		for(long long i=n; i>=1; i--){
			while(true){
				if(j==n){
					break;				
				}else if((X[i-1]+Y[j])>=a_m){
					break;
				}else{
					j = j+1;
				}
			}
			P[i-1] = j;
		} 
		//compute Q
		j = 0;
		for(long long i=1; i<=n; i++){
			while(true){
				if(j==n){
					break;				
				}else if((X[i-1]+Y[n-j-1])<=a_m){
					break;
				}else{
					j = j+1;
				}
			}
			Q[i-1] = n+1-j;
		} 
		if(K<=std::accumulate(P.begin(),P.end(),0LL)){
			Rb = P;
		}else if(K>std::accumulate(Q.begin(),Q.end(),-n)){
			Lb = Q;
		}else{
			return(a_m);
		}
		L = std::accumulate(Lb.begin(),Lb.end(),0LL)-n;
		R = std::accumulate(Rb.begin(),Rb.end(),0LL);
		
	}

	std::vector<double> A;
	for(long long i=2; i<=n; i++){
		if(Lb[i-1]<=Rb[i-1]){
			for(long long j=Lb[i-1];j<=Rb[i-1];j++){
				A.push_back(X[i-1]+Y[j-1]);
			}
		}
	}
	std::sort(A.begin(),A.end());

	return std::pow(2.2191*A[(K-L)-1],2);
	
}


