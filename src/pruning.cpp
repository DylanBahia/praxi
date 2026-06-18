#include <vector>
#include <list>
#include <tuple>
#include <utility>
#include <iterator>
#include <algorithm>
#include <limits>
#include "cm.h"
#include "quad_compare.h"

void pruning(std::list<cm>& cms,
			const int& p){

if(cms.size()-(p-1)-p>1){

auto benchmark = std::prev(cms.end(),p-1+p);

for(auto pos = std::make_reverse_iterator(benchmark); pos!=cms.rend();){
	std::pair<double,double> bounds = quad_compare(benchmark->coeffs,pos->coeffs);
	if(bounds.first==bounds.second){
		auto temp = pos.base();
		--temp;
		pos = std::list<cm>::reverse_iterator(cms.erase(temp));
	}else{
		std::pair<double,double> new_range(std::max(bounds.first,pos->mu_range.first),
							std::min(bounds.second,pos->mu_range.second));
		if(new_range.first>=new_range.second){
			auto temp = pos.base();
			--temp;
			pos = std::list<cm>::reverse_iterator(cms.erase(temp));
		}else{
			pos->mu_range = new_range;
			++pos;
		}

	}
}

}	
}