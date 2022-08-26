#include <iostream>
#include <boost/heap/binomial_heap.hpp>
#include "bloom_filter.hpp"

int main(){
    bloom_parameters param;
    param.projected_element_count = 500;
    param.false_positive_probability = 0.001;
    param.compute_optimal_parameters();
    bloom_filter filt = bloom_filter(param);
    int64_t num = 5;
    filt.insert(num);
    std::cout << filt.contains(num)<<std::endl;

    return 0;
}
