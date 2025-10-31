#include "reduce_sum.hpp"
#include <omp.h>

// reduce_sum_omp
void reduce_sum_omp(unsigned nthreads) {

  omp_set_num_threads(nthreads);

  double sum = 0.0;

  #pragma omp parallel for reduction(+: sum)
  for(size_t i=0; i<vec.size(); ++i) {
    sum += vec[i];
  }
}

std::chrono::microseconds measure_time_omp(unsigned num_threads) {
  auto beg = std::chrono::high_resolution_clock::now();
  reduce_sum_omp(num_threads);
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
}
