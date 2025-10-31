// This program demonstrates loop-based parallelism using:
//   + STL-styled iterators
//   + plain integral indices

#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>

// Procedure: for_each
void for_each(int N) {

  tf::Executor executor;
  tf::Taskflow taskflow;

  std::vector<int> range(N);
  std::iota(range.begin(), range.end(), 0);

  taskflow.for_each(range.begin(), range.end(), [&] (int i) {
    printf("for_each on container item: %d\n", i);
  }, tf::StaticPartitioner());

  executor.run(taskflow).get();
}

// Procedure: for_each_by_index
void for_each_by_index(int N) {

  tf::Executor executor;
  tf::Taskflow taskflow;

  // [0, N) with a step size of 2
  taskflow.for_each_index(0, N, 2, [] (int i) {
    printf("for_each_index on index: %d\n", i);
  });

  executor.run(taskflow).wait();

  // [0, N) with a step size of 2 using tf::IndexRange
  tf::IndexRange<int> range(0, N, 2);
  
  taskflow.for_each_by_index(range, [](tf::IndexRange<int> subrange) {
    for(int i=subrange.begin(); i<subrange.end(); i+=subrange.step_size()) {
      printf("for_each_by_index on index (subrange): %d\n", i);
    }
  });
  
  executor.run(taskflow).wait();
}

// ----------------------------------------------------------------------------

// Function: main
int main(int argc, char* argv[]) {

  if(argc != 2) {
    std::cerr << "Usage: ./parallel_for num_iterations" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  for_each(std::atoi(argv[1]));
  for_each_by_index(std::atoi(argv[1]));

  return 0;
}






