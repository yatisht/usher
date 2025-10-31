#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/sort.hpp>

// ----------------------------------------------------------------------------
// Data Type
// ----------------------------------------------------------------------------

struct MoveOnly1{

  int a {-1234};
  
  MoveOnly1() = default;

  MoveOnly1(const MoveOnly1&) = delete;
  MoveOnly1(MoveOnly1&&) = default;

  MoveOnly1& operator = (const MoveOnly1& rhs) = delete;
  MoveOnly1& operator = (MoveOnly1&& rhs) = default;

};

// ----------------------------------------------------------------------------
// parallel sort
// ----------------------------------------------------------------------------

template <typename T>
void ps_pod(size_t W, size_t N) {

  std::srand(static_cast<unsigned int>(time(NULL)));

  std::vector<T> data(N);

  for(auto& d : data) {
    d = ::rand() % 1000 - 500;
  }

  tf::Taskflow taskflow;
  tf::Executor executor(W);

  taskflow.sort(data.begin(), data.end());

  executor.run(taskflow).wait();

  REQUIRE(std::is_sorted(data.begin(), data.end()));
}

TEST_CASE("ParallelSort.int.1.100000" * doctest::timeout(300)) {
  ps_pod<int>(1, 100000);
}

TEST_CASE("ParallelSort.int.2.100000" * doctest::timeout(300)) {
  ps_pod<int>(2, 100000);
}

TEST_CASE("ParallelSort.int.3.100000" * doctest::timeout(300)) {
  ps_pod<int>(3, 100000);
}

TEST_CASE("ParallelSort.int.4.100000" * doctest::timeout(300)) {
  ps_pod<int>(4, 100000);
}

TEST_CASE("ParallelSort.ldouble.1.100000" * doctest::timeout(300)) {
  ps_pod<long double>(1, 100000);
}

TEST_CASE("ParallelSort.ldouble.2.100000" * doctest::timeout(300)) {
  ps_pod<long double>(2, 100000);
}

TEST_CASE("ParallelSort.ldouble.3.100000" * doctest::timeout(300)) {
  ps_pod<long double>(3, 100000);
}

TEST_CASE("ParallelSort.ldouble.4.100000" * doctest::timeout(300)) {
  ps_pod<long double>(4, 100000);
}

struct Object {

  std::array<int, 10> integers;

  int sum() const {
    int s = 0;
    for(const auto i : integers) {
      s += i;
    }
    return s;
  }
};

void ps_object(size_t W, size_t N) {

  std::srand(static_cast<unsigned int>(time(NULL)));

  std::vector<Object> data(N);

  for(auto& d : data) {
    for(auto& i : d.integers) {
      i = ::rand();
    }
  }

  tf::Taskflow taskflow;
  tf::Executor executor(W);

  taskflow.sort(data.begin(), data.end(), [](const auto& l, const auto& r){
    return l.sum() < r.sum();
  });

  executor.run(taskflow).wait();

  REQUIRE(std::is_sorted(data.begin(), data.end(),
    [](const auto& l, const auto& r){ return l.sum() < r.sum(); }
  ));
}

TEST_CASE("ParallelSort.object.1.100000" * doctest::timeout(300)) {
  ps_object(1, 100000);
}

TEST_CASE("ParallelSort.object.2.100000" * doctest::timeout(300)) {
  ps_object(2, 100000);
}

TEST_CASE("ParallelSort.object.3.100000" * doctest::timeout(300)) {
  ps_object(3, 100000);
}

TEST_CASE("ParallelSort.object.4.100000" * doctest::timeout(300)) {
  ps_object(4, 100000);
}

void move_only_ps(unsigned W) {
  
  std::vector<MoveOnly1> vec(1000000);
  for(auto& i : vec) {
    i.a = rand()%100;
  }

  tf::Taskflow taskflow;
  tf::Executor executor(W);

  taskflow.sort(vec.begin(), vec.end(),
    [](const MoveOnly1& m1, const MoveOnly1&m2) {
      return m1.a < m2.a;
    }
  );

  executor.run(taskflow).wait();

  for(size_t i=1; i<vec.size(); i++) {
    REQUIRE(vec[i-1].a <= vec[i].a);
  }

}

TEST_CASE("ParallelSort.MoveOnlyObject.1thread" * doctest::timeout(300)) {
  move_only_ps(1);
}

TEST_CASE("ParallelSort.MoveOnlyObject.2threads" * doctest::timeout(300)) {
  move_only_ps(2);
}

TEST_CASE("ParallelSort.MoveOnlyObject.3threads" * doctest::timeout(300)) {
  move_only_ps(3);
}

TEST_CASE("ParallelSort.MoveOnlyObject.4threads" * doctest::timeout(300)) {
  move_only_ps(4);
}

// ----------------------------------------------------------------------------
// Parallel Sort with  Async Tasks
// ----------------------------------------------------------------------------

void async(size_t W) {

  std::srand(static_cast<unsigned int>(time(NULL)));
  
  tf::Executor executor(W);
  std::vector<int> data;

  for(size_t n=0; n < 100000; n = (n ? n*10 : 1)) {
    
    data.resize(n);

    for(auto& d : data) {
      d = ::rand() % 1000 - 500;
    }
  
    executor.async(tf::make_sort_task(data.begin(), data.end()));
    executor.wait_for_all();
    REQUIRE(std::is_sorted(data.begin(), data.end()));
  }
}

TEST_CASE("ParallelSort.Async.1thread" * doctest::timeout(300)) {
  async(1);
}

TEST_CASE("ParallelSort.Async.2threads" * doctest::timeout(300)) {
  async(2);
}

TEST_CASE("ParallelSort.Async.3threads" * doctest::timeout(300)) {
  async(3);
}

TEST_CASE("ParallelSort.Async.4threads" * doctest::timeout(300)) {
  async(4);
}

// ----------------------------------------------------------------------------
// Parallel Sort with Dependent Async Tasks
// ----------------------------------------------------------------------------

void dependent_async(size_t W) {

  std::srand(static_cast<unsigned int>(time(NULL)));
  
  tf::Executor executor(W);
  std::vector<int> data;

  for(size_t n=0; n < 100000; n = (n ? n*10 : 1)) {
    
    data.resize(n);

    for(auto& d : data) {
      d = ::rand() % 1000 - 500;
    }
  
    executor.dependent_async(tf::make_sort_task(data.begin(), data.end()));
    executor.wait_for_all();
    REQUIRE(std::is_sorted(data.begin(), data.end()));
  }
}

TEST_CASE("ParallelSort.DependentAsync.1thread" * doctest::timeout(300)) {
  dependent_async(1);
}

TEST_CASE("ParallelSort.DependentAsync.2threads" * doctest::timeout(300)) {
  dependent_async(2);
}

TEST_CASE("ParallelSort.DependentAsync.3threads" * doctest::timeout(300)) {
  dependent_async(3);
}

TEST_CASE("ParallelSort.DependentAsync.4threads" * doctest::timeout(300)) {
  dependent_async(4);
}

// ----------------------------------------------------------------------------
// Parallel Sort with Silent Async Tasks
// ----------------------------------------------------------------------------

void silent_async(size_t W) {

  std::srand(static_cast<unsigned int>(time(NULL)));
  
  tf::Executor executor(W);
  std::vector<int> data;

  for(size_t n=0; n < 100000; n = (n ? n*10 : 1)) {
    
    data.resize(n);

    for(auto& d : data) {
      d = ::rand() % 1000 - 500;
    }
  
    executor.silent_async(tf::make_sort_task(data.begin(), data.end()));
    executor.wait_for_all();
    REQUIRE(std::is_sorted(data.begin(), data.end()));
  }
}

TEST_CASE("ParallelSort.SilentAsync.1thread" * doctest::timeout(300)) {
  silent_async(1);
}

TEST_CASE("ParallelSort.SilentAsync.2threads" * doctest::timeout(300)) {
  silent_async(2);
}

TEST_CASE("ParallelSort.SilentAsync.3threads" * doctest::timeout(300)) {
  silent_async(3);
}

TEST_CASE("ParallelSort.SilentAsync.4threads" * doctest::timeout(300)) {
  silent_async(4);
}

// ----------------------------------------------------------------------------
// Parallel Sort with Silent Dependent Async Tasks
// ----------------------------------------------------------------------------

void silent_dependent_async(size_t W) {

  std::srand(static_cast<unsigned int>(time(NULL)));
  
  tf::Executor executor(W);
  std::vector<int> data;

  for(size_t n=0; n < 100000; n = (n ? n*10 : 1)) {
    
    data.resize(n);

    for(auto& d : data) {
      d = ::rand() % 1000 - 500;
    }
  
    executor.silent_dependent_async(tf::make_sort_task(data.begin(), data.end()));
    executor.wait_for_all();
    REQUIRE(std::is_sorted(data.begin(), data.end()));
  }
}

TEST_CASE("ParallelSort.SilentDependentAsync.1thread" * doctest::timeout(300)) {
  silent_dependent_async(1);
}

TEST_CASE("ParallelSort.SilentDependentAsync.2threads" * doctest::timeout(300)) {
  silent_dependent_async(2);
}

TEST_CASE("ParallelSort.SilentDependentAsync.3threads" * doctest::timeout(300)) {
  silent_dependent_async(3);
}

TEST_CASE("ParallelSort.SilentDependentAsync.4threads" * doctest::timeout(300)) {
  silent_dependent_async(4);
}


// --------------------------------------------------------
// Testcase: BubbleSort
// --------------------------------------------------------
void bubble_sort(unsigned W) {

  tf::Executor executor(W);
  tf::Taskflow taskflow;
    
  std::vector<int> data;
 
  for(int end=1; end <= 1000; end *= 10) {
 
    taskflow.clear();
    data.resize(end);
 
    for(auto& d : data) d = ::rand()%100;
 
    auto gold = data;
    std::sort(gold.begin(), gold.end());
 
    std::atomic<bool> swapped;
 
    // init task
    auto init = taskflow.emplace([&swapped](){ swapped = false; });
    auto cond = taskflow.emplace([&swapped](){
      if(swapped) {
        swapped = false;
        return 0;
      }
      return 1;
    });
    auto stop = taskflow.emplace([](){});
 
    auto even_phase = taskflow.emplace([&](tf::Subflow& sf){
      for(size_t i=0; i<data.size(); i+=2) {
        sf.emplace([&, i](){
          if(i+1 < data.size() && data[i] > data[i+1]) {
            std::swap(data[i], data[i+1]);
            swapped = true;
          }
        });
      }
    });
 
    auto odd_phase = taskflow.emplace([&](tf::Subflow& sf) {
      for(size_t i=1; i<data.size(); i+=2) {
        sf.emplace([&, i](){
          if(i+1 < data.size() && data[i] > data[i+1]) {
            std::swap(data[i], data[i+1]);
            swapped = true;
          }
        });
      }
    });
 
    init.precede(even_phase).name("init");
    even_phase.precede(odd_phase).name("even-swap");
    odd_phase.precede(cond).name("odd-swap");
    cond.precede(even_phase, stop).name("cond");
 
    executor.run(taskflow).wait();
 
    REQUIRE(gold == data);
  }
}

TEST_CASE("BubbleSort.1thread" * doctest::timeout(300)) {
  bubble_sort(1);
}

TEST_CASE("BubbleSort.2threads" * doctest::timeout(300)) {
  bubble_sort(2);
}

TEST_CASE("BubbleSort.3threads" * doctest::timeout(300)) {
  bubble_sort(3);
}

TEST_CASE("BubbleSort.4threads" * doctest::timeout(300)) {
  bubble_sort(4);
}

TEST_CASE("BubbleSort.5threads" * doctest::timeout(300)) {
  bubble_sort(5);
}

TEST_CASE("BubbleSort.6threads" * doctest::timeout(300)) {
  bubble_sort(6);
}

TEST_CASE("BubbleSort.7threads" * doctest::timeout(300)) {
  bubble_sort(7);
}

TEST_CASE("BubbleSort.8threads" * doctest::timeout(300)) {
  bubble_sort(8);
}


// --------------------------------------------------------
// Testcase: SelectionSort
// --------------------------------------------------------

void selection_sort(unsigned W) {

  std::function<
    void(tf::Subflow& sf, std::vector<int>&, int, int, int&)
  > spawn;

  spawn = [&] (
    tf::Subflow& sf,
    std::vector<int>& data,
    int beg,
    int end,
    int& min
  ) mutable {

    if(!(beg < end)) {
      min = -1;
      return;
    }

    if(end - beg == 1) {
      min = beg;
      return;
    }

    int m = (beg + end + 1) / 2;

    auto minl = new int(-1);
    auto minr = new int(-1);

    auto SL = sf.emplace(
      [&spawn, &data, beg, m, l=minl] (tf::Subflow& sf2) mutable {
      spawn(sf2, data, beg, m, *l);
    }).name(std::string("[")
          + std::to_string(beg)
          + ':'
          + std::to_string(m)
          + ')');

    auto SR = sf.emplace(
      [&spawn, &data, m, end, r=minr] (tf::Subflow& sf2) mutable {
      spawn(sf2, data, m, end, *r);
    }).name(std::string("[")
          + std::to_string(m)
          + ':'
          + std::to_string(end)
          + ')');

    auto SM = sf.emplace([&data, &min, minl, minr] () {
      if(*minl == -1) {
        min = *minr;
      }
      else if(*minr == -1) {
        min = *minl;
        return;
      }
      else {
        min = data[*minl] < data[*minr] ? *minl : *minr;
      }
      delete minl;
      delete minr;
    }).name(std::string("merge [")
          + std::to_string(beg)
          + ':'
          + std::to_string(end) + ')');

    SM.succeed(SL, SR);
  };

  tf::Executor executor(W);
  tf::Taskflow taskflow;
  std::vector<int> data;

  for(int end=1; end <= 256; end <<= 1) {

    taskflow.clear();
    data.resize(end);

    for(auto& d : data) d = ::rand()%100;

    auto gold = data;
    std::sort(gold.begin(), gold.end());

    int beg = 0;
    int min = -1;

    auto start = taskflow.emplace([](){});

    auto argmin = taskflow.emplace(
      [&spawn, &data, &beg, end, &min](tf::Subflow& sf) mutable {
      spawn(sf, data, beg, end, min);
    }).name(std::string("[0")
          + ":"
          + std::to_string(end) + ")");

    auto putmin = taskflow.emplace([&](){
      std::swap(data[beg], data[min]);
      //std::cout << "select " << data[beg] << '\n';
      beg++;
      if(beg < end) {
        min = -1;
        return 0;
      }
      else return 1;
    });

    start.precede(argmin);
    argmin.precede(putmin);
    putmin.precede(argmin);

    executor.run(taskflow).wait();

    REQUIRE(gold == data);
    //std::exit(1);
  }
}

TEST_CASE("SelectionSort.1thread" * doctest::timeout(300)) {
  selection_sort(1);
}

TEST_CASE("SelectionSort.2threads" * doctest::timeout(300)) {
  selection_sort(2);
}

TEST_CASE("SelectionSort.3threads" * doctest::timeout(300)) {
  selection_sort(3);
}

TEST_CASE("SelectionSort.4threads" * doctest::timeout(300)) {
  selection_sort(4);
}

TEST_CASE("SelectionSort.5threads" * doctest::timeout(300)) {
  selection_sort(5);
}

TEST_CASE("SelectionSort.6threads" * doctest::timeout(300)) {
  selection_sort(6);
}

TEST_CASE("SelectionSort.7threads" * doctest::timeout(300)) {
  selection_sort(7);
}

TEST_CASE("SelectionSort.8threads" * doctest::timeout(300)) {
  selection_sort(8);
}

// --------------------------------------------------------
// Testcase: MergeSort
// --------------------------------------------------------

void merge_sort(unsigned W) {

  std::function<void(tf::Subflow& sf, std::vector<int>&, int, int)> spawn;

  spawn = [&] (tf::Subflow& sf, std::vector<int>& data, int beg, int end) mutable {

    if(!(beg < end) || end - beg == 1) {
      return;
    }

    if(end - beg <= 5) {
      std::sort(data.begin() + beg, data.begin() + end);
      return;
    }

    int m = (beg + end + 1) / 2;

    auto SL = sf.emplace([&spawn, &data, beg, m] (tf::Subflow& sf2) {
      spawn(sf2, data, beg, m);
    }).name(std::string("[")
          + std::to_string(beg)
          + ':'
          + std::to_string(m)
          + ')');

    auto SR = sf.emplace([&spawn, &data, m, end] (tf::Subflow& sf2) {
      spawn(sf2, data, m, end);
    }).name(std::string("[")
          + std::to_string(m)
          + ':'
          + std::to_string(end)
          + ')');

    auto SM = sf.emplace([&data, beg, end, m] () {
      std::vector<int> tmpl, tmpr;
      for(int i=beg; i<m; ++i) tmpl.push_back(data[i]);
      for(int i=m; i<end; ++i) tmpr.push_back(data[i]);

      // merge to data
      size_t i=0, j=0, k=beg;
      while(i<tmpl.size() && j<tmpr.size()) {
        data[k++] = (tmpl[i] < tmpr[j] ? tmpl[i++] : tmpr[j++]);
      }

      // remaining SL
      for(; i<tmpl.size(); ++i) data[k++] = tmpl[i];

      // remaining SR
      for(; j<tmpr.size(); ++j) data[k++] = tmpr[j];
    }).name(std::string("merge [")
          + std::to_string(beg)
          + ':'
          + std::to_string(end) + ')');

    SM.succeed(SL, SR);
  };

  tf::Executor executor(W);
  tf::Taskflow taskflow;
  std::vector<int> data;

  for(int end=10; end <= 10000; end *= 10) {

    taskflow.clear();
    data.resize(end);

    for(auto& d : data) d = ::rand()%100;

    auto gold = data;

    taskflow.emplace([&spawn, &data, end](tf::Subflow& sf){
      spawn(sf, data, 0, end);
    }).name(std::string("[0")
          + ":"
          + std::to_string(end) + ")");

    executor.run(taskflow).wait();

    std::sort(gold.begin(), gold.end());

    REQUIRE(gold == data);
  }
}

TEST_CASE("MergeSort.1thread" * doctest::timeout(300)) {
  merge_sort(1);
}

TEST_CASE("MergeSort.2threads" * doctest::timeout(300)) {
  merge_sort(2);
}

TEST_CASE("MergeSort.3threads" * doctest::timeout(300)) {
  merge_sort(3);
}

TEST_CASE("MergeSort.4threads" * doctest::timeout(300)) {
  merge_sort(4);
}

TEST_CASE("MergeSort.5threads" * doctest::timeout(300)) {
  merge_sort(5);
}

TEST_CASE("MergeSort.6threads" * doctest::timeout(300)) {
  merge_sort(6);
}

TEST_CASE("MergeSort.7threads" * doctest::timeout(300)) {
  merge_sort(7);
}

TEST_CASE("MergeSort.8threads" * doctest::timeout(300)) {
  merge_sort(8);
}

// --------------------------------------------------------
// Testcase: QuickSort
// --------------------------------------------------------
void quick_sort(unsigned W) {

  using itr_t = std::vector<int>::iterator;

  std::function<void(tf::Subflow& sf, std::vector<int>&, itr_t, itr_t)> spawn;

  spawn = [&] (
    tf::Subflow& sf,
    std::vector<int>& data,
    itr_t beg,
    itr_t end
  ) mutable {

    if(!(beg < end) || std::distance(beg, end) == 1) {
      return;
    }

    if(std::distance(beg, end) <= 5) {
      std::sort(beg, end);
      return;
    }

    auto pvt = beg + std::distance(beg, end) / 2;

    std::iter_swap(pvt, end-1);

    pvt = std::partition(beg, end-1, [end] (int item) {
      return item < *(end - 1);
    });

    std::iter_swap(pvt, end-1);

    sf.emplace([&spawn, &data, beg, pvt] (tf::Subflow& sf2) {
      spawn(sf2, data, beg, pvt);
    }).name(std::string("[")
          + std::to_string(beg-data.begin())
          + ':'
          + std::to_string(pvt-data.begin())
          + ')');

    sf.emplace([&spawn, &data, pvt, end] (tf::Subflow& sf2) {
      spawn(sf2, data, pvt+1, end);
    }).name(std::string("[")
          + std::to_string(pvt-data.begin())
          + ':'
          + std::to_string(end-data.begin())
          + ')');
  };

  tf::Executor executor(W);
  tf::Taskflow taskflow;
  std::vector<int> data;

  for(size_t end=1; end <= 10000; end *= 10) {

    taskflow.clear();
    data.resize(end);

    for(auto& d : data) d = ::rand()%100;

    auto gold = data;

    taskflow.emplace([&spawn, &data](tf::Subflow& sf){
      spawn(sf, data, data.begin(), data.end());
    }).name(std::string("[0")
          + ":"
          + std::to_string(end) + ")");

    executor.run(taskflow).wait();

    std::sort(gold.begin(), gold.end());

    REQUIRE(gold == data);
  }
  
}

TEST_CASE("QuickSort.1thread" * doctest::timeout(300)) {
  quick_sort(1);
}

TEST_CASE("QuickSort.2threads" * doctest::timeout(300)) {
  quick_sort(2);
}

TEST_CASE("QuickSort.3threads" * doctest::timeout(300)) {
  quick_sort(3);
}

TEST_CASE("QuickSort.4threads" * doctest::timeout(300)) {
  quick_sort(4);
}

TEST_CASE("QuickSort.5threads" * doctest::timeout(300)) {
  quick_sort(5);
}

TEST_CASE("QuickSort.6threads" * doctest::timeout(300)) {
  quick_sort(6);
}

TEST_CASE("QuickSort.7threads" * doctest::timeout(300)) {
  quick_sort(7);
}

TEST_CASE("QuickSort.8threads" * doctest::timeout(300)) {
  quick_sort(8);
}



//// ----------------------------------------------------------------------------
//// Exception
//// ----------------------------------------------------------------------------
//
//void parallel_sort_exception(unsigned W) {
//
//  tf::Taskflow taskflow;
//  tf::Executor executor(W);
//
//  std::vector<int> data(1000000);
//
//  // for_each
//  taskflow.sort(data.begin(), data.end(), [](int a, int b){
//    throw std::runtime_error("x");
//    return a < b;
//  });
//  REQUIRE_THROWS_WITH_AS(executor.run(taskflow).get(), "x", std::runtime_error);
//}
//
//TEST_CASE("ParallelSort.Exception.1thread") {
//  parallel_sort_exception(1);
//}
//
//TEST_CASE("ParallelSort.Exception.2threads") {
//  parallel_sort_exception(2);
//}
//
//TEST_CASE("ParallelSort.Exception.3threads") {
//  parallel_sort_exception(3);
//}
//
//TEST_CASE("ParallelSort.Exception.4threads") {
//  parallel_sort_exception(4);
//}


