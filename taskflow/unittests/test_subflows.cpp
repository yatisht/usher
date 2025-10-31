#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/reduce.hpp>

// --------------------------------------------------------
// Testcase: JoinedSubflow
// --------------------------------------------------------

void joined_subflow(unsigned W) {

  using namespace std::literals::chrono_literals;

  SUBCASE("Trivial") {
    tf::Executor executor(W);
    tf::Taskflow tf;

    // empty flow with future
    tf::Task subflow3, subflow3_;
    //std::future<int> fu3, fu3_;
    std::atomic<int> fu3v{0}, fu3v_{0};

    // empty flow
    auto subflow1 = tf.emplace([&] (tf::Subflow& fb) {
      fu3v++;
      fb.join();
    }).name("subflow1");

    // nested empty flow
    auto subflow2 = tf.emplace([&] (tf::Subflow& fb) {
      fu3v++;
      fb.emplace([&] (tf::Subflow& fb2) {
        fu3v++;
        fb2.emplace( [&] (tf::Subflow& fb3) {
          fu3v++;
          fb3.join();
        }).name("subflow2_1_1");
      }).name("subflow2_1");
    }).name("subflow2");

    subflow3 = tf.emplace([&] (tf::Subflow& fb) {

      REQUIRE(fu3v == 4);

      fu3v++;
      fu3v_++;

      subflow3_ = fb.emplace([&] (tf::Subflow& fb2) {
        REQUIRE(fu3v_ == 3);
        fu3v++;
        fu3v_++;
        //return 200;
        fb2.join();
      });
      subflow3_.name("subflow3_");

      // hereafter we use 100us to avoid dangling reference ...
      auto s1 = fb.emplace([&] () {
        fu3v_++;
        fu3v++;
      }).name("s1");

      auto s2 = fb.emplace([&] () {
        fu3v_++;
        fu3v++;
      }).name("s2");

      auto s3 = fb.emplace([&] () {
        fu3v++;
        REQUIRE(fu3v_ == 4);
      }).name("s3");

      s1.precede(subflow3_);
      s2.precede(subflow3_);
      subflow3_.precede(s3);

      REQUIRE(fu3v_ == 1);

      //return 100;
    });
    subflow3.name("subflow3");

    // empty flow to test future
    auto subflow4 = tf.emplace([&] () {
      fu3v++;
    }).name("subflow4");

    subflow1.precede(subflow2);
    subflow2.precede(subflow3);
    subflow3.precede(subflow4);

    executor.run(tf).get();
    // End of for loop
  }

  // Mixed intra- and inter- operations
  SUBCASE("Complex") {
    tf::Executor executor(W);
    tf::Taskflow tf;

    std::vector<int> data;
    int sum {0};

    auto A = tf.emplace([&data] () {
      for(int i=0; i<10; ++i) {
        data.push_back(1);
      }
    });

    std::atomic<size_t> count {0};

    auto B = tf.emplace([&count, &data, &sum](tf::Subflow& fb){

      //auto [src, tgt] = fb.reduce(data.begin(), data.end(), sum, std::plus<int>());
      auto task = fb.reduce(data.begin(), data.end(), sum, std::plus<int>());

      fb.emplace([&sum] () { REQUIRE(sum == 0); }).precede(task);

      task.precede(fb.emplace([&sum] () { REQUIRE(sum == 10); }));

      for(size_t i=0; i<10; i ++){
        ++count;
      }

      auto n = fb.emplace([&count](tf::Subflow& fb2){

        REQUIRE(count == 20);
        ++count;

        auto prev = fb2.emplace([&count](){
          REQUIRE(count == 21);
          ++count;
        });

        for(size_t i=0; i<10; i++){
          auto next = fb2.emplace([&count, i](){
            REQUIRE(count == 22+i);
            ++count;
          });
          prev.precede(next);
          prev = next;
        }
      });

      for(size_t i=0; i<10; i++){
        fb.emplace([&count](){ ++count; }).precede(n);
      }
    });

    A.precede(B);

    executor.run(tf).get();
    REQUIRE(count == 32);
    REQUIRE(sum == 10);

  }
}

TEST_CASE("JoinedSubflow.1thread" * doctest::timeout(300)){
  joined_subflow(1);
}

TEST_CASE("JoinedSubflow.2threads" * doctest::timeout(300)){
  joined_subflow(2);
}

TEST_CASE("JoinedSubflow.3threads" * doctest::timeout(300)){
  joined_subflow(3);
}

TEST_CASE("JoinedSubflow.4threads" * doctest::timeout(300)){
  joined_subflow(4);
}

TEST_CASE("JoinedSubflow.5threads" * doctest::timeout(300)){
  joined_subflow(5);
}

TEST_CASE("JoinedSubflow.6threads" * doctest::timeout(300)){
  joined_subflow(6);
}

TEST_CASE("JoinedSubflow.7threads" * doctest::timeout(300)){
  joined_subflow(7);
}

TEST_CASE("JoinedSubflow.8threads" * doctest::timeout(300)){
  joined_subflow(8);
}

//// --------------------------------------------------------
//// Testcase: DetachedSubflow
//// --------------------------------------------------------
//
//void detached_subflow(unsigned W) {
//
//  using namespace std::literals::chrono_literals;
//
//  SUBCASE("Trivial") {
//    tf::Executor executor(W);
//    tf::Taskflow tf;
//
//    // empty flow with future
//    tf::Task subflow3, subflow3_;
//    std::atomic<int> fu3v{0}, fu3v_{0};
//
//    // empty flow
//    auto subflow1 = tf.emplace([&] (tf::Subflow& fb) {
//      fu3v++;
//      fb.detach();
//    }).name("subflow1");
//
//    // nested empty flow
//    auto subflow2 = tf.emplace([&] (tf::Subflow& fb) {
//      fu3v++;
//      fb.emplace([&] (tf::Subflow& fb2) {
//        fu3v++;
//        fb2.emplace( [&] (tf::Subflow& fb3) {
//          fu3v++;
//          fb3.join();
//        }).name("subflow2_1_1");
//        fb2.detach();
//      }).name("subflow2_1");
//      fb.detach();
//    }).name("subflow2");
//
//    subflow3 = tf.emplace([&] (tf::Subflow& fb) {
//
//      REQUIRE((fu3v >= 2 && fu3v <= 4));
//
//      fu3v++;
//      fu3v_++;
//
//      subflow3_ = fb.emplace([&] (tf::Subflow& fb2) {
//        REQUIRE(fu3v_ == 3);
//        fu3v++;
//        fu3v_++;
//        fb2.join();
//      });
//      subflow3_.name("subflow3_");
//
//      // hereafter we use 100us to avoid dangling reference ...
//      auto s1 = fb.emplace([&] () {
//        fu3v_++;
//        fu3v++;
//      }).name("s1");
//
//      auto s2 = fb.emplace([&] () {
//        fu3v_++;
//        fu3v++;
//      }).name("s2");
//
//      auto s3 = fb.emplace([&] () {
//        fu3v++;
//        REQUIRE(fu3v_ == 4);
//      }).name("s3");
//
//      s1.precede(subflow3_);
//      s2.precede(subflow3_);
//      subflow3_.precede(s3);
//
//      REQUIRE(fu3v_ == 1);
//
//      fb.detach();
//
//      //return 100;
//    });
//    subflow3.name("subflow3");
//
//    // empty flow to test future
//    auto subflow4 = tf.emplace([&] () {
//      REQUIRE((fu3v >= 3 && fu3v <= 9));
//      fu3v++;
//    }).name("subflow4");
//
//    subflow1.precede(subflow2);
//    subflow2.precede(subflow3);
//    subflow3.precede(subflow4);
//
//    executor.run(tf).get();
//
//    REQUIRE(fu3v  == 10);
//    REQUIRE(fu3v_ == 4);
//
//  }
//}
//
//TEST_CASE("DetachedSubflow.1thread" * doctest::timeout(300)) {
//  detached_subflow(1);
//}
//
//TEST_CASE("DetachedSubflow.2threads" * doctest::timeout(300)) {
//  detached_subflow(2);
//}
//
//TEST_CASE("DetachedSubflow.3threads" * doctest::timeout(300)) {
//  detached_subflow(3);
//}
//
//TEST_CASE("DetachedSubflow.4threads" * doctest::timeout(300)) {
//  detached_subflow(4);
//}
//
//TEST_CASE("DetachedSubflow.5threads" * doctest::timeout(300)) {
//  detached_subflow(5);
//}
//
//TEST_CASE("DetachedSubflow.6threads" * doctest::timeout(300)) {
//  detached_subflow(6);
//}
//
//TEST_CASE("DetachedSubflow.7threads" * doctest::timeout(300)) {
//  detached_subflow(7);
//}
//
//TEST_CASE("DetachedSubflow.8threads" * doctest::timeout(300)) {
//  detached_subflow(8);
//}
//
//
//// --------------------------------------------------------
//// Testcase: TreeSubflow
//// --------------------------------------------------------
//void detach_spawn(const int max_depth, std::atomic<int>& counter, int depth, tf::Subflow& subflow)  {
//  if(depth < max_depth) {
//    counter.fetch_add(1, std::memory_order_relaxed);
//    subflow.emplace([&, max_depth, depth=depth+1](tf::Subflow& sfl){
//      detach_spawn(max_depth, counter, depth, sfl); }
//    );
//    subflow.emplace([&, max_depth, depth=depth+1](tf::Subflow& sfr){
//      detach_spawn(max_depth, counter, depth, sfr); }
//    );
//    subflow.detach();
//  }
//}
//
//void join_spawn(const int max_depth, std::atomic<int>& counter, int depth, tf::Subflow& subflow)  {
//  if(depth < max_depth) {
//    counter.fetch_add(1, std::memory_order_relaxed);
//    subflow.emplace([&, max_depth, depth=depth+1](tf::Subflow& sfl){
//      join_spawn(max_depth, counter, depth, sfl); }
//    );
//    subflow.emplace([&, max_depth, depth=depth+1](tf::Subflow& sfr){
//      join_spawn(max_depth, counter, depth, sfr); }
//    );
//  }
//}
//
//void mix_spawn(
//  const int max_depth, std::atomic<int>& counter, int depth, tf::Subflow& subflow
//) {
//
//  if(depth < max_depth) {
//    auto ret = counter.fetch_add(1, std::memory_order_relaxed);
//    subflow.emplace([&, max_depth, depth=depth+1](tf::Subflow& sfl){
//      mix_spawn(max_depth, counter, depth, sfl); }
//    ).name(std::string("left") + std::to_string(ret%2));
//    subflow.emplace([&, max_depth, depth=depth+1](tf::Subflow& sfr){
//      mix_spawn(max_depth, counter, depth, sfr); }
//    ).name(std::string("right") + std::to_string(ret%2));
//    if(ret % 2) {
//      subflow.detach();
//    }
//  }
//}
//
//TEST_CASE("TreeSubflow" * doctest::timeout(300)) {
//
//  SUBCASE("AllDetach") {
//    constexpr int max_depth {10};
//    for(int W=1; W<=4; W++) {
//      std::atomic<int> counter {0};
//      tf::Taskflow tf;
//      tf.emplace([&](tf::Subflow& subflow){
//        detach_spawn(max_depth, counter, 0, subflow);
//      });
//
//      tf::Executor executor(W);
//      executor.run(tf).get();
//      REQUIRE(counter == (1<<max_depth) - 1);
//    }
//  }
//
//
//  SUBCASE("AllJoin") {
//    constexpr int max_depth {10};
//    for(int W=1; W<=4; W++) {
//      std::atomic<int> counter {0};
//      tf::Taskflow tf;
//      tf.emplace([&](tf::Subflow& subflow){
//        join_spawn(max_depth, counter, 0, subflow);
//      });
//      tf::Executor executor(W);
//      executor.run(tf).get();
//      REQUIRE(counter == (1<<max_depth) - 1);
//    }
//  }
//
//  SUBCASE("Mix") {
//    constexpr int max_depth {10};
//    for(int W=1; W<=4; W++) {
//      std::atomic<int> counter {0};
//      tf::Taskflow tf;
//      tf.emplace([&](tf::Subflow& subflow){
//        mix_spawn(max_depth, counter, 0, subflow);
//      }).name("top task");
//
//      tf::Executor executor(W);
//      executor.run(tf).get();
//      REQUIRE(counter == (1<<max_depth) - 1);
//    }
//  }
//}

// --------------------------------------------------------
// Testcase: FibSubflow
// --------------------------------------------------------
int fibonacci_spawn(int n, tf::Subflow& sbf) {
  if (n < 2) return n;
  int res1, res2;
  sbf.emplace([&res1, n] (tf::Subflow& sbfl) { res1 = fibonacci_spawn(n - 1, sbfl); } );
  sbf.emplace([&res2, n] (tf::Subflow& sbfr) { res2 = fibonacci_spawn(n - 2, sbfr); } );
  REQUIRE(sbf.joinable() == true);
  sbf.join();
  REQUIRE(sbf.joinable() == false);
  return res1 + res2;
}

void fibonacci(size_t W) {

  int N = 20;
  int res = -1;  // result

  tf::Executor executor(W);
  tf::Taskflow taskflow;

  taskflow.emplace([&res, N] (tf::Subflow& sbf) {
    res = fibonacci_spawn(N, sbf);
  });

  executor.run(taskflow).wait();

  REQUIRE(res == 6765);
}

TEST_CASE("FibSubflow.1thread" * doctest::timeout(300)) {
  fibonacci(1);
}

TEST_CASE("FibSubflow.2threads" * doctest::timeout(300)) {
  fibonacci(2);
}

TEST_CASE("FibSubflow.4threads" * doctest::timeout(300)) {
  fibonacci(4);
}

TEST_CASE("FibSubflow.5threads" * doctest::timeout(300)) {
  fibonacci(5);
}

TEST_CASE("FibSubflow.6threads" * doctest::timeout(300)) {
  fibonacci(6);
}

TEST_CASE("FibSubflow.7threads" * doctest::timeout(300)) {
  fibonacci(7);
}

TEST_CASE("FibSubflow.8threads" * doctest::timeout(300)) {
  fibonacci(8);
}

// ----------------------------------------------------------------------------
// multiple subflow runs
// ----------------------------------------------------------------------------
void multiple_subflow_runs(unsigned W) {

  tf::Executor executor(W);
  tf::Taskflow taskflow;

  std::atomic<size_t> count {0};

  auto A = taskflow.emplace([&](){ count ++; });
  auto B = taskflow.emplace([&](tf::Subflow& subflow){
    count ++;
    auto B1 = subflow.emplace([&](){ count++; });
    auto B2 = subflow.emplace([&](){ count++; });
    auto B3 = subflow.emplace([&](){ count++; });
    B1.precede(B3); B2.precede(B3);
  });
  auto C = taskflow.emplace([&](){ count ++; });
  auto D = taskflow.emplace([&](){ count ++; });

  A.precede(B, C);
  B.precede(D);
  C.precede(D);

  std::list<tf::Future<void>> fu_list;
  for(size_t i=0; i<500; i++) {
    if(i == 499) {
      executor.run(taskflow).get();   // Synchronize the first 500 runs
      executor.run_n(taskflow, 500);  // Run 500 times more
    }
    else if(i % 2) {
      fu_list.push_back(executor.run(taskflow));
    }
    else {
      fu_list.push_back(executor.run(taskflow, [&, i=i](){
        REQUIRE(count == (i+1)*7); })
      );
    }
  }

  executor.wait_for_all();

  for(auto& fu: fu_list) {
    REQUIRE(fu.valid());
    REQUIRE(fu.wait_for(std::chrono::seconds(1)) == std::future_status::ready);
  }

  REQUIRE(count == 7000);
}

TEST_CASE("MultipleSubflowRuns.1thread" * doctest::timeout(300)) {
  multiple_subflow_runs(1);
}

TEST_CASE("MultipleSubflowRuns.2threads" * doctest::timeout(300)) {
  multiple_subflow_runs(2);
}

TEST_CASE("MultipleSubflowRuns.3threads" * doctest::timeout(300)) {
  multiple_subflow_runs(3);
}

TEST_CASE("MultipleSubflowRuns.4threads" * doctest::timeout(300)) {
  multiple_subflow_runs(4);
}

TEST_CASE("MultipleSubflowRuns.4threads" * doctest::timeout(300)) {
  multiple_subflow_runs(4);
}

TEST_CASE("MultipleSubflowRuns.5threads" * doctest::timeout(300)) {
  multiple_subflow_runs(5);
}

TEST_CASE("MultipleSubflowRuns.6threads" * doctest::timeout(300)) {
  multiple_subflow_runs(6);
}

TEST_CASE("MultipleSubflowRuns.7threads" * doctest::timeout(300)) {
  multiple_subflow_runs(7);
}

TEST_CASE("MultipleSubflowRuns.8threads" * doctest::timeout(300)) {
  multiple_subflow_runs(8);
}

// ----------------------------------------------------------------------------
// Multiple subflow runs with change
// ----------------------------------------------------------------------------

void multiple_subflow_runs_with_changed_taskflow(unsigned W) {

  tf::Executor executor(W);
  tf::Taskflow taskflow;

  std::atomic<size_t> count {0};

  auto A = taskflow.emplace([&](){ count ++; });
  auto B = taskflow.emplace([&](tf::Subflow& subflow){
    count ++;
    auto B1 = subflow.emplace([&](){ count++; });
    auto B2 = subflow.emplace([&](){ count++; });
    auto B3 = subflow.emplace([&](){ count++; });
    B1.precede(B3); B2.precede(B3);
  });
  auto C = taskflow.emplace([&](){ count ++; });
  auto D = taskflow.emplace([&](){ count ++; });

  A.precede(B, C);
  B.precede(D);
  C.precede(D);

  executor.run_n(taskflow, 10).get();
  REQUIRE(count == 70);

  auto E = taskflow.emplace([](){});
  D.precede(E);
  executor.run_n(taskflow, 10).get();
  REQUIRE(count == 140);

  auto F = taskflow.emplace([](){});
  E.precede(F);
  executor.run_n(taskflow, 10);
  executor.wait_for_all();
  REQUIRE(count == 210);

}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.1thread" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(1);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.2threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(2);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.3threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(3);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.4threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(4);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.4threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(4);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.5threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(5);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.6threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(6);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.7threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(7);
}

TEST_CASE("MultipleSubflowRuns.ChangedTaskflow.8threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_changed_taskflow(8);
}

// ----------------------------------------------------------------------------
// multiple_subflow_runs_with_predicate
// ----------------------------------------------------------------------------

void multiple_subflow_runs_with_predicate(unsigned W) {

  tf::Executor executor(W);
  tf::Taskflow taskflow;
  
  std::atomic<size_t> count {0};
  auto A = taskflow.emplace([&](){ count ++; });
  auto B = taskflow.emplace([&](tf::Subflow& subflow){
    count ++;
    auto B1 = subflow.emplace([&](){ count++; });
    auto B2 = subflow.emplace([&](){ count++; });
    auto B3 = subflow.emplace([&](){ count++; });
    B1.precede(B3); B2.precede(B3);
  });
  auto C = taskflow.emplace([&](){ count ++; });
  auto D = taskflow.emplace([&](){ count ++; });

  A.precede(B, C);
  B.precede(D);
  C.precede(D);

  executor.run_until(taskflow, [run=10]() mutable { return run-- == 0; },
    [&](){
      REQUIRE(count == 70);
      count = 0;
    }
  ).get();


  executor.run_until(taskflow, [run=10]() mutable { return run-- == 0; },
    [&](){
      REQUIRE(count == 70);
      count = 0;
  });

  executor.run_until(taskflow, [run=10]() mutable { return run-- == 0; },
    [&](){
      REQUIRE(count == 70);
      count = 0;
    }
  ).get();
}

TEST_CASE("MultipleSubflowRuns.Predicate.1thread" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(1);
}

TEST_CASE("MultipleSubflowRuns.Predicate.2threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(2);
}

TEST_CASE("MultipleSubflowRuns.Predicate.3threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(3);
}

TEST_CASE("MultipleSubflowRuns.Predicate.4threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(4);
}

TEST_CASE("MultipleSubflowRuns.Predicate.4threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(4);
}

TEST_CASE("MultipleSubflowRuns.Predicate.5threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(5);
}

TEST_CASE("MultipleSubflowRuns.Predicate.6threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(6);
}

TEST_CASE("MultipleSubflowRuns.Predicate.7threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(7);
}

TEST_CASE("MultipleSubflowRuns.Predicate.8threads" * doctest::timeout(300)) {
  multiple_subflow_runs_with_predicate(8);
}

// ----------------------------------------------------------------------------
// subflow state test
// ----------------------------------------------------------------------------

void bit_state(unsigned W) {
  tf::Executor executor(W);
  tf::Taskflow taskflow;

  auto init = taskflow.emplace([](){});

  auto task = taskflow.emplace([](tf::Subflow& sf){
    // each newly spawned subflow should have clean status
    REQUIRE(sf.joinable());
    REQUIRE(sf.retain() == false);
    sf.join();
    sf.retain(true);
  });

  auto cond = taskflow.emplace([i=0]() mutable {
    return (i++ < 100) ? 0 : 1;
  });

  init.precede(task);
  task.precede(cond);
  cond.precede(task);

  executor.run(taskflow).wait();
}

TEST_CASE("Subflow.BitState.1thread") {
  bit_state(1);
}

TEST_CASE("Subflow.BitState.2threads") {
  bit_state(2);
}

TEST_CASE("Subflow.BitState.3threads") {
  bit_state(3);
}

TEST_CASE("Subflow.BitState.4threads") {
  bit_state(4);
}







