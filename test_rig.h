#ifndef TEST_RIG_H_
#define TEST_RIG_H_

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

class TestRig {
public:
  TestRig(std::shared_ptr<void> mem,
          std::function<void(std::shared_ptr<void> mem)> run,
          std::function<void(std::shared_ptr<void> mem)> init)
      : mem_(mem), test_init_(init), test_run_(run), mp_good_(CheckOpenMP()){};
  void Init(int num_threads);
  void Run(double sz);

  double MaxMegaMults();
  double MinMegaMults();

private:
  bool CheckOpenMP();
  std::shared_ptr<void> mem_;
  std::function<void(std::shared_ptr<void> mem)> test_init_;
  std::function<void(std::shared_ptr<void> mem)> test_run_;
  bool mp_good_;
  std::vector<double> mega_mults_;
};

#endif // TEST_RIG_H_