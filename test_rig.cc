#include "test_rig.h"

#include <algorithm>
#include <cmath>
#include <omp.h>
#include <stdio.h>

void TestRig::Init(int num_threads) {
  if (mp_good_) {
    // Set the number of threads we want to use.
    omp_set_num_threads(num_threads);
    fprintf(stderr, "Using %d threads\n", num_threads);

    // Initialize test memory.
    if(test_init_) {
      test_init_(mem_);
    }

    // Clear our records.
    mega_mults_.clear();
  } else {
    printf("Init error: OpenMP not supported!\n");
  }
}

void TestRig::Run(double sz) {
  if (mp_good_) {
    // Get the starting time for our test.
    double start = omp_get_wtime();
    // Run. That. Test!
    test_run_(mem_);
    // Get the ending time for our test.
    double stop = omp_get_wtime();
    // Calculate the multiplications per second we accomplished.
    double mults = sz / (stop - start);
    // Convert into megamults.
    double mega_mults = mults / 1000000.00;
    // Add results to our records.
    mega_mults_.push_back(mega_mults);
  } else {
    printf("Run error: OpenMP not supported!\n");
  }
}

double TestRig::MaxMegaMults() {
  return *std::max_element(mega_mults_.begin(), mega_mults_.end());
}

double TestRig::MinMegaMults() {
  return *std::min_element(mega_mults_.begin(), mega_mults_.end());
}

bool TestRig::CheckOpenMP() {
#ifndef _OPENMP
  fprintf(stderr, "OpenMP is not supported!\n");
  return false;
#endif
  return true;
}