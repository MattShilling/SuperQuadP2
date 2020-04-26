#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>

#include "test_rig.h"

// These will be default if not specified in the command line.
#define NUMNODES 256
#define NUMT 8
#define N 4

// This defines the boundaries of our superquadric.
#define XMIN -1.0
#define XMAX 1.0
#define YMIN -1.0
#define YMAX 1.0

double CalcSpeedup(double test_two, double test_one) {
  return test_two / test_one;
}

float CalcFp(double speedup, double threads) {
  return (threads / (threads - 1.0)) * (1. - (1. / speedup));
}

float Height(int, int, int);

struct test_mem {
  int numnodes;
  float volume;
};

void calc_superquad(std::shared_ptr<void> mem) {

  // Get the number of nodes to be used from memory
  test_mem *data = static_cast<test_mem *>(mem.get());
  const int numnodes = data->numnodes;

  // The area of a single full-sized tile: dA. NOTE: There should be no
  // "numnodes-1" term because the total number of nodes is numnodes^2, not
  // (numnodes-1)^2.
  const float dA =
      (XMAX - XMIN) * (YMAX - YMIN) / pow(static_cast<float>(numnodes), 2.0);

  // Use an OpenMP for loop and a reduction to find weighted heights.
  float total_height = 0.0;

#pragma omp parallel for collapse(2), default(none), reduction(+ : total_height)
  for (int iv = 0; iv < numnodes; iv++) {
    for (int iu = 0; iu < numnodes; iu++) {
      float weight = 0.0;
      // Logic to decide if we are on a side of the superquadric or not. If both
      // are true then we are on a corner.
      const bool side_iv = (iv == numnodes - 1) || (iv == 0);
      const bool side_iu = (iu == numnodes - 1) || (iu == 0);
      // XOR
      if (!side_iv != !side_iu) {
        // We are on a side so weight height by 0.5.
        weight = 0.5;
      } else if (side_iv && side_iu) {
        // We are on a corner so weight height by 0.25.
        weight = 0.25;
      } else {
        // We have a regular cell so the normal height * dA will apply.
        weight = 1.0;
      }

      total_height += Height(iu, iv, numnodes);
    }
  }

  // Multiply weighted heights by the dA for each cell, multiply by two to
  // approximate total volume.
  data->volume = total_height * dA * 2.0;

  std::cout << "Volume = " << data->volume << std::endl;
}

int main(int argc, char *argv[]) {
  // Taking in some command line arguments to control the program.
  int threads = NUMT;
  int numnodes = NUMNODES;
  bool just_volume = false;
  std::string file;

  if (argc >= 2) {
    threads = std::stoi(std::string(argv[1]));
  }

  if (argc >= 3) {
    numnodes = std::stoi(std::string(argv[2]));
  }

  if (argc >= 4) {
    std::string temp(argv[3]);
    if (temp == "-vol") {
      if (argc >= 5) {
        file = std::string(argv[4]);
        just_volume = true;
      }
    } else {
      file = std::string(argv[3]);
    }
  }

  std::shared_ptr<test_mem> mem = std::make_shared<test_mem>();

  mem->numnodes = numnodes;

  TestRig superquad(mem, calc_superquad,
                    std::function<void(std::shared_ptr<void> mem)>());

  std::cout << 1.00 / 32768.00 << std::endl;

  superquad.Init(threads);
  superquad.Run(static_cast<double>(numnodes * numnodes));

  if (!file.empty()) {
    float superquad_two = 0.0;
    float superquad_one = 0.0;
    double speedup = 0.0;
    float fp = 0.0;

    if (!just_volume) {
      // We only need to run a test with one thread if we want to find the
      // speedup and parallel fraction.
      superquad_two = superquad.MaxMegaMults();

      superquad.Init(1);
      superquad.Run(static_cast<double>(numnodes * numnodes));
      superquad_one = superquad.MaxMegaMults();

      speedup = CalcSpeedup(superquad_two, superquad_one);
      fp = CalcFp(speedup, threads);
    }

    std::ofstream outfile;
    outfile.open(file, std::ios_base::app);

    // Setting the proper precision.
    outfile << std::fixed;
    outfile << std::setprecision(3);
    std::cout << std::fixed;
    std::cout << std::setprecision(3);

    if (!just_volume) {
      // We care about all the data.
      std::cout << "Writing: \n" << threads << '\t' << numnodes << '\t'
                << superquad_two << '\t' << speedup << '\t' << fp << "\nto "
                << file << std::endl;

      outfile << threads << '\t' << numnodes << '\t' << mem->volume << '\t'
              << superquad_two << '\t' << speedup << '\t' << fp << std::endl;
    } else {
      // We only care about the volume results.
      std::cout << "Writing: \n" << threads << '\t' << numnodes << '\t'
                << mem->volume << "\nto " << file << std::endl;

      outfile << threads << '\t' << numnodes << '\t' << mem->volume
              << std::endl;
    }
    outfile.close();
  }
}

// iu,iv = 0 .. NUMNODES-1
float Height(int iu, int iv, int numnodes) {

  double x = -1.00 + 2.00 * ((double)iu / (double)(numnodes - 1)); // -1. to +1.
  double y = -1.00 + 2.00 * ((double)iv / (double)(numnodes - 1)); // -1. to +1.

  double xn = pow(fabs(x), (double)N);
  double yn = pow(fabs(y), (double)N);
  double r = 1.00 - xn - yn;
  if (r < 0.00)
    return 0.00;
  float height = pow(r, 1.00 / (float)N);
  return height;
}