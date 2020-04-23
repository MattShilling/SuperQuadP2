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

// Setting the number of threads.
#ifndef NUMT
#define NUMT 8
#endif

// Setting the number of trials in the monte carlo simulation.
#ifndef NUMTRIALS
#define NUMTRIALS 1000000
#endif

// How many tries to discover the maximum performance.
#ifndef NUMTRIES
#define NUMTRIES 10
#endif

// ranges for the random numbers:
const float XCMIN = -1.0;
const float XCMAX = 1.0;
const float YCMIN = 0.0;
const float YCMAX = 2.0;
const float RMIN = 0.5;
const float RMAX = 2.0;

const float k_tn = tan((M_PI / 180.) * 30.);

// Test memory.
struct test_mem {
  float *xcs;
  float *ycs;
  float *rs;
  int trials;
  bool initialized;
  float prob;
};

void mc_run(std::shared_ptr<void> mem) {
  // Get our test memory.
  test_mem *data = static_cast<test_mem *>(mem.get());
  int hits = 0;
#pragma omp parallel for default(none) shared(data) reduction(+ : hits)
  for (int n = 0; n < data->trials; n++) {
    // randomize the location and radius of the circle:
    float xc = data->xcs[n];
    float yc = data->ycs[n];
    float r = data->rs[n];

    // solve for the intersection using the quadratic formula:
    float a = 1. + k_tn * k_tn;
    float b = -2. * (xc + yc * k_tn);
    float c = xc * xc + yc * yc - r * r;
    float d = b * b - 4. * a * c;

    // We missed the circle.
    if (d < 0.0) {
      continue;
    }

    // hits the circle:
    // get the first intersection:
    d = sqrt(d);
    float t1 = (-b + d) / (2. * a); // time to intersect the circle
    float t2 = (-b - d) / (2. * a); // time to intersect the circle
    float tmin = t1 < t2 ? t1 : t2; // only care about the first intersection

    // Laser pointer is inside circle.
    if (tmin < 0.0) {
      continue;
    }

    // where does it intersect the circle?
    float xcir = tmin;
    float ycir = tmin * k_tn;

    // get the unitized normal vector at the point of intersection:
    float nx = xcir - xc;
    float ny = ycir - yc;
    float nxy = sqrt(nx * nx + ny * ny);
    nx /= nxy; // unit vector
    ny /= nxy; // unit vector

    // get the unitized incoming vector:
    float inx = xcir - 0.;
    float iny = ycir - 0.;
    float in = sqrt(inx * inx + iny * iny);
    inx /= in; // unit vector
    iny /= in; // unit vector

    // get the outgoing (bounced) vector:
    float dot = inx * nx + iny * ny;
    float outx =
        inx - 2. * nx * dot; // angle of reflection = angle of incidence`
    float outy =
        iny - 2. * ny * dot; // angle of reflection = angle of incidence`

    // find out if it hits the infinite plate:
    float tt = (0. - ycir) / outy;

    // The beam hits the infinite plate.
    if (tt >= 0.0) {
      hits++;
    }
  }
  // std::cout << numHits << std::endl;
  data->prob = (float)hits / (float)data->trials;
}

void mc_init(std::shared_ptr<void> mem) {
  test_mem *data = static_cast<test_mem *>(mem.get());

  if(!data->initialized) {
    data->xcs = new float[data->trials];
    data->ycs = new float[data->trials];
    data->rs = new float[data->trials];
    data->initialized = true;
  }

  std::default_random_engine generator;
  // Most Monte Carlo sampling or integration techniques assume a “random number
  // generator,” which generates uniform statistically independent values
  // - MONTE CARLO TECHNIQUES, 2009 by G. Cowan
  std::uniform_real_distribution<float> xcs_dist(XCMIN, XCMAX);
  std::uniform_real_distribution<float> ycs_dist(YCMIN, YCMAX);
  std::uniform_real_distribution<float> rs_dist(RMIN, RMAX);

  for (int n = 0; n < data->trials; n++) {
    data->xcs[n] = xcs_dist(generator);
    data->ycs[n] = ycs_dist(generator);
    data->rs[n] = rs_dist(generator);
  }
}

double CalcSpeedup(double test_two, double test_one) {
  return test_two / test_one;
}

float CalcFp(double speedup, double threads) {
  return (threads / (threads - 1.0)) * (1. - (1. / speedup));
}

int main(int argc, char *argv[]) {
  std::shared_ptr<test_mem> mc_mem = std::make_shared<test_mem>();
  TestRig monte_carlo(mc_mem, mc_run, mc_init);

  // Taking in some command line arguments to control the program.
  float threads = NUMT;
  int trials = NUMTRIALS;

  if(argc >= 2) {
    threads = std::stof(std::string(argv[1]));
  }

  if(argc >= 3) {
    trials = std::stoi(std::string(argv[2]));
  }

  // Some initialization of variables. 
  mc_mem->trials = trials;
  mc_mem->initialized = false;
  double megatrials_one, megatrials_two;
  
  // Setting the precision for output.
  std::cout << std::fixed;
  std::cout << std::setprecision(3);

  if (argc >= 4 && std::string(argv[3]) == "-s") {
    // Test with one thread.
    monte_carlo.Init(1);
    for (int t = 0; t < NUMTRIES; t++) {
      monte_carlo.Run(static_cast<double>(trials));
    }

    megatrials_one = monte_carlo.MaxMegaMults();
  }

  // Test with one thread.
  monte_carlo.Init(threads);
  for (int t = 0; t < NUMTRIES; t++) {
    monte_carlo.Run(static_cast<double>(trials));
  }

  std::cout << "Freeing memory..." << std::endl;
  delete [] mc_mem->xcs;
  delete [] mc_mem->ycs;
  delete [] mc_mem->rs;

  megatrials_two = monte_carlo.MaxMegaMults();

  std::cout << threads << '\t' << trials << '\t' << mc_mem->prob << '\t'
            << monte_carlo.MaxMegaMults() << std::endl;

  if (argc >= 4 && std::string(argv[3]) == "-s") {
    double speedup = CalcSpeedup(megatrials_two, megatrials_one);
    float fp = CalcFp(speedup, threads);

    std::cout << "Speedup = " << speedup << std::endl;
    std::cout << "Parallel fraction = " << fp << std::endl;

    std::cout << "Writing to records.csv..." << std::endl;

    std::ofstream outfile;
    outfile.open("records.csv", std::ios_base::app);

    // Setting the precision for output.
    outfile << std::fixed;
    outfile << std::setprecision(3);
    outfile << threads << '\t' << trials << '\t' << mc_mem->prob << '\t'
              << monte_carlo.MaxMegaMults() << '\t' << fp << std::endl;
    outfile.close();
  }

  return 0;
}
