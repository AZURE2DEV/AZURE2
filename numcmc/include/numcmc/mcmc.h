#ifndef MCMC_H
#define MCMC_H

#include <math.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "numcmc/ensemble.h"
#include "numcmc/walker.h"

namespace nu {
class Mcmc {
   private:
    Ensemble sample;
    int nwalkers;
    int ndim;
    bool should_stop;

    // Cumulative Inverse of g(z)
    double F(double z) { return 0.5 * (z * z + 2 * z + 1); };
    [[maybe_unused]] bool constrained = true;

    std::vector<std::vector<double>> load_state(char *file_name);

   public:
    int run(double (*func)(std::vector<double> &), int steps);
    int run_parallel(double (*func)(std::vector<double> &), int steps,
                     int threads);
    int run_notlog(double (*func)(std::vector<double> &), int steps);
    
    // Enhanced run method with progress callback and stop functionality
    int run_with_callback(double (*func)(std::vector<double> &), int steps,
                         void (*progress_callback)(int current_step, int total_steps, double log_prob) = nullptr,
                         bool (*should_stop_callback)() = nullptr,
                         void (*iteration_callback)(int current_step, int total_steps) = nullptr);
    
    // Enhanced parallel run method with progress callback and stop functionality
    int run_parallel_with_callback(double (*func)(std::vector<double> &), int steps,
                                   int threads,
                                   void (*progress_callback)(int current_step, int total_steps, double log_prob) = nullptr,
                                   bool (*should_stop_callback)() = nullptr,
                                   void (*iteration_callback)(int current_step, int total_steps) = nullptr);
    
    int reset();
    int save_state(char *file_name);
    void request_stop() { should_stop = true; }
    bool is_stop_requested() const { return should_stop; }

    std::vector<std::vector<double>> get_chain();
    std::vector<std::vector<double>> get_chain_walkers();
    int save_chain(char *file_name, char *header);
    int save_chain_walker(char *file_name, char *header);

    Mcmc(int K, int N, std::vector<std::vector<double>> init_sample);
    Mcmc(int K, int N, char *file_name);
    ~Mcmc(){};
};

}  // namespace nu

#endif