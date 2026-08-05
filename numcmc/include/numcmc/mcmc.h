#ifndef MCMC_H
#define MCMC_H

#include <math.h>
#include <time.h>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "numcmc/ensemble.h"
#include "numcmc/walker.h"

namespace nu {

/*!
 * Affine-invariant ensemble sampler (Goodman & Weare 2010, stretch move),
 * following the parallel two-half scheme of Foreman-Mackey et al. 2013 --
 * i.e. the same algorithm as emcee's ``StretchMove``.
 *
 * The log-probability of a walker's current position is cached on the walker,
 * so each step costs one log-probability evaluation per walker rather than
 * two.  A proposal that lands outside the support (non-finite log-probability)
 * is simply rejected: redrawing until the proposal is valid would condition
 * the proposal on the support without the compensating ratio and so would bias
 * the stationary distribution.
 */
class Mcmc {
   private:
    Ensemble sample;
    int nwalkers;
    int ndim;
    bool should_stop;

    // Stretch-move scale parameter `a` (emcee's default is 2).
    double stretch_a = 2.0;

    // One generator per walker.  A walker is only ever advanced by one thread
    // at a time, so per-walker state needs no locking and -- unlike thread-
    // local state -- makes a run reproducible for a given seed regardless of
    // the thread count.
    std::vector<std::mt19937_64> rngs;
    std::uint64_t seed = 12345u;

    // Acceptance bookkeeping, over the lifetime of the sampler.
    long naccepted = 0;
    long nproposed = 0;

    // Walkers whose initial position had a non-finite log-probability. Such a
    // walker can never move, so the caller should warn about it.
    int ninvalid_initial = 0;

    // Inverse CDF of g(z) ~ 1/sqrt(z) on [1/a, a]; identical to emcee's
    // ``((a - 1) * u + 1)**2 / a``.
    double F(double u) const {
        double t = (stretch_a - 1.0) * u + 1.0;
        return t * t / stretch_a;
    };

    void seed_rngs();

    std::vector<std::vector<double>> load_state(char *file_name);

   public:
    /*! Called once per completed step with every walker's current position and
     *  cached log-probability.  This -- not the log-probability function -- is
     *  where a caller should record the chain: it fires exactly once per step
     *  per walker and reports accepted states, including the repeats that a
     *  rejection must contribute. */
    typedef void (*SampleCallback)(
        int current_step, const std::vector<std::vector<double>> &positions,
        const std::vector<double> &logps);

    int run(double (*func)(std::vector<double> &), int steps);
    int run_parallel(double (*func)(std::vector<double> &), int steps,
                     int threads);

    // Enhanced run method with progress callback and stop functionality
    int run_with_callback(double (*func)(std::vector<double> &), int steps,
                         void (*progress_callback)(int current_step, int total_steps, double log_prob) = nullptr,
                         bool (*should_stop_callback)() = nullptr,
                         void (*iteration_callback)(int current_step, int total_steps) = nullptr,
                         SampleCallback sample_callback = nullptr);

    // Enhanced parallel run method with progress callback and stop functionality
    int run_parallel_with_callback(double (*func)(std::vector<double> &), int steps,
                                   int threads,
                                   void (*progress_callback)(int current_step, int total_steps, double log_prob) = nullptr,
                                   bool (*should_stop_callback)() = nullptr,
                                   void (*iteration_callback)(int current_step, int total_steps) = nullptr,
                                   SampleCallback sample_callback = nullptr);

    int reset();
    int save_state(char *file_name);
    void request_stop() { should_stop = true; }
    bool is_stop_requested() const { return should_stop; }

    /*! Seed for the per-walker generators.  Runs are reproducible for a given
     *  (seed, nwalkers) pair whatever the thread count. */
    void set_seed(std::uint64_t new_seed) {
        this->seed = new_seed;
        seed_rngs();
    };
    /*! Stretch-move scale; must be > 1.  Larger means bolder proposals and a
     *  lower acceptance fraction. */
    void set_stretch_scale(double a) {
        if (a > 1.0) this->stretch_a = a;
    };
    /*! Fraction of proposals accepted so far.  Healthy ensemble sampling sits
     *  roughly in 0.2-0.5; near 0 or near 1 means the chain is not exploring. */
    double acceptance_fraction() const {
        return (nproposed > 0) ? double(naccepted) / double(nproposed) : 0.0;
    };
    int num_invalid_initial() const { return ninvalid_initial; };

    /*! The walkers in their current state, for persisting and later resuming
     *  the ensemble. */
    std::vector<Walker> get_walkers() { return sample.getWalkers(); };

    std::vector<std::vector<double>> get_chain();
    std::vector<std::vector<double>> get_chain_walkers();
    int save_chain(char *file_name, char *header);
    int save_chain_walker(char *file_name, char *header);

    Mcmc(int K, int N, std::vector<std::vector<double>> init_sample);
    Mcmc(int K, int N, char *file_name);
    ~Mcmc(){};

   private:
    int run_impl(double (*func)(std::vector<double> &), int steps, int threads,
                 void (*progress_callback)(int, int, double),
                 bool (*should_stop_callback)(),
                 void (*iteration_callback)(int, int),
                 SampleCallback sample_callback);
};

}  // namespace nu

#endif
