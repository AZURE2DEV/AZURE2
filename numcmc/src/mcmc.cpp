#include "numcmc/mcmc.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>

namespace nu {
namespace {

// SplitMix64, used only to turn one user seed into well-separated per-walker
// seeds.  Seeding N generators with seed, seed+1, ... would leave their output
// streams correlated.
std::uint64_t splitmix64(std::uint64_t &state) {
    std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

}  // namespace

void Mcmc::seed_rngs() {
    std::uint64_t state = this->seed;
    this->rngs.clear();
    this->rngs.reserve(this->nwalkers);
    for (int k = 0; k < this->nwalkers; k++) {
        this->rngs.push_back(std::mt19937_64(splitmix64(state)));
    };
}

Mcmc::Mcmc(int K, int N, std::vector<std::vector<double>> init_sample) {
    std::vector<Walker> walkers;
    for (int k_walk = 0; k_walk < K; k_walk++) {
        Walker single_walker(N, init_sample[k_walk]);
        walkers.push_back(single_walker);
    };
    Ensemble create_sample(K, N, walkers);
    this->sample = create_sample;
    this->nwalkers = K;
    this->ndim = N;
    this->should_stop = false;
    seed_rngs();
}

Mcmc::Mcmc(int K, int N, char *file_name) {
    std::vector<std::vector<double>> init_sample = this->load_state(file_name);
    std::vector<Walker> walkers;
    for (int k_walk = 0; k_walk < K; k_walk++) {
        Walker single_walker(N, init_sample[k_walk]);
        walkers.push_back(single_walker);
    };
    Ensemble create_sample(K, N, walkers);
    this->sample = create_sample;
    this->nwalkers = K;
    this->ndim = N;
    this->should_stop = false;
    seed_rngs();
}

std::vector<std::vector<double>> Mcmc::get_chain() {
    std::vector<std::vector<double>> samples_chain;

    for (int k = 1; k <= this->nwalkers; k++) {
        // By reference: getHistory() copies the whole history, which would make
        // this scan quadratic in the number of steps.
        const std::vector<std::vector<double>> &history =
            this->sample.getWalkerRef(k).historyRef();

        for (std::size_t i = 0; i < history.size(); i++) {
            samples_chain.push_back(history[i]);
        };
    };

    return samples_chain;
}

std::vector<std::vector<double>> Mcmc::get_chain_walkers() {
    std::vector<std::vector<double>> samples_chain;

    for (int k = 1; k <= this->nwalkers; k++) {
        const std::vector<std::vector<double>> &history =
            this->sample.getWalkerRef(k).historyRef();

        for (std::size_t i = 0; i < history.size(); i++) {
            std::vector<double> pos_aux = history[i];
            pos_aux.push_back(k);
            samples_chain.push_back(pos_aux);
        };
    };

    return samples_chain;
}

int Mcmc::save_state(char *file_name) {
    std::ofstream file;
    file.open(file_name);

    if (!(file.is_open())) {
        std::cout << "\n\n ****You should pass an opened file!****\n\n";
        return 0;
    };

    /*Writing samples*/

    for (int i = 1; i <= this->nwalkers; i++) {
        Walker aux_walker;
        aux_walker = this->sample.getWalkerCopy(i);
        for (int j = 0; j < this->ndim; j++) {
            file << aux_walker.getPos()[j];
            if (j != ndim - 1) file << ",";
        };
        if (i != nwalkers) file << std::endl;
    };
    file.close();
    return 1;
}

std::vector<std::vector<double>> Mcmc::load_state(char *file_name) {
    std::ifstream file;
    std::vector<std::vector<double>> ensemble_init;
    file.open(file_name);

    if (!(file.is_open())) {
        std::cout << "\n\n ****You should pass an opened file!****\n\n";
        return ensemble_init;
    };

    this->reset();

    /*Load samples*/
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // Create a stringstream from line
        std::stringstream ss(line);
        std::string col;
        std::vector<double> walker_pos;
        // Extract each column name
        while (std::getline(ss, col, ',')) {
            double xx = std::stod(col);
            walker_pos.push_back(xx);
        };
        if (!walker_pos.empty()) ensemble_init.push_back(walker_pos);
    };
    return ensemble_init;
}

int Mcmc::save_chain_walker(char *file_name, char *header) {
    std::ofstream sampling_file;
    sampling_file.open(file_name);

    std::vector<std::vector<double>> sampling_points;
    sampling_points = this->get_chain_walkers();

    sampling_file << header;

    /*Writing samples*/
    std::size_t nrows = sampling_points.size();
    for (std::size_t i = 0; i < nrows; i++) {
        for (int j = 0; j < ndim + 1; j++) {
            sampling_file << sampling_points[i][j];
            if (j != ndim) sampling_file << ",";
        };
        if (i != nrows - 1) sampling_file << std::endl;
    };
    return 1;
}

int Mcmc::save_chain(char *file_name, char *header) {
    std::ofstream sampling_file;
    sampling_file.open(file_name);

    std::vector<std::vector<double>> sampling_points;
    sampling_points = this->get_chain();

    sampling_file << header;

    /*Writing samples*/
    std::size_t nrows = sampling_points.size();
    for (std::size_t i = 0; i < nrows; i++) {
        for (int j = 0; j < ndim; j++) {
            sampling_file << sampling_points[i][j];
            if (j != ndim - 1) sampling_file << ",";
        };
        if (i != nrows - 1) sampling_file << std::endl;
    };
    return 1;
}

int Mcmc::run(double (*func)(std::vector<double> &), int steps) {
    return run_impl(func, steps, 1, nullptr, nullptr, nullptr, nullptr);
}

int Mcmc::run_parallel(double (*func)(std::vector<double>&), int steps,
                       int threads) {
    return run_impl(func, steps, threads, nullptr, nullptr, nullptr, nullptr);
}

int Mcmc::run_with_callback(double (*func)(std::vector<double> &), int steps,
                           void (*progress_callback)(int current_step, int total_steps, double log_prob),
                           bool (*should_stop_callback)(),
                           void (*iteration_callback)(int current_step, int total_steps),
                           SampleCallback sample_callback) {
    return run_impl(func, steps, 1, progress_callback, should_stop_callback,
                    iteration_callback, sample_callback);
}

int Mcmc::run_parallel_with_callback(double (*func)(std::vector<double> &), int steps,
                                     int threads,
                                     void (*progress_callback)(int current_step, int total_steps, double log_prob),
                                     bool (*should_stop_callback)(),
                                     void (*iteration_callback)(int current_step, int total_steps),
                                     SampleCallback sample_callback) {
    return run_impl(func, steps, threads, progress_callback,
                    should_stop_callback, iteration_callback, sample_callback);
}

int Mcmc::run_impl(double (*func)(std::vector<double> &), int steps,
                   int threads,
                   void (*progress_callback)(int, int, double),
                   bool (*should_stop_callback)(),
                   void (*iteration_callback)(int, int),
                   SampleCallback sample_callback) {
    this->should_stop = false;

    // The stretch move draws the complementary walker from the other half, so
    // each half needs at least two members for the move to be well defined.
    if (this->nwalkers < 4) return 3;

    // Address the walkers by index in one flat vector. The previous version
    // pushed and erased walkers on two shared Ensembles while several threads
    // were reading them, which reallocated the underlying storage under other
    // threads' feet and corrupted the heap.
    std::vector<Walker> walkers = this->sample.getWalkers();
    if ((int)walkers.size() != this->nwalkers) return 2;
    if ((int)this->rngs.size() != this->nwalkers) seed_rngs();

    if (threads < 1) threads = 1;
    const int half = nwalkers / 2;

    // One log-probability evaluation per walker up front; from here on a
    // walker's current log-probability is whatever is cached on it, so a step
    // costs one evaluation per walker instead of two.
#pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (int k = 0; k < nwalkers; k++) {
        if (!walkers[k].hasLogP()) {
            std::vector<double> pos = walkers[k].getPos();
            walkers[k].setLogP(func(pos));
        }
    };

    this->ninvalid_initial = 0;
    for (int k = 0; k < nwalkers; k++) {
        if (!std::isfinite(walkers[k].getLogP())) this->ninvalid_initial++;
    };

    for (int i = 1; i <= steps; i++) {
        double current_log_prob = -std::numeric_limits<double>::infinity();

        // Check stop condition at beginning of each step
        if (this->should_stop || (should_stop_callback && should_stop_callback())) {
            this->sample.setWalkers(walkers);
            return 1; // Stopped early
        }

        // Goodman & Weare stretch move: the two halves are updated one after
        // the other, never simultaneously. While half [lo,hi) moves, the
        // complementary half [clo,chi) is strictly read-only, so every walker
        // in the moving half can be updated concurrently without sharing state.
        for (int ip = 0; ip <= 1; ip++) {
            const int lo = (ip == 0) ? 0 : half;
            const int hi = (ip == 0) ? half : nwalkers;
            const int clo = (ip == 0) ? half : 0;
            const int chi = (ip == 0) ? nwalkers : half;
            const int csize = chi - clo;
            if (csize <= 0 || hi <= lo) continue;

            long accepted_here = 0;

#pragma omp parallel for num_threads(threads) schedule(dynamic) \
    reduction(+ : accepted_here)
            for (int k = lo; k < hi; k++) {
                std::mt19937_64 &rng = this->rngs[k];
                std::uniform_real_distribution<double> unif(0.0, 1.0);

                std::vector<double> Xk = walkers[k].getPos();
                const double logPX = walkers[k].getLogP();

                // Complementary walker, drawn from the other half only.
                std::uniform_int_distribution<int> pick(0, csize - 1);
                const std::vector<double> Xj = walkers[clo + pick(rng)].getPos();

                const double z = F(unif(rng));
                std::vector<double> Y(this->ndim);
                for (int l = 0; l < this->ndim; l++) {
                    Y[l] = Xj[l] + z * (Xk[l] - Xj[l]);
                };

                const double logPY = func(Y);

                // Metropolis-Hastings in log space: exp() of a large positive
                // difference overflows to inf, and a proposal outside the
                // support must simply be rejected (redrawing until it lands
                // inside would bias the stationary distribution).
                bool accept = false;
                if (std::isfinite(logPY)) {
                    const double lndiff =
                        (this->ndim - 1) * std::log(z) + logPY - logPX;
                    accept = (lndiff >= 0.0) || (std::log(unif(rng)) <= lndiff);
                }

                if (accept) {
                    walkers[k].setPos(Y, logPY);
                    accepted_here++;
                } else {
                    // Record the repeat: a rejected step still contributes the
                    // current position to the chain.
                    walkers[k].setPos(Xk, logPX);
                }

                const double kept_log_prob = walkers[k].getLogP();
#pragma omp critical
                {
                    current_log_prob = std::max(current_log_prob, kept_log_prob);
                }
            };

            this->naccepted += accepted_here;
            this->nproposed += (hi - lo);
        };

        // Hand the completed step to the caller: every walker's accepted state,
        // exactly once.
        if (sample_callback) {
            std::vector<std::vector<double>> positions(nwalkers);
            std::vector<double> logps(nwalkers);
            for (int k = 0; k < nwalkers; k++) {
                positions[k] = walkers[k].getPos();
                logps[k] = walkers[k].getLogP();
            };
            sample_callback(i, positions, logps);
        }

        // Call iteration callback every step
        if (iteration_callback) {
            iteration_callback(i, steps);
        }

        // Call progress callback every step or on last step
        if (progress_callback) {
            progress_callback(i, steps, current_log_prob);
        }
    };

    this->sample.setWalkers(walkers);

    return 0;
}

int Mcmc::reset() {
    this->sample.cleanEnsembleHistory();
    this->should_stop = false;
    this->naccepted = 0;
    this->nproposed = 0;
    this->ninvalid_initial = 0;
    return 1;
}
}  // namespace nu
