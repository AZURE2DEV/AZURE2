#include "numcmc/mcmc.h"

#include <algorithm>
#include <cstdint>
#include <random>

namespace nu {
namespace {

// Maximum proposals drawn for one walker before giving up and rejecting. A
// walker that cannot reach the support of the posterior would otherwise spin
// in the redraw loop forever.
const int kMaxProposalAttempts = 1000;

// rand() is not reentrant. Sharing it across OpenMP threads serialises the
// sampler at best and corrupts its internal state at worst, so give every
// thread its own generator.
std::mt19937& tls_rng() {
    static thread_local std::mt19937 gen(
        std::random_device{}() +
        static_cast<unsigned>(reinterpret_cast<std::uintptr_t>(&gen)));
    return gen;
}

// Same 1/5001 discretisation the serial path uses.
double tls_uniform() { return double(tls_rng()() % 5001) / 5001.; }

int tls_index(int n) { return int(tls_rng()() % static_cast<unsigned>(n)); }

}  // namespace
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
}

std::vector<std::vector<double>> Mcmc::get_chain() {
    std::vector<std::vector<double>> samples_chain;

    for (int k = 1; k <= this->nwalkers; k++) {
        Walker walker;
        walker = this->sample.getWalkerCopy(k);
        int steps = walker.getSteps();

        for (int i = 1; i <= steps; i++) {
            samples_chain.push_back(walker.getHistory()[i - 1]);
        };
    };

    return samples_chain;
}

std::vector<std::vector<double>> Mcmc::get_chain_walkers() {
    std::vector<std::vector<double>> samples_chain;

    for (int k = 1; k <= this->nwalkers; k++) {
        Walker walker;
        walker = this->sample.getWalkerCopy(k);
        int steps = walker.getSteps();

        for (int i = 1; i <= steps; i++) {
            std::vector<double> pos_aux = walker.getHistory()[i - 1];
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
    do {
        std::string line, col;
        std::getline(file, line);

        // Create a stringstream from line
        std::stringstream ss(line);
        std::vector<double> walker_pos;
        // Extract each column name
        while (std::getline(ss, col, ',')) {
            double xx = std::stold(col);
            walker_pos.push_back(xx);
        };
        ensemble_init.push_back(walker_pos);
    } while (!file.eof());
    return ensemble_init;
}

int Mcmc::save_chain_walker(char *file_name, char *header) {
    std::ofstream sampling_file;
    sampling_file.open(file_name);

    std::vector<std::vector<double>> sampling_points;
    sampling_points = this->get_chain_walkers();

    sampling_file << header;

    /*Writing samples*/
    double nsteps = sampling_points.size();
    for (int i = 0; i < nsteps; i++) {
        for (int j = 0; j < ndim + 1; j++) {
            sampling_file << sampling_points[i][j];
            if (j != ndim) sampling_file << ",";
        };
        if (i != nsteps * nwalkers - 1) sampling_file << std::endl;
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
    double nsteps = sampling_points.size();
    for (int i = 0; i < nsteps; i++) {
        for (int j = 0; j < ndim; j++) {
            sampling_file << sampling_points[i][j];
            if (j != ndim - 1) sampling_file << ",";
        };
        if (i != nsteps * nwalkers - 1) sampling_file << std::endl;
    };
    return 1;
}

int Mcmc::run(double (*func)(std::vector<double> &), int steps) {
    for (int i = 1; i <= steps; i++) {
        // 1:
        for (int k = 1; k <= nwalkers; k++) {
            // 0.0001 s
            Walker walker, walker_aux;
            double logPY, logPX, q, z, r;
            std::vector<double> Y;
            std::vector<double> Xj;
            std::vector<double> Xk;

            walker = sample.getWalker(1);
            // 2:
            // The loop constraint some region
            do {
                Xj.clear();
                Xk.clear();
                Y.clear();
                walker_aux = sample.getRandomWalkerCopy();
                // 3:
                double uniform_rand = double((rand() % 5001)) / 5001.;
                z = F(uniform_rand);
                // 4:
                for (int l = 0; l < this->ndim; l++) {
                    Xj.push_back(walker_aux.getPos()[l]);
                    Xk.push_back(walker.getPos()[l]);
                    Y.push_back(Xj[l] + z * (Xk[l] - Xj[l]));
                };
                // 5:
                logPY = func(Y);
            } while (logPY == -std::numeric_limits<double>::infinity());

            logPX = func(Xk);
            q = pow(z, ndim - 1) * exp(logPY - logPX);
            // 6:
            r = double((rand() % 5001)) / 5001.;
            // 7:"
            if (r <= q) {  // 8:
                walker.setPos(Y);
            }  // 9:
            else {
                // 10:
                walker.setPos(Xk);
            }
            sample.includeWalker(walker);
        };
    };

    return 0;
}

int Mcmc::run_parallel(double (*func)(std::vector<double>&), int steps,
                       int threads) {
    return run_parallel_with_callback(func, steps, threads, nullptr, nullptr,
                                      nullptr);
}

int Mcmc::run_notlog(double (*func)(std::vector<double>&), int steps) {
    for (int i = 1; i <= steps; i++) {
        // 1:
        for (int k = 1; k <= nwalkers; k++) {
            // 0.0001 s
            Walker walker, walker_aux;
            double logPY, logPX, q, z, r;
            std::vector<double> Y;
            std::vector<double> Xj;
            std::vector<double> Xk;

            walker = sample.getWalker(1);
            // 2:
            // The loop constraint some region
            do {
                Xj.clear();
                Xk.clear();
                Y.clear();
                walker_aux = sample.getRandomWalkerCopy();
                // 3:
                double uniform_rand = double((rand() % 5001)) / 5001.;
                z = F(uniform_rand);
                // 4:
                for (int l = 0; l < this->ndim; l++) {
                    Xj.push_back(walker_aux.getPos()[l]);
                    Xk.push_back(walker.getPos()[l]);
                    Y.push_back(Xj[l] + z * (Xk[l] - Xj[l]));
                };
                // 5:
                logPY = func(Y);
            } while (logPY == -std::numeric_limits<double>::infinity());

            logPX = func(Xk);
            q = pow(z, ndim - 1) * logPY / logPX;
            // 6:
            r = double((rand() % 5001)) / 5001.;
            // 7:"
            if (r <= q) {  // 8:
                walker.setPos(Y);
            }  // 9:
            else {
                // 10:
                walker.setPos(Xk);
            }
            sample.includeWalker(walker);
        };
    };

    return 0;
}

int Mcmc::run_with_callback(double (*func)(std::vector<double> &), int steps,
                           void (*progress_callback)(int current_step, int total_steps, double log_prob),
                           bool (*should_stop_callback)(),
                           void (*iteration_callback)(int current_step, int total_steps)) {
    this->should_stop = false;
    
    for (int i = 1; i <= steps; i++) {
        double current_log_prob = -std::numeric_limits<double>::infinity();
        
        // Run one step for all walkers
        for (int k = 1; k <= nwalkers; k++) {
            // Check stop condition before each walker
            if (this->should_stop || (should_stop_callback && should_stop_callback())) {
                return 1; // Stopped early
            }
            
            Walker walker, walker_aux;
            double logPY, logPX, q, z, r;
            std::vector<double> Y;
            std::vector<double> Xj;
            std::vector<double> Xk;

            walker = sample.getWalker(1);
            
            do {
                Xj.clear();
                Xk.clear();
                Y.clear();
                walker_aux = sample.getRandomWalkerCopy();
                double uniform_rand = double((rand() % 5001)) / 5001.;
                z = F(uniform_rand);
                
                for (int l = 0; l < this->ndim; l++) {
                    Xj.push_back(walker_aux.getPos()[l]);
                    Xk.push_back(walker.getPos()[l]);
                    Y.push_back(Xj[l] + z * (Xk[l] - Xj[l]));
                };
                
                logPY = func(Y);
            } while (logPY == -std::numeric_limits<double>::infinity());

            logPX = func(Xk);
            q = pow(z, ndim - 1) * exp(logPY - logPX);
            r = double((rand() % 5001)) / 5001.;
            
            if (r <= q) {
                walker.setPos(Y);
                current_log_prob = std::max(current_log_prob, logPY);
            } else {
                walker.setPos(Xk);
                current_log_prob = std::max(current_log_prob, logPX);
            }
            sample.includeWalker(walker);
        }
        
        // Call iteration callback every step
        if (iteration_callback) {
            iteration_callback(i, steps);
        }
        
        // Call progress callback every step or on last step
        if (progress_callback && (i % 1 == 0 || i == steps)) {
            progress_callback(i, steps, current_log_prob);
        }
    }
    
    return 0;
}

int Mcmc::run_parallel_with_callback(double (*func)(std::vector<double> &), int steps, 
                                     int threads,
                                     void (*progress_callback)(int current_step, int total_steps, double log_prob),
                                     bool (*should_stop_callback)(),
                                     void (*iteration_callback)(int current_step, int total_steps)) {
    this->should_stop = false;

    // Address the walkers by index in one flat vector. The previous version
    // pushed and erased walkers on two shared Ensembles while several threads
    // were reading them, which reallocated the underlying storage under other
    // threads' feet and corrupted the heap.
    std::vector<Walker> walkers = this->sample.getWalkers();
    const int half = nwalkers / 2;

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

#pragma omp parallel for num_threads(threads) schedule(dynamic)
            for (int k = lo; k < hi; k++) {
                std::vector<double> Xk = walkers[k].getPos();
                std::vector<double> Y(this->ndim);
                double logPY = -std::numeric_limits<double>::infinity();
                double z = 1.0;

                // Redraw while the proposal lands outside the support, but
                // give up eventually rather than looping forever.
                for (int attempt = 0; attempt < kMaxProposalAttempts; attempt++) {
                    std::vector<double> Xj = walkers[clo + tls_index(csize)].getPos();
                    z = F(tls_uniform());

                    for (int l = 0; l < this->ndim; l++) {
                        Y[l] = Xj[l] + z * (Xk[l] - Xj[l]);
                    };

                    logPY = func(Y);
                    if (logPY != -std::numeric_limits<double>::infinity()) break;
                }

                double logPX = func(Xk);
                double q = pow(z, ndim - 1) * exp(logPY - logPX);
                double r = tls_uniform();

                // The logPY guard matters when every proposal was rejected:
                // q is then 0, and r can legitimately draw 0 as well.
                double accepted_log_prob;
                if (r <= q && logPY != -std::numeric_limits<double>::infinity()) {
                    walkers[k].setPos(Y);
                    accepted_log_prob = logPY;
                } else {
                    walkers[k].setPos(Xk);
                    accepted_log_prob = logPX;
                }

                // Thread-safe update of current_log_prob
#pragma omp critical
                {
                    current_log_prob = std::max(current_log_prob, accepted_log_prob);
                }
            };
        };

        // Call iteration callback every step
        if (iteration_callback) {
            iteration_callback(i, steps);
        }

        // Call progress callback every step or on last step
        if (progress_callback && (i % 1 == 0 || i == steps)) {
            progress_callback(i, steps, current_log_prob);
        }
    };

    this->sample.setWalkers(walkers);

    return 0;
}

int Mcmc::reset() {
    this->sample.cleanEnsembleHistory();
    this->should_stop = false;
    return 1;
}
}  // namespace nu
