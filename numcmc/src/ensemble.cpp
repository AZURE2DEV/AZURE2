#include "numcmc/ensemble.h"

namespace nu {

Ensemble::Ensemble(int K, int N, std::vector<Walker> initial_walkers) {
    this->walkers = initial_walkers;
    this->nwalkers = K;
    this->ndim = N;
}

int Ensemble::includeWalker(Walker walker) {
    this->walkers.push_back(walker);
    this->nwalkers++;
    return 1;
};

int Ensemble::substituteWalker(int k, Walker walker) {
    this->walkers[k - 1] = walker;
    return 1;
};

Walker Ensemble::getWalker(int k) {
    Walker walker = walkers[k - 1];
    this->nwalkers--;
    walkers.erase(walkers.begin() + k - 1);
    return walker;
};

Walker Ensemble::getWalkerCopy(int k) { return this->walkers[k - 1]; };

Walker Ensemble::getRandomWalker() {
    int k = rand() % nwalkers + 1;
    Walker walker = walkers[k - 1];
    walkers.erase(walkers.begin() + k - 1);
    return walker;
};

Walker Ensemble::getRandomWalkerCopy() {
    int k = rand() % nwalkers + 1;
    return this->walkers[k - 1];
};

int Ensemble::cleanEnsembleHistory() {
    for (int k = 0; k < nwalkers; k++) {
        this->walkers[k].clearHistory();
    };
    return 1;
};

std::vector<Ensemble> Ensemble::divideEnsemble() {
    std::size_t const half_size = this->walkers.size() / 2;
    std::vector<Walker> split_lo(this->walkers.begin(),
                                 this->walkers.begin() + half_size);
    std::vector<Walker> split_hi(this->walkers.begin() + half_size,
                                 this->walkers.end());

    Ensemble sample1((int)this->walkers.size() / 2, ndim, split_lo);
    Ensemble sample2((int)this->walkers.size() / 2, ndim, split_hi);

    std::vector<Ensemble> samples;
    samples.push_back(sample1);
    samples.push_back(sample2);

    return samples;
}

std::vector<Walker> Ensemble::getWalkers() { return this->walkers; }
}  // namespace nu