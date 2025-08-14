#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <stdlib.h>

#include <vector>

#include "numcmc/walker.h"

/************************************************
 *  The Ensemble was tested and properly working.

 *  The order of the walkers does not matter,
    they don't need to have the same history size.
 *  Ensemble constructor K is the number of walkers
    N is the number o dimensions and you should.
************************************************/

namespace nu {

class Ensemble {
   private:
    std::vector<Walker> walkers;
    int nwalkers = 0;
    int ndim;

   public:
    Walker getRandomWalker();
    Walker getRandomWalkerCopy();
    Walker getWalker(int k);
    Walker getWalkerCopy(int k);
    std::vector<Walker> getWalkers();
    void setWalkers(std::vector<Walker> walkers) { this->walkers = walkers; };

    int includeWalker(Walker walker);
    int substituteWalker(int k, Walker walker);
    int cleanEnsembleHistory();
    std::vector<Ensemble> divideEnsemble();

    Ensemble(int K, int N, std::vector<Walker> intial_walkers);
    Ensemble(){};
    ~Ensemble(){};
};
}  // namespace nu

#endif