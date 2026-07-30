#ifndef WALKER_H
#define WALKER_H

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>


/************************************************

 *  The Walker was tested and properly working

 *  We should point that the first hystory position
    is refered with the number 1 and the last
    position is refered as steps

 *  Each recorded position carries the log-probability
    that was evaluated for it.  Keeping the two together
    is what lets the sampler reuse the current position's
    log-probability instead of recomputing it every step.

************************************************/

namespace nu {
// Walkers, the history begins with 1 and ends with steps
class Walker {
   private:
    std::vector<std::vector<double> > history;
    std::vector<double> logpHistory;
    std::vector<double> pos;

    // Log-probability of `pos`.  NaN means "not evaluated yet"; the sampler
    // tests hasLogP() before trusting it.
    double logp = std::numeric_limits<double>::quiet_NaN();

    int steps = 0;
    int ndim;

   public:
    std::vector<double> getPos();
    std::vector<double> getPos(int n);
    int setPos(std::vector<double> new_pos);
    int setPos(std::vector<double> new_pos, double new_logp);
    int getSteps();
    std::vector<std::vector<double> > getHistory() { return this->history; };
    // Reference access, for callers that only read.  getHistory() copies the
    // whole history, which turns a sequential scan into O(steps^2).
    const std::vector<std::vector<double> > &historyRef() const {
        return this->history;
    };
    const std::vector<double> &logpHistoryRef() const {
        return this->logpHistory;
    };
    int clearHistory();

    double getLogP() const { return this->logp; };
    void setLogP(double new_logp) { this->logp = new_logp; };
    bool hasLogP() const { return !std::isnan(this->logp); };

    Walker(){};
    Walker(int ndim);
    Walker(int ndim, std::vector<double> initial_position);
    ~Walker();
};

}  // namespace nu

#endif
