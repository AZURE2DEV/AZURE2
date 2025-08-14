#include "numcmc/walker.h"

namespace nu {
Walker::Walker(int ndim) { this->ndim = ndim; };
Walker::Walker(int ndim, std::vector<double> initial_position) {
    this->ndim = ndim;
    this->setPos(initial_position);
};
Walker::~Walker(){};

std::vector<double> Walker::getPos(int n) { return this->history[n - 1]; }

std::vector<double> Walker::getPos() { return this->pos; }

int Walker::getSteps() { return this->steps; }

int Walker::setPos(std::vector<double> new_pos) {
    this->pos = new_pos;
    this->history.push_back(new_pos);
    this->steps++;

    return 1;
}

int Walker::clearHistory() {
    this->history.clear();
    this->steps = 1;
    this->history.push_back(pos);
    return 1;
}

}  // namespace nu