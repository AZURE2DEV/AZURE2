#ifndef WALKER_H
#define WALKER_H

#include <iostream>
#include <vector>


/************************************************

 *  The Walker was tested and properly working

 *  We should point that the first hystory position
    is refered with the number 1 and the last
    position is refered as steps

************************************************/

namespace nu {
// Walkers, the history begins with 1 and ends with steps
class Walker {
   private:
    std::vector<std::vector<double> > history;
    std::vector<double> pos;

    int steps = 0;
    int ndim;

   public:
    std::vector<double> getPos();
    std::vector<double> getPos(int n);
    int setPos(std::vector<double> new_pos);
    int getSteps();
    std::vector<std::vector<double> > getHistory() { return this->history; };
    int clearHistory();

    Walker(){};
    Walker(int ndim);
    Walker(int ndim, std::vector<double> initial_position);
    ~Walker();
};

}  // namespace nu

#endif