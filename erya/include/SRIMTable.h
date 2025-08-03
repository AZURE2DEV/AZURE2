#ifndef SRIMTABLE_H
#define SRIMTABLE_H

#include <vector>

class SRIMTable
{
  private:
    std::vector<double> Energies;
    std::vector<double> StoppingPower;

  public:
    SRIMTable( int index );
    ~SRIMTable( );

    // Get entire vectors
    int    getTableSize(){ return Energies.size(); }
    double getEnergy(int index){ return Energies[index]; }
    double getStoppingPower(int index){ return StoppingPower[index]; }

};

#endif