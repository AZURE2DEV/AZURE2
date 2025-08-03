#ifndef SRIM_H
#define SRIM_H

#include <string>
#include <vector>

#include "include/SRIMTable.h"

class SRIM
{
  private:
    std::vector<std::string> Elements;
    std::vector<double> Values1;
    std::vector<double> Values2;
    std::vector<double> Values3;
    std::vector<double> Values4;
    std::vector<double> Values5;
    std::vector<double> Values6;
    std::vector<double> Values7;
    std::vector<double> Values8;
    std::vector<double> Values9;
    std::vector<double> Values10;
    std::vector<double> Values11;
    std::vector<double> Values12;
    std::vector<double> Mass;
    std::vector<double> Density;
    std::vector<double> Bloch;
    std::vector<SRIMTable> Tables;

  public:
    SRIM( );
    ~SRIM( );

    // Get values for given index
    std::string getElementName(int index){ return Elements[index-1]; }
    double getElementValue1(int index){ return Values1[index-1]; }
    double getElementValue2(int index){ return Values2[index-1]; }
    double getElementValue3(int index){ return Values3[index-1]; }
    double getElementValue4(int index){ return Values4[index-1]; }
    double getElementValue5(int index){ return Values5[index-1]; }
    double getElementValue6(int index){ return Values6[index-1]; }
    double getElementValue7(int index){ return Values7[index-1]; }
    double getElementValue8(int index){ return Values8[index-1]; }
    double getElementValue9(int index){ return Values9[index-1]; }
    double getElementValue10(int index){ return Values10[index-1]; }
    double getElementValue11(int index){ return Values11[index-1]; }
    double getElementValue12(int index){ return Values12[index-1]; }
    double getElementMass(int index){ return Mass[index-1]; }
    double getElementDensity(int index){ return Density[index-1]; }
    double getElementBloch(int index){ return Bloch[index-1]; }
    double getElementIonization(int index){ return Bloch[index-1] * (index-1) * 0.001; }
    SRIMTable getElementTable(int index){ return Tables[index-1]; }
    
};

#endif