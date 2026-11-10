#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
#include <vector>

#include "azure2/physics/Constants.h"

// Physical constants: CODATA 2018 (since 2026-11-05; before, pi was rounded to 3.141592650 -- 1.1e-9 low -- and
// hbarc, uconv, fstruc, the Avogadro number and the nuclear magneton were CODATA 2006 or rounded).
const double pi = 3.14159265358979323846;
// hbar c, the atomic mass unit and alpha are physics/'s (azure2/physics/Constants.h), so that the channel kinematics
// of physics/ and of the engine use one set of values.
const double hbarc = azure2::physics::kHbarC_MeV_fm;  // MeV fm (exact in the 2019 SI)
const double uconv = azure2::physics::kAtomicMassUnit_MeV;  // MeV per u
const double fstruc = azure2::physics::kFineStructure;
const double boltzConst = 8.617333262e-2;  // MeV/GK (CODATA 2018); was 8.6171e-2, 2.7e-5 low
const double lightSpeedInCmPerS = 29979245800.;
const double avagadroNum = 6.02214076e23;  // exact in the 2019 SI
const double nuclearMagneton = 0.1051544551;  // e fm: hbarc / (2 m_p c^2), m_p c^2 = 938.27208816 MeV
// Two masses kept at the values the code has always used, so that no result moves (maintainers' decision
// 2026-11-10).  Each is used for one purpose only; to adopt the current value, change it here and nowhere else.
// m_e c^2 (MeV) in the argument of the beta-decay Fermi function: CODATA 2002.  With CODATA 2018
// (azure2::physics::kElectronMass_MeV = 0.51099895) the beta-delayed outputs move by about 1e-7 relative.
const double electronMassOfFermiFunction = 0.510998903;
// Proton mass (u) in the M1 external-capture amplitude (effective charge and magnetic-moment term).  CODATA 2018:
// 1.007276466621 (8.6e-8 relative above this value).
const double protonMassOfM1Capture = 1.00727638;
const unsigned char isE1 = 1 << 0;
const unsigned char isM1 = 1 << 1;
const unsigned char isE2 = 1 << 2;
const int maxECMult = 2;

typedef std::complex<double> complex;
typedef std::vector<double> vector_r;
typedef std::vector<std::complex<double>> vector_c;
typedef std::vector<std::vector<double>> matrix_r;
typedef std::vector<std::vector<std::complex<double>>> matrix_c;
typedef std::vector<std::vector<std::vector<double>>> vector_matrix_r;
typedef std::vector<std::vector<std::vector<std::complex<double>>>> vector_matrix_c;

#endif
