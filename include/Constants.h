#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
#include <vector>

// Physical constants: CODATA 2018 (since 2026-11-05; before, pi was rounded to 3.141592650 -- 1.1e-9 low -- and
// hbarc, uconv, fstruc, the Avogadro number and the nuclear magneton were CODATA 2006 or rounded).
const double pi = 3.14159265358979323846;
const double hbarc = 197.3269804;           // MeV fm (exact in the 2019 SI)
const double uconv = 931.49410242;          // MeV per u
const double fstruc = 1.00 / 137.035999084;
const double boltzConst = 8.617333262e-2;  // MeV/GK (CODATA 2018); was 8.6171e-2, 2.7e-5 low
const double lightSpeedInCmPerS = 29979245800.;
const double avagadroNum = 6.02214076e23;  // exact in the 2019 SI
const double nuclearMagneton = 0.1051544551;  // e fm: hbarc / (2 m_p c^2), m_p c^2 = 938.27208816 MeV
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
