#ifndef THMMATRIXFUNC_H
#define THMMATRIXFUNC_H

#include "AMatrixFunc.h"

/// Modified R-matrix (Trojan Horse Method) cross-section calculator.

/*!
 * THMMatrixFunc reuses the shared R-matrix interior of AMatrixFunc
 * (ClearMatrices / FillMatrices / InvertMatrices give the level matrix A,
 * unchanged) and replaces only the channel-surface assembly with the
 * half-off-energy-shell (HOES) cross section of the modified R-matrix
 * formalism used for the Trojan Horse Method:
 *
 *   sigma(E) = sum_J (2J+1) (k_f/mu_f) sum_{c_exit} 2 P_{c_exit}
 *                sum_{s_in} | sum_{la,lap} gamma_{la,c_exit} A_{la,lap}
 *                                sum_{c_in in s_in} gamma_{lap,c_in} M_l |^2
 *
 * coherent over levels (through A) and over entrance partial waves sharing a
 * channel spin (each carrying its own transfer form factor M_l, stored on the
 * point), incoherent over distinct entrance channel spins, exit channels and
 * J^pi groups. The entrance penetrability of conventional R-matrix is replaced
 * by M_l; the exit vertex keeps sqrt(2 P). The absolute scale is arbitrary
 * (folded into the segment normalization), as in mrmpy. See
 * docs/THM_IMPLEMENTATION.md and mrmpy MRMModel.unsmeared.
 */

class THMMatrixFunc : public AMatrixFunc {
 public:
  THMMatrixFunc(CNuc *compound, const Config &configure);
  /// Computes the HOES cross section for the point and stores it as the fit XS.
  void CalculateTHMCrossSection(EPoint *point);
};

#endif
