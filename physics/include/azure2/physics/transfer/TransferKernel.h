#ifndef AZURE2_PHYSICS_TRANSFER_TRANSFERKERNEL_H
#define AZURE2_PHYSICS_TRANSFER_TRANSFERKERNEL_H

// The finite-range radial integrals of Eq. (6.4), for Transfer.cpp: an internal header of physics/transfer (not an
// interface for callers; use Transfer.h).

#include <complex>
#include <cstddef>
#include <vector>

#include "azure2/physics/transfer/Transfer.h"

namespace azure2 {
namespace physics {
namespace transfer {

/// I^L_{l_f l_i} of Eq. (6.4) for L = L_lo ... L_lo + nL - 1, l_f <= lf_max, l_i <= li_max (fm^3 MeV in the unit of
/// the overlaps' u).
struct KernelIntegrals {
  int L_lo = 0, nL = 0, lf_max = 0, li_max = 0;
  std::vector<std::complex<double>> I;
  std::size_t Index(int iL, int lf, int li) const {
    return (static_cast<std::size_t>(iL) * (lf_max + 1) + lf) * (li_max + 1) + li;
  }
};

/// The kernel of Eq. (6.4) with the interaction of reaction.interaction (Eq. (17.1)) on the grid.
KernelIntegrals FiniteRangeIntegrals(const TransferReaction &reaction, const DistortedWaves &entrance,
                                     const DistortedWaves &exit, const TransferGrid &grid);

/// i^n, the highest l of a set of waves, u_l of a spinless set, Y_lm(theta, 0) from a table of stride m_max + 1.
std::complex<double> IPower(int n);
int MaxL(const DistortedWaves &w);
const std::vector<std::complex<double>> &Wave(const DistortedWaves &w, int l);
double Ylm(const std::vector<double> &table, int stride, int l, int m);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
