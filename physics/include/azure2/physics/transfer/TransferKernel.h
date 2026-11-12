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

/// I^L of Eq. (6.4) for L = L_lo ... L_lo + nL - 1 and every pair of partial waves (exit wave wf, entrance wave wi:
/// the positions in DistortedWaves::waves, so (l, j) with a spin-orbit term and l alone without), in fm^3 MeV in the
/// unit of the overlaps' u.  The zero range fills the same table (Transfer.cpp).
struct KernelIntegrals {
  int L_lo = 0, nL = 0, n_exit = 0, n_entrance = 0;
  std::vector<std::complex<double>> I;
  std::size_t Index(int iL, int wf, int wi) const {
    return (static_cast<std::size_t>(iL) * n_exit + wf) * n_entrance + wi;
  }
};

/// The kernel of Eq. (6.4) with the interaction of reaction.interaction (Eq. (17.1)) on the grid.
KernelIntegrals FiniteRangeIntegrals(const TransferReaction &reaction, const DistortedWaves &entrance,
                                     const DistortedWaves &exit, const TransferGrid &grid);

/// i^n, the highest l of a set of waves, u of the wave at position k of a set (its u kept), Y_lm(theta, 0) from a table
/// of stride m_max + 1.
std::complex<double> IPower(int n);
int MaxL(const DistortedWaves &w);
const std::vector<std::complex<double>> &Wave(const DistortedWaves &w, std::size_t k);
double Ylm(const std::vector<double> &table, int stride, int l, int m);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
