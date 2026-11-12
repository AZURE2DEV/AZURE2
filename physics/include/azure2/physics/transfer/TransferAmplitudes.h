#ifndef AZURE2_PHYSICS_TRANSFER_TRANSFERAMPLITUDES_H
#define AZURE2_PHYSICS_TRANSFER_TRANSFERAMPLITUDES_H

// The angular and spin algebra of the transfer amplitude (docs/upgrade/transfer.md section 19), for Transfer.cpp: an
// internal header of physics/transfer (not an interface for callers; use Transfer.h).

#include <complex>
#include <vector>

#include "azure2/physics/transfer/Transfer.h"
#include "azure2/physics/transfer/TransferKernel.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * T(mu, M_F; m, M_A) of Eq. (19.6) at one angle: m the projection of the projectile a (spin J_a) along k_i, M_A of
 * the target, mu of the ejectile s (spin J_s) and M_F of the residual along z, k_f at (theta, 0) (MeV fm^3, the
 * spectroscopic amplitudes included).  With waves that have no spin-orbit term the sums over j of the distorted
 * waves are identities and T is that of Eqs. (6.3)-(6.5).
 */
struct AmplitudeTable {
  int two_Ja = 0, two_JA = 0, two_Js = 0, two_JF = 0;
  std::vector<std::complex<double>> T;
  std::size_t Index(int two_mu, int two_MF, int two_m, int two_MA) const {
    return ((static_cast<std::size_t>((two_mu + two_Js) / 2) * (two_JF + 1) + (two_MF + two_JF) / 2) * (two_Ja + 1) +
            (two_m + two_Ja) / 2) * (two_JA + 1) + (two_MA + two_JA) / 2;
  }
};

/*!
 * The theta-independent part of Eq. (19.5): for each projection set (m_F, m_a, m', m, mu', mu) and each (l_f, n) the
 * coefficient C of Y_{l_f n}(theta, 0), from the radial integrals K (Eq. (6.4), per pair of waves) and the
 * prefactor (4 pi)^2 / (k_i k_f (1 - alpha beta)^3).  Computed once per reaction; AmplitudesAt evaluates it.
 */
struct AngularCoefficients {
  int la = 0, lF = 0, two_Ja = 0, two_Js = 0, lf_max = 0, n_max = 0;
  std::vector<std::complex<double>> C;  ///< by projection set, l_f, n + n_max
};
AngularCoefficients AngularCoefficientsOf(const TransferReaction &reaction, const DistortedWaves &entrance,
                                          const DistortedWaves &exit, const KernelIntegrals &K, double prefactor);

/// T of Eq. (19.6) at theta (deg).
AmplitudeTable AmplitudesAt(const TransferReaction &reaction, const AngularCoefficients &c, double theta_deg);

/// dsigma/dOmega of Eq. (6.6) (mb/sr in the unit of the residual overlap) from the amplitudes.
double CrossSectionOf(const TransferReaction &reaction, const DistortedWaves &entrance, const DistortedWaves &exit,
                      const AmplitudeTable &t);

/// The vector analysing power iT11 of the projectile (Madison frame, Eq. (19.8)); 0 for a spinless projectile.
double VectorAnalysingPowerOf(const AmplitudeTable &t);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
