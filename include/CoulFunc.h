#ifndef COULFUNC_H
#define COULFUNC_H

#include <cmath>
#include <complex>
#include <memory>
#include <optional>
#include <unordered_map>
#include "Kinematics.h"
#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/HybridPotential.h"

class PPair;
class CoulombTable;

/// The return structure of the CoulFunc function class.

/*!
 * The CoulWaves structure contains both the irregular and regular solutions to the
 * Coulomb equation, as well as their derivatives with respect to \f$ \rho \f$.
 */

struct CoulWaves {
  /// Regular solution
  double F;
  /// Derivative of regular solution with respect to \f$ \rho \f$
  double dF;
  /// Irregular solution
  double G;
  /// Derivative of irregular solution with respect to \f$ \rho \f$
  double dG;
};

/// F, F', G, G' at (l, eta, rho) from physics/coulomb (docs/upgrade/coulomb.md) by the engine's algorithm, coulsm;
/// throws std::runtime_error when the evaluation fails (the library returns NaN with a status, never a number).
CoulWaves CoulombWavesAt(int l, double eta, double rho);

/*!
 * exp(i phi) of the hard-sphere phase, (G - i F)/sqrt(F^2 + G^2).  Just above a
 * charged-particle threshold (eta >~ 110) G^2 overflows and the plain formula
 * returns 0 instead of a unit phase; the scaled form (hypot) is used only then,
 * so every value the plain formula gets right is unchanged to the bit.
 */
/// The engine's Coulomb waves as physics/'s (same numbers; physics/ formulas take them).
inline azure2::physics::coulomb::CoulombWaves PhysicsWaves(const CoulWaves &c) {
  return azure2::physics::coulomb::CoulombWaves{c.F, c.dF, c.G, c.dG, 0.0, azure2::physics::coulomb::Status::ok};
}
inline std::complex<double> CoulombHardSpherePhase(const CoulWaves &c) { return azure2::physics::exp_minus_i_phi(PhysicsWaves(c)); }

/// A function class to calculate Coulomb functions for positive energy channels

/*!
 * The CoulFunc function class calculates the solutions to the Coulomb equation,
 * as well as other useful quantities such as shift functions and their energy
 * derivative and penetrabilities.
 */

/*!
 * The Coulomb waves of each (pair, l, 2J, energy) block at the radii asked for, keyed by the exact radius: what the
 * external-capture integrals of one entrance l and energy share (the adaptive quadrature asks for the same radii for
 * every term).  Everything CoulFunc::operator() depends on besides (l, radius, energy) is in the block's key, so a
 * stored value is the value operator() computes.  Owned by the caller (the run's per-thread scratch, RunScratch.h),
 * handed to a CoulFunc with SetRadialMemo; at most kMaxRadii radii, emptied when full.
 */
struct CoulRadialMemo {
  struct Key {
    int z1, z2, l, two_j;
    double redmass, mtot, energy;
    std::optional<azure2::physics::HybridPotential> potential;
    bool operator==(const Key &o) const {
      return z1 == o.z1 && z2 == o.z2 && l == o.l && two_j == o.two_j && redmass == o.redmass && mtot == o.mtot &&
             energy == o.energy && potential == o.potential;
    }
  };
  struct Hash {
    std::size_t operator()(const Key &k) const {
      std::size_t h = std::hash<double>()(k.energy);
      for (std::size_t v : {std::hash<int>()(k.z1), std::hash<int>()(k.z2), std::hash<int>()(k.l), std::hash<int>()(k.two_j),
                            std::hash<double>()(k.redmass), std::hash<double>()(k.mtot),
                            std::hash<double>()(k.potential ? k.potential->central_depth_MeV : -1.0)})
        h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      return h;
    }
  };
  std::unordered_map<Key, std::unordered_map<double, CoulWaves>, Hash> blocks;
  std::size_t size = 0;
  static constexpr std::size_t kMaxRadii = std::size_t(1) << 20;
};

class CoulFunc {
 public:
  /// Build for a particle pair, with that pair's hybrid nuclear potential if it has one.
  explicit CoulFunc(const PPair *pPair);
  /// Charge number of the first particle.
  int z1() const;
  /// Charge number of the second particle.
  int z2() const;
  /// Reduced mass of the pair, in u.
  double redmass() const;
  /// Orbital angular momentum of the last evaluation.
  int lLast() const;
  /// The 2J key (PPair::SpinOrbitTwoJ) of the last evaluation.
  int twoJLast() const { return twoJLast_; }
  /// Radius of the last evaluation, fm.
  double radiusLast() const;
  /// Centre-of-mass energy of the last evaluation, MeV.
  double energyLast() const;
  /// The last functions computed; a one-deep memo for repeated calls at the same point.
  struct CoulWaves coulLast() const;
  /// Record an evaluation as the memoized one.
  void setLast(int l, int two_j, double radius, double energy, CoulWaves);
  /// Every function below is of a channel (l, two_j): two_j is the key PPair::SpinOrbitTwoJ makes from the channel's
  /// J (kNoSpinOrbit for a pair without a spin-orbit term, whose functions do not depend on J).
  /// Coulomb functions at (l, two_j, radius fm, centre-of-mass energy MeV).
  CoulWaves operator()(int l, int two_j, double radius, double energy);
  /// Penetrability \f$P_l\f$ at (l, two_j, radius, energy).
  double Penetrability(int l, int two_j, double radius, double energy);
  /// P of the waves outside the radius alone: Penetrability without the interior factor of a potential that also acts
  /// inside the radius (equal to it for any other pair).
  double OuterPenetrability(int l, int two_j, double radius, double energy);
  /// Shift function \f$S_l\f$ at (l, two_j, radius, energy), positive energies.
  double PEShift(int l, int two_j, double radius, double energy);
  /// Hard-sphere phase factor exp(-i phi_l) = (G - iF)/|G + iF| at (l, two_j, radius, energy):
  /// from the Coulomb table when one is in use and covers the point, else from the waves.
  std::complex<double> HardSpherePhase(int l, int two_j, double radius, double energy);
  /// Energy derivative \f$dS_l/dE\f$ -- the term that makes the observed-width transformation singular when it grows too large.
  double PEShift_dE(int l, int two_j, double radius, double energy, double step = 1.0e-6);
  /// Second energy derivative of the shift function (numerical).
  double PEShift_d2E(int l, int two_j, double radius, double energy);

  /// Every value of the pair's potentials (PPair::HybridPotentialKey), for the keys of tables of the whole pair.
  std::vector<double> PotentialKeyOfPair() const;
  /// Whether P, S and the phase come from the universal Coulomb table.
  bool UsesTable() const { return table_ != nullptr; }
  /// Does the pair's potential also act inside the channel radius (HybridPotential::interior_boundary)?  Then P, S
  /// and the phase are the primed functions of docs/upgrade/interior_potential.md Eqs. (4), (5).
  bool Interior(int l) const { return Pot(l) && Pot(l)->interior_boundary; }  // (Pot: CoulFunc.cpp)
  /// The potential of the channels of orbital angular momentum l (PPair::GetHybridPotential(l)).
  const std::optional<azure2::physics::HybridPotential> &Pot(int l) const;

  /// The external-capture integrals ask for F, G at the radii the adaptive quadrature chooses, the same radii for
  /// every term of one entrance l and energy.  With a memo (the caller's, used by one thread at a time), operator()
  /// keeps the waves there, keyed by the block and the exact radius: values that recomputing gives
  /// (docs/upgrade/channel_table.md section 3).  Null: none.
  void SetRadialMemo(CoulRadialMemo *memo) { radialMemo_ = memo; }

 private:
  /// The waves of the radial equation with the pair's potential (spin_orbit_potential.md section 1.1): the
  /// solution equal to H+ = G + iF beyond the matching radius, integrated inward to the radius (physics/
  /// IntegrateInward); the plain Coulomb functions at a radius beyond the matching radius.
  CoulWaves computeHybrid(int l, int two_j, double radius, double energy);
  /// The point-Coulomb functions (physics/coulomb, CoulombWavesAt).
  CoulWaves PlainWaves(int l, double radius, double energy) const;
  /// S', P' and exp(-i phi') of an interior potential at (l, two_j, radius, energy > 0): the outer functions of the
  /// waves and the regular solution's logarithmic derivative (physics/ RegularSolution, InteriorChannelFunctions).
  /// A one-deep memo: a point asks for the three in turn.
  const azure2::physics::ChannelTriple &InteriorFunctions(int l, int two_j, double radius, double energy);
  int iL_ = -1, iTwoJ_ = 0;
  double iR_ = 0.0, iE_ = 0.0;
  azure2::physics::ChannelTriple iTriple_;

  CoulRadialMemo *radialMemo_ = nullptr;
  int z1_;
  int z2_;
  int lLast_;
  int twoJLast_ = azure2::physics::kNoSpinOrbit;
  /// The pair (its spins and potential give <2 l.i> of a two_j key: PPair::SpinOrbitFactor).
  const PPair *pair_;
  double redmass_;
  /// k and eta of this pair (relativistic or not, Kinematics.h).
  Kinematics::Channel kin_;
  double radiusLast_;
  double energyLast_;
  struct CoulWaves coulLast_;

  // The universal Coulomb table (charged pairs without a nuclear potential, when
  // Config::useCoulombTable is set), with a one-deep memo of its last lookup: a
  // point asks for P, S and the phase of the same (l, radius, energy) in turn.
  const CoulombTable *table_ = nullptr;
  /// dS/dE by a 1 eV difference of the exact S (AZURE_SHIFT_DE_1EV, PairFunctionSettings).
  bool shiftDerivative1eV_ = false;
  /// The run's count of table lookups (PairFunctionSettings::tableUse).
  CoulombTableUse *tableUse_ = nullptr;
  bool TableLookup(int l, double radius, double energy);
  /// PEShift by the exact routine (never the table): one route for every point of a difference stencil.
  double PEShiftExact(int l, int two_j, double radius, double energy);
  int tLastL_ = -1;
  double tLastR_ = 0.0, tLastE_ = 0.0, tP_ = 0.0, tS_ = 0.0, tPhi_ = 0.0;
  bool tOk_ = false;

};

#endif
