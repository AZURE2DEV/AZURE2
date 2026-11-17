#ifndef PPAIR_H
#define PPAIR_H

#include <optional>
#include <string>

#include "azure2/physics/HybridPotential.h"

namespace azure2 {
namespace io {
class LevelRow;
}
}  // namespace azure2

struct CoulombTableUse;

/*!
 * What the run's Config says about how a pair's functions are computed: the same for every pair of a run, given when
 * the pair is built (CompoundSystem::Fill).  The readers deep in the channel functions (Kinematics::Channel, CoulFunc, ShftFunc,
 * WhitFunc) take it from their pair, never from a process-wide pointer, so two runs in one process (two pyazr sessions,
 * the web UI's server) with different settings never read each other's.
 */
struct PairFunctionSettings {
  /// Relativistic kinematics ("relativistic=1" in the <config> block, Config::relativistic; Kinematics.h).
  bool relativistic = false;
  /// The universal Coulomb table where it covers the pair (Config::useCoulombTable; --no-coulomb-table turns it off).
  bool useCoulombTable = false;
  /// The table's file (Config::coulombTablePath, "" = coulombTableDefault), loaded on first use (CoulombTable::Shared).
  std::string coulombTablePath;
  /// The table's file when none is named (Config::environment.coulombTableDefault: AZURE2_COULOMB_TABLE, else
  /// ~/.azure2/coulomb_table.bin).
  std::string coulombTableDefault;
  /// The older routes of the shift function, from the run's environment (Config::environment): the bound-state S
  /// from the numerical derivative (AZURE_BOUND_SHIFT_NUMERICAL), the neutral bound S without its closed form
  /// (AZURE_NO_NEUTRAL_CLOSED_FORM), dS/dE by a 1 eV difference (AZURE_SHIFT_DE_1EV).
  bool boundShiftNumerical = false;
  bool noNeutralClosedForm = false;
  bool shiftDerivative1eV = false;
  /// The run's count of Coulomb-table lookups (RunDiagnostics::coulombTable); null: not counted.
  CoulombTableUse *tableUse = nullptr;
};

/// An AZURE Particle Pair

/*!
 * In R-Matrix theory, the configuration space in the external region is decomposed into combinations of particle
 * pairs, traditionally given by the symbol \f$ \alpha \f$.  In AZURE, these particle pair are represented by a PPair object.
 * The reaction pathways between pairs are in the pathway table (PathwayTable.h).
 */

class PPair {
 public:
  /// Build from a line of the nuclear input file, with the run's settings of its functions.
  PPair(const azure2::io::LevelRow &, const PairFunctionSettings &);
  /// Is this the entrance pair the calculation is driven from?
  bool IsEntrance() const;
  /// Charge number of particle 1 or 2.
  int GetZ(int) const;
  /// Parity of particle 1 or 2, \f$\pm1\f$.
  int GetPi(int) const;
  /// Pair type: 0 particle+particle, 10 particle+gamma.
  int GetPType() const;
  /// The pair's key as written in the input file -- its 1-based position, and
  /// what an AChannel's pair number refers to.
  int GetPairKey() const;
  /// Charge number of the compound nucleus this pair forms, inferred from the
  /// pair columns: Z1+Z2 for a particle pair, Z2 for a gamma pair, Z2-Z1 for a
  /// beta pair (Z1 = -1 for beta-minus).
  int GetNucleusZ() const;
  /// Mass number of the compound nucleus: round(M1+M2) for a particle pair,
  /// round(M2) otherwise.
  int GetNucleusA() const;
  /// 1-based index of the compound nucleus in the CompoundSystem's nucleus list, set by
  /// CompoundSystem::AddPair; 0 until the pair is added.
  int GetNucleus() const;
  void SetNucleus(int);
  /// Mass of particle 1 or 2, in u.
  double GetM(int) const;
  /// g-factor of particle 1 or 2.
  double GetG(int) const;
  /// Intrinsic spin of particle 1 or 2.
  double GetJ(int) const;
  /// Total excitation energy of the pair (MeV), e1 + e2 -- light plus heavy
  /// particle.  The channel threshold is GetSepE() + GetExE().
  double GetExE() const;
  /// Excitation energy of the light particle (MeV), the optional 32nd column
  /// of a <levels> line (e1; 0 when absent).  Always 0 for gamma / beta pairs.
  double GetLightExE() const;
  /// Excitation energy of the heavy particle (MeV), column 17 (e2); for a
  /// gamma pair the energy of the final level.
  double GetHeavyExE() const;
  /// Separation energy of the pair (MeV). Subtract it, and the excitation
  /// energy, from a level's energy to get that channel's centre-of-mass energy.
  double GetSepE() const;
  /// The pair's radius (fm), column 28 of a level line: the radius a channel of the pair takes at set-up
  /// (AChannel::SetRadius).  A channel quantity -- penetrability, shift function, hard-sphere phase, boundary
  /// condition, Wigner limit, a bound-state integral -- reads AChannel::GetRadius(), not this.  This one is for
  /// what belongs to the pair and to no channel: the beta-decay Fermi function, the THM vertex's form factor, the
  /// signatures of the tables, the reports' "Channel Radius" line.
  double GetPairRadius() const;
  /// Change the channel radius. Everything derived from it -- penetrabilities,
  /// shift functions, boundary conditions, external-capture integrals -- must
  /// be recomputed afterwards.
  void SetChRad(double);
  /// THM binding energy B_xs (MeV) of the transferred particle in the
  /// Trojan-Horse nucleus; enters the half-off-shell momentum of the
  /// entrance transfer form factor. 0 for a conventional pair.
  double GetBindingEnergy() const;
  /// Reduced mass of the pair, in u.
  double GetRedMass() const;
  /// \f$1/[(2I_1+1)(2I_2+1)]\f$, the spin-statistics denominator of the entrance channel.
  double GetI1I2Factor() const;
  /// Returns true if the two particles in the pair are identical
  /// (same Z, mass, intrinsic spin, parity, g-factor and excitation energy,
  /// e1 = e2), for a particle pair only.
  /// Detected automatically in the constructor.
  bool IsIdentical() const;
  /// Boson/fermion sign: +1 for identical bosons (2j even), -1 for
  /// identical fermions. Returns +1 when IsIdentical() is false.
  int GetIdenticalSign() const;
  /// True for an identical pair whose members carry spin (j != 0).  For these
  /// the exchange term of the elastic amplitude depends on the channel spin s
  /// through (-1)^s, so no single sign symmetrizes it -- the spin-0 shortcut
  /// (one Coulomb amplitude f_C(theta) + f_C(pi - theta) for every channel)
  /// does not apply.
  bool HasSpinDependentExchange() const;
  /// \f$1+\delta_{12}\f$: 2 for an identical pair, 1 otherwise.  An identical
  /// pair in the ENTRANCE channel multiplies every cross section out of it --
  /// reaction, capture or elastic, angle-integrated or differential -- by this
  /// factor, sigma = (pi/k^2) sum_J g_J (1+delta_12) sum |T|^2 over the channels
  /// symmetry allows -- for spin-0 bosons the familiar
  /// sigma = (pi/k^2) sum_l (2l+1) [1+(-1)^l] T_l (e.g. Hagino, Rowley and Yao,
  /// EPJ Web Conf. 122 (2016) 07002).  It is the convention of the reciprocity
  /// theorem w_12 k_12^2 sigma_12->34 / (1+delta_12) =
  /// w_34 k_34^2 sigma_34->12 / (1+delta_34), w = (2i_1+1)(2i_2+1)
  /// (Blatt and Weisskopf 1952; Fowler,
  /// Caughlan and Zimmerman, ARA&A 5 (1967) 525) with a symmetric U.  An
  /// identical pair in the EXIT channel contributes no factor: the cross section
  /// counts each reaction once, not each outgoing particle.  The elastic
  /// differential cross section, which by convention counts particles, carries
  /// a further 2 of its own (GenMatrixFunc).
  double GetEntranceSymmetryFactor() const;
  /// The nuclear potential of the hybrid model for this pair (spin_orbit_potential.md section 1.1), none for a pair
  /// with plain Coulomb functions.  Set by CompoundSystem::Fill from Config::potentials.
  const std::optional<azure2::physics::HybridPotential> &GetHybridPotential() const { return hybrid_potential_; }
  /// Throws std::invalid_argument (with the code of experiment/SpinOrbitChecks.h) for a spin-orbit term the pair cannot
  /// have: both particles spinning, or the term on a particle of spin 0.
  void SetHybridPotential(const std::optional<azure2::physics::HybridPotential> &p);
  /// The 2J of the keys of this pair's channel functions (channel table, radial scratch, shift-derivative memo, EC
  /// scratch) for a channel of total J: 2J when the pair's potential has a spin-orbit term (its functions depend on
  /// J), azure2::physics::kNoSpinOrbit otherwise.  Every call site forms the key with this function.
  int SpinOrbitTwoJ(double J) const;
  /// <2 l.i> of the pair's spin-orbit term for a channel of orbital l and key two_j (SpinOrbitTwoJ): 0 for
  /// kNoSpinOrbit, else j(j+1) - l(l+1) - i(i+1) with i the spin of the particle the term acts on.
  double SpinOrbitFactor(int l, int two_j) const;
  void SetEntrance();
  /// The run's settings of this pair's functions (relativistic kinematics, the Coulomb table).
  const PairFunctionSettings &FunctionSettings() const { return function_settings_; }

 private:
  bool entrance_;
  bool is_identical_;
  int identical_sign_;
  int pair_z_[2];
  int pair_pi_[2];
  int pair_ptype_;
  int pair_key_;
  int nucleus_z_;
  int nucleus_a_;
  int nucleus_;
  double pair_m_[2];
  double pair_g_[2];
  double pair_j_[2];
  double pair_light_ex_e_;
  double pair_heavy_ex_e_;
  double pair_sep_e_;
  double pair_ch_rad_;
  double binding_energy_;
  double red_mass_;
  double i1i2factor_;
  std::optional<azure2::physics::HybridPotential> hybrid_potential_;
  PairFunctionSettings function_settings_;
};

#endif
