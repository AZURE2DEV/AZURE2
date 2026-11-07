#ifndef AZURE2_PHYSICS_CONVENTION_H
#define AZURE2_PHYSICS_CONVENTION_H

namespace azure2 {
namespace physics {

/// What the level energies E_lambda and amplitudes gamma_{lambda c} of a project are.
enum class LevelBasis {
  Formal,  ///< formal R-matrix parameters (Lane and Thomas), with constant boundary values B_c
  Brune,   ///< Brune's alternative parameters (PRC 66, 044611 (2002)): observed energies, B_c = S_c(E_lambda)
  Park     ///< Park's observed-amplitude parameters (PRC 104, 064612 (2021)), on Brune's level energies
};

/// How the channel spin and the orbital angular momentum of a particle channel are coupled.
enum class CouplingScheme {
  LS,  ///< channel spin: s = j_a + j_b, then J from s and l (the internal scheme)
  JJ   ///< j = l + j_x of one particle x, then J from j and the other particle's spin j_y
};

/// Which particle's spin comes first.  LS: s = j1 + j2 (LightFirst, internal) or s = j2 + j1 (HeavyFirst).
/// JJ: the particle whose spin couples with l first, x = light (LightFirst: j = l + j1, J = j + j2) or heavy.
enum class SpinOrder { LightFirst, HeavyFirst };

/// Order of the orbital coupling.  LS: J = s + l (SpinFirst, internal; (s nu l m | J M)) or J = l + s (OrbitalFirst).
/// JJ: j = j_x + l (SpinFirst) or j = l + j_x (OrbitalFirst); J = j + j_y in both.
enum class OrbitalOrder { SpinFirst, OrbitalFirst };

/// Phase of the orbital part of the channel surface functions.  None: Y_l^m (internal); IPowL: i^l Y_l^m;
/// IPowMinusL: i^-l Y_l^m.  For real amplitudes the two i^(+-l) conventions give the same real signs
/// (level_parameters.md section 5), so both convert alike.
enum class OrbitalPhase { None, IPowL, IPowMinusL };

/// Coulomb-phase reference of the collision matrix: U relative to e^{2 i sigma_0} (SigmaLMinusSigma0, internal,
/// omega_c = sigma_l - sigma_0) or the nuclear S-matrix with the full Coulomb phase (SigmaL).  The reference changes
/// U by a phase common to each partition and no reduced width amplitude; the conversion of gamma is the identity.
enum class CoulombPhaseReference { SigmaLMinusSigma0, SigmaL };

/// Where the constant boundary value B_c of a J^pi block of the formal basis is taken.
enum class BoundaryRule {
  FirstListedLevel,  ///< B_c = S_c(E of the first level listed in the block): the engine's rule (default)
  LowestLevel        ///< B_c = S_c(E of the lowest level of the block): proposal, conventions.rst; OFF by default
};

/// Sign of each formal level's amplitude vector, which the eigenproblem of the Brune transform leaves free.
enum class FormalSignRule {
  AsComputed,      ///< the eigensolver's choice (SymmetricEigen.h: largest eigenvector component positive) (default)
  LargestPositive  ///< the amplitude of largest |gamma| of each formal level is positive: proposal; OFF by default
};

/*!
 * What a project's level parameters mean (docs/upgrade/ARCHITECTURE.md section 5; docs/source/theory/conventions.rst,
 * "Convention record").  The engine computes in one internal convention, the default-constructed record: channel-spin
 * coupling s = j1 + j2 with the light particle first, J = s + l, no i^l, the engine's photon multipole phases,
 * omega_c = sigma_l - sigma_0, the Brune basis, B_c from the first listed level, formal signs as computed.
 * level_parameters (LevelParameters.h) converts a project's parameters from its convention to the internal one; that
 * is the only place a conversion happens (docs/upgrade/RULES.md section 1.6).
 */
struct Convention {
  CouplingScheme coupling = CouplingScheme::LS;
  SpinOrder spin_order = SpinOrder::LightFirst;
  OrbitalOrder orbital_order = OrbitalOrder::SpinFirst;
  OrbitalPhase orbital_phase = OrbitalPhase::None;
  /// Sign by which the project's amplitudes of electric (E) and magnetic (M) photon channels differ from the
  /// engine's: gamma_internal = sign * gamma_project.  The engine's own multipole phases are those of conventions.rst
  /// Eqs. conv-ec-e and conv-ec-m; a source whose E or M amplitudes carry the opposite real sign sets -1 here.
  int electric_sign = 1;
  int magnetic_sign = 1;
  CoulombPhaseReference coulomb_phase = CoulombPhaseReference::SigmaLMinusSigma0;
  LevelBasis basis = LevelBasis::Brune;
  BoundaryRule boundary = BoundaryRule::FirstListedLevel;
  FormalSignRule formal_sign = FormalSignRule::AsComputed;
};

/// True when the record differs from the internal convention only in the parameter basis and the two proposals
/// (boundary, formal_sign): the gamma columns need no recoupling or sign change.
inline bool is_internal_coupling(const Convention &c) {
  return c.coupling == CouplingScheme::LS && c.spin_order == SpinOrder::LightFirst &&
      c.orbital_order == OrbitalOrder::SpinFirst && c.orbital_phase == OrbitalPhase::None && c.electric_sign == 1 &&
      c.magnetic_sign == 1;
}

}  // namespace physics
}  // namespace azure2

#endif
