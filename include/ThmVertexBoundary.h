#ifndef THMVERTEXBOUNDARY_H
#define THMVERTEXBOUNDARY_H

#include "ChannelFunc.h"
#include "Constants.h"
#include "LevelState.h"
#include "azure2/physics/HybridPotential.h"

#include <optional>

class ALevel;
class CompoundSystem;
class Config;
class Segment;
class JGroup;
class PPair;

/*!
 * The last kSize values of S_c at level energies, keyed by every input of S_c -- charges, reduced mass, channel
 * radius, the hybrid-potential state of the pair, l, 2J, the resonance energy and the Coulomb routine -- not by the
 * PPair address (a later compound nucleus can get the same address back, and a pair's radius can change in place:
 * keyed by the address, a session at 6.1 fm read S_c of an earlier 4.1 fm session, 18O(p,a)).  Within one evaluation
 * every THM point asks for the same few (pair, l, energy).  Owned by the run (RunScratch.h: one per thread), whose
 * Coulomb-table setting is fixed.
 */
struct ThmShiftMemo {
  struct Entry {
    int z1, z2;
    double redmass, radius;
    std::optional<azure2::physics::HybridPotential> potential;
    int l, two_j;
    double e;
    double s;
  };
  static constexpr int kSize = 32;
  Entry entries[kSize];
  int filled = 0, next = 0;
};

/*!
 * The boundary value B of the THM entrance vertex, for one entrance channel
 * of one J group.
 *
 * Physics.  The HOES vertex is the surface term the prior-form transfer
 * amplitude leaves on the channel surface r = a after Green's theorem:
 *
 *   M_l = (B - 1) j_l(rho) - rho j_l'(rho) [+ C_l]   plane waves (Unit),
 *   M_l = (B - 1) a_l - d_l                          DW vertex (ThmDwVertex),
 *
 * where B multiplies the value of the transferred wave and stands for the
 * logarithmic derivative a O_l'(ka)/O_l(ka) of the outgoing wave of the
 * entrance pair.  R-matrix theory replaces that derivative by a boundary
 * value; Config::ThmOptions::vertex says which:
 *
 *  - constant: the shift S_c(E_min) at the lowest level of the J group, in
 *    either formalism, whatever the order of the levels in the file (the
 *    channel constant of the R matrix is tied to the first level read);
 *  - perlevel: the shift S_c(E_lambda) of each level under Brune (mrmpy
 *    vertex_boundary="per_level"; Tumino et al. 2021 eq. 51) -- the value
 *    CompoundSystem::CalcShiftFunctions stores in the level state -- and the channel
 *    boundary constant B_c without Brune;
 *  - onshell: the log-derivative of the outgoing wave itself,
 *    L_c(E) = S_c(E) + i P_c(E) at the energy of the point (Tribble et al.
 *    2014 eq. 2.76): complex, the same for every level.
 *
 * A term without a level (the coherent background c(E) M_l, cbackground=)
 * takes the level-free value: B_c under perlevel, else as above.
 *
 * Every user of the vertex goes through this class: the model
 * (THMMatrixFunc), the vertex report (WorkingCopy::ThmVertexTable, pyazr
 * thm_vertex) and the GUI diagnostics.  The model passes in its own on-shell
 * value, L_o + B_c from the point; the report and the GUI, which work on an
 * energy grid of their own, take OnShellAt(E) from the Coulomb functions.
 */
class ThmVertexBoundary {
 public:
  enum Rule { CHANNEL_CONSTANT, LEVEL_SHIFT, LOWEST_LEVEL_SHIFT, ON_SHELL };

  /// The rule of a configuration (perlevel without Brune is CHANNEL_CONSTANT).
  static Rule RuleOf(const Config &configure);

  /// Channel `channel` (1-based) of J group `jgroup` (number `j`, 1-based); `pair` is its particle pair; `state`
  /// the levels' fit values (the evaluation's LevelState), which the object reads while it is used.
  ThmVertexBoundary(const Config &configure, const PPair *pair, const JGroup *jgroup, const LevelState &state, int j,
                    int channel);

  Rule rule() const { return rule_; }
  bool OnShell() const { return rule_ == ON_SHELL; }

  /// Index (1-based) of the lowest in-R-matrix level of the J group, 0 if none.
  int LowestLevel() const { return lowest_; }

  /// The real B of level `la` (1-based) of this J group; not for ON_SHELL.
  double Level(int la) const;
  /// The real B of a term without a level; not for ON_SHELL.
  double LevelFree() const { return fixed_; }

  /// B for level `la` (1-based; 0: a term without a level), given the on-shell value.
  complex At(int la, complex onShell) const {
    if (rule_ == ON_SHELL) return onShell;
    return complex(la > 0 ? Level(la) : LevelFree(), 0.0);
  }
  /// L_c(E) = S_c(E) + i P_c(E) at c.m. energy E above the pair's threshold (the channel functions keep their last
  /// values, so the call is not const).
  complex OnShellAt(double e);

  /*! S_c of `pair` at a compound-system level energy (MeV); with a memo (the run's per-thread scratch), through it. */
  static double ShiftAtLevelEnergy(const PPair *pair, int l, int two_j, double levelEnergy,
                                   ThmShiftMemo *memo = nullptr);

 private:
  Rule rule_;
  const LevelState *state_;
  int j_;
  int channel_;
  int l_;
  int two_j_;  ///< the channel's 2J key (PPair::SpinOrbitTwoJ)
  int lowest_;
  double fixed_;  // B_c, or S_c(E_min) for LOWEST_LEVEL_SHIFT
  ChannelFunc channelFunc_;
};

/*!
 * The vertex uses the point-Coulomb wave of the entrance pair (ThmFunc.h); a hybrid nuclear potential on that pair
 * would bend the line shape's channel functions and not the vertex, a model that is neither.  False, with the
 * ERROR E_THM_HYBRID_POTENTIAL written, for a THM segment (or test segment, isTest) #key whose entrance pair has a
 * potential (docs/upgrade/spin_orbit_potential.md section 3.4); true otherwise.
 */
bool ThmVertexAllowsPair(const CompoundSystem *compound, const Segment &segment, int key, bool isTest, const Config &configure);

#endif
