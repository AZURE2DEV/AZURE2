#ifndef THMVERTEXBOUNDARY_H
#define THMVERTEXBOUNDARY_H

#include "ChannelFunc.h"
#include "Constants.h"

class ALevel;
class Config;
class JGroup;
class PPair;

/*!
 * The boundary value B of the THM entrance vertex, for one entrance channel
 * of one J group.
 *
 * Physics.  The HOES vertex is the surface term the prior-form transfer
 * amplitude leaves on the channel surface r = a after Green's theorem:
 *
 *   M_l = (B - 1) j_l(rho) - rho j_l'(rho) [+ C_l]   plane waves (EPoint),
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
 *    CNuc::CalcShiftFunctions stores on the level -- and the channel
 *    boundary constant B_c without Brune;
 *  - onshell: the log-derivative of the outgoing wave itself,
 *    L_c(E) = S_c(E) + i P_c(E) at the energy of the point (Tribble et al.
 *    2014 eq. 2.76): complex, the same for every level.
 *
 * A term without a level (the coherent background c(E) M_l, cbackground=)
 * takes the level-free value: B_c under perlevel, else as above.
 *
 * Every user of the vertex goes through this class: the model
 * (THMMatrixFunc), the vertex report (EData::ThmVertexTable, pyazr
 * thm_vertex) and the GUI diagnostics.  The model passes in its own on-shell
 * value, L_o + B_c from the point; the report and the GUI, which work on an
 * energy grid of their own, take OnShellAt(E) from the Coulomb functions.
 */
class ThmVertexBoundary {
 public:
  enum Rule { CHANNEL_CONSTANT, LEVEL_SHIFT, LOWEST_LEVEL_SHIFT, ON_SHELL };

  /// The rule of a configuration (perlevel without Brune is CHANNEL_CONSTANT).
  static Rule RuleOf(const Config &configure);

  /// Channel `channel` (1-based) of J group `jgroup`; `pair` is its particle pair.
  ThmVertexBoundary(const Config &configure, PPair *pair, JGroup *jgroup, int channel);

  Rule rule() const { return rule_; }
  bool OnShell() const { return rule_ == ON_SHELL; }

  /// Index (1-based) of the lowest in-R-matrix level of the J group, 0 if none.
  int LowestLevel() const { return lowest_; }

  /// The real B of `level` (a level of this J group); not for ON_SHELL.
  double Level(const ALevel *level) const;
  /// The real B of a term without a level; not for ON_SHELL.
  double LevelFree() const { return fixed_; }

  /// B for `level` (nullptr: a term without a level), given the on-shell value.
  complex At(const ALevel *level, complex onShell) const {
    if (rule_ == ON_SHELL) return onShell;
    return complex(level ? Level(level) : LevelFree(), 0.0);
  }
  /// L_c(E) = S_c(E) + i P_c(E) at c.m. energy E above the pair's threshold.
  complex OnShellAt(double e) const;

  /*! S_c of `pair` at a compound-system level energy (MeV), through a small
      per-thread memo keyed by every input of S_c (charges, reduced mass,
      radius, hybrid-potential state, l, energy, Coulomb routine). */
  static double ShiftAtLevelEnergy(PPair *pair, int l, double levelEnergy, bool useGSL);

 private:
  Rule rule_;
  int channel_;
  int l_;
  int lowest_;
  double fixed_;  // B_c, or S_c(E_min) for LOWEST_LEVEL_SHIFT
  mutable ChannelFunc channelFunc_;
};

#endif
