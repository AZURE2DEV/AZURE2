#ifndef THM_OPTIONS_H
#define THM_OPTIONS_H

#include <atomic>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ThmExperiment.h"

/*!
 * An energy-dependent weight w(E) multiplying the THM (HOES) model cross
 * section of one segment, read from a two-column file named by the
 * `weight[<segment key>]=` key of the <thm> block (Config::ReadThmBlock).
 * It is the hook for a correction AZURE2 does not compute itself, typically
 * Mukhamedzhanov's Coulomb-distortion factor R(E) = |M(E)|^2/|M(E_norm)|^2
 * (docs/source/theory/thm_implementation.rst).  Interpolated log-linearly
 * (linear in E, linear in ln w); constant beyond the ends of the table.
 */
struct ThmWeightTable {
  /// The file as named in the <thm> block, and as resolved against the .azr.
  std::string name, path;
  /// Strictly increasing E_cm (MeV) of the THM entrance pair, w(E) and ln w(E).
  std::vector<double> e, w, lnw;
  /// Set once the first evaluation outside [e.front(), e.back()] has been reported.
  mutable std::atomic<bool> warned{false};
  /// Reads the file; returns an empty string on success, else the reason.
  std::string Read(const std::string &file);
  /// w(E); *outside (if given) is set when E lies beyond the table.
  double operator()(double energy, bool *outside = nullptr) const;
  bool Covers(double energy) const { return energy >= e.front() && energy <= e.back(); }
};

/// Choices in the THM (modified R-matrix, HOES) observable, read from the
/// optional <thm> block.  The defaults reproduce the Typel--Baur / La Cognata
/// working formula (Tumino et al. 2021, eqs. 50-51) that the code has always
/// used, except that entrance partial waves are summed incoherently.
struct ThmOptions {
  /// Boundary in the entrance vertex (B - 1) j_l - rho j_l' (key `vertex=`):
  ///  perlevel  S_c(E_lambda) of each level under Brune, as mrmpy (alias `real`);
  ///  constant  B_c = S_c(E_1) at the lowest-energy level of the J group,
  ///            applied after the level sum.  gamma^T A gamma is
  ///            the same in the Brune and formal representations, so this is
  ///            the formal-R-matrix vertex at that B (La Cognata et al.,
  ///            ApJ 723 (2010) 1512 use B = S(E_1));
  ///  onshell   L_c(E) = S_c(E) + i P_c(E) (Tribble et al. 2014 eq. 2.76),
  ///            independent of any boundary constant.
  /// perlevel and constant coincide for a J group with a single level.
  enum Vertex { PER_LEVEL, CONSTANT, ON_SHELL };
  Vertex vertex = CONSTANT;
  /// Which kinematic factors the model carries, set by how the HOES data
  /// were extracted from the triple cross section (key `kinematics=`):
  ///  lacognata  exit k_f/mu_f, as La Cognata's working formula (Tumino et
  ///             al. 2021 eq. 50) -- the default, what the code always did;
  ///  triple     exit Gamma_f alone: raw d3sigma / |phi|^2 (Mukhamedzhanov
  ///             2017 eq. 34);
  ///  kf3body    exit Gamma_f / (mu_f k_f): data divided by the full
  ///             three-body phase-space KF (Typel & Baur 2003 eq. 16);
  ///  lambda32   exit Gamma_f, entrance 1/k_i on shell: data divided by
  ///             KF = lambda3/lambda2 (Pizzone et al. 2011; Tumino eq. 22).
  enum Kinematics { LA_COGNATA, TRIPLE, KF_THREE_BODY, LAMBDA32 };
  Kinematics kinematics = LA_COGNATA;
  /// Sum entrance partial waves of one channel spin coherently, as mrmpy and
  /// AZURE2 before Sep 2026 did.  For an observable integrated over the exit
  /// direction the l cross terms vanish, so incoherent is the default.
  /// Key `entranceL=coherent|incoherent`.
  bool coherentL = false;
  /// Add the external Coulomb term 2 eta k \int_a^inf O_l(kr)/O_l(ka) j_l(pr) dr to
  /// the vertex (Tribble eq. 2.79; Typel & Baur eq. A.4).  Key `coulombIntegral=0|1`.
  bool coulombIntegral = false;
  /// Mean spectator kinetic energy <p_sx^2>/2mu_sx (MeV) added to E + B in the
  /// half-off-shell momentum (Typel & Baur eq. 11).  Key `spectatorEnergy=` for
  /// every THM pair, `spectatorEnergy[<pair key>]=` for one.
  double spectatorEnergy = 0.0;
  std::map<int, double> spectatorEnergyByPair;
  double SpectatorEnergy(int pairKey) const {
    std::map<int, double>::const_iterator it = spectatorEnergyByPair.find(pairKey);
    return it == spectatorEnergyByPair.end() ? spectatorEnergy : it->second;
  }
  /// Energy-dependent weight of the THM model cross section, per segment:
  /// key `weight[<k>]=<file>` for the k-th line of <segmentsData> (counting
  /// inactive lines, the numbering of segment_<k>_norm and the output
  /// files), `weightTest[<k>]=<file>` for the k-th line of <segmentsTest>.
  /// A relative path is taken from the directory of the .azr file.
  std::map<int, std::shared_ptr<const ThmWeightTable>> weightBySegment, weightByTestSegment;
  /// THM experiments, `experiment[<name>] key=value ...` (ThmExperiment.h):
  /// segments sharing one profiled norm and an optional background.  Empty
  /// unless the block has such lines; then nothing changes.
  std::vector<ThmExperiment> experiments;
  /// Set once a THM point below E = -B has been reported (once per session:
  /// copies of this Config share it).
  std::shared_ptr<std::atomic<bool>> belowBWarned = std::make_shared<std::atomic<bool>>(false);
};

#endif
