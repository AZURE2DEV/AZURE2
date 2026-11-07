#ifndef AZURE2_PHYSICS_LEVELPARAMETERS_H
#define AZURE2_PHYSICS_LEVELPARAMETERS_H

#include <complex>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "azure2/physics/Convention.h"
#include "azure2/physics/System.h"

namespace azure2 {
namespace physics {

/*!
 * level_parameters: the conversions between a project's level parameters and the internal ones (formal R-matrix
 * parameters, Lane-Thomas, channel-spin coupling s = j1 + j2), docs/upgrade/ARCHITECTURE.md section 5 and
 * docs/upgrade/level_parameters.md.  Every function works on one J^pi block b of a System; the System gives the
 * channels' partitions, spins and thresholds and is not modified.
 *
 * Basis conversions (reproducing CNuc::TransformIn, CNuc::CalcBoundaryConditions, CNuc::TransformOut and
 * CNuc::ConvertAmplitudeBasis):
 *   amplitudes_from_observed   observed widths (eV), ANCs (fm^-1/2), amplitudes -> Brune or Park amplitudes
 *   observed_from_amplitudes   the inverse
 *   boundary_values            B_c of the formal basis
 *   subtract_external_capture  total photon amplitude -> internal one (external capture given by the caller)
 *   formal_from_brune          generalized eigenproblem N a = E M a (Brune, PRC 66, 044611 (2002))
 *   brune_from_formal          the inverse, by the engine's fixed-point iteration
 *   park_from_brune, brune_from_park   gamma_Park = sqrt(J) gamma_Brune (Park, PRC 104, 064612 (2021))
 * Convention conversions:
 *   project_channels, recoupling_gaps, convention_transform: an orthogonal matrix on the gamma columns of a block
 *   ((l, s) <-> jj by 6j symbols, coupling-order and i^l signs, photon E/M signs) and its Jacobian.
 * Composite:
 *   level_parameters           project parameters in any convention -> internal formal E, gamma, B.
 *
 * The channel functions are an input: physics/ has no Coulomb functions yet (channel_functions waits for
 * pycoulomb), so the caller supplies S_c, dS_c/dE and the width factor P_c at any energy through ChannelFunctions.
 * The engine's adapter (src/CNucLevelParameters.cpp) supplies them from ChannelFunc and WhitFunc.
 */

/// Channel functions of one block channel at one level energy, supplied by the caller.
struct ChannelValues {
  /// Width factor P with Gamma = 2 P gamma^2 (in Brune's basis, with the level's normalization):
  ///  - particle channel above threshold (E_c > 0): the penetrability P_c(E_c) (dimensionless);
  ///  - particle channel below threshold (E_c < 0): mu a_c / (hbar^2 W_l(2 kappa a_c)^2) (MeV^-1 fm^-1), so that
  ///    an ANC C (fm^-1/2) gives gamma^2 = C^2 / (2 P) at first order (conventions.rst Eq. conv-brune);
  ///  - photon channel: (E_gamma / hbar c)^(2L+1) (fm^-(2L+1)), or 1 with --use-rmc; for a static moment
  ///    (static_moment true) the engine's factor (CNuc.cpp, "Static moments" in conventions.rst);
  ///  - beta channel: 1 (not used).
  double P = 0.0;
  double S = 0.0;      ///< shift function S_c(E_c) (dimensionless), particle channels
  double dS_dE = 0.0;  ///< dS_c/dE (MeV^-1), particle channels
  /// Multiplies the input value before it is squared: 1, except an E2 static moment (100 (alpha hbar c)^1/2).
  double input_scale = 1.0;
  /// Photon channel to the level's own state (|E - E_f| < 1 keV, same J^pi): a static M1 or E2 moment.
  bool static_moment = false;
};

/// (block b, block channel k, level energy E_x in MeV, excitation energy of the compound nucleus) -> values.
using ChannelFunctions = std::function<ChannelValues(std::size_t b, std::size_t k, double E_MeV)>;

/// A gamma column as the project labels it.  In the (l, s) scheme two_x = 2s; in the jj scheme, for particle
/// channels, two_x = 2j with j the coupling of l and the spin of particle x (Convention.h).  Photon and beta channels
/// are labelled as in the System (two_x = 2s).
struct ProjectChannel {
  int partition = 0;
  ChannelKind kind = ChannelKind::Particle;
  int l = 0;
  int two_x = 0;
};

/// The parameters of one block as the project gives them.
struct ProjectLevels {
  enum class Form {
    Observed,   ///< as in <levels>: Gamma (eV), ANC C (fm^-1/2), photon Gamma (eV), beta amplitude; gammaIsRWA
                ///< entries are Brune amplitudes (MeV^1/2)
    Amplitudes  ///< reduced width amplitudes of the convention's basis (MeV^1/2), as with --no-transform
  };
  Form form = Form::Observed;
  std::vector<double> E;                    ///< E_lambda (MeV, E_x), in the order listed
  std::vector<double> value;                ///< row-major, value[lambda * n_columns + m]
  std::vector<unsigned char> is_amplitude;  ///< gammaIsRWA flags, same layout (Observed form only; may be empty)
};

/// External-capture amplitude of one (level, photon channel), supplied by the caller (CNuc::CalcExternalWidth).
struct ExternalWidth {
  bool present = false;
  std::complex<double> value{0.0, 0.0};
};

/// The block's own channels as ProjectChannel labels (the internal convention).
std::vector<ProjectChannel> internal_channels(const System &system, std::size_t b);

// ---------------------------------------------------------------------------------------------------- bases

struct BasisAmplitudes {
  std::vector<double> gamma;      ///< row-major [lambda * n_columns + m], MeV^1/2 (Brune or Park)
  std::vector<double> park_norm;  ///< Park's J_lambda (Eq. 28 of Park 2021) per level
  std::vector<std::string> warnings;
};

/*!
 * Observed widths and ANCs -> Brune (basis Brune) or Park (basis Park) amplitudes (CNuc.cpp, TransformIn, first
 * loop).  columns label the value columns; each must have a block channel with the same partition, kind and l from
 * which the channel functions are taken (any channel spin: P_c, S_c depend on (alpha, l) only).
 */
BasisAmplitudes amplitudes_from_observed(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                                         const ProjectLevels &levels, LevelBasis basis, const ChannelFunctions &cf);

/// Brune or Park amplitudes -> the observed form (the inverse of amplitudes_from_observed; is_amplitude marks
/// the columns to report as Brune amplitudes).  Returns a ProjectLevels in Form::Observed.
ProjectLevels observed_from_amplitudes(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                                       const std::vector<double> &E, const std::vector<double> &gamma,
                                       const std::vector<unsigned char> &is_amplitude, LevelBasis basis,
                                       const ChannelFunctions &cf);

/// B_c of the block's channels for the formal basis: S_c at the first listed level (rule FirstListedLevel, the
/// engine's) or at the lowest level (LowestLevel).  Non-particle channels take the value of block channel 0 (0 when
/// channel 0 is not a particle channel), as the engine does; they do not enter the level matrix.  No levels: 0.
std::vector<double> boundary_values(const System &system, std::size_t b, const std::vector<double> &E,
                                    BoundaryRule rule, const ChannelFunctions &cf);

/// Photon channels: total amplitude -> internal one, gamma_int = sign(gamma) (gamma^2 - Im^2)^1/2 - Re of the
/// external amplitude (CNuc.cpp, TransformIn, second loop).  ext has the layout of gamma; entries not present are
/// left alone, and so are |gamma| < 1e-8 when ignore_zero_widths.
void subtract_external_capture(const System &system, std::size_t b, std::vector<double> &gamma,
                               const std::vector<ExternalWidth> &ext, bool ignore_zero_widths,
                               std::vector<std::string> &warnings);

struct LevelSet {
  std::vector<double> E;      ///< MeV
  std::vector<double> gamma;  ///< row-major, block channels
  std::vector<int> iterations;  ///< brune_from_formal: iterations per level (0 for the others)
  std::vector<std::string> warnings;
};

/*!
 * Brune -> formal with constant B (CNuc.cpp, TransformIn, second loop).  Level lambda keeps its position: it
 * receives the eigenvalue of the rank of its Brune energy.  Throws std::domain_error when M is not positive definite
 * (no formal parameters correspond) and std::invalid_argument for two levels at the same energy.
 */
LevelSet formal_from_brune(const System &system, std::size_t b, const std::vector<double> &E,
                           const std::vector<double> &gamma, const std::vector<double> &B, const ChannelFunctions &cf,
                           FormalSignRule sign = FormalSignRule::AsComputed);

struct FormalToBruneOptions {
  double energy_tolerance_MeV = 1e-6;  ///< the engine's (CNuc::TransformOut)
  int max_iterations = 1000;           ///< the engine's
};

/// Formal -> Brune by the engine's iteration (CNuc::TransformOut): per level, B_c moved to S_c(E) and the matrix
/// E - sum_c (S_c - B_c) gamma gamma rediagonalized until the level's energy moves by <= the tolerance.  A level
/// that does not converge keeps its formal values, with a warning, as in the engine.
LevelSet brune_from_formal(const System &system, std::size_t b, const std::vector<double> &E,
                           const std::vector<double> &gamma, const std::vector<double> &B, const ChannelFunctions &cf,
                           const FormalToBruneOptions &options = FormalToBruneOptions());

/// Brune -> Park: every amplitude of a level times sqrt(J), J = 1 / (1 + sum_c gamma_c^2 dS_c/dE).
std::vector<double> park_from_brune(const System &system, std::size_t b, const std::vector<double> &E,
                                    const std::vector<double> &gamma, const ChannelFunctions &cf,
                                    std::vector<double> *park_norm = nullptr);
/// Park -> Brune: divided by sqrt(J), J = 1 - sum_c gamma_c^2 dS_c/dE; throws std::domain_error for J <= 0.
std::vector<double> brune_from_park(const System &system, std::size_t b, const std::vector<double> &E,
                                    const std::vector<double> &gamma, const ChannelFunctions &cf,
                                    std::vector<double> *park_norm = nullptr);

// ---------------------------------------------------------------------------------------------- conventions

/// <internal (l, s) channel | project channel> for one particle channel of a partition with spins j1 (light),
/// j2 (heavy), orbital l and total J: a 6j recoupling (jj) or a sign (LS).  two_x is 2s (LS) or 2j (jj).
double recoupling_coefficient(int two_j1, int two_j2, int l, int two_s, int two_J, const Convention &convention,
                              int two_x);

/// The block's channels relabeled in the project's convention: LS keeps the labels; jj replaces each (alpha, l)
/// group of particle channels by every j allowed with J (ascending).  Photon and beta channels are kept.
std::vector<ProjectChannel> project_channels(const System &system, std::size_t b, const Convention &convention);

/// What a recoupling between the block's channels and the project's columns would lose: one line per channel
/// that has a non-zero coefficient with a channel the other side does not have, or a project column with no
/// counterpart in the block.  Empty: the transformation is orthogonal on the columns given.
std::vector<std::string> recoupling_gaps(const System &system, std::size_t b, const Convention &convention,
                                         const std::vector<ProjectChannel> &columns);

/*!
 * gamma_internal[lambda, k] = sum_m T[k, m] gamma_project[lambda, m].  T is block-orthogonal (one block per (alpha, l)
 * group of particle channels, a +-1 for photon and beta channels), so the inverse is T^T and
 * d gamma_internal / d gamma_project = T (jacobian()).  T does not depend on the parameters.
 */
class ConventionTransform {
 public:
  ConventionTransform(std::size_t n_internal, std::size_t n_project, std::vector<double> T) :
    n_internal_(n_internal), n_project_(n_project), T_(std::move(T)) {}
  std::size_t num_internal() const { return n_internal_; }
  std::size_t num_project() const { return n_project_; }
  double operator()(std::size_t k, std::size_t m) const { return T_[k * n_project_ + m]; }
  /// The Jacobian d gamma_internal[k] / d gamma_project[m] (row-major n_internal x n_project), = T.
  const std::vector<double> &jacobian() const { return T_; }
  /// gamma rows (levels) from project columns to internal columns.
  std::vector<double> to_internal(const std::vector<double> &gamma_project, std::size_t n_levels) const;
  /// gamma rows from internal columns to project columns (T^T).
  std::vector<double> to_project(const std::vector<double> &gamma_internal, std::size_t n_levels) const;

 private:
  std::size_t n_internal_, n_project_;
  std::vector<double> T_;
};

/// Throws std::invalid_argument, naming every gap (recoupling_gaps), when the transformation is not orthogonal.
ConventionTransform convention_transform(const System &system, std::size_t b, const Convention &convention,
                                         const std::vector<ProjectChannel> &columns);

// ------------------------------------------------------------------------------------------------ composite

struct InternalParameters {
  std::vector<double> E;            ///< formal E_lambda (MeV), in the listed order
  std::vector<double> gamma;        ///< formal gamma_{lambda c} (MeV^1/2), block channels
  std::vector<double> B;            ///< B_c of the block channels
  LevelBasis basis = LevelBasis::Brune;  ///< basis of basis_gamma (the project's, Brune for an observed formal one)
  std::vector<double> basis_E;      ///< level energies in that basis (MeV)
  std::vector<double> basis_gamma;  ///< amplitudes in that basis, internal coupling, after external capture
  std::vector<double> park_norm;    ///< Park's J per level (Observed form)
  std::vector<std::string> warnings;
};

struct LevelParameterOptions {
  const std::vector<ExternalWidth> *external = nullptr;  ///< photon external-capture amplitudes (internal columns)
  bool ignore_zero_widths = false;                       ///< the engine's IGNORE_ZERO_WIDTHS for that subtraction
  bool formal = true;  ///< compute the formal parameters (false: basis_* only, E/gamma left empty)
};

/*!
 * A project's parameters of block b, in the convention of the record and with the given columns, -> the internal
 * formal parameters.  Observed form: widths -> Brune (or Park) amplitudes in the project's columns; recoupling to
 * the block's channels; external capture subtracted; B_c by the record's rule; Park -> Brune; Brune -> formal with
 * the record's sign rule.  Amplitudes form in the formal basis: recoupling only.
 */
InternalParameters level_parameters(const System &system, std::size_t b, const Convention &convention,
                                    const std::vector<ProjectChannel> &columns, const ProjectLevels &levels,
                                    const ChannelFunctions &cf,
                                    const LevelParameterOptions &options = LevelParameterOptions());

}  // namespace physics
}  // namespace azure2

#endif
