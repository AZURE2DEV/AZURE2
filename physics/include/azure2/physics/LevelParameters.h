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
 * Basis conversions (called by CompoundSystem::TransformIn, CompoundSystem::CalcBoundaryConditions, CompoundSystem::TransformOut; reproducing
 * CompoundSystem::ConvertAmplitudeBasis):
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
 * The engine's adapter (src/EngineChannelFunctions.cpp) supplies them from ChannelFunc and WhitFunc.
 */

/// Channel functions of one block channel at one level energy, supplied by the caller.
struct ChannelValues {
  /// Width factor P with Gamma = 2 P gamma^2 (in Brune's basis, with the level's normalization):
  ///  - particle channel above threshold (E_c > 0): the penetrability P_c(E_c) (dimensionless);
  ///  - particle channel below threshold (E_c < 0): mu a_c / (hbar^2 W_l(2 kappa a_c)^2) (MeV^-1 fm^-1), so that
  ///    an ANC C (fm^-1/2) gives gamma^2 = C^2 / (2 P) at first order (conventions.rst Eq. conv-brune);
  ///  - photon channel: (E_gamma / hbar c)^(2L+1) (fm^-(2L+1)), or 1 with --use-rmc; for a static moment
  ///    (static_moment true) the engine's factor (LevelTransforms.cpp, "Static moments" in conventions.rst);
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


/// External-capture amplitude of one (level, photon channel), supplied by the caller (ExternalWidthsOfInputs).
struct ExternalWidth {
  bool present = false;
  std::complex<double> value{0.0, 0.0};
};

/// The block's own channels as ProjectChannel labels (the internal convention).
std::vector<ProjectChannel> internal_channels(const System &system, std::size_t b);

// ---------------------------------------------------------------------------------------------------- bases

/// Why amplitudes_from_observed did not convert a column's input as a width (one per level and column).
enum class ColumnNote : unsigned char {
  None,
  ZeroPenetrability,     ///< particle channel just above a Coulomb threshold, P = 0 in double precision: gamma = 0
  AtThreshold,           ///< particle channel exactly at threshold, where an ANC defines no width: gamma = 0
  AmplitudeFlagIgnored   ///< the amplitude flag (gammaIsRWA) on a photon or beta column: read as for no flag
};

/// One column's share of Brune's denominator D_lambda (Eq. 2.1 of level_parameters.md): the term |Gamma_c| S'_c / P_c
/// (|Gamma_c| in MeV; C_c^2 and P^b_c for a bound channel) with the P_c and S'_c (MeV^-1) it was formed from.
struct DenominatorTerm {
  bool present = false;  ///< false for columns that add nothing to D (amplitude input, photon, beta, P = 0)
  double term = 0.0;
  double P = 0.0;
  double dS_dE = 0.0;
};

struct BasisAmplitudes {
  std::vector<double> gamma;      ///< row-major [lambda * n_columns + m], MeV^1/2 (Brune or Park)
  std::vector<double> park_norm;  ///< Park's J_lambda (Eq. 28 of Park 2021) per level
  std::vector<std::string> warnings;
  /// Diagnostics, so that a caller can word its own messages: D_lambda / 2 per level (negative: no reduced width
  /// reproduces the widths), and per (level, column), row-major, each column's term of D and its ColumnNote.  A note
  /// ZeroPenetrability or AtThreshold is set only when the input was not 0 (an input of 0 gives gamma = 0 exactly).
  std::vector<double> half_denominator;
  std::vector<DenominatorTerm> denominator_terms;
  std::vector<ColumnNote> notes;
};

/*!
 * Observed widths and ANCs -> Brune (basis Brune) or Park (basis Park) amplitudes (LevelTransforms.cpp, TransformIn, first
 * loop).  columns label the value columns; each must have a block channel with the same partition, kind and l from
 * which the channel functions are taken (any channel spin: P_c, S_c depend on (alpha, l) only).
 */
BasisAmplitudes amplitudes_from_observed(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                                         const ProjectLevels &levels, LevelBasis basis, const ChannelFunctions &cf);

/// What the observed value of a (level, column) is (cnuc_split.md section 8, the table).
enum class ObservedKind : unsigned char {
  Width,             ///< W = Gamma (MeV): open particle channel (0 at threshold or for P = 0), photon channel
  ANC,               ///< W = C^2 (fm^-1): particle channel below threshold
  MagneticMoment,    ///< W = mu^2 (mu_N^2): static M1 moment (photon channel to the level's own state)
  QuadrupoleMoment,  ///< W = (Q 100 (alpha hbar c)^1/2)^2: static E2 moment
  Beta               ///< W = B: beta-decay feeding amplitude (Brune normalization)
};

/// The observed quantities of the levels of a block (observed_levels).  Row-major [lambda * n_columns + m].
struct ObservedLevels {
  std::vector<ObservedKind> kind;
  /// W_lambda,m = s 2 P |gamma_tot|^2 / N, s = sign Re gamma_tot (reversed for a static moment of half-integer J);
  /// for a beta column Re gamma_tot (/ J^1/2 in Park's basis).  Units by kind (above).
  std::vector<double> W;
  /// The value as a <levels> line gives it: Gamma (eV), C (fm^-1/2), mu (mu_N), Q (b), B, all signed; the Brune
  /// amplitude (MeV^1/2) of a column flagged is_amplitude.  amplitudes_from_observed reads it back.
  std::vector<double> value;
  /// Per level: N_lambda = 1 + sum_c gamma_c^2 dS_c/dE (Brune), or Park's J_lambda = 1 - sum_c gamma_c^2 dS_c/dE.
  std::vector<double> norm;
};

/*!
 * The observed quantities of the levels of block b from their total amplitudes gamma_int + gamma_ext (cnuc_split.md
 * section 8, Eq. 8.1): E are the Brune (Park) energies (MeV, E_x), gamma the internal Brune (Park) amplitudes
 * (MeV^1/2, row-major over columns), ext the external-capture amplitudes in the same layout (empty: none; entries not
 * present count 0).  The normalization N sums the internal particle amplitudes only.  The one definition of
 * Gamma, C, mu, Q and B from amplitudes; observed_from_amplitudes is its value with ext empty.
 */
ObservedLevels observed_levels(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                               const std::vector<double> &E, const std::vector<double> &gamma,
                               const std::vector<ExternalWidth> &ext, const std::vector<unsigned char> &is_amplitude,
                               LevelBasis basis, const ChannelFunctions &cf);

/// Brune or Park amplitudes -> the observed form (the inverse of amplitudes_from_observed; is_amplitude marks
/// the columns to report as Brune amplitudes): observed_levels(...).value without external capture.  Returns a
/// ProjectLevels in Form::Observed.
ProjectLevels observed_from_amplitudes(const System &system, std::size_t b, const std::vector<ProjectChannel> &columns,
                                       const std::vector<double> &E, const std::vector<double> &gamma,
                                       const std::vector<unsigned char> &is_amplitude, LevelBasis basis,
                                       const ChannelFunctions &cf);

/// B_c of the block's channels for the formal basis: S_c at the lowest of the level energies E (MeV, E_x; the first of
/// equal ones), so that the formal parameters do not depend on the order in which the levels are listed (HANDOFF
/// decision 10).  Non-particle channels take the value of block channel 0 (0 when channel 0 is not a particle
/// channel), as the engine does; they do not enter the level matrix.  No levels: 0.
std::vector<double> boundary_values(const System &system, std::size_t b, const std::vector<double> &E,
                                    const ChannelFunctions &cf);

/// Photon channels: total amplitude -> internal one, gamma_int = sign(gamma) (gamma^2 - Im^2)^1/2 - Re of the
/// external amplitude (LevelTransforms.cpp, TransformIn, second loop).  ext has the layout of gamma; entries not present are
/// left alone: the --ignore-externals mask is applied by the caller, which leaves a masked channel's entry absent
/// (external_capture.md section 3).
/// exceeded, if given, receives the index (lambda * n_c + k) of every entry where |Im| exceeded the total amplitude.
void subtract_external_capture(const System &system, std::size_t b, std::vector<double> &gamma,
                               const std::vector<ExternalWidth> &ext, std::vector<std::string> &warnings, std::vector<std::size_t> *exceeded = nullptr);

struct LevelSet {
  std::vector<double> E;      ///< MeV
  std::vector<double> gamma;  ///< row-major, block channels
  std::vector<int> iterations;  ///< brune_from_formal: iterations per level (0 for the others)
  std::vector<std::string> warnings;
  /// brune_from_formal only, so that a caller can word its own messages.  Per level: status (Converged, or the
  /// formal values kept because the eigenproblem failed or the iteration did not converge), the last energy step
  /// |e_k - E_k| (MeV) when it did not converge; per (level, block channel), row-major, the boundary-value shift
  /// S_c(E) - B_c of the last iteration (0 for non-particle channels).
  /// RoundingLimited (FormalToBruneOptions::evaluation only): the tolerance was not met in the iterations, but the
  /// last step is within rounding_bound_MeV, and the last iterate is returned.
  /// Bracketed (FormalToBruneOptions::bracket): the iteration did not converge (the fixed point repels, e' < -1), and
  /// the level is the root of E - e(E) found by bisection (level_parameters.md section 10); residual_MeV is |E - e(E)|
  /// there and boundary_shift S_c(E) - B_c.
  enum class Status : unsigned char { Converged, EigenFailed, NotConverged, RoundingLimited, Bracketed };
  std::vector<Status> status;
  std::vector<double> residual_MeV;
  std::vector<double> boundary_shift;
  /// brune_from_formal with FormalToBruneOptions::want_transform: T (n x n, row-major, input level order), the Brune
  /// amplitudes as combinations of the formal ones, gamma~_lambda c = sum_mu T_lambda mu gamma_mu c; row lambda is the
  /// unit eigenvector of lambda's last problem (a unit vector for a level whose formal values are kept).  Empty
  /// otherwise.  The masked channel capture of the formal basis uses it (external_capture.md section 3).
  std::vector<double> transform;
};

/*!
 * Brune -> formal with constant B (LevelTransforms.cpp, TransformIn, second loop), by generalized_symmetric_eigen.  Level lambda
 * keeps its position: it receives the eigenvalue of the rank of its Brune energy, and the orientation of its Brune
 * level: of the two signs of the eigenvector, the one whose component on the level's own Brune state is positive.  So
 * the result does not depend on the eigensolver, the formal amplitudes keep the input signs where levels mix weakly,
 * and an external-capture final state keeps the sign the photon amplitudes of other levels refer to
 * (level_parameters.md section 4).  Throws std::domain_error when M is not positive definite (no formal parameters
 * correspond) and std::invalid_argument for two levels at the same energy.
 */
LevelSet formal_from_brune(const System &system, std::size_t b, const std::vector<double> &E,
                           const std::vector<double> &gamma, const std::vector<double> &B, const ChannelFunctions &cf);

struct FormalToBruneOptions {
  double energy_tolerance_MeV = 1e-6;  ///< the engine's (CompoundSystem::TransformOut)
  /// Converged when the step is <= max(energy_tolerance_MeV, relative_tolerance |E|); 0: the absolute test alone.
  double relative_tolerance = 0.0;
  int max_iterations = 1000;           ///< the engine's
  bool want_transform = false;         ///< fill LevelSet::transform
  /// Inside an evaluation (level_parameters.md section 9): a level that does not meet the tolerance takes its last
  /// iterate if the last step is <= rounding_bound_MeV (status RoundingLimited), else NaN energy, amplitudes and row
  /// of T; never its formal values.  false: the formal values are kept, with a warning (the report iteration).
  bool evaluation = false;
  double rounding_bound_MeV = 1e-9;
  /// A level that would end NotConverged is solved as the root of E - e(E) by bisection (status Bracketed,
  /// level_parameters.md section 10); false: the fixed-point outcome alone.
  bool bracket = true;
};

/// Formal -> Brune by the engine's iteration (CompoundSystem::TransformOut): per level, B_c moved to S_c(E) and the matrix
/// E - sum_c (S_c - B_c) gamma gamma rediagonalized (symmetric_eigen) until the level's energy moves by <= the
/// tolerance.  A level that does not converge is the root of E - e(E), found by bisection (level_parameters.md
/// section 10); if no root is bracketed, it keeps its formal values, with a warning, as in the engine (with
/// FormalToBruneOptions::evaluation: its last iterate or NaN, level_parameters.md section 9).  Each
/// eigenvector is oriented as in formal_from_brune (its component on its own level positive), so formal -> Brune ->
/// formal and Brune -> formal -> Brune return the signs they started from.
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
  bool formal = true;  ///< compute the formal parameters (false: basis_* only, E/gamma left empty)
};

/*!
 * A project's parameters of block b, in the convention of the record and with the given columns, -> the internal
 * formal parameters.  Observed form: widths -> Brune (or Park) amplitudes in the project's columns; recoupling to
 * the block's channels; external capture subtracted; B_c (boundary_values); Park -> Brune; Brune -> formal
 * (formal_from_brune).  Amplitudes form in the formal basis: recoupling only.
 */
InternalParameters level_parameters(const System &system, std::size_t b, const Convention &convention,
                                    const std::vector<ProjectChannel> &columns, const ProjectLevels &levels,
                                    const ChannelFunctions &cf,
                                    const LevelParameterOptions &options = LevelParameterOptions());

}  // namespace physics
}  // namespace azure2

#endif
