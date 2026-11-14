// The messages of the level-parameter transformations (CNucTransformMessages.h).  Reads the nucleus through const
// access only.
#include "CNucTransformMessages.h"

#include <cmath>
#include <iostream>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "Config.h"
#include "JGroup.h"
#include "PPair.h"
#include "ParameterLabel.h"

namespace phys = azure2::physics;

namespace {
// AZURELabel takes non-const pointers (include/ParameterLabel.h); it only reads them.
std::string LevelLabel(const CNuc &compound, int j, int la) {
  const JGroup *group = const_cast<CNuc &>(compound).GetJGroup(j);
  return AZURELabel::Level(group, group->GetLevel(la), j, la);
}
std::string ChannelLabel(const CNuc &compound, int j, int ch) {
  CNuc &c = const_cast<CNuc &>(compound);
  return AZURELabel::Channel(&c, c.GetJGroup(j), ch);
}
}  // namespace

void WarnInputColumn(const CNuc &compound, const Config &configure, int j, int la, int ch, phys::ColumnNote note,
                     double input) {
  const JGroup *group = compound.GetJGroup(j);
  const ALevel *level = group->GetLevel(la);
  switch (note) {
    case phys::ColumnNote::None: break;
    case phys::ColumnNote::AmplitudeFlagIgnored:
      configure.outStream << "WARNING: Reduced width amplitude input flag is only "
                          << "supported for particle channels; ignored for J = "
                          << group->GetJ() << " E = " << level->GetE()
                          << " MeV channel " << ch << "." << std::endl;
      break;
    case phys::ColumnNote::ZeroPenetrability: {
      const PPair *pair = compound.GetPair(group->GetChannel(ch)->GetPairNum());
      configure.outStream << "**WARNING: " << LevelLabel(compound, j, la)
                          << ", " << ChannelLabel(compound, j, ch)
                          << ": the level is " << level->GetE() - pair->GetSepE() - pair->GetExE()
                          << " MeV above the channel threshold, where P = 0 in double precision, "
                          << "so its observed width of " << input
                          << " eV cannot be converted; the reduced width is set to 0. "
                          << "Give a reduced-width amplitude (gammaIsRWA) for a level at threshold."
                          << std::endl;
      break;
    }
    case phys::ColumnNote::AtThreshold:
      configure.outStream << "**WARNING: " << LevelLabel(compound, j, la)
                          << ", " << ChannelLabel(compound, j, ch)
                          << ": the level sits exactly at the channel threshold, where an ANC "
                          << "does not define a reduced width; the reduced width is set to 0. "
                          << "Give a reduced-width amplitude (gammaIsRWA) for a level at threshold."
                          << std::endl;
      break;
  }
}

void WarnObservedOutOfReach(const CNuc &compound, const Config &configure, int j, int la,
                            const phys::BasisAmplitudes &amp, const phys::ProjectLevels &in, std::size_t lb,
                            std::size_t nc) {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  const JGroup *group = compound.GetJGroup(j);
  const ALevel *level = group->GetLevel(la);
  const double parkNorm = amp.park_norm[lb];
  const double halfDenom = amp.half_denominator[lb];
  if (park && parkNorm <= 0.) {
    // In Brune's parametrization this limit is only reached for an infinite reduced width; in Park's it is an
    // ordinary point of the parameter space, with nothing in the level matrix to stop it.
    configure.outStream << "**WARNING: Park norm J = 1 - sum gamma^2 dS/dE = " << parkNorm
                        << " is not positive" << std::endl
                        << "    " << LevelLabel(compound, j, la) << std::endl
                        << "  The widths of this level exceed what the channel radii allow; no"
                        << std::endl
                        << "  standard R-matrix (real reduced widths) corresponds to it."
                        << std::endl;
  } else if (!park && halfDenom < 0.) {
    // An observed width Gamma_c = 2 P_c gamma_c^2 / (1 + sum gamma^2 dS/dE) saturates at 2 P_c / (dS_c/dE) as
    // gamma_c -> infinity, so the input is out of reach of any reduced width when
    // 1 - sum_c Gamma_c (dS_c/dE) / (2 P_c) < 0.  The usual cause is input in the wrong units (a reduced-width
    // amplitude read as a width), so name the channel and the numbers.
    configure.outStream << "**WARNING: Denominator less than zero while transforming"
                        << std::endl
                        << "    " << LevelLabel(compound, j, la) << std::endl
                        << "  No reduced width reproduces these observed widths: "
                        << "1 - sum_c Gamma_c (dS_c/dE)/(2 P_c) = " << halfDenom << " < 0."
                        << std::endl;
    int worst = 0;  // channel number of the largest term, first of equal ones
    for (std::size_t k = 0; k < nc; k++) {
      const phys::DenominatorTerm &t = amp.denominator_terms[lb * nc + k];
      if (t.present && (worst == 0 || t.term > amp.denominator_terms[lb * nc + worst - 1].term)) worst = (int)k + 1;
    }
    if (worst > 0) {
      const phys::DenominatorTerm &t = amp.denominator_terms[lb * nc + worst - 1];
      const PPair *pair = compound.GetPair(group->GetChannel(worst)->GetPairNum());
      bool bound = level->GetE() - pair->GetSepE() - pair->GetExE() <= 0.0;
      configure.outStream << "  Largest term: " << ChannelLabel(compound, j, worst)
                          << (bound ? " (bound; ANC input)" : "")
                          << ": input " << in.value[lb * nc + worst - 1]
                          << (bound ? " fm^(-1/2)" : " eV")
                          << ", P = " << t.P << ", dS/dE = " << t.dS_dE
                          << " MeV^-1, Gamma (dS/dE)/(2P) = " << t.term / 2.0 << "." << std::endl;
      if (!bound && t.dS_dE > 0.0)
        configure.outStream << "  Even alone, this channel's observed width saturates at 2P/(dS/dE) = "
                            << 2.0 * t.P / t.dS_dE * 1e6 << " eV." << std::endl;
    }
    // Orders of magnitude over the limit is the signature of a unit mistake; a sum just past 1 is a level broader
    // than its channels allow at these radii.
    if (-halfDenom > 100.0)
      configure.outStream << "  So large an excess usually means wrong units: widths in eV (ANCs in fm^-1/2), "
                          << "not reduced-width amplitudes -- pyazr's calculate() and calculate_energies() "
                          << "take physical values, calculate_rwa() amplitudes." << std::endl;
    else
      configure.outStream << "  The level is broader than its channels allow at these radii: "
                          << "reduce the widths or enlarge the channel radii." << std::endl;
    configure.outStream << "  The transformation may not have been successful for this level."
                        << std::endl;
  }
}

void WarnExternalExceedsTotal(const CNuc &compound, const Config &configure, int j, int la, int ch) {
  CNuc &c = const_cast<CNuc &>(compound);  // AZURELabel only reads it
  const JGroup *group = c.GetJGroup(j);
  configure.outStream << "**WARNING: Imaginary portion of the external width is greater "
                      << "than the total width" << std::endl
                      << "    "
                      << AZURELabel::LevelAndChannel(&c, group, group->GetLevel(la), j, la, ch)
                      << std::endl
                      << "    [j=" << j << " la=" << la << " ch=" << ch << "]"
                      << std::endl;
}

void WarnParkNormNotPositive(const CNuc &compound, const Config &configure, int j, int la) {
  configure.outStream << "**WARNING: Park norm J = " << compound.State().ParkNorm(j, la)
                      << " is not positive" << std::endl
                      << "    " << LevelLabel(compound, j, la) << std::endl
                      << "  The widths of this level exceed what the channel radii allow; no"
                      << std::endl
                      << "  standard R-matrix (real reduced widths) corresponds to it."
                      << std::endl;
}

void WarnFormalToBrune(const CNuc &compound, const Config &configure, int j, int la, const phys::LevelSet &out,
                       std::size_t lb, std::size_t nc, const phys::FormalToBruneOptions &options) {
  if (out.status[lb] == phys::LevelSet::Status::EigenFailed) {
    configure.outStream << "**WARNING: Could not transform level (eigenvalue problem failed)" << std::endl
                        << "    " << LevelLabel(compound, j, la) << std::endl
                        << "  The input energy and widths are kept for this level; its "
                        << "transformed values are unreliable." << std::endl;
  } else if (out.status[lb] == phys::LevelSet::Status::Bracketed) {
    // The fixed point repels (e' < -1): the root of E - e(E) by bisection (level_parameters.md section 10).
    configure.outStream << "**NOTE: Level transformed by a bracketed root (the iteration does not converge here)"
                        << std::endl
                        << "    " << LevelLabel(compound, j, la) << std::endl
                        << "  Energy residual " << out.residual_MeV[lb] << " MeV (" << out.iterations[lb]
                        << " iterations and evaluations of E - e(E))." << std::endl;
  } else if (out.status[lb] == phys::LevelSet::Status::NotConverged) {
    // The iteration chases the energy eigenvalue while each channel's boundary condition moves with it, so the
    // channel whose boundary condition is still shifting most is the one holding up convergence -- name it, since
    // that is the width to look at.
    int worstChannel = 0;
    double worstDiff = 0.0;
    for (std::size_t k = 0; k < nc; k++)
      if (std::fabs(out.boundary_shift[lb * nc + k]) >= worstDiff) {
        worstDiff = std::fabs(out.boundary_shift[lb * nc + k]);
        worstChannel = (int)k + 1;
      }
    configure.outStream << "**WARNING: Could not transform level after "
                        << options.max_iterations << " iterations" << std::endl
                        << "    " << LevelLabel(compound, j, la) << std::endl
                        << "  Energy residual " << out.residual_MeV[lb]
                        << " MeV, tolerance " << options.energy_tolerance_MeV << " MeV." << std::endl;
    if (worstChannel > 0)
      configure.outStream << "  Least converged channel (largest boundary-condition shift, "
                          << worstDiff << "):" << std::endl
                          << "    " << ChannelLabel(compound, j, worstChannel) << std::endl;
    configure.outStream << "  The input energy and widths are kept for this level; its "
                        << "transformed values are unreliable." << std::endl;
  }
}
