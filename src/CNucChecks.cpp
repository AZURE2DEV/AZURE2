// The checks of a compound nucleus and their messages (CNucChecks.h).  Reads the nucleus through const access only.
#include "CNucChecks.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "Config.h"
#include "Constants.h"
#include "JGroup.h"
#include "PPair.h"
#include "ParameterLabel.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/LevelParameters.h"

/*!
 * An identical pair with spin (p+p, d+d, 3He+3He) may not have a channel forbidden by Bose/Fermi symmetry, one with
 * L+S odd.  With spin the elastic cross section symmetrizes each pathway by 1 + (-1)^(L'+S'), which the cross-section
 * route applies as a factor 4 on every channel it is given: a forbidden channel would be counted four times instead of
 * not at all, and the result would be wrong without any sign of it.  Such channels were also what the previous fermion
 * rule, (-1)^(L+S) = -1, asked for (3S1 instead of 1S0 in p+p), so an old model is stopped here with the reason rather
 * than silently re-interpreted.  For spin-0 pairs such a channel only draws a warning (WarnForbiddenSpinlessChannels).
 */
int RefuseForbiddenIdenticalChannels(const CNuc &compound, const Config &configure) {
  int bad = 0;
  for (int j = 1; j <= compound.NumJGroups(); j++) {
    const JGroup *jg = compound.GetJGroup(j);
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      const AChannel *channel = jg->GetChannel(ch);
      if (channel->GetRadType() != 'P') continue;
      const PPair *pp = compound.GetPair(channel->GetPairNum());
      if (!pp->HasSpinDependentExchange()) continue;
      const int twoS = (int)std::lround(2.0 * channel->GetS());
      if ((channel->GetL() + twoS / 2) % 2 == 0) continue;
      configure.outStream << "ERROR: Identical-particle pair (Z=" << pp->GetZ(1)
                          << ", A=" << pp->GetM(1) << ", j=" << pp->GetJ(1)
                          << ") has a channel with J=" << jg->GetJ()
                          << ", L=" << channel->GetL() << ", S=" << channel->GetS()
                          << ", which Bose/Fermi symmetry forbids (L+S must be even)."
                          << std::endl;
      bad++;
    }
  }
  if (bad) {
    configure.outStream << "       Remove them from the model: allowed channels have L+S even (1S0, 3P0,1,2, 1D2, ... for p+p)."
                        << std::endl;
    return -1;
  }
  return 0;
}

/*!
 * Bose/Fermi symmetry requires L+S even (for bosons the same as (-1)^(L+S) == +1).  Such channels still contribute to
 * the calculation, but the result does not respect identical-particle symmetry.  Only spin-0 pairs get this far: for a
 * pair with spin CNuc::Fill has refused the model (RefuseForbiddenIdenticalChannels).
 */
void WarnForbiddenSpinlessChannels(const CNuc &compound, const Config &configure) {
  for (int j = 1; j <= compound.NumJGroups(); j++) {
    const JGroup *jg = compound.GetJGroup(j);
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      const AChannel *channel = jg->GetChannel(ch);
      if (channel->GetRadType() != 'P') continue;
      const PPair *pp = compound.GetPair(channel->GetPairNum());
      if (!pp->IsIdentical()) continue;
      int chL = channel->GetL();
      int chS_twice = (int)(2.0 * channel->GetS() + 0.5);
      // Half-integer S in an identical pair would only arise for half-integer-spin particles in an asymmetric
      // coupling; the even/odd test is skipped there to avoid false positives.
      if (chS_twice % 2 != 0) continue;
      if ((chL + chS_twice / 2) % 2 != 0) {
        configure.outStream << "**WARNING: Identical-particle pair (Z="
                            << pp->GetZ(1) << ", A=" << pp->GetM(1)
                            << ") has channel with L=" << chL
                            << ", S=" << channel->GetS()
                            << " violating Bose/Fermi symmetry "
                            << "(L+S must be even). This channel should be removed from the input."
                            << std::endl;
      }
    }
  }
}

/*!
 * Legendre orders above Config::maxLOrder are never computed, so a particle channel with l > maxLOrder/2 contributes
 * terms (l + l' > maxLOrder) that are silently dropped from every differential cross section (angle-integrated values
 * are complete).  Said once per such channel (checklist t3-24).
 */
void WarnTruncatedLegendreOrders(const CNuc &compound, const Config &configure) {
  for (int j = 1; j <= compound.NumJGroups(); j++) {
    const JGroup *jg = compound.GetJGroup(j);
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      const AChannel *channel = jg->GetChannel(ch);
      if (channel->GetRadType() != 'P' || 2 * channel->GetL() <= Config::maxLOrder) continue;
      configure.outStream << "**WARNING: channel l = " << channel->GetL() << " (J = " << jg->GetJ() << ", pair "
                          << compound.GetPair(channel->GetPairNum())->GetPairKey() << ") exceeds " << Config::maxLOrder / 2
                          << ": Legendre terms above order " << Config::maxLOrder
                          << " are not computed, so differential cross sections involving it are truncated." << std::endl;
    }
  }
}

/*!
 * Separation energies against masses, per compound nucleus (DeBoer 2026-10-29): S_b - S_a should equal
 * (m1b + m2b - m1a - m2a) c^2.  Thresholds come from the separation energies and the kinematics from the masses, so a
 * mismatch is never an error in the calculation, but it is how typos and mixed mass conventions show up.  Every
 * particle pair (pType 0) of a nucleus is compared with the first one.  Gamma and beta pairs are left out -- their
 * heavy mass is only a label (the photon energies come from excitation energies), and many projects enter it as a
 * mass number.  A mismatch within 5 keV of a multiple of the electron mass gets a hint: atomic and nuclear (or
 * bare-proton) masses mixed.  Nothing in the calculation changes.
 */
int WarnMassSeparationMismatch(const CNuc &compound, const Config &configure) {
  const double me = azure2::physics::kElectronMass_MeV;
  const int kMaxListed = 3;
  // A pair whose masses are mass numbers (both integers, other than 12C's exact 12) carries no Q value: it is
  // reported once, not compared.
  auto massNumbers = [](const PPair *p) {
    auto integral = [](double m) { return std::fabs(m - std::round(m)) < 1e-6; };
    return integral(p->GetM(1)) && integral(p->GetM(2)) && !(std::round(p->GetM(1)) == 12 && std::round(p->GetM(2)) == 12);
  };
  int warnings = 0, unlisted = 0;
  std::map<int, int> reference;  // nucleus -> pair number of its first particle pair with real masses
  std::set<int> massNumberPairs;
  for (int i = 1; i <= compound.NumPairs(); i++) {
    const PPair *pb = compound.GetPair(i);
    if (pb->GetPType() != 0) continue;
    if (massNumbers(pb)) {
      massNumberPairs.insert(pb->GetPairKey());
      continue;
    }
    auto it = reference.find(pb->GetNucleus());
    if (it == reference.end()) {
      reference[pb->GetNucleus()] = i;
      continue;
    }
    const PPair *pa = compound.GetPair(it->second);
    const double fromS = pb->GetSepE() - pa->GetSepE();
    const double fromM = (pb->GetM(1) + pb->GetM(2) - pa->GetM(1) - pa->GetM(2)) * uconv;
    const double diff = fromS - fromM;
    if (!(std::fabs(diff) > kMassCheckTolerance)) continue;
    if (++warnings > kMaxListed) {
      ++unlisted;
      continue;
    }
    std::ostringstream msg;
    msg << std::setprecision(6) << "**WARNING: separation energies and masses disagree";
    if (compound.NumNuclei() > 1) msg << " in " << compound.NucleusLabel(pb->GetNucleus());
    msg << ": pair " << pb->GetPairKey() << " (S = " << pb->GetSepE() << " MeV, m1 + m2 = " << std::setprecision(10)
        << pb->GetM(1) + pb->GetM(2) << " u) and pair " << pa->GetPairKey() << " (S = " << std::setprecision(6)
        << pa->GetSepE() << " MeV, m1 + m2 = " << std::setprecision(10) << pa->GetM(1) + pa->GetM(2)
        << " u): the separation energies differ by " << std::fixed << std::setprecision(4) << fromS
        << " MeV, the masses by " << fromM << " MeV (" << std::showpos << std::setprecision(1) << diff * 1e3
        << std::noshowpos << " keV).";
    const double ne = diff / me;
    if (std::fabs(ne) > 0.5 && std::fabs(ne - std::round(ne)) * me < 5.0e-3)
      msg << "  That is " << std::lround(std::fabs(ne)) << " electron mass" << (std::lround(std::fabs(ne)) > 1 ? "es" : "")
          << ": atomic and nuclear masses mixed (a bare-proton 1.00728 u among atomic masses, for one).";
    msg << "  Thresholds use the separation energies and the kinematics the masses; take both from one mass table "
           "(S = m1 + m2 - M of the compound nucleus).";
    configure.outStream << msg.str() << std::endl;
  }
  if (unlisted > 0)
    configure.outStream << "**WARNING: " << unlisted << " more pair" << (unlisted > 1 ? "s" : "")
                        << " disagree in the same way (separation energies against masses)." << std::endl;
  if (!massNumberPairs.empty()) {
    std::ostringstream msg;
    msg << "**NOTE: pair" << (massNumberPairs.size() > 1 ? "s" : "");
    for (int k : massNumberPairs) msg << " " << k;
    msg << " use" << (massNumberPairs.size() > 1 ? "" : "s") << " mass numbers as masses: the kinematics are good to "
        << "about 0.1 % and the separation energies cannot be checked against them.";
    configure.outStream << msg.str() << std::endl;
  }
  return warnings;
}

/*!
 * The R-matrix description of radiative capture treats the photon channels perturbatively (the gamma widths are
 * assumed to make a negligible contribution to the total width entering the level matrix), so a level with
 * Gamma_gamma >~ 0.1 Gamma_particle is outside the range in which the calculated cross section can be trusted.  The
 * widths are the observed ones the levels hold (parameters.out's: open particle channels and photon channels, with the
 * level's normalization and the external part of the photon amplitude).  Bound channels and static moments carry no
 * width and are skipped; levels with no open particle channel are not checked.  One warning lists every flagged
 * level, worst first.
 */
void WarnLargeRadiativeWidths(const CNuc &compound, const ObservedLevelTable &table, const Config &configure) {
  namespace phys = azure2::physics;
  const double warningRatio = 0.1;
  struct FlaggedLevel {
    double ratio;
    double radWidth;
    double particleWidth;
    std::string label;
  };
  std::vector<FlaggedLevel> flagged;

  for (int j = 1; j <= compound.NumJGroups(); j++) {
    const JGroup *theJGroup = compound.GetJGroup(j);
    if (!theJGroup->IsInRMatrix()) continue;
    for (int la = 1; la <= theJGroup->NumLevels(); la++) {
      const ALevel *theLevel = theJGroup->GetLevel(la);
      if (!theLevel->IsInRMatrix()) continue;
      const ObservedLevel &observed = table.Level(j, la);
      double particleWidth = 0.0;
      double radiativeWidth = 0.0;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        if (static_cast<phys::ObservedKind>(observed.kind[ch - 1]) != phys::ObservedKind::Width) continue;
        const char type = theJGroup->GetChannel(ch)->GetRadType();
        if (type == 'P')
          particleWidth += fabs(observed.W[ch - 1]);
        else if (type == 'M' || type == 'E')
          radiativeWidth += fabs(observed.W[ch - 1]);
      }
      if (particleWidth <= 0.0 || radiativeWidth <= 0.0) continue;

      double ratio = radiativeWidth / particleWidth;
      if (ratio > warningRatio) {
        FlaggedLevel entry;
        entry.ratio = ratio;
        entry.radWidth = radiativeWidth;
        entry.particleWidth = particleWidth;
        entry.label = AZURELabel::Level(theJGroup, compound.Inputs().E(j, la), j, la);
        flagged.push_back(entry);
      }
    }
  }
  if (flagged.empty()) return;

  std::sort(flagged.begin(), flagged.end(), [](const FlaggedLevel &a, const FlaggedLevel &b) { return a.ratio > b.ratio; });
  const int numFlagged = (int)flagged.size();
  configure.outStream << std::endl
                      << "**WARNING: " << numFlagged << " level"
                      << ((numFlagged == 1) ? " has" : "s have")
                      << " a radiative width larger than " << warningRatio * 100.
                      << "% of the particle width." << std::endl
                      << "  R-Matrix capture assumes G_gamma << G_particle, so the calculated"
                      << " cross section may be incorrect for "
                      << ((numFlagged == 1) ? "this level." : "these levels.") << std::endl;
  // Cap the list so a badly-configured model cannot bury the rest of the log.
  const int maxListed = 10;
  const int listed = std::min(numFlagged, maxListed);
  for (int i = 0; i < listed; i++) {
    configure.outStream << "    " << flagged[i].label << std::endl
                        << "      G_gamma/G_particle = "
                        << std::scientific << std::setprecision(2) << flagged[i].ratio
                        << " (G_gamma = " << flagged[i].radWidth * 1e6
                        << " eV, G_particle = " << flagged[i].particleWidth * 1e6 << " eV)"
                        << std::endl;
    configure.outStream.unsetf(std::ios::scientific);
    configure.outStream.precision(6);
  }
  if (numFlagged > listed)
    configure.outStream << "    ...and " << (numFlagged - listed)
                        << " more level" << (((numFlagged - listed) == 1) ? "" : "s")
                        << " with a smaller ratio." << std::endl;
  configure.outStream << std::endl;
}
