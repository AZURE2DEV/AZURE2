#include "CNucLevelParameters.h"

#include <cmath>

#include "CNuc.h"
#include "CNucExternalCapture.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "Constants.h"
#include "WhitFunc.h"
#include "azure2/physics/Widths.h"

namespace phys = azure2::physics;

phys::ChannelFunctions ChannelFunctionsFromCNuc(CNuc &compound, const CNucSystemMap &map, const Config &configure) {
  return ChannelFunctionsFromCNuc(compound, map, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC),
                                  !!(configure.paramMask & Config::USE_RMC_FORMALISM));
}

phys::ChannelFunctions ChannelFunctionsFromCNuc(CNuc &compound, const CNucSystemMap &map, bool gsl, bool rmc) {
  CNuc *cn = &compound;
  const CNucSystemMap m = map;
  return [cn, gsl, rmc, m](std::size_t b, std::size_t k, double E) {
    const JGroup *group = cn->GetJGroup(m.groupOfBlock.at(b));
    const AChannel *channel = group->GetChannel(m.channelOfBlockChannel.at(b).at(k));
    const PPair *pair = cn->GetPair(channel->GetPairNum());
    const double Ec = E - pair->GetSepE() - pair->GetExE();  // (E_x - S) - e, as at a data point
    const int l = channel->GetL();
    const int twoJ = pair->SpinOrbitTwoJ(group->GetJ());
    phys::ChannelValues v;
    const char type = channel->GetRadType();
    if (type == 'P') {
      ChannelFunc cf(pair, gsl);
      v.S = cf.Shift(l, twoJ, Ec);
      if (Ec > 0.0) {
        v.P = cf.Penetrability(l, twoJ, Ec);
      } else if (Ec < 0.0) {
        // The ANC conversion factor mu a / (hbar^2 W_l(2 kappa a)^2) (level_parameters.md section 1).
        const double radius = pair->GetChRad();
        WhitFunc whit(pair);
        const double W = whit(l, twoJ, radius, std::fabs(Ec));
        v.P = pair->GetRedMass() * radius * uconv / std::pow(hbarc, 2.0) / std::pow(W, 2.0);
      }
      v.dS_dE = cf.ShiftDerivative(l, twoJ, Ec);
    } else if (type == 'E' || type == 'M') {
      if (std::fabs(E - pair->GetHeavyExE()) < 1.e-3 && group->GetJ() == pair->GetJ(2) && group->GetPi() == pair->GetPi(2)) {
        // A static moment (CNuc.cpp TransformIn, "Static moments" in conventions.rst).
        const double J = group->GetJ();
        v.static_moment = true;
        v.P = 1e-10;
        if (type == 'M' && l == 1) {
          v.P = 3.0 * J / 4.0 / (J + 1.) / nuclearMagneton / nuclearMagneton;
        } else if (type == 'E' && l == 2) {
          v.P = 60.0 * J * (2. * J - 1.) / (J + 1.) / (2. * J + 3.);
          v.input_scale = 100 * std::sqrt(fstruc * hbarc);
        }
      } else {
        v.P = rmc ? 1.0 : azure2::physics::photon_penetrability(Ec, l);  // physics/Widths.h
      }
    } else {
      v.P = 1.0;
    }
    return v;
  };
}

phys::ProjectLevels ProjectLevelsFromCNuc(CNuc &filled, const CNucSystemMap &map, std::size_t b) {
  const JGroup *group = filled.GetJGroup(map.groupOfBlock.at(b));
  const std::vector<int> &levels = map.levelOfBlockLevel.at(b);
  const std::vector<int> &channels = map.channelOfBlockChannel.at(b);
  phys::ProjectLevels out;
  out.form = phys::ProjectLevels::Form::Observed;
  for (int la : levels) {
    const ALevel *level = group->GetLevel(la);
    out.E.push_back(level->GetE());
    for (int ch : channels) {
      out.value.push_back(level->GetGamma(ch));
      out.is_amplitude.push_back(level->GammaIsRWA(ch) ? 1 : 0);
    }
  }
  return out;
}

std::vector<phys::ExternalWidth> ExternalWidthsFromCNuc(CNuc &compound, const CNucSystemMap &map, std::size_t b,
                                                        const Config &configure) {
  const JGroup *group = compound.GetJGroup(map.groupOfBlock.at(b));
  const std::vector<int> &levels = map.levelOfBlockLevel.at(b);
  const std::vector<int> &channels = map.channelOfBlockChannel.at(b);
  std::vector<phys::ExternalWidth> out(levels.size() * channels.size());
  if (!(configure.paramMask & Config::USE_EXTERNAL_CAPTURE)) return out;
  for (std::size_t la = 0; la < levels.size(); la++)
    for (std::size_t k = 0; k < channels.size(); k++) {
      const AChannel *channel = group->GetChannel(channels[k]);
      if (channel->GetRadType() != 'E' && channel->GetRadType() != 'M') continue;
      // --ignore-externals: a channel whose input width is zero has no external part (external_capture.md
      // section 3); gamma^ext is not evaluated.
      if ((configure.paramMask & Config::IGNORE_ZERO_WIDTHS) && group->GetLevel(levels[la])->ZeroInputWidth(channels[k]))
        continue;
      out[la * channels.size() + k].present = true;
      out[la * channels.size() + k].value =
          ExternalReducedWidth(compound, group, group->GetLevel(levels[la]), channel, true, configure);
    }
  return out;
}
