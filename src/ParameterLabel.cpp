#include "ParameterLabel.h"

#include <cmath>
#include <sstream>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "JGroup.h"
#include "PPair.h"
#include "azure2/physics/Elements.h"

namespace AZURELabel {

std::string Spin(double j) {
  std::ostringstream out;
  const double twice = 2.0 * j;
  const long twiceRounded = std::lround(twice);

  if (std::fabs(twice - double(twiceRounded)) > 1.e-6) {
    // Not a multiple of 1/2; print it as-is rather than inventing a fraction.
    out << j;
    return out.str();
  }

  if (twiceRounded % 2 == 0) {
    out << (twiceRounded / 2);
  } else {
    out << twiceRounded << "/2";
  }
  return out.str();
}

std::string JPi(double j, int parity) {
  return Spin(j) + (parity < 0 ? "-" : "+");
}

std::string Nuclide(int z, double mass) {
  const int a = (int)std::lround(mass);

  // Conventional shorthands for the light projectiles. The alpha keeps its
  // "4He" form: a bare "a" is easy to misread in a pair like "n+a".
  if (z == 0 && a == 1) return "n";
  if (z == 1 && a == 1) return "p";
  if (z == 1 && a == 2) return "d";
  if (z == 1 && a == 3) return "t";

  std::ostringstream out;
  out << a;
  if (z >= 0 && z <= azure2::physics::kMaxElementZ) {
    out << azure2::physics::ElementSymbol(z);
  } else {
    out << "(Z=" << z << ")";
  }
  return out.str();
}

std::string Pair(CNuc *compound, int pairNum) {
  if (!compound || pairNum < 1 || pairNum > compound->NumPairs()) {
    std::ostringstream out;
    out << "pair " << pairNum;
    return out.str();
  }

  const PPair *pair = compound->GetPair(pairNum);

  // PType 10 marks a photon pair; naming its two "particles" is meaningless.
  if (pair->GetPType() == 10) return "gamma";

  return Nuclide(pair->GetZ(1), pair->GetM(1)) + "+" +
      Nuclide(pair->GetZ(2), pair->GetM(2));
}

std::string Level(const JGroup *jgroup, double energy, int jgroupIndex, int levelIndex) {
  std::ostringstream out;
  if (!jgroup) {
    out << "level " << levelIndex << " of J-group " << jgroupIndex;
    return out.str();
  }

  out.precision(4);
  out << "J^pi=" << JPi(jgroup->GetJ(), jgroup->GetPi()) << ", E=" << std::fixed
      << energy << " MeV (J-group " << jgroupIndex << ", level "
      << levelIndex << ")";
  // A project with several compound nuclei says which one; a single-nucleus
  // label is unchanged.
  if (jgroup->NucleusLabelled())
    out << " [Z=" << jgroup->GetNucleusZ() << " A=" << jgroup->GetNucleusA() << "]";
  return out.str();
}

std::string Channel(CNuc *compound, const JGroup *jgroup, int channelIndex) {
  std::ostringstream out;
  out << "channel " << channelIndex;

  if (!jgroup || channelIndex < 1 || channelIndex > jgroup->NumChannels()) {
    return out.str();
  }

  const AChannel *channel = jgroup->GetChannel(channelIndex);
  const char radType = channel->GetRadType();

  out << ": ";
  if (radType == 'E' || radType == 'M') {
    // Radiative channel: the multipolarity is the informative part.
    out << radType << channel->GetL() << " gamma to "
        << Pair(compound, channel->GetPairNum());
  } else {
    out << Pair(compound, channel->GetPairNum()) << ", L=" << channel->GetL()
        << ", S=" << Spin(channel->GetS());
    if (radType != 'P') out << " (type " << radType << ")";
  }
  return out.str();
}

std::string LevelAndChannel(CNuc *compound, const JGroup *jgroup, int jgroupIndex, int levelIndex, int channelIndex) {
  return Level(jgroup, compound->Inputs().E(jgroupIndex, levelIndex), jgroupIndex, levelIndex) + ", " +
      Channel(compound, jgroup, channelIndex);
}

int NumParameters(CNuc *compound, WorkingCopy *data) {
  int n = 0;
  if (compound) n += compound->ParameterLayout().Size();
  if (data) {
    n += data->Data().Abundances().NumParams();  // isotopic abundances, right after the level block
    const std::vector<Segment> &segments = data->Data().GetSegments();
    for (size_t s = 0; s < segments.size(); s++)
      if (segments[s].IsVaryNorm()) n++;
    n += 2 * (int)segments.size();  // energy shift and sqrt(E) coefficient, every segment
    for (size_t s = 0; s < segments.size(); s++)  // shift-function parameters of the segments declaring one
      if (segments[s].OwnsShiftFn() && segments[s].HasShiftFn()) n += segments[s].NumShiftFnParams();
    n += data->Data().NumScreeningParams();  // electron-screening potentials
    n += data->Data().NumDensityParams();    // fitted target densities
    n += data->Data().NumDopplerParams();    // fitted Doppler temperatures
    n += data->Data().NumThicknessParams();  // fitted sample thicknesses
    n += data->Data().NumApparatusParams();   // apparatus parameters (stopping scales, kernel widths)
  }
  return n;
}

std::string Parameter(CNuc *compound, WorkingCopy *data, int minuitIndex) {
  if (minuitIndex < 0) return std::string();

  // The level block first, as the compound's layout places it.
  int index = 0;
  if (compound) {
    const LevelParameterLayout &layout = compound->ParameterLayout();
    if (minuitIndex < layout.Size()) {
      const LevelParameterLayout::Slot &slot = layout.Slots()[minuitIndex];
      const JGroup *jgroup = compound->GetJGroup(slot.group);
      if (slot.kind == LevelParameterLayout::Slot::Kind::kEnergy)
        return "level energy: " + Level(jgroup, compound->Inputs().E(slot.group, slot.level), slot.group, slot.level);
      return "width: " + LevelAndChannel(compound, jgroup, slot.group, slot.level, slot.channel);
    }
    index = layout.Size();
    // then the isotopic abundances, one per (sample, nucleus) (AbundanceSamples::AddParameters)
    const int numAbundances = data ? data->Data().Abundances().NumParams() : 0;
    for (int k = 1; k <= numAbundances; k++) {
      if (index == minuitIndex) {
        int sample, nucleus;
        data->Data().Abundances().Of(k, sample, nucleus);
        const std::string &name = data->Data().Abundances().Sample(sample).name;
        return "abundance of nucleus " + compound->NucleusLabel(nucleus) + (name.empty() ? "" : " in sample " + name);
      }
      index++;
    }
  }

  if (!data) return std::string();

  // Then WorkingCopy::FillMnParams: norms for segments that vary one (a free THM
  // norm is profiled, not a parameter), then an energy shift and a sqrt(E)
  // coefficient for every segment, then the THM coherent backgrounds.
  const std::vector<Segment> &segments = data->Data().GetSegments();
  for (size_t s = 0; s < segments.size(); s++) {
    if (!segments[s].IsVaryNorm() || segments[s].IsProfiledNorm()) continue;
    if (index == minuitIndex) {
      std::ostringstream out;
      out << "normalization of segment " << segments[s].GetSegmentKey()
          << " (" << segments[s].GetDataFile() << ")";
      return out.str();
    }
    index++;
  }
  for (size_t s = 0; s < segments.size(); s++) {
    if (index == minuitIndex) {
      std::ostringstream out;
      out << "energy shift of segment " << segments[s].GetSegmentKey()
          << " (" << segments[s].GetDataFile() << ")";
      return out.str();
    }
    index++;
  }
  for (size_t s = 0; s < segments.size(); s++) {
    if (index == minuitIndex) {
      std::ostringstream out;
      out << "sqrt(E) energy-shift coefficient of segment " << segments[s].GetSegmentKey()
          << " (" << segments[s].GetDataFile() << ")";
      return out.str();
    }
    index++;
  }
  for (size_t s = 0; s < segments.size(); s++) {
    if (!segments[s].OwnsShiftFn() || !segments[s].HasShiftFn()) continue;
    for (int j = 0; j < segments[s].NumShiftFnParams(); j++) {
      if (index == minuitIndex) {
        std::ostringstream out;
        out << "parameter a" << j << " of the energy-shift function " << segments[s].GetShiftFnExpression()
            << " of segment " << segments[s].GetSegmentKey() << " (" << segments[s].GetDataFile() << ")";
        return out.str();
      }
      index++;
    }
  }
  for (int k = 1; k <= data->Data().NumScreeningParams(); k++) {
    if (index == minuitIndex) {
      std::ostringstream out;
      out << "electron-screening potential U_e (eV) of experimental effect " << data->Data().GetScreeningEffect(k);
      return out.str();
    }
    index++;
  }
  for (int k = 1; k <= data->Data().NumDensityParams(); k++) {
    if (index == minuitIndex) {
      std::ostringstream out;
      const TargetEffect *te = data->Data().GetDensityTargetEffect(k);
      out << "target density (atoms/cm^2) of experimental effect " << data->Data().GetDensityEffect(k);
      if (!te->GetDensityGroup().empty()) out << " (density group " << te->GetDensityGroup() << ")";
      return out.str();
    }
    index++;
  }
  for (int k = 1; k <= data->Data().NumDopplerParams(); k++) {
    if (index == minuitIndex) {
      std::ostringstream out;
      out << "free-gas Doppler temperature (K) of experimental effect " << data->Data().GetDopplerEffect(k);
      return out.str();
    }
    index++;
  }
  for (int k = 1; k <= data->Data().NumThicknessParams(); k++) {
    if (index == minuitIndex) {
      std::ostringstream out;
      const Segment *seg = data->Data().GetSegment(data->Data().GetThicknessSegment(k));
      out << "sample thickness (atoms/barn) of the " << (seg->IsTransmission() ? "transmission" : "yield")
          << " of segment " << seg->GetSegmentKey();
      if (!seg->GetThicknessGroup().empty()) out << " (thickness group " << seg->GetThicknessGroup() << ")";
      return out.str();
    }
    index++;
  }
  for (int k = 0; k < data->Data().NumThmCoherentParams(); k++) {
    if (index == minuitIndex)
      return "THM coherent background " + data->Data().GetThmCoherentParam(k).name;
    index++;
  }
  for (int k = 1; k <= data->Data().NumApparatusParams(); k++) {
    if (index == minuitIndex) return ApparatusLabel(data->Data(), k);
    index++;
  }

  return std::string();
}

}  // namespace AZURELabel
