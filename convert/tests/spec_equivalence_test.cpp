// Spec equivalence (docs/upgrade/project_format3.md step 3b, G2): for every old project of the repository that the
// converter converts, the engine's set-up specs made from the old lines (the temporary adaptor, SpecsFromOldLines.h)
// equal the specs made from the conversion (ConvertToFormat3 with exact grids -> ReadRecords -> Validate -> Build),
// member by member, doubles by their bits.
//
// The rules of the comparison (each is a meaning the old engine gave a field, or a difference that section 8 records):
// - Data sets and extrapolations by key (file order, inactive ones counted).  Of an inactive one only the key, the
//   flag and the data file: the old engine read nothing else of it.
// - Shared values (normalization, energy shift, sample thickness; WorkingCopy::Resolve*Groups): the groups as partitions of
//   the active keys (a group of one is no group; labels differ: the converter names the records), the values on the
//   first active member of each group and on ungrouped data sets only.  A follower's own old values are the recorded
//   difference of meaning of section 8.1 (the old engine never used them; the converter drops them with a note).
// - A prior centre: the effective centre (the prior_centre row, else the norm; the shift's, else 0).  A THM data set
//   with a profiled norm: the norm value is not compared (profiled, Segment::IsProfiledNorm).
// - The angle window and grid of an angle-integrated observable are not compared (Segment::IsInSegment tests the
//   window of differential data only; an old test line on several angles repeated its energy grid, written once).
// - A target effect: per active key, the effect the engine applies to it (a later <targetInt> line takes over a key),
//   its members with their defaults made explicit (straggling coefficient 0 = Bohr, width multiplier 20, points per
//   width 50); and as partitions: the keys sharing one effect, and the keys sharing one fitted density.  Old lines
//   that apply no effect to an active key are not effects.
// - Abundances: per sample (old name mapped to the target the data sets name) and compound nucleus.
// - THM: the experiment of each active key, and the weight files.
// - Parameter names (B, the director's decision on G2 item 4): the converter's map (ConvertResult::parameterNames)
//   sends each old name to one new name and no two old names to the same one (injective, every file); and the name
//   the old-line adaptor gives each parameter the engine makes for an active data set (norm when it varies, shift,
//   sqrt shift, the function terms it owns, a fitted thickness) maps to the name Build gives it.
// - Model data, in the engine's row order (src/ModelData.cpp, ModelData::Read): the old <scatteringLengths> lines that
//   give a residual (active, an uncertainty > 0) as the engine's constraint (pair key, s, Re a, sigma, Im a, sigma),
//   then the datums of an old <modelData> section (BuildModelData against the tables of the old <levels> rows); the
//   format-3 side: the datums Build made, a scattering length as the same constraint.  Doubles by their bits.
// - Hybrid potentials: for every particle pair of <levels>, the potential the old <potential> block gives it (the
//   engine's parser, convert/PotentialBlock.h: its own section, else the default) and the one its `potential` record
//   gives, as physics::PotentialKey (shape, V0, R, a, r0, R_C, r_m, step), doubles by their bits.  Photon and beta
//   pairs are not compared: the converter writes no record for them (no channel function of a potential is computed).
// AZURE2_SPEC_INJECT=1 moves one value of the format-3 side by one ulp, to show that the test fails on a difference
// (=2: the value of its first model datum; =3: the name of its first data set's energy-shift parameter; =4: the
// Coulomb radius of its first potential).
// The old projects are those of tests/ and examples/, an original kept in tests/convert/fixtures/ read in place of the
// committed conversion of the same path (docs/upgrade/format3_tree.md).
//   spec_equivalence_test <repository root>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/convert/Convert.h"
#include "azure2/convert/PotentialBlock.h"
#include "azure2/convert/SpecsFromOldLines.h"
#include "azure2/convert/ScatteringLengthLine.h"
#include "azure2/convert/ThmLines.h"
#include "azure2/experiment/ModelDatum.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/Convert.h"
#include "azure2/io/Document.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/RecordReader.h"

using namespace azure2;
using namespace azure2::experiment;

namespace {

std::vector<std::string> problems;  // of the current file

// Files whose specs differ, with the reason (each a construct the old engine refused or skipped, which a format-3
// record cannot express; or a converter defect left for the director, named so).
const std::map<std::string, const char *> kRecorded = {
    // none since natural_target's bad_key and bad_sample: the converter refuses both, as the old engine did
};
int injected = 0;                   // AZURE2_SPEC_INJECT: change one value of the format-3 side (shows the test fails)

template <class T>
void Eq(const std::string &what, const T &a, const T &b) {
  if (!(a == b)) problems.push_back(what);
}
void Bits(const std::string &what, double a, double b) {
  if (std::memcmp(&a, &b, sizeof a) != 0) {
    char buf[96];
    std::snprintf(buf, sizeof buf, " (%.17g vs %.17g)", a, b);
    problems.push_back(what + buf);
  }
}
void Bits(const std::string &what, const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size()) {
    problems.push_back(what + ": " + std::to_string(a.size()) + " vs " + std::to_string(b.size()) + " values");
    return;
  }
  for (std::size_t i = 0; i < a.size(); i++) Bits(what + "[" + std::to_string(i) + "]", a[i], b[i]);
}

bool Slurp(const std::string &path, std::string &out) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  std::ostringstream s;
  s << in.rdbuf();
  out = s.str();
  return true;
}

void Find(const std::string &dir, const std::string &rel, std::vector<std::string> &out) {
  namespace fs = std::filesystem;
  std::error_code ec;
  for (const fs::directory_entry &e : fs::directory_iterator(dir, ec)) {
    const std::string name = e.path().filename().string();
    if (name == "output" || name == "checks" || rel + name == "tests/convert") continue;
    if (e.is_directory(ec)) Find(e.path().string(), rel + name + "/", out);
    else if (name.size() > 4 && name.substr(name.size() - 4) == ".azr" && name.find(".f3.") == std::string::npos)
      out.push_back(rel + name);
  }
}

std::unique_ptr<std::istream> Section(const io::Document &doc, const std::string &name) {
  const io::Document::Section *s = doc.Find(name);
  if (!s || s->close == io::Document::Section::npos) return nullptr;
  std::string text;
  for (std::size_t i = s->open; i <= s->close; i++) text += doc.Line(i) + "\n";
  return std::unique_ptr<std::istream>(new std::istringstream(text));
}

bool Differential(Observable o) {
  return o == Observable::DifferentialLab || o == Observable::DifferentialCM || o == Observable::AnalyzingPower;
}

// A partition of keys by a label ("" = none); groups of one dropped.
std::set<std::set<int>> Partition(const std::map<int, std::string> &labelOfKey) {
  std::map<std::string, std::set<int>> groups;
  for (const auto &kv : labelOfKey)
    if (!kv.second.empty()) groups[kv.second].insert(kv.first);
  std::set<std::set<int>> out;
  for (const auto &g : groups)
    if (g.second.size() > 1) out.insert(g.second);
  return out;
}

void CompareReaction(const std::string &w, const ReactionSpec &a, const ReactionSpec &b) {
  Eq(w + " entrance", a.entranceKey, b.entranceKey);
  // total capture: the engine replaces the exit by the capture pairs of the entrance nucleus (WorkingCopy::MakePoints, Fill)
  if (a.observable != Observable::TotalCapture) Eq(w + " exit", a.exitKey, b.exitKey);
  Eq(w + " observable", (int)a.observable, (int)b.observable);
  Eq(w + " thm", a.thm, b.thm);
  if (a.observable == Observable::PhaseShift) {
    Bits(w + " J", a.phaseJ, b.phaseJ);
    Eq(w + " l", a.phaseL, b.phaseL);
    if (a.phaseS >= 0) Bits(w + " s", a.phaseS, b.phaseS);  // absent: the one spin of (J, l), which the converter writes
  }
  // a composite block that lists no component was read as a plain reaction (with a warning; ComponentList.h)
  const auto plain = [](const ReactionSpec &r) {
    return r.components.empty() && r.combination == ReactionSpec::Combination::Sum ? ReactionSpec::Combination::None : r.combination;
  };
  Eq(w + " combination", (int)plain(a), (int)plain(b));
  Eq(w + " components", a.components.size(), b.components.size());
  for (std::size_t i = 0; i < std::min(a.components.size(), b.components.size()); i++) {
    const ComponentSpec &x = a.components[i], &y = b.components[i];
    const std::string c = w + " component " + std::to_string(i + 1);
    Eq(c + " entrance", x.entranceKey, y.entranceKey);
    Eq(c + " exit", x.exitKey, y.exitKey);
    Eq(c + " hasAngle", x.hasAngle, y.hasAngle);
    Bits(c + " angle", x.angle, y.angle);
    Bits(c + " scaling", x.scaling, y.scaling);
  }
  Eq(w + " sample given", a.sample != "0", b.sample != "0");
  Eq(w + " total", a.total, b.total);
  Bits(w + " transmission", a.transmission, b.transmission);
  Bits(w + " yield thickness", a.yieldThickness, b.yieldThickness);
  Eq(w + " Coulomb reference", (int)a.coulombReference, (int)b.coulombReference);
  Bits(w + " model factor", a.scale.factor, b.scale.factor);
  Eq(w + " S factor", a.scale.sFactor, b.scale.sFactor);
  Eq(w + " S factor shape g", a.scale.sFactorShape_perMeV, b.scale.sFactorShape_perMeV);
  Eq(w + " Legendre expansion", static_cast<int>(a.legendre), static_cast<int>(b.legendre));
  Eq(w + " reference table rows", a.reference.energyLab_MeV.size(), b.reference.energyLab_MeV.size());
}

// The first active member of each group of `label`.
std::set<int> Leaders(const std::vector<DatasetSpec> &ds, std::string DatasetSpec::*label) {
  std::set<int> leaders;
  std::set<std::string> seen;
  for (const DatasetSpec &d : ds) {
    if (!d.active) continue;
    if ((d.*label).empty() || seen.insert(d.*label).second) leaders.insert(d.key);
  }
  return leaders;
}

void CompareDatasets(const std::vector<DatasetSpec> &a, const std::vector<DatasetSpec> &b) {
  Eq(std::string("number of data sets"), a.size(), b.size());
  if (a.size() != b.size()) return;
  std::map<int, std::string> na, nb, sa, sb, fa, fb, ta, tb;
  for (std::size_t i = 0; i < a.size(); i++) {
    const DatasetSpec &x = a[i], &y = b[i];
    const std::string w = "data set " + std::to_string(x.key);
    Eq(w + " key", x.key, y.key);
    Eq(w + " active", x.active, y.active);
    Eq(w + " file", x.dataFile, y.dataFile);
    if (!x.active || !y.active) continue;
    Bits(w + " energy column factor", x.columns.energy, y.columns.energy);
    Eq(w + " energy column c.m.", x.columns.energyCM, y.columns.energyCM);
    Bits(w + " angle column factor", x.columns.angle, y.columns.angle);
    Bits(w + " value column factor", x.columns.value, y.columns.value);
    CompareReaction(w, x.reaction, y.reaction);
    Bits(w + " minE", x.minE, y.minE);
    Bits(w + " maxE", x.maxE, y.maxE);
    if (Differential(x.reaction.observable)) {
      Bits(w + " minA", x.minA, y.minA);
      Bits(w + " maxA", x.maxA, y.maxA);
    }
    if (x.reaction.observable == Observable::AngularCoefficients) Eq(w + " absolute", x.absolute, y.absolute);
    Eq(w + " covariance", x.covarianceFile, y.covarianceFile);
    Eq(w + " upos", x.upos, y.upos);
    if (x.upos) {
      Eq(w + " upos L", x.secondaryDecayL, y.secondaryDecayL);
      Bits(w + " upos Ic", x.secondaryFinalSpin, y.secondaryFinalSpin);
      Bits(w + " upos delta", x.secondaryMixingRatio, y.secondaryMixingRatio);
    }
    // sqrt term and function: each data set's own (the function on its owner, the first active member of a group)
    Bits(w + " sqrt", x.hasSqrtShift ? x.energyShiftSqrt : 0.0, y.hasSqrtShift ? y.energyShiftSqrt : 0.0);
    Eq(w + " sqrt vary", x.hasSqrtShift && x.varyEnergyShiftSqrt, y.hasSqrtShift && y.varyEnergyShiftSqrt);
    Bits(w + " sqrt error", x.hasSqrtShift ? x.energyShiftSqrtError : 0.0, y.hasSqrtShift ? y.energyShiftSqrtError : 0.0);
    Eq(w + " function", x.shiftFunction, y.shiftFunction);
    Bits(w + " function values", x.shiftFunctionValues, y.shiftFunctionValues);
    Bits(w + " function errors", x.shiftFunctionErrors, y.shiftFunctionErrors);
    Eq(w + " function vary", x.shiftFunctionVary, y.shiftFunctionVary);
    na[x.key] = x.normGroup;
    nb[y.key] = y.normGroup;
    sa[x.key] = x.shiftGroup;
    sb[y.key] = y.shiftGroup;
    fa[x.key] = x.shiftFunctionGroup;
    fb[y.key] = y.shiftFunctionGroup;
    if (x.thicknessFit) ta[x.key] = x.thicknessGroup;
    if (y.thicknessFit) tb[y.key] = y.thicknessGroup;
    Eq(w + " thickness fit", x.thicknessFit, y.thicknessFit);
  }
  Eq(std::string("norm groups"), Partition(na), Partition(nb));
  Eq(std::string("shift groups"), Partition(sa), Partition(sb));
  Eq(std::string("shift-function groups"), Partition(fa), Partition(fb));
  Eq(std::string("thickness groups"), Partition(ta), Partition(tb));
  const std::set<int> normLeaders = Leaders(a, &DatasetSpec::normGroup), shiftLeaders = Leaders(a, &DatasetSpec::shiftGroup);
  const std::set<int> thicknessLeaders = Leaders(a, &DatasetSpec::thicknessGroup);
  for (std::size_t i = 0; i < a.size(); i++) {
    const DatasetSpec &x = a[i], &y = b[i];
    if (!x.active || !y.active) continue;
    const std::string w = "data set " + std::to_string(x.key);
    if (normLeaders.count(x.key)) {
      Eq(w + " vary norm", x.varyNorm, y.varyNorm);
      if (!(x.reaction.thm && x.varyNorm)) {
        Bits(w + " norm", x.norm, y.norm);
        Bits(w + " norm error (%)", x.normError, y.normError);
        Bits(w + " norm prior centre", x.hasNormCentre ? x.normCentre : x.norm, y.hasNormCentre ? y.normCentre : y.norm);
      }
    }
    if (shiftLeaders.count(x.key)) {
      Bits(w + " shift", x.energyShift, y.energyShift);
      Eq(w + " vary shift", x.varyEnergyShift, y.varyEnergyShift);
      if (x.varyEnergyShift) Bits(w + " shift error", x.energyShiftError, y.energyShiftError);  // a fixed shift pays no prior
      Bits(w + " shift prior centre", x.hasShiftCentre ? x.shiftCentre : 0.0, y.hasShiftCentre ? y.shiftCentre : 0.0);
    }
    if (x.thicknessFit && thicknessLeaders.count(x.key)) {
      Eq(w + " thickness vary", x.thicknessVary, y.thicknessVary);
      Bits(w + " thickness sigma", x.thicknessSigma, y.thicknessSigma);
    }
  }
}

void CompareExtrapolations(const std::vector<ExtrapolationSpec> &a, const std::vector<ExtrapolationSpec> &b) {
  Eq(std::string("number of extrapolations"), a.size(), b.size());
  if (a.size() != b.size()) return;
  for (std::size_t i = 0; i < a.size(); i++) {
    const ExtrapolationSpec &x = a[i], &y = b[i];
    const std::string w = "extrapolation " + std::to_string(x.key);
    Eq(w + " key", x.key, y.key);
    Eq(w + " active", x.active, y.active);
    if (!x.active || !y.active) continue;
    CompareReaction(w, x.reaction, y.reaction);
    Bits(w + " energies", x.energies, y.energies);
    if (Differential(x.reaction.observable)) Bits(w + " angles", x.angles, y.angles);
    if (x.reaction.observable == Observable::AngularDistribution) Eq(w + " max order", x.maxAngDistOrder, y.maxAngDistOrder);
  }
}

// The effect the engine applies to each active key (a later line takes over a key), and the line's index.
std::map<int, int> EffectOfKey(const SpecList<TargetEffectSpec> &effects, const std::set<int> &activeKeys) {
  std::map<int, int> out;
  for (std::size_t i = 0; i < effects.items.size(); i++) {
    const TargetEffectSpec &e = effects.items[i].spec;
    if (!e.active) continue;
    for (int k : e.datasetKeys)
      if (activeKeys.count(k)) out[k] = (int)i;
  }
  return out;
}

bool FoldsAnything(const TargetEffectSpec &e) {
  return e.convolution || e.convolutionCoefficients || e.targetIntegration || e.doppler || e.beamProfile || e.udr;
}

void CompareEffect(const std::string &w, const TargetEffectSpec &x, const TargetEffectSpec &y, bool densityLeader) {
  Eq(w + " integration", x.targetIntegration, y.targetIntegration);
  if (x.targetIntegration && densityLeader) {  // a density-group member takes its first line's (section 8.1)
    Bits(w + " density", x.density, y.density);
  }
  if (x.targetIntegration) {
    Eq(w + " stopping", x.stoppingPowerEquation, y.stoppingPowerEquation);
    Bits(w + " stopping parameters", x.stoppingPowerParameters, y.stoppingPowerParameters);
    const bool sa = x.hasStragglingFlag && x.straggling, sb = y.hasStragglingFlag && y.straggling;
    Eq(w + " straggling", sa, sb);
    if (sa) Bits(w + " straggling coefficient", x.hasStragglingCoefficient ? x.stragglingCoefficient : 0.0,
                 y.hasStragglingCoefficient ? y.stragglingCoefficient : 0.0);
  }
  if (x.targetIntegration && densityLeader) {
    Eq(w + " density fit", x.densityFit, y.densityFit);
    if (x.densityFit) {
      Eq(w + " density vary", x.densityVary, y.densityVary);
      Eq(w + " density sigma", x.densitySigmaText.empty() ? std::string("0") : x.densitySigmaText,
         y.densitySigmaText.empty() ? std::string("0") : y.densitySigmaText);
      Bits(w + " density max", x.densityMax, y.densityMax);
    }
  }
  const bool convolution = x.convolution && !x.convolutionCoefficients;  // the engine uses the expression when both
  Eq(w + " convolution", convolution, y.convolution);
  if (convolution) Bits(w + " sigma", x.sigma, y.sigma);
  Eq(w + " sigma(E)", x.convolutionCoefficients, y.convolutionCoefficients);
  if (x.convolutionCoefficients) {
    Eq(w + " sigma(E) equation", x.convolutionEquation, y.convolutionEquation);
    Bits(w + " sigma(E) coefficients", x.convolutionCoefficientValues, y.convolutionCoefficientValues);
    Eq(w + " sigmalab", x.sigmaLab, y.sigmaLab);
  }
  Eq(w + " Q", x.qCoefficients, y.qCoefficients);
  if (x.qCoefficients) Bits(w + " Q values", x.qCoefficientValues, y.qCoefficientValues);
  if (FoldsAnything(x)) {
    Eq(w + " points", x.numIntegrationPoints, y.numIntegrationPoints);
    // given or left out matters (TargetEffect::GetResonanceWidthMultiplier: a pure convolution's core when left out)
    Eq(w + " width multiplier given", x.hasResonanceWidthMultiplier, y.hasResonanceWidthMultiplier);
    Bits(w + " width multiplier", x.hasResonanceWidthMultiplier ? x.resonanceWidthMultiplier : 20.0,
         y.hasResonanceWidthMultiplier ? y.resonanceWidthMultiplier : 20.0);
    Bits(w + " points per width", x.hasPointsPerWidth ? x.pointsPerWidth : 50.0, y.hasPointsPerWidth ? y.pointsPerWidth : 50.0);
    const bool ra = x.hasRanges && (x.autoWindows || !x.ranges.empty()), rb = y.hasRanges && (y.autoWindows || !y.ranges.empty());
    Eq(w + " windows", ra, rb);
    if (ra) {
      Eq(w + " auto windows", x.autoWindows, y.autoWindows);
      if (!x.autoWindows) Eq(w + " ranges", x.ranges, y.ranges);
      Eq(w + " auto tolerance given", x.hasAutoWindowTolerance, y.hasAutoWindowTolerance);
      Bits(w + " auto tolerance", x.autoWindowTolerance, y.autoWindowTolerance);
      Eq(w + " auto margin given", x.hasAutoWindowMargin, y.hasAutoWindowMargin);
      Bits(w + " auto margin", x.autoWindowMargin, y.autoWindowMargin);
      Bits(w + " blend", x.transitionWidth, y.transitionWidth);
      Bits(w + " tolerance", x.autoTolerance, y.autoTolerance);
    }
  }
  Eq(w + " profile", x.beamProfile, y.beamProfile);
  if (x.beamProfile) {
    Eq(w + " profile components", x.beamProfileComponents.size(), y.beamProfileComponents.size());
    for (std::size_t i = 0; i < std::min(x.beamProfileComponents.size(), y.beamProfileComponents.size()); i++) {
      const auto &p = x.beamProfileComponents[i], &q = y.beamProfileComponents[i];
      Bits(w + " xi", p.xi, q.xi);
      Bits(w + " omega", p.omega, q.omega);
      Bits(w + " alpha", p.alpha, q.alpha);
      Bits(w + " weight", p.weight, q.weight);
    }
    Bits(w + " tpc sigma", x.beamTpcSigma, y.beamTpcSigma);
    Bits(w + " truncation", x.beamTruncation, y.beamTruncation);
    Eq(w + " detailed balance", x.beamDetailedBalance, y.beamDetailedBalance);
  }
  Eq(w + " screening", x.screening, y.screening);
  if (x.screening) {
    Bits(w + " Ue", x.screeningUe, y.screeningUe);
    Eq(w + " Ue vary", x.screeningVary, y.screeningVary);
    Bits(w + " Ue sigma", x.screeningSigma, y.screeningSigma);
  }
  Eq(w + " doppler", x.doppler, y.doppler);
  if (x.doppler) {
    Bits(w + " temperature", x.dopplerTemperature, y.dopplerTemperature);
    Eq(w + " doppler fit", x.dopplerFit, y.dopplerFit);
    if (x.dopplerFit) {
      Bits(w + " doppler vary", std::strtod(x.dopplerVaryText.c_str(), nullptr), std::strtod(y.dopplerVaryText.c_str(), nullptr));
      Bits(w + " doppler sigma", std::strtod(x.dopplerSigmaText.c_str(), nullptr), std::strtod(y.dopplerSigmaText.c_str(), nullptr));
      Bits(w + " doppler max", std::strtod(x.dopplerMaxText.c_str(), nullptr), std::strtod(y.dopplerMaxText.c_str(), nullptr));
    }
  }
  Eq(w + " udr", x.udr, y.udr);
  if (x.udr) {
    Eq(w + " udr file", x.udrFile, y.udrFile);
    Bits(w + " flight path", x.udrFlightPath, y.udrFlightPath);
    Bits(w + " burst", x.udrBurstFwhm, y.udrBurstFwhm);
    Bits(w + " channel", x.udrChannelWidth, y.udrChannelWidth);
    Eq(w + " centred", x.udrCentred, y.udrCentred);
  }
  Eq(w + " kernel", x.kernel, y.kernel);
  Eq(w + " kernel parameters", x.kernelParameters.size(), y.kernelParameters.size());
  for (size_t i = 0; i < x.kernelParameters.size() && i < y.kernelParameters.size(); i++)
    Bits(w + " kernel parameter", x.kernelParameters[i], y.kernelParameters[i]);
}

void CompareEffects(const std::string &what, const SpecList<TargetEffectSpec> &a, const SpecList<TargetEffectSpec> &b,
                    const std::set<int> &activeKeys) {
  std::map<int, int> ea = EffectOfKey(a, activeKeys), eb = EffectOfKey(b, activeKeys);
  // old lines that apply nothing (no fold, no Q, no screening) are no effect in format 3
  for (auto it = ea.begin(); it != ea.end();) {
    const TargetEffectSpec &e = a.items[it->second].spec;
    if (!FoldsAnything(e) && !e.qCoefficients && !e.screening) it = ea.erase(it);
    else ++it;
  }
  std::set<int> ka, kb;
  for (const auto &kv : ea) ka.insert(kv.first);
  for (const auto &kv : eb) kb.insert(kv.first);
  Eq(what + ": keys with an effect", ka, kb);
  // the first line of each density group (WorkingCopy::SetUpDensityParams: the group's density is its first line's)
  std::set<std::string> groupSeen;
  std::set<int> densityLeaderLines;
  for (std::size_t i = 0; i < a.items.size(); i++) {
    const TargetEffectSpec &e = a.items[i].spec;
    if (!e.active) continue;
    if (e.densityGroup.empty() || groupSeen.insert(e.densityGroup).second) densityLeaderLines.insert((int)i);
  }
  // effects whose copies share a parameter of their own (screening_<n>, doppler_<n>): which keys share one; effects
  // without one differ only in their number, which names nothing (the converter merges identical ones)
  std::map<int, std::string> sharedA, sharedB, densityA, densityB;
  for (const auto &kv : ea) {
    if (!eb.count(kv.first)) continue;
    const TargetEffectSpec &x = a.items[kv.second].spec, &y = b.items[eb[kv.first]].spec;
    CompareEffect(what + " of key " + std::to_string(kv.first), x, y, densityLeaderLines.count(kv.second) != 0);
    if (x.screening || x.dopplerFit) sharedA[kv.first] = "line" + std::to_string(kv.second);
    if (y.screening || y.dopplerFit) sharedB[kv.first] = "effect" + std::to_string(eb[kv.first]);
    const std::string la = "line" + std::to_string(kv.second), lb = "effect" + std::to_string(eb[kv.first]);
    const TargetEffectSpec &lead = x;
    if (lead.densityFit || !lead.densityGroup.empty()) densityA[kv.first] = x.densityGroup.empty() ? la : "group " + x.densityGroup;
    if (y.densityFit) densityB[kv.first] = y.densityGroup.empty() ? lb : "group " + y.densityGroup;
  }
  Eq(what + ": keys sharing one effect parameter", Partition(sharedA), Partition(sharedB));
  Eq(what + ": keys sharing one fitted density", Partition(densityA), Partition(densityB));
}

}  // namespace

// A model datum as the engine's row: a scattering length as its constraint (ModelData::Read), else its fields.
std::string DatumRow(const ModelDatum &d) {
  char buf[512];
  if (d.quantity == ModelQuantity::ScatteringLength)
    std::snprintf(buf, sizeof buf, "scattering length pair %d s %a Re %a +- %a Im %a +- %a", d.pair,
                  d.coherent ? -1.0 : 0.5 * d.twoS, d.value, d.error, d.hasImag ? d.imag : 0.0, d.hasImag ? d.imagError : 0.0);
  else
    std::snprintf(buf, sizeof buf, "%s level %d (%d, %d) %s pairs %d %d %d frame '%s' l %d 2s %d value %a +- %a unit %s",
                  d.quantityName.c_str(), d.levelId, d.group, d.level, d.jpi.c_str(), d.entrancePair, d.exitPair, d.pair,
                  d.frame.c_str(), d.l, d.twoS, d.value, d.error, d.unit.c_str());
  return buf;
}

// The model-data rows of an old project in the engine's order: <scatteringLengths> lines with a residual, then the
// datums of its <modelData> section.
std::vector<std::string> OldModelDataRows(const io::Document &doc) {
  std::vector<std::string> rows;
  if (const io::Document::Section *s = doc.Find("scatteringLengths"))
    for (const std::string &line : doc.Body(*s)) {
      const convert::ScatteringLengthLine l = convert::ReadScatteringLengthLine(line);
      if (l.kind != convert::ScatteringLengthLine::Constraint || l.active == 0) continue;
      if (!(l.sigma > 0.0) && !(l.imSigma > 0.0)) continue;
      char buf[256];
      std::snprintf(buf, sizeof buf, "scattering length pair %d s %a Re %a +- %a Im %a +- %a", l.pairKey, l.s, l.value,
                    l.sigma, l.hasImaginary ? l.imValue : 0.0, l.hasImaginary ? l.imSigma : 0.0);
      rows.push_back(buf);
    }
  if (const io::Document::Section *s = doc.Find("modelData")) {
    const io::RecordReadResult read = io::ReadRecords(doc, *s, Schema::KindsOfSection("modelData"));
    const ValidationResult v = Schema::Validate({{"modelData", read.set.records}});
    bool ok = false;
    const LevelTables levels = LevelTablesOfRows(LevelRowsOf(doc, ok));
    for (const ModelDatum &d : BuildModelData(v, levels.levels, levels.pairs).data) rows.push_back(DatumRow(d));
  }
  return rows;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: spec_equivalence_test <repository root>\n");
    return 2;
  }
  const std::string root = argv[1];
  if (const char *inject = std::getenv("AZURE2_SPEC_INJECT")) injected = std::atoi(inject);
  std::vector<std::string> files;
  // The tree's projects and the originals of the converted ones (tests/convert/fixtures/<path>, the original of
  // <path>; docs/upgrade/format3_tree.md): an original is read in place of its conversion, with the data files of
  // its own folder.
  const std::string fixtures = root + "/tests/convert/fixtures";
  Find(root + "/tests", "tests/", files);
  Find(root + "/examples", "examples/", files);
  Find(fixtures + "/tests", "tests/", files);
  Find(fixtures + "/examples", "examples/", files);
  std::sort(files.begin(), files.end());
  files.erase(std::unique(files.begin(), files.end()), files.end());
  int passed = 0, failed = 0, notConverted = 0, filesNotAtHand = 0, recorded = 0;
  for (const std::string &rel : files) {
    std::string text;
    if (!Slurp(fixtures + "/" + rel, text) && !Slurp(root + "/" + rel, text)) continue;
    const std::string dir = root + "/" + rel.substr(0, rel.find_last_of('/'));
    convert::ConvertOptions o;
    o.projectDir = dir;
    o.exactGrids = true;
    const convert::ConvertResult r = convert::ConvertToFormat3(text, o);
    if (!r.ok) {
      notConverted++;
      std::printf("  not converted  %s: %s\n", rel.c_str(), r.refusal.c_str());
      continue;
    }
    problems.clear();
    // old lines -> specs (the engine's reading of a format-2 document)
    io::Document doc = io::Document::Parse(text);
    io::ConvertToFormat2(doc);
    const bool relativistic = RelativisticFlag(doc);
    DataSpecs oldData;
    ExtrapolationSpecs oldExtrap;
    if (auto s = Section(doc, "segmentsData")) {
      oldData.datasets = convert::DatasetSpecsFromOldLines(*s, relativistic);
      auto settings = Section(doc, "parameterSettings");
      if (!oldData.datasets.failed) oldData.priorCentreError = convert::PriorCentresFromOldRows(settings.get(), oldData.datasets);
    }
    if (auto s = Section(doc, "segmentsTest")) oldExtrap.extrapolations = convert::ExtrapolationSpecsFromOldLines(*s, relativistic);
    if (auto s = Section(doc, "targetInt")) oldData.targetEffects = oldExtrap.targetEffects = convert::TargetEffectSpecsFromOldLines(*s);
    auto abundances = Section(doc, "abundances");
    const AbundanceSpecs oldAbundances = convert::AbundanceSpecsFromOldLines(abundances.get());
    // the old THM membership: experiment[<name>] segments=..., weight[k]=, weightTest[k]= (the old grammar's parsers)
    std::map<int, std::string> oldExperimentOf, oldWeightOfData, oldWeightOfTest;
    if (const io::Document::Section *s = doc.Find("thm"))
      for (std::string line : doc.Body(*s)) {
        line = line.substr(0, line.find('#'));
        std::istringstream words(line);
        std::string first, token;
        words >> first;
        if (first.compare(0, 11, "experiment[") == 0) {
          const std::string name = first.substr(11, first.find(']') - 11);
          while (words >> token)
            if (token.compare(0, 9, "segments=") == 0) {
              std::vector<int> keys;
              convert::ParseThmSegmentList(token.substr(9), keys);
              for (int k : keys) oldExperimentOf[k] = name;
            }
          continue;
        }
        const std::size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq), value = line.substr(eq + 1);
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t\r") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t\r") + 1);
        convert::ThmWeightKey w;
        bool valid = false;
        if (convert::ReadThmWeightKey(key, value, w, valid) && valid) (w.test ? oldWeightOfTest : oldWeightOfData)[w.segment] = w.file;
      }
    // the conversion -> records -> Validate -> Build
    // the engine opens data files relative to its working directory, the test's top-level folder (tools/run_one.sh);
    // files that a test's check.sh writes are not at hand: then the build runs without reading data files
    Format3ReadOptions read;
    read.fileName = rel;
    read.projectDir = root + "/" + rel.substr(0, rel.find('/', rel.find('/') + 1));
    read.setUpChecks = false;  // the specs the engine builds: its set-up makes the external-capture checks itself
    Format3Project p = ReadFormat3Text(r.text, read);
    bool onlyFiles = !p.Ok();
    for (const Diagnostic &d : p.diagnostics) onlyFiles = onlyFiles && (d.severity != Severity::Error || d.code == codes::kDataFile);
    if (onlyFiles) {
      read.readDataFiles = false;
      p = ReadFormat3Text(r.text, read);
      filesNotAtHand++;
    }
    if (!p.Ok()) {
      for (const Diagnostic &d : p.diagnostics) problems.push_back("format 3: " + d.ToText());
    } else {
      BuiltProject b = p.built;
      if (injected == 1 && !b.data.datasets.items.empty()) b.data.datasets.items[0].spec.maxE = std::nextafter(b.data.datasets.items[0].spec.maxE, 1e300);
      std::vector<DatasetSpec> da, db;
      for (const auto &i : oldData.datasets.items) da.push_back(i.spec);
      for (const auto &i : b.data.datasets.items) db.push_back(i.spec);
      CompareDatasets(da, db);
      std::vector<ExtrapolationSpec> xa, xb;
      for (const auto &i : oldExtrap.extrapolations.items) xa.push_back(i.spec);
      for (const auto &i : b.extrapolations.extrapolations.items) xb.push_back(i.spec);
      CompareExtrapolations(xa, xb);
      std::set<int> activeData, activeExtrap;
      for (const DatasetSpec &d : da)
        if (d.active) activeData.insert(d.key);
      for (const ExtrapolationSpec &x : xa)
        if (x.active) activeExtrap.insert(x.key);
      CompareEffects("data effect", oldData.targetEffects, b.data.targetEffects, activeData);
      CompareEffects("extrapolation effect", oldExtrap.targetEffects, b.extrapolations.targetEffects, activeExtrap);
      // abundances: per (sample, compound nucleus); the old sample name maps to the format-3 one of its data sets
      std::map<std::string, std::string> sampleName;
      for (std::size_t i = 0; i < std::min(da.size(), db.size()); i++)
        if (da[i].active && da[i].reaction.sample != "0") sampleName[da[i].reaction.sample] = db[i].reaction.sample;
      std::map<std::string, std::string> rowsA, rowsB;
      auto key = [](const std::string &s, const AbundanceRowSpec &r) { return s + " " + std::to_string(r.z) + "_" + std::to_string(r.a); };
      auto row = [](const AbundanceRowSpec &r) {
        char buf[128];
        std::snprintf(buf, sizeof buf, "%a %d %a", r.value, (int)r.vary, r.sigma);
        return std::string(buf);
      };
      std::set<std::string> usedNew;
      for (const auto &kv : sampleName) usedNew.insert(kv.second);
      for (const AbundanceRowSpec &r : oldAbundances.abundances.rows) {
        if (r.isSum) continue;
        std::string s = r.sample.empty() ? "1" : r.sample;
        if (!sampleName.count(s) && sampleName.count("1") && r.sample.empty()) s = "1";
        if (!sampleName.count(s)) {
          // "1" names the unnamed sample or the only one (CNuc::DefaultSample)
          if (sampleName.count("1")) s = "1";
          else continue;  // a sample no active data set uses: dropped by the converter with a note
        }
        rowsA[key(sampleName[s], r)] = row(r);
      }
      for (const AbundanceRowSpec &r : b.abundances.abundances.rows)
        if (!r.isSum && usedNew.count(r.sample)) rowsB[key(r.sample, r)] = row(r);
      for (auto &kv : rowsA)
        if (rowsB.count(kv.first)) Eq("abundance " + kv.first, kv.second, rowsB[kv.first]);
        else problems.push_back("abundance " + kv.first + " missing in format 3");
      // THM
      std::map<int, std::string> expA, expB;
      for (const auto &kv : oldExperimentOf)
        if (activeData.count(kv.first)) expA[kv.first] = kv.second;
      for (const auto &x : b.thm.experiments)
        for (int k : x.datasetKeys) expB[k] = x.name;
      Eq(std::string("THM experiments"), expA, expB);
      std::map<int, std::string> wa = oldWeightOfData, wb;
      for (const auto &w : b.thm.datasetWeights) wb[w.first] = w.second;
      Eq(std::string("THM data weights"), wa, wb);
      std::map<int, std::string> twa = oldWeightOfTest, twb;
      for (const auto &w : b.thm.extrapolationWeights) twb[w.first] = w.second;
      for (auto it = twa.begin(); it != twa.end();) it = activeExtrap.count(it->first) ? std::next(it) : twa.erase(it);
      Eq(std::string("THM extrapolation weights"), twa, twb);
      std::vector<std::string> modelB;
      for (std::size_t k = 0; k < b.modelData.size(); k++) {
        ModelDatum d = b.modelData[k];
        if (injected == 2 && k == 0) d.value = std::nextafter(d.value, 1e300);
        modelB.push_back(DatumRow(d));
      }
      // parameter names: the converter's map is injective and takes the adaptor's names to Build's
      std::map<std::string, std::string> map;
      std::set<std::string> targets;
      for (const auto &pn : r.parameterNames) {
        if (map.count(pn.first) && map[pn.first] != pn.second) problems.push_back("parameter map: " + pn.first + " maps twice");
        if (!map.count(pn.first) && !targets.insert(pn.second).second) problems.push_back("parameter map: two names map to " + pn.second);
        map[pn.first] = pn.second;
      }
      if (injected == 3 && !db.empty()) db[0].names.energyShift += "x";
      auto name = [&](const std::string &oldName, const std::string &newName) {
        const auto it = map.find(oldName);
        if (it == map.end()) problems.push_back("parameter map: no new name for " + oldName + " (Build: " + newName + ")");
        else if (it->second != newName) problems.push_back("parameter name " + oldName + ": map " + it->second + ", Build " + newName);
      };
      for (std::size_t i = 0; i < std::min(da.size(), db.size()); i++) {
        if (!da[i].active) continue;
        const ParameterNames &a = da[i].names, &n = db[i].names;
        if (da[i].varyNorm && !da[i].reaction.thm) name(a.norm, n.norm);  // a THM norm is profiled (Segment::IsProfiledNorm)
        name(a.energyShift, n.energyShift);
        name(a.energyShiftSqrt, n.energyShiftSqrt);
        for (std::size_t j = 0; j < std::min(a.shiftFunction.size(), n.shiftFunction.size()); j++)
          name(a.shiftFunction[j], n.shiftFunction[j]);
        if (da[i].thicknessFit) name(a.thickness, n.thickness);
      }
      // hybrid potentials, per particle pair
      experiment::HybridPotentialSettings oldPotentials;
      if (auto s = Section(doc, "potential")) oldPotentials = convert::ReadPotentialBlock(*s).settings;
      HybridPotentialSettings newPotentials = b.potentials;
      if (injected == 4 && !newPotentials.pairs.empty())
        newPotentials.pairs.begin()->second.coulomb_radius_fm = std::nextafter(newPotentials.pairs.begin()->second.coulomb_radius_fm, 1e300);
      for (const PairInfo &pair : p.levels.pairs) {
        if (pair.kind != 0) continue;
        auto keyText = [](const std::optional<physics::HybridPotential> &v) {
          std::string t;
          char buf[64];
          for (double x : physics::PotentialKey(v)) {
            std::snprintf(buf, sizeof buf, " %a", x);
            t += buf;
          }
          return t;
        };
        const std::string a = keyText(oldPotentials.PotentialFor(pair.key)), n = keyText(newPotentials.PotentialFor(pair.key));
        if (a != n) problems.push_back("potential of pair " + std::to_string(pair.key) + ":" + a + " vs" + n);
      }
      const std::vector<std::string> modelA = OldModelDataRows(doc);
      for (std::size_t k = 0; k < std::max(modelA.size(), modelB.size()); k++)
        if (k >= modelA.size() || k >= modelB.size() || modelA[k] != modelB[k])
          problems.push_back("model data row " + std::to_string(k + 1) + ": " + (k < modelA.size() ? modelA[k] : "-") +
                             " vs " + (k < modelB.size() ? modelB[k] : "-"));
    }
    const auto known = kRecorded.find(rel);
    if (known != kRecorded.end()) {
      if (problems.empty()) {
        failed++;
        std::printf("  FAIL  %s: listed as a recorded difference, but the specs are equal now\n", rel.c_str());
      } else {
        recorded++;
        std::printf("  recorded  %s: %s\n", rel.c_str(), known->second);
      }
    } else if (problems.empty()) {
      passed++;
    } else {
      failed++;
      std::printf("  FAIL  %s (%zu)\n", rel.c_str(), problems.size());
      for (std::size_t i = 0; i < problems.size() && i < 12; i++) std::printf("        %s\n", problems[i].c_str());
    }
  }
  std::printf("spec_equivalence: %zu files, %d converted: %d equal, %d recorded differences, %d differ; %d not converted; "
              "%d built without their data files (written by a check.sh)\n",
              files.size(), passed + failed + recorded, passed, recorded, failed, notConverted, filesNotAtHand);
  return failed == 0 && passed > 0 ? 0 : 1;
}
