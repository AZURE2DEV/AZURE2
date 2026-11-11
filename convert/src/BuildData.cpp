// Data sets, extrapolations and the targets of neutron-sample thicknesses (Builder.h; section 8.1).
#include "Builder.h"

namespace azure2 {
namespace convert {

// ================================================================ neutron samples
std::string Builder::SampleThicknessTarget(const std::string &record, int entrance, double n, const std::string &line, bool fit,
                                           int vary, double sigma, const std::string &group, int key, bool isData) {
  Count(Construct::kThickness);
  if (!group.empty() && isData) {
    Count(Construct::kThicknessGroup);
    const auto it = thicknessTarget_.find(group);
    if (it != thicknessTarget_.end()) {
      const OldDataLine &lead = *activeData_[thicknessLeader_[group]];
      const double leadN = lead.parsed->transmission() > 0 ? lead.parsed->transmission() : lead.parsed->yieldThickness();
      if (n != leadN)
        Note(Construct::kThicknessGroup, record, "its own sample thickness " + NumberText(n, line) +
                                                     " was not used by the old engine (the group's is its first line's): dropped");
      // a later member's own parameter, fixed by the engine at the group's value
      if (fit) Param("segment_" + std::to_string(key) + "_thickness", it->second + ".areal_density[" + record + "]");
      return it->second;
    }
  }
  io::Record t;
  t.kind = "target";
  t.name = Unique(group.empty() || !isData ? "tgt_" + record : Sanitized(group));
  Add(t, "kind", "sample");
  const OldPair *p = Pair(entrance);
  Add(t, "active", p ? Nuclide(p->z2, p->A2()) : "?");
  const std::string nt = NumberText(n, line);
  Add(t, "areal_density", Nuisance(nt, "atoms/b", fit && vary == 1 && isData, fit && sigma > 0 && isData, nt, NumberText(sigma, line)));
  out_.targets.records.push_back(t);
  if (fit && isData) Param("segment_" + std::to_string(key) + "_thickness", t.name + ".areal_density");
  if (!group.empty() && isData) {
    thicknessTarget_[group] = t.name;
    thicknessLeader_[group] = key;
  }
  return t.name;
}

// ================================================================ data sets
void Builder::BuildData() {
  // Groups: the leader is the first active line with the label (file order); its values are the group's.
  std::map<std::string, int> normLeader, shiftLeader, fnLeader;
  std::map<std::string, std::string> normName, scaleName, fnName;
  for (const OldDataLine &d : old_.data) {
    if (!activeData_.count(d.key)) continue;
    const SegLine &s = *d.parsed;
    if (!s.normGroup().empty() && !normLeader.count(s.normGroup())) normLeader[s.normGroup()] = d.key;
    if (!s.shiftGroup().empty() && !shiftLeader.count(s.shiftGroup())) shiftLeader[s.shiftGroup()] = d.key;
    if (!s.shiftFnGroup().empty() && !fnLeader.count(s.shiftFnGroup())) fnLeader[s.shiftFnGroup()] = d.key;
  }
  // Group records first: their names before the data sets reference them.
  for (const auto &g : normLeader) normName[g.first] = Unique(Sanitized(g.first));
  for (const auto &g : shiftLeader) scaleName[g.first] = Unique(Sanitized(g.first));
  for (const auto &g : fnLeader) fnName[g.first] = Unique(Sanitized(g.first));

  // The function of a shift-function owner as an energy_scale's fields.
  auto addFunction = [&](io::Record &r, const OldDataLine &owner, const std::string &scale) {
    const SegLine &s = *owner.parsed;
    std::vector<io::Token> fn = {Quoted(s.shiftFn()), Tok("MeV")};
    Add(r, "function", fn);
    for (int j = 0; j < (int)s.shiftFnValues().size(); j++) {
      const std::string v = NumberText(s.shiftFnValues()[j], owner.text);
      const double err = s.shiftFnErrors()[j];
      std::vector<io::Token> p = {Tok("a" + std::to_string(j))};
      for (const io::Token &t : Nuisance(v, "", s.shiftFnVary()[j] == 1, err != 0, "0", Abs(NumberText(err, owner.text)))) p.push_back(t);
      Add(r, "param", p);
      Param("segment_" + std::to_string(owner.key) + "_energy_shift_fn" + std::to_string(j), scale + ".a" + std::to_string(j));
    }
  };
  auto offsetNuisance = [&](const OldDataLine &d) {
    const SegLine &s = *d.parsed;
    const std::string centre = old_.shiftCentreText.count(d.key) ? old_.shiftCentreText.at(d.key) : "0";
    return Nuisance(NumberText(s.energyShift(), d.text), "MeV", s.varyEnergyShift() == 1, s.energyShiftError() != 0, centre,
                    Abs(NumberText(s.energyShiftError(), d.text)));
  };

  for (const auto &g : normLeader) {
    Count(Construct::kNormGroup);
    const OldDataLine &d = *activeData_[g.second];
    const SegLine &s = *d.parsed;
    io::Record r;
    r.kind = "normalization";
    r.name = normName[g.first];
    const std::string n = NumberText(s.dataNorm(), d.text);
    const std::string centre = old_.normCentreText.count(d.key) ? old_.normCentreText.at(d.key) : n;
    Add(r, "value", Nuisance(n, "", s.varyNorm() == 1, s.dataNormError() != 0, centre, Abs(NumberText(s.dataNormError(), d.text)) + "%"));
    groupRecords_.push_back(r);
  }
  for (const auto &g : shiftLeader) {
    Count(Construct::kShiftGroup);
    const OldDataLine &d = *activeData_[g.second];
    io::Record r;
    r.kind = "energy_scale";
    r.name = scaleName[g.first];
    Add(r, "offset", offsetNuisance(d));
    if (!d.parsed->shiftFn().empty()) addFunction(r, d, r.name);
    groupRecords_.push_back(r);
  }
  for (const auto &g : fnLeader) {
    Count(Construct::kShiftFnGroup);
    const OldDataLine &d = *activeData_[g.second];
    io::Record r;
    r.kind = "energy_scale";
    r.name = fnName[g.first];
    if (d.parsed->shiftFn().empty()) {
      res_.refusal = "the first line of shift-function group " + g.first + " declares no shiftfn, which the old engine refused";
      return;
    }
    addFunction(r, d, r.name);
    groupRecords_.push_back(r);
  }

  for (const OldDataLine &d : old_.data) {
    Count(Construct::kDataLine);
    const SegLine &s = *d.parsed;
    const std::string &name = dataName_[d.key];
    const bool active = activeData_.count(d.key) != 0;
    io::Record r;
    r.kind = "dataset";
    r.name = name;
    r.comments.push_back("# <segmentsData> line " + std::to_string(d.key));
    const bool thm = s.isDiff() >= 10;
    const int code = thm ? s.isDiff() - 10 : s.isDiff();
    Count(Construct::kDataCode);
    if (thm) Count(Construct::kThm);
    const Apparatus *app = active && apparatusOfData_.count(d.key) ? &apparatus_[apparatusOfData_[d.key]] : nullptr;
    const OldPair *in = Pair(s.entranceKey());
    const bool beta = in && in->pType == 20;
    const bool thickTarget = app && app->hasTarget && app->line->line.isTargetIntegration == 1 && !(app->line->line.density > 0.0);
    const bool thinTarget = app && app->hasTarget && app->line->line.isTargetIntegration == 1 && app->line->line.density > 0.0;
    const bool ratio = s.isAdvanced() == 1 && s.operationType() == 1 && !s.components().empty();
    // The old engine refused the exit key -1 with any data type but total capture (src/EData.cpp, TotalCaptureCodeOk):
    // format 3 writes both halves as one token ([total]), so converting the line would make a refused file run.
    if (active && in && s.exitKey() == -1 && code != 3) {
      res_.refusal = "data line " + std::to_string(d.key) +
                     ": the exit key -1 (total capture) needs the total-capture data type (3 on a data line); the old engine refused it";
      return;
    }

    // reaction
    const std::string angleFrame = (code == 4 || code == 7) ? "cm" : "lab";
    std::vector<io::Token> reaction;
    std::vector<std::pair<std::string, std::string>> pairKeys;
    std::string why;
    if (!ReactionTokens(name, s.entranceKey(), (code == 3 || s.exitKey() == -1) ? -1 : s.exitKey(), s.isAdvanced(), s.operationType(), s.components(),
                        angleFrame, reaction, why, pairKeys)) {
      if (!active) {
        Note(Construct::kDataLine, name, "inactive data line " + std::to_string(d.key) + ": " + why + "; dropped");
        continue;
      }
      Note(Construct::kDataLine, name, "data line " + std::to_string(d.key) + ": " + why + "; the old engine skipped the line with a warning; dropped");
      continue;
    }
    Add(r, "reaction", reaction);
    for (const auto &k : pairKeys) Add(r, k.first, k.second);

    // quantity, unit, columns
    std::string quantity, unit, angleRole, quantityNote;
    std::vector<std::pair<std::string, std::string>> extra;  // quantity-specific fields
    switch (code) {
      case 0: quantity = "sigma"; unit = "b"; break;
      case 1: quantity = "dsigma_domega"; unit = "b/sr"; angleRole = "theta_lab"; break;
      case 2: quantity = "phase_shift"; unit = "deg"; break;
      case 3: quantity = "sigma"; unit = "b"; break;
      case 4: quantity = "dsigma_domega"; unit = "b/sr"; angleRole = "theta_cm"; break;
      case 5: quantity = "sigma"; unit = "b"; extra.push_back({"multipole", "E1"}); break;
      case 6: quantity = "sigma"; unit = "b"; extra.push_back({"multipole", "E2"}); break;
      case 7: quantity = "analyzing_power"; unit = "1"; angleRole = "theta_cm"; break;
      case 9:
        quantity = "legendre";
        unit = s.absolute() == 1 ? "b/sr" : "1";
        extra.push_back({"expansion", s.absolute() == 1 ? "absolute" : "a_over_a0"});
        break;
      default:
        res_.refusal = "data line " + std::to_string(d.key) + " has isDiff " + std::to_string(s.isDiff()) + ", which the old engine did not read";
        return;
    }
    if (s.total() == 1) {
      Count(Construct::kTotal);
      quantity = "sigma_total";
    }
    const bool sampleThickness = s.transmission() > 0 || s.yieldThickness() > 0;
    if (sampleThickness) {
      quantity = s.transmission() > 0 ? "transmission" : "neutron_yield";
      unit = "1";
    }
    if (ratio) unit = "1";
    if (thickTarget) {
      if (code == 0 || code == 3) {
        quantity = "yield";
        unit = Implicit(Construct::kThick, name, "unit", "1/particle",
                        "a thick target (density 0): the old engine compared the numbers as reactions per incident particle (yield in 1/particle)");
      } else if (code == 1) {
        quantity = "yield";
        unit = Implicit(Construct::kThick, name, "unit", "1/particle/sr",
                        "a thick target (density 0): the old engine compared the numbers as reactions per incident particle and steradian");
      } else {
        quantity = Ask(Construct::kThick, name + ".quantity", "a thick target on isDiff " + std::to_string(code) + " data: format 3 has no such yield", "");
      }
    } else if (thinTarget && (code == 0 || code == 3 || code == 5 || code == 6) && !sampleThickness) {
      quantityNote = "a thin target: the old engine compared the numbers as the target-averaged cross section in b";
    }
    std::string quantityComment;
    if (quantity == "sigma") {
      int terms = 0;
      const int charged = ChargedElasticTerms(s.entranceKey(), s.exitKey(), s.components(), terms);
      if (charged == terms) {
        quantity = Implicit(Construct::kDataCode, name, "quantity", "sigma_nuclear", kChargedElastic);
        quantityComment = TakeComment();
      } else if (charged > 0) {
        quantity = Ask(Construct::kDataCode, name + ".quantity",
                       "a sum or ratio of an elastic integral of a charged pair and other terms: no quantity holds it", "");
      }
    }
    if (beta) {
      // A beta-delayed spectrum: the old engine compared the numbers in the line's normalization, any units.
      quantity = "beta_delayed_spectrum";
      unit = "arbitrary";
    }
    const bool profiled = thm && s.varyNorm() == 1;  // ESegment::IsProfiledNorm
    if (profiled) unit = "arbitrary";
    Add(r, "quantity", quantity, quantityComment);
    if (!quantityNote.empty()) unit = Implicit(Construct::kDataCode, name, "unit", unit, quantityNote);
    Add(r, "unit", unit, TakeComment());
    for (const auto &x : extra) Add(r, x.first, x.second);
    if (code == 2) {
      Add(r, "J", SpinText(s.phaseJ()));
      Add(r, "l", std::to_string(s.phaseL()));
      std::string spin;
      if (s.phaseS() >= 0) {
        Count(Construct::kSpin);
        spin = SpinText(s.phaseS());
      } else {
        spin = PhaseSpin(name, s.entranceKey(), s.phaseJ(), s.phaseL());
        if (spin.empty()) return;
      }
      Add(r, "s", spin);
    }
    if (s.isUPOS() == 1) {
      Count(Construct::kUpos);
      Add(r, "secondary_gamma", "L " + std::to_string(s.secondaryDecayL()) + " Ic " + SpinText(s.Ic()) + " delta " + NumberText(s.delta(), d.text));
    }

    // file and columns
    Count(Construct::kFile);
    Add(r, "file", std::vector<io::Token>{Tok(s.dataFile())});
    const std::string path = (!s.dataFile().empty() && (s.dataFile()[0] == '/' || opt_.projectDir.empty())) ? s.dataFile() : opt_.projectDir + "/" + s.dataFile();
    int columns = DataColumns(path);
    if (columns < 4) {
      if (active) Note(Construct::kFile, name, "the data file " + s.dataFile() + " cannot be read: four columns assumed");
      columns = 4;
    }
    std::vector<std::string> roles = {beta ? "E_decay[MeV]" : "E_lab[MeV]",
                                      code == 9 ? "order" : angleRole.empty() ? "_" : angleRole + "[deg]", "value", "error"};
    for (int c = 4; c < columns; c++) roles.push_back("_");
    if (app && app->profileWindows) {
      if (columns >= 6) {
        roles[4] = "window_lo";
        roles[5] = "window_hi";
      } else {
        Note(Construct::kProfile, name, "a beam-profile effect, but the data file has no window columns 5-6");
      }
    }
    std::string rolesText;
    for (const std::string &x : roles) rolesText += (rolesText.empty() ? "" : " ") + x;
    Add(r, "columns", rolesText);

    // windows
    Count(Construct::kEnergyWindow);
    Add(r, "energy_window", NumberText(s.minE(), d.text) + " " + NumberText(s.maxE(), d.text) + " MeV");
    if (!angleRole.empty()) {
      Count(Construct::kAngleWindow);
      Add(r, "angle_window", NumberText(s.minA(), d.text) + " " + NumberText(s.maxA(), d.text) + " deg");
    }

    // norm (ESegment: the prior is about the nominal norm, the prior_centre row's value when there is one)
    Count(Construct::kNorm);
    if (!s.normGroup().empty() && normName.count(s.normGroup())) {
      if (normLeader[s.normGroup()] != d.key && (s.dataNorm() != activeData_[normLeader[s.normGroup()]]->parsed->dataNorm() ||
                                                  s.varyNorm() != activeData_[normLeader[s.normGroup()]]->parsed->varyNorm()))
        Note(Construct::kNormGroup, name, "its own norm " + NumberText(s.dataNorm(), d.text) + " (vary " + std::to_string(s.varyNorm()) +
                                             ") was not used by the old engine (the group's is its first line's): dropped");
      Add(r, "norm", normName[s.normGroup()]);
      // the group's parameter is its first active member's; a later member keeps its own, fixed
      if (active && s.varyNorm() == 1)
        Param("segment_" + std::to_string(d.key) + "_norm",
              normName[s.normGroup()] + ".value" + (normLeader[s.normGroup()] == d.key ? "" : "[" + name + "]"));
    } else if (profiled) {
      Add(r, "norm", "profiled");
    } else {
      const std::string n = NumberText(s.dataNorm(), d.text);
      const std::string centre = old_.normCentreText.count(d.key) ? old_.normCentreText.at(d.key) : n;
      Add(r, "norm", Nuisance(n, "", s.varyNorm() == 1, s.dataNormError() != 0, centre, Abs(NumberText(s.dataNormError(), d.text)) + "%"));
      if (active && s.varyNorm() == 1) Param("segment_" + std::to_string(d.key) + "_norm", name + ".norm");
    }

    // energy shift: offset, sqrt, function (summed by the engine in this order)
    const std::string sk = "segment_" + std::to_string(d.key);
    if (!s.shiftGroup().empty() && scaleName.count(s.shiftGroup())) {
      Count(Construct::kShift);
      Add(r, "shift", scaleName[s.shiftGroup()]);
      const OldDataLine &lead = *activeData_[shiftLeader[s.shiftGroup()]];
      if (lead.key != d.key && (s.energyShift() != 0 || s.varyEnergyShift() == 1))
        Note(Construct::kShiftGroup, name, "its own shift " + NumberText(s.energyShift(), d.text) +
                                               " was not used by the old engine (the group's is its first line's): dropped");
      if (lead.key != d.key && !s.shiftFn().empty())
        Note(Construct::kShiftGroup, name, "a shiftfn on a follower of a shift group: the old engine refused it");
      if (active) Param(sk + "_energy_shift", scaleName[s.shiftGroup()] + ".offset" + (lead.key == d.key ? "" : "[" + name + "]"));
    } else {
      if (s.energyShift() != 0 || s.varyEnergyShift() == 1) {
        Count(Construct::kShift);
        std::vector<io::Token> t = {Tok("offset")};
        for (const io::Token &x : offsetNuisance(d)) t.push_back(x);
        Add(r, "shift", t);
      }
      if (active) Param(sk + "_energy_shift", name + ".shift.offset");  // every data set has it (fixed 0 without a shift)
    }
    if (s.hasSqrtShift() && (s.energyShiftSqrt() != 0 || s.varyEnergyShiftSqrt() == 1)) {
      Count(Construct::kSqrtShift);
      std::vector<io::Token> t = {Tok("sqrt")};
      for (const io::Token &x : Nuisance(NumberText(s.energyShiftSqrt(), d.text), "MeV^1/2", s.varyEnergyShiftSqrt() == 1,
                                         s.energyShiftSqrtError() != 0, "0", Abs(NumberText(s.energyShiftSqrtError(), d.text))))
        t.push_back(x);
      Add(r, "shift", t);
    }
    if (active) Param(sk + "_energy_shift_sqrt", name + ".shift.sqrt");  // every data set has it (fixed 0 without one)
    if (!s.shiftFnGroup().empty() && fnName.count(s.shiftFnGroup())) {
      Add(r, "shift", fnName[s.shiftFnGroup()]);
    } else if (!s.shiftFn().empty() && s.shiftGroup().empty()) {
      Count(Construct::kShiftFn);
      if (active) {
        io::Record g;
        g.kind = "energy_scale";
        g.name = Unique(name + "_scale");
        addFunction(g, d, g.name);
        groupRecords_.push_back(g);
        Add(r, "shift", g.name);
      }
    } else if (!s.shiftFn().empty()) {
      Count(Construct::kShiftFn);
    }

    if (!s.covarianceFile().empty()) {
      Count(Construct::kCovariance);
      Add(r, "covariance", std::vector<io::Token>{Tok(s.covarianceFile())});
    }

    // apparatus
    if (app) {
      if (thm && app->line->line.isTargetIntegration == 1) {
        res_.refusal = "data line " + std::to_string(d.key) + " is a THM segment with target integration, which the old engine refused";
        return;
      }
      if (app->hasTarget) Add(r, "target", app->targetName);
      if (app->hasBeam) Add(r, "beam", app->beamName);
      if (app->hasDetector) Add(r, "detector", app->detectorName);
      if (app->folded) Add(r, "quadrature", app->quadrature);
      if (!app->foldWindows.empty()) Add(r, "fold_windows", app->foldWindows);
    }
    if (sampleTargetOfData_.count(d.key) && !(app && app->hasTarget)) Add(r, "target", sampleTargetOfData_[d.key]);
    if (sampleThickness && active) {
      if ((app && app->hasTarget) || sampleTargetOfData_.count(d.key)) {
        Ask(Construct::kThickness, name + ".target", "a sample thickness and a <targetInt> or <abundances> target: two areal densities", "");
      } else {
        const double n = s.transmission() > 0 ? s.transmission() : s.yieldThickness();
        Add(r, "target", SampleThicknessTarget(name, s.entranceKey(), n, d.text, s.thicknessFit(), s.thicknessVary(),
                                               s.thicknessSigma(), s.thicknessGroup(), d.key, true));
      }
    }

    if (thm) {
      const auto e = old_.thm.experimentOfSegment.find(d.key);
      Add(r, "thm", e != old_.thm.experimentOfSegment.end() ? e->second : "global");
      const auto w = old_.thm.weightOfData.find(d.key);
      if (w != old_.thm.weightOfData.end()) Add(r, "thm_weight", std::vector<io::Token>{Tok(w->second)});
    }
    if (!active) {
      Count(Construct::kInactive);
      Add(r, "active", "no");
    }
    out_.data.records.push_back(r);
  }
}

// ================================================================ extrapolations
void Builder::BuildExtrapolations() {
  for (const OldTestLine &t : old_.tests) {
    Count(Construct::kTestLine);
    const ExtrapLine &s = *t.parsed;
    const std::string &name = testName_[t.key];
    const bool active = activeTest_.count(t.key) != 0;
    io::Record r;
    r.kind = "extrapolation";
    r.name = name;
    r.comments.push_back("# <segmentsTest> line " + std::to_string(t.key));
    const bool thm = s.isDiff() >= 10;
    const int code = thm ? s.isDiff() - 10 : s.isDiff();
    Count(Construct::kTestCode);
    if (thm) Count(Construct::kThm);
    const OldPair *in = Pair(s.entranceKey());
    const bool ratio = s.isAdvanced() == 1 && s.operationType() == 1 && !s.components().empty();
    // As for a data line: exit key -1 needs the total-capture type (4 on a test line), or the old engine refused it.
    if (active && in && s.exitKey() == -1 && code != 4) {
      res_.refusal = "test line " + std::to_string(t.key) +
                     ": the exit key -1 (total capture) needs the total-capture data type (4 on a test line); the old engine refused it";
      return;
    }
    const std::string angleFrame = (code == 5 || code == 7) ? "cm" : "lab";
    std::vector<io::Token> reaction;
    std::vector<std::pair<std::string, std::string>> pairKeys;
    std::string why;
    if (!ReactionTokens(name, s.entranceKey(), (code == 4 || s.exitKey() == -1) ? -1 : s.exitKey(), s.isAdvanced(), s.operationType(), s.components(),
                        angleFrame, reaction, why, pairKeys)) {
      Note(Construct::kTestLine, name, (active ? "" : "inactive ") + std::string("test line ") + std::to_string(t.key) + ": " + why + "; dropped");
      continue;
    }
    Add(r, "reaction", reaction);
    for (const auto &k : pairKeys) Add(r, k.first, k.second);
    std::string quantity, unit;
    bool angles = false;
    switch (code) {
      case 0: quantity = "sigma"; unit = "b"; break;
      case 1: quantity = "dsigma_domega"; unit = "b/sr"; angles = true; break;
      case 2: quantity = "phase_shift"; unit = "deg"; break;
      case 3: quantity = "legendre"; unit = "1"; break;
      case 4: quantity = "sigma"; unit = "b"; break;
      case 5: quantity = "dsigma_domega"; unit = "b/sr"; angles = true; break;
      case 7: quantity = "analyzing_power"; unit = "1"; angles = true; break;
      default:
        res_.refusal = "test line " + std::to_string(t.key) + " has isDiff " + std::to_string(s.isDiff()) + ", which the old engine did not compute";
        return;
    }
    if (s.total() == 1) {
      Count(Construct::kTotal);
      quantity = "sigma_total";
    }
    if (s.transmission() > 0 || s.yieldThickness() > 0) {
      quantity = s.transmission() > 0 ? "transmission" : "neutron_yield";
      unit = "1";
    }
    if (ratio) unit = "1";
    const Apparatus *app = active && apparatusOfTest_.count(t.key) ? &apparatus_[apparatusOfTest_[t.key]] : nullptr;
    const bool thickTarget = app && app->hasTarget && app->line->line.isTargetIntegration == 1 && !(app->line->line.density > 0.0);
    if (thickTarget && (code == 0 || code == 4)) {
      quantity = "yield";
      unit = "1/particle";
    } else if (thickTarget && code == 1) {
      quantity = "yield";
      unit = "1/particle/sr";
    }
    if (in && in->pType == 20) {
      quantity = "beta_delayed_spectrum";
      unit = "arbitrary";
    }
    std::string quantityComment;
    if (quantity == "sigma") {
      int terms = 0;
      const int charged = ChargedElasticTerms(s.entranceKey(), s.exitKey(), s.components(), terms);
      if (charged == terms) {
        quantity = Implicit(Construct::kTestCode, name, "quantity", "sigma_nuclear", kChargedElastic);
        quantityComment = TakeComment();
      } else if (charged > 0) {
        quantity = Ask(Construct::kTestCode, name + ".quantity",
                       "a sum or ratio of an elastic integral of a charged pair and other terms: no quantity holds it", "");
      }
    }
    Add(r, "quantity", quantity, quantityComment);
    Add(r, "unit", unit);
    if (code == 3) {
      Add(r, "expansion", "a_over_a0");
      Add(r, "orders", "0 " + std::to_string(s.maxAngDistOrder()));
    }
    if (code == 2) {
      Add(r, "J", SpinText(s.phaseJ()));
      Add(r, "l", std::to_string(s.phaseL()));
      std::string spin;
      if (s.phaseS() >= 0) {
        Count(Construct::kSpin);
        spin = SpinText(s.phaseS());
      } else {
        spin = PhaseSpin(name, s.entranceKey(), s.phaseJ(), s.phaseL());
        if (spin.empty()) return;
      }
      Add(r, "s", spin);
    }
    // the energy grid: the old loop accumulates energy += step from Emin while energy <= Emax (EData::MakePoints)
    // (one pass when the step is 0)
    const std::string lo = NumberText(s.minE(), t.text), hi = NumberText(s.maxE(), t.text), step = NumberText(s.eStep(), t.text);
    std::vector<double> energyValues;
    for (double e = s.minE(); e <= s.maxE(); e += s.eStep()) {
      energyValues.push_back(e);
      if (s.eStep() == 0.0 || energyValues.size() > 1000000) break;
    }
    if (opt_.exactGrids || energyValues.size() <= 1 || !(s.eStep() > 0)) {
      std::string list = "list";
      for (double e : energyValues) list += " " + ShortestText(e);
      if (energyValues.empty()) {
        list += " " + lo;
        Note(Construct::kTestLine, name, "the energy range " + lo + " to " + hi + " holds no point (the old engine computed none); one point written");
      }
      Add(r, "energies", list + " MeV lab");
    } else {
      Add(r, "energies", "range " + lo + " " + hi + " MeV lab step " + step);
    }
    // the angle grid: an outer loop angle += aStep from Amin while angle <= Amax, one pass when aStep is 0
    std::vector<double> angleValues;
    for (double a = s.minA(); a <= s.maxA(); a += s.aStep()) {
      angleValues.push_back(a);
      if (s.aStep() == 0.0 || angleValues.size() > 1000000) break;
    }
    if (angles) {
      if (opt_.exactGrids || angleValues.size() <= 1 || !(s.aStep() > 0)) {
        std::string list = "list";
        for (double a : angleValues) list += " " + ShortestText(a);
        if (angleValues.empty()) list += " " + NumberText(s.minA(), t.text);
        Add(r, "angles", list + " deg " + angleFrame);
      } else {
        Add(r, "angles", "range " + NumberText(s.minA(), t.text) + " " + NumberText(s.maxA(), t.text) + " deg " + angleFrame +
                             " step " + NumberText(s.aStep(), t.text));
      }
    } else if (angleValues.size() > 1) {
      Note(Construct::kTestLine, name, "an angle-integrated quantity on " + std::to_string(angleValues.size()) +
                                           " angles: the old engine repeated the energy grid once per angle; written once");
    }
    if (app && Implicit(Construct::kEffectLine, name, "apparatus", "yes",
                        "the <targetInt> line that lists line " + std::to_string(t.key) +
                            " also folded this extrapolation (the same index, EData::MakePoints); written explicitly") != "yes")
      app = nullptr;
    if (app) {
      const std::string comment = TakeComment();
      if (app->hasTarget) Add(r, "target", app->targetName, comment);
      if (app->hasBeam) Add(r, "beam", app->beamName, app->hasTarget ? "" : comment);
      if (app->hasDetector) Add(r, "detector", app->detectorName, app->hasTarget || app->hasBeam ? "" : comment);
      if (app->folded) Add(r, "quadrature", app->quadrature);
      if (!app->foldWindows.empty()) Add(r, "fold_windows", app->foldWindows);
    }
    pendingComment_.clear();
    if (sampleTargetOfTest_.count(t.key) && !(app && app->hasTarget)) Add(r, "target", sampleTargetOfTest_[t.key]);
    if ((s.transmission() > 0 || s.yieldThickness() > 0) && active) {
      if ((app && app->hasTarget) || sampleTargetOfTest_.count(t.key)) {
        Ask(Construct::kThickness, name + ".target", "a sample thickness and a <targetInt> or <abundances> target: two areal densities", "");
      } else {
        if (s.thicknessFitIgnored())
          Note(Construct::kThickness, name, "the fit form or thicknessgroup of a test line was ignored by the old engine: the thickness is fixed");
        Add(r, "target", SampleThicknessTarget(name, s.entranceKey(), s.transmission() > 0 ? s.transmission() : s.yieldThickness(),
                                               t.text, false, 0, 0.0, "", t.key, false));
      }
    }
    if (thm) {
      Add(r, "thm", "global");
      const auto w = old_.thm.weightOfTest.find(t.key);
      if (w != old_.thm.weightOfTest.end()) Add(r, "thm_weight", std::vector<io::Token>{Tok(w->second)});
    }
    if (!active) {
      Count(Construct::kInactive);
      Add(r, "active", "no");
    }
    out_.extrapolations.records.push_back(r);
  }
}

}  // namespace convert
}  // namespace azure2
