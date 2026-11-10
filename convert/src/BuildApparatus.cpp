// The <targetInt> lines as targets, beams and detectors (F7), the samples of <abundances>, and the names of the
// shared records (Builder.h).
#include "Builder.h"

namespace azure2 {
namespace convert {

// ================================================================ <targetInt>
void Builder::BuildApparatus() {
  for (const OldEffectLine &e : old_.effects) {
    Count(Construct::kEffectLine);
    if (e.activeIndex == 0) {
      Note(Construct::kEffectLine, "", "<targetInt> line " + std::to_string(e.index) + " is inactive: dropped (the old engine did not read it)");
      continue;
    }
    const TargetIntLine &l = e.line;
    Apparatus a;
    a.line = &e;
    const std::string n = std::to_string(e.activeIndex);
    // The data and test lines it applies to (EData::ReadTargetEffectsFile: a listed key that is an active segment).
    for (int k : e.segments) {
      const bool data = activeData_.count(k) != 0, test = activeTest_.count(k) != 0;
      if (!data && !test) {
        Note(Construct::kNotASegment, "", "<targetInt> line " + std::to_string(e.index) + " lists " + std::to_string(k) +
                                              ", which is no active data or test line: dropped (the old engine skipped it)");
        Count(Construct::kNotASegment);
      }
      if (data) a.dataKeys.push_back(k);
      if (test) a.testKeys.push_back(k);
    }
    // target: integration, screening, Doppler
    const bool integrates = l.isTargetIntegration == 1;
    if (integrates || l.hasScreening || l.isDoppler) {
      a.hasTarget = true;
      a.target.kind = "target";
      // The active nucleus: the target nucleus of the first data (or test) line the effect applies to.
      std::string active;
      int firstKey = !a.dataKeys.empty() ? a.dataKeys[0] : !a.testKeys.empty() ? a.testKeys[0] : 0;
      int entrance = 0;
      if (firstKey && activeData_.count(firstKey)) entrance = activeData_[firstKey]->parsed->entranceKey();
      else if (firstKey && activeTest_.count(firstKey)) entrance = activeTest_[firstKey]->parsed->entranceKey();
      if (const OldPair *p = Pair(entrance)) active = Nuclide(p->z2, p->A2());
      if (integrates) {
        const std::string text = e.text;
        const bool thick = !(l.density > 0.0);
        Add(a.target, "kind", thick ? "infinite" : "layer");
        if (thick) {
          Count(Construct::kThick);
        }
        if (active.empty()) active = Ask(Construct::kEffectLine, "line" + n + ".active", "the target nucleus of <targetInt> line " + std::to_string(e.index) + " (it applies to no line)", "");
        Add(a.target, "active", active);
        if (!thick) {
          const std::string density = NumberText(l.density, text);
          if (l.densityFitRead) {
            Count(Construct::kDensityFit);
            const std::string &tok = l.densitySigmaToken;
            const bool percent = !tok.empty() && tok.back() == '%';
            const double s = std::strtod(percent ? tok.substr(0, tok.size() - 1).c_str() : tok.c_str(), nullptr);
            std::string g = density + " atoms/cm2 " + (l.densityVary ? "fit" : "fixed");
            if (s > 0) g += " prior " + density + " " + tok;
            if (l.densityMaxRead > 0.0) g += " grid_max " + NumberText(l.densityMaxRead, text);
            Add(a.target, "areal_density", g);
            Param("density_" + n, "@target.areal_density");
          } else {
            Add(a.target, "areal_density", density + " atoms/cm2 fixed");
          }
        }
        std::string stopping = "formula";
        std::vector<io::Token> st = {Tok("formula"), Quoted(TargetIntUnquoted(l.stoppingPowerEq))};
        for (int i = 0; i < (int)l.parameters.size(); i++)
          st.push_back(Tok("a" + std::to_string(i) + "=" + NumberText(l.parameters[i], text)));
        st.push_back(Tok("MeV*cm2/atom"));
        Add(a.target, "stopping", st);
        Count(Construct::kStraggling);
        if (l.isStraggling != 1)
          Add(a.target, "straggling", "none");
        else if (l.stragglingCoefficient == 0.0)
          Add(a.target, "straggling", "bohr");
        else
          Add(a.target, "straggling", "coefficient " + NumberText(l.stragglingCoefficient, text) + " keV^1/2 cm");
      } else {
        Add(a.target, "kind", "none");
      }
      if (l.isDoppler) {
        Count(Construct::kDoppler);
        const std::string T = NumberText(l.dopplerTemperature, e.text);
        std::string d = "free_gas " + T + " K";
        if (l.dopplerFitRead) {
          const double sigma = std::strtod(l.dopplerSigmaToken.c_str(), nullptr);
          d += std::strtod(l.dopplerVaryToken.c_str(), nullptr) == 1.0 ? " fit" : " fixed";
          if (sigma > 0) d += " prior " + T + " " + l.dopplerSigmaToken;
          if (!l.dopplerMaxToken.empty()) d += " grid_max " + l.dopplerMaxToken;
          Param("doppler_" + n, "@target.doppler");
        } else {
          d += " fixed";
        }
        Add(a.target, "doppler", d);
      }
      if (l.hasScreening) {
        Count(Construct::kScreening);
        const std::string ue = NumberText(l.screeningNominal, e.text);
        Add(a.target, "screening",
            Nuisance(ue, "eV", l.screeningVary != 0, l.screeningSigma > 0, ue, NumberText(l.screeningSigma, e.text)));
        Param("screening_" + n, "@target.screening");
      }
      if (!l.densityGroup.empty()) Count(Construct::kDensityGroup);
    }
    // beam: a Gaussian spread or a beam profile
    const bool constant = l.isConvolution == 1, formula = l.isConvCoefficients == 1;
    if (constant || formula || l.isBeamProfile) {
      a.hasBeam = true;
      a.beam.kind = "beam";
      if (constant && formula)
        Note(Construct::kSpread, "", "<targetInt> line " + std::to_string(e.index) + " has a constant sigma and an expression; the engine uses the expression");
      const std::string spreadNote = "a Gaussian convolution is written as the beam's energy spread (the kernel is the same; "
                                     "the old line does not say whether it was a detector resolution)";
      if (formula) {
        Count(Construct::kSpread);
        std::vector<io::Token> sp = {Tok("formula"), Quoted(TargetIntUnquoted(l.convolutionEq))};
        for (int i = 0; i < (int)l.convCoefficients.size(); i++)
          sp.push_back(Tok("a" + std::to_string(i) + "=" + NumberText(l.convCoefficients[i], e.text)));
        sp.push_back(Tok("MeV"));
        sp.push_back(Tok(l.sigmaLab ? "lab" : "cm"));
        std::string text;
        for (const io::Token &t : sp) text += (text.empty() ? "" : " ") + (t.quoted ? "\"" + t.text + "\"" : t.text);
        if (Implicit(Construct::kSpread, "targetInt" + std::to_string(e.index), "spread", text, spreadNote) != text) sp = Words("?");
        Add(a.beam, "spread", sp, TakeComment());
      } else if (constant) {
        Count(Construct::kSpread);
        const std::string spread = Implicit(Construct::kSpread, "targetInt" + std::to_string(e.index), "spread",
                                            NumberText(l.sigma, e.text) + " MeV lab", spreadNote);
        Add(a.beam, "spread", spread, TakeComment());
      }
      if (l.isBeamProfile) {
        Count(Construct::kProfile);
        for (const TargetIntLine::ProfileComponent &c : l.beamProfile)
          Add(a.beam, "profile", "skewed_gaussian xi " + NumberText(c.xi, e.text) + " MeV omega " + NumberText(c.omega, e.text) +
                                     " MeV alpha " + NumberText(c.alpha, e.text) + " weight " + NumberText(c.weight, e.text));
        Add(a.beam, "truncation", NumberText(l.beamTruncation, e.text));
        Add(a.beam, "detailed_balance", l.beamDetailedBalance == 1 ? "yes" : "no");
        a.profileWindows = true;
      }
    }
    // detector: Q coefficients, the beam profile's resolution, UDR
    if (l.isQCoefficients == 1 || l.isBeamProfile || l.isUdr) {
      a.hasDetector = true;
      a.detector.kind = "detector";
      if (l.isQCoefficients == 1) {
        Count(Construct::kQ);
        std::string q;
        for (double v : l.qCoefficients) q += (q.empty() ? "" : " ") + NumberText(v, e.text);
        if (q.empty()) q = "?";
        Add(a.detector, "attenuation", q);
      }
      if (l.isBeamProfile) Add(a.detector, "energy_resolution", NumberText(l.beamTpcSigma, e.text) + " MeV lab");
      if (l.isUdr) {
        Count(Construct::kUdr);
        std::vector<io::Token> u = {Tok("file"), Tok(l.udrFile), Tok("flight_path"), Tok(NumberText(l.udrFlightPath, e.text)),
                                    Tok("m"), Tok("burst"), Tok(NumberText(l.udrBurstFwhm, e.text)), Tok("ns"), Tok("channel"),
                                    Tok(NumberText(l.udrChannelWidth, e.text)), Tok("ns"), Tok("centred"),
                                    Tok(l.udrCentred != 0 ? "yes" : "no")};
        Add(a.detector, "udr", u);
      }
    }
    // the fold's numerical settings, on the data sets
    a.folded = a.hasBeam || integrates || l.isDoppler || l.isBeamProfile || l.isUdr;
    if (a.folded) {
      const std::string rwm = l.hasResonanceWidthMultiplier ? NumberText(l.resonanceWidthMultiplier, e.text) : "20";
      const std::string ppw = l.hasPointsPerWidth ? NumberText(l.pointsPerWidth, e.text) : "50";
      a.quadrature = Words("points " + std::to_string(l.numIntegrationPoints) + " width_multiplier " + rwm + " per_width " + ppw);
      if (l.hasRangesToken && (l.autoWindows || !l.ranges.empty())) {
        Count(Construct::kWindows);
        std::string w;
        if (l.autoWindows) {
          w = "auto";
          if (l.hasAutoWindowTolerance) w += ":" + ShortestText(l.autoWindowTolerance);
          if (l.hasAutoWindowMargin) w += ":" + ShortestText(l.autoWindowMargin);
          if (!l.ranges.empty())
            Note(Construct::kWindows, "", "<targetInt> line " + std::to_string(e.index) + ": automatic windows and explicit ranges; format 3 writes 'auto' (the ranges are dropped)");
        } else {
          for (const auto &r : l.ranges) w += (w.empty() ? "" : " ") + ShortestText(r.first) + "-" + ShortestText(r.second);
          w += " MeV";
        }
        if (l.transitionWidth > 0) w += " blend " + ShortestText(l.transitionWidth) + " MeV";
        if (l.autoTolerance > 0) w += " auto_tolerance " + ShortestText(l.autoTolerance);
        a.foldWindows = Words(w);
      }
    } else if (l.hasRangesToken && (l.autoWindows || !l.ranges.empty())) {
      Note(Construct::kWindows, "", "<targetInt> line " + std::to_string(e.index) + ": energy windows on an effect that folds nothing: dropped");
    }
    if (!a.hasTarget && !a.hasBeam && !a.hasDetector) {
      Note(Construct::kEffectLine, "", "<targetInt> line " + std::to_string(e.index) + " applies no effect: dropped");
      continue;
    }
    // Later lines take over a key (EData::ReadTargetEffectsFile overwrites the segment's effect number).
    const int index = (int)apparatus_.size();
    for (int k : a.dataKeys) {
      if (apparatusOfData_.count(k))
        Note(Construct::kEffectLine, dataName_[k], "data line " + std::to_string(k) + " is listed by two <targetInt> lines; the later one applies (as in the old engine)");
      apparatusOfData_[k] = index;
    }
    for (int k : a.testKeys) apparatusOfTest_[k] = index;
    apparatus_.push_back(a);
  }
}

// ================================================================ <abundances>
void Builder::BuildSamples() {
  for (const OldAbundance &a : old_.abundances) {
    Count(Construct::kAbundances);
    Sample *s = nullptr;
    for (Sample &x : samples_)
      if (x.name == a.line.sample) s = &x;
    if (!s) {
      samples_.push_back(Sample());
      samples_.back().name = a.line.sample;
      s = &samples_.back();
    }
    if (a.line.kind == AbundanceLine::Sum) {
      s->hasSum = true;
      s->sum = a.line.sumSigma;
      s->sumText = NumberText(a.line.sumSigma, a.text);
    } else {
      s->entries[{a.line.z, a.line.a}] = &a;
    }
  }
  if (doc_.Find("abundances") && samples_.empty()) {
    samples_.push_back(Sample());  // an empty block: one unnamed sample at the defaults (CNuc::ReadAbundances)
  }
  // The sample of each line: "1" is the unnamed sample, or the only one (CNuc::DefaultSample).
  auto sampleIndex = [&](const std::string &token) -> int {
    if (token == "0") return -1;
    if (token == "1") {
      for (int i = 0; i < (int)samples_.size(); i++)
        if (samples_[i].name.empty()) return i;
      return samples_.size() == 1 ? 0 : -2;
    }
    for (int i = 0; i < (int)samples_.size(); i++)
      if (samples_[i].name == token) return i;
    return -2;
  };
  // data and test lines alike: a sample <abundances> lacks is dropped with a note naming the line
  auto sampleOfLine = [&](const std::string &record, const std::string &natural) {
    const int i = sampleIndex(natural);
    if (i == -2)
      Note(Construct::kSample, record, "sample '" + natural + "' is not in <abundances>: the old engine refused the line; the record has no isotopes");
    return i;
  };
  for (const auto &kv : activeData_) {
    const int i = sampleOfLine(dataName_[kv.first], kv.second->parsed->natural());
    if (i >= 0) sampleOfData_[kv.first] = i;
  }
  for (const auto &kv : activeTest_) {
    const int i = sampleOfLine(testName_[kv.first], kv.second->parsed->natural());
    if (i >= 0) sampleOfTest_[kv.first] = i;
  }
}

// ================================================================ names and isotopes of the apparatus
void Builder::NameApparatus() {
  // Samples: the isotopes of each sample, on the target of its data sets.
  for (int si = 0; si < (int)samples_.size(); si++) {
    const Sample &smp = samples_[si];
    std::map<std::pair<int, int>, std::string> nucleusOfCompound;  // compound (Z, A) -> target nucleus label
    std::vector<int> users;
    auto addTerms = [&](int entrance, const std::vector<ComponentEntry> &components) {
      std::vector<int> keys = {entrance};
      for (const ComponentEntry &c : components) keys.push_back(c.entrance);
      for (int k : keys)
        if (const OldPair *p = Pair(k)) {
          const std::pair<int, int> compound = {p->z1 + p->z2, p->A1() + p->A2()};
          const std::string nucleus = Nuclide(p->z2, p->A2());
          auto it = nucleusOfCompound.find(compound);
          if (it == nucleusOfCompound.end()) nucleusOfCompound[compound] = nucleus;
          else if (it->second != nucleus)
            it->second = Ask(Construct::kSample, (smp.name.empty() ? "natural" : smp.name) + ".isotope." + std::to_string(compound.first) + "_" + std::to_string(compound.second),
                             "the compound nucleus (" + std::to_string(compound.first) + ", " + std::to_string(compound.second) +
                                 ") of the sample is reached from " + it->second + " and from " + nucleus, it->second);
        }
    };
    for (const auto &kv : sampleOfData_)
      if (kv.second == si) {
        users.push_back(kv.first);
        addTerms(activeData_[kv.first]->parsed->entranceKey(), activeData_[kv.first]->parsed->components());
      }
    for (const auto &kv : sampleOfTest_)
      if (kv.second == si) addTerms(activeTest_[kv.first]->parsed->entranceKey(), activeTest_[kv.first]->parsed->components());
    if (users.empty() && nucleusOfCompound.empty()) {
      Note(Construct::kAbundances, "", "sample '" + smp.name + "' is used by no active line: its abundances are dropped");
      continue;
    }
    Count(Construct::kSample);
    // isotope fields, in compound-nucleus order
    std::vector<io::Field> isotopes;
    std::vector<std::pair<std::string, std::string>> params;  // old name, isotope label
    for (const auto &kv : nucleusOfCompound) {
      const auto e = smp.entries.find(kv.first);
      io::Field f;
      f.key = "isotope";
      f.tokens.push_back(Tok(kv.second));
      std::string v = "1 fixed";
      if (e != smp.entries.end()) {
        const AbundanceLine &a = e->second->line;
        const std::string value = NumberText(a.value, e->second->text);
        for (const io::Token &t : Nuisance(value, "", a.vary == 1, a.sigma > 0, value, NumberText(a.sigma, e->second->text))) f.tokens.push_back(t);
      } else {
        for (const io::Token &t : Words(v)) f.tokens.push_back(t);
      }
      isotopes.push_back(f);
      params.push_back({"abundance_" + (smp.name.empty() ? "" : smp.name + "_") + std::to_string(kv.first.first) + "_" +
                            std::to_string(kv.first.second),
                        kv.second});
    }
    for (const auto &kv : smp.entries)
      if (!nucleusOfCompound.count(kv.first))
        Note(Construct::kAbundances, "", "sample '" + smp.name + "': compound nucleus (" + std::to_string(kv.first.first) + ", " +
                                             std::to_string(kv.first.second) + ") is reached by no data set of the sample: dropped");
    io::Field sum;
    sum.key = "isotope_sum";
    const double sigma = smp.hasSum && smp.sum >= 0 ? smp.sum : old_.abundanceSumSigma;
    const std::string sigmaText = smp.hasSum && smp.sum >= 0 ? smp.sumText : old_.abundanceSumSigmaText;
    sum.tokens = sigma > 0 ? Words("prior 1 " + sigmaText) : Words("none");
    // Where the isotopes go: the data set's <targetInt> target, else one target per sample.
    // Data sets and extrapolations alike (one path): a member on a <targetInt> target puts the isotopes there, and a
    // target shared with a line of another sample (or none) is a question.
    std::set<int> targetApparatus;
    bool ownTarget = false;
    struct Member {
      bool isData;
      int key;
    };
    std::vector<Member> members;
    for (int k : users) members.push_back({true, k});
    for (const auto &kv : sampleOfTest_)
      if (kv.second == si) members.push_back({false, kv.first});
    auto sampleOf = [&](bool isData, int key) {
      const std::map<int, int> &samples = isData ? sampleOfData_ : sampleOfTest_;
      return samples.count(key) ? samples.at(key) : -1;
    };
    for (const Member &m : members) {
      std::map<int, int> &apparatusOf = m.isData ? apparatusOfData_ : apparatusOfTest_;
      if (!(apparatusOf.count(m.key) && apparatus_[apparatusOf[m.key]].hasTarget)) {
        ownTarget = true;
        continue;
      }
      const int ai = apparatusOf[m.key];
      targetApparatus.insert(ai);
      const std::string &name = (m.isData ? dataName_ : testName_)[m.key];
      for (const bool otherIsData : {true, false})
        for (int other : otherIsData ? apparatus_[ai].dataKeys : apparatus_[ai].testKeys)
          if (sampleOf(otherIsData, other) != si)
            Ask(Construct::kSample, name + ".target",
                "its target (<targetInt> line " + std::to_string(apparatus_[ai].line->index) + ") is shared with " +
                    (otherIsData ? dataName_ : testName_)[other] + ", which has not the same sample",
                "");
    }
    if (targetApparatus.size() + (ownTarget ? 1 : 0) > 1)
      Ask(Construct::kSample, (smp.name.empty() ? "natural" : smp.name) + ".target",
          "the data sets of the sample have different targets; the old engine shared one set of abundances", "");
    for (int ai : targetApparatus) {
      for (const io::Field &f : isotopes) apparatus_[ai].target.fields.push_back(f);
      apparatus_[ai].target.fields.push_back(sum);
    }
    std::string sampleTarget;
    if (ownTarget) {
      io::Record t;
      t.kind = "target";
      t.name = sampleTarget = Unique(Sanitized(smp.name.empty() ? "natural" : smp.name));
      t.comments.push_back("# <abundances> sample " + (smp.name.empty() ? std::string("(unnamed)") : smp.name));
      Add(t, "kind", "none");
      for (const io::Field &f : isotopes) t.fields.push_back(f);
      t.fields.push_back(sum);
      sampleTargets_.push_back(t);
      for (int k : users)
        if (!(apparatusOfData_.count(k) && apparatus_[apparatusOfData_[k]].hasTarget)) sampleTargetOfData_[k] = sampleTarget;
      for (const auto &kv : sampleOfTest_)
        if (kv.second == si && !(apparatusOfTest_.count(kv.first) && apparatus_[apparatusOfTest_[kv.first]].hasTarget))
          sampleTargetOfTest_[kv.first] = sampleTarget;
    }
    for (const auto &p : params)
      Param(p.first, (ownTarget ? sampleTarget : "@sample" + std::to_string(*targetApparatus.begin())) + ".isotope." + p.second);
  }

  // Density groups: the lines of a group are one target (one areal density, the first line's).
  std::map<std::string, std::vector<int>> groupLines;
  for (int i = 0; i < (int)apparatus_.size(); i++)
    if (!apparatus_[i].line->line.densityGroup.empty() && apparatus_[i].hasTarget &&
        (!apparatus_[i].dataKeys.empty() || !apparatus_[i].testKeys.empty()))
      groupLines[apparatus_[i].line->line.densityGroup].push_back(i);
  std::map<int, std::string> groupTargetName;
  for (const auto &g : groupLines) {
    const Apparatus &lead = apparatus_[g.second[0]];
    if (!lead.line->line.densityFitRead) {
      res_.refusal = "density group " + g.first + ": its first <targetInt> line has no densityfit block, which the old engine refused";
      return;
    }
    const std::string name = Unique(Sanitized(g.first));
    for (int i : g.second) {
      Apparatus &a = apparatus_[i];
      if (i != g.second[0]) {
        // Equal apart from the areal density?  Compare with the leader's fields.
        io::Record x = a.target, y = lead.target;
        auto strip = [](io::Record &r) {
          r.fields.erase(std::remove_if(r.fields.begin(), r.fields.end(), [](const io::Field &f) { return f.key == "areal_density"; }),
                         r.fields.end());
        };
        strip(x);
        strip(y);
        if (Signature(x) != Signature(y))
          Ask(Construct::kDensityGroup, g.first + ".split",
              "the lines of density group " + g.first + " differ in stopping, straggling, Doppler or screening; one target cannot hold them", "");
        if (a.line->line.density != lead.line->line.density)
          Note(Construct::kDensityGroup, name, "<targetInt> line " + std::to_string(a.line->index) + "'s own density " +
                                                   NumberText(a.line->line.density, a.line->text) +
                                                   " was not used by the old engine (the group's is its first line's): dropped");
        a.target = lead.target;
      }
      groupTargetName[i] = name;
    }
  }

  // Merge identical records without free parameters (F7), then name each by its users.
  auto nameFor = [&](const std::string &prefix, const std::vector<int> &dataUsers, const std::vector<int> &testUsers, int lineIndex) {
    std::set<std::string> exps;
    bool allThm = !dataUsers.empty();
    for (int k : dataUsers) {
      const auto e = old_.thm.experimentOfSegment.find(k);
      if (e == old_.thm.experimentOfSegment.end()) allThm = false;
      else exps.insert(e->second);
    }
    if (dataUsers.size() == 1) return Unique(prefix + dataName_[dataUsers[0]]);
    if (allThm && exps.size() == 1) return Unique(prefix + Sanitized(*exps.begin()));
    if (dataUsers.empty() && testUsers.size() == 1) return Unique(prefix + testName_[testUsers[0]]);
    return Unique(prefix + "targetInt" + std::to_string(lineIndex));
  };
  struct Group {
    std::vector<int> members;
    std::vector<int> dataUsers, testUsers;
  };
  auto mergeAndName = [&](bool Apparatus::*has, io::Record Apparatus::*rec, std::string Apparatus::*name, const std::string &prefix,
                          bool allowParameters) {
    std::vector<Group> groups;
    std::map<std::string, int> bySignature;
    for (int i = 0; i < (int)apparatus_.size(); i++) {
      Apparatus &a = apparatus_[i];
      if (!(a.*has)) continue;
      if (a.dataKeys.empty() && a.testKeys.empty()) continue;
      int gi = -1;
      const auto dg = groupTargetName.find(i);
      if (dg != groupTargetName.end() && prefix == "tgt_") {
        a.*name = dg->second;
        continue;
      }
      const std::string sig = Signature(a.*rec);
      if (!HasFreeParameter(a.*rec) || allowParameters) {
        const auto it = bySignature.find(sig);
        if (it != bySignature.end()) gi = it->second;
      }
      if (gi < 0) {
        gi = (int)groups.size();
        groups.push_back(Group());
        if (!HasFreeParameter(a.*rec)) bySignature[sig] = gi;
      }
      groups[gi].members.push_back(i);
      for (int k : a.dataKeys)
        if (apparatusOfData_.count(k) && apparatusOfData_[k] == i) groups[gi].dataUsers.push_back(k);
      for (int k : a.testKeys)
        if (apparatusOfTest_.count(k) && apparatusOfTest_[k] == i) groups[gi].testUsers.push_back(k);
    }
    for (Group &g : groups) {
      std::sort(g.dataUsers.begin(), g.dataUsers.end());
      g.dataUsers.erase(std::unique(g.dataUsers.begin(), g.dataUsers.end()), g.dataUsers.end());
      const std::string n = nameFor(prefix, g.dataUsers, g.testUsers, apparatus_[g.members[0]].line->index);
      for (int i : g.members) apparatus_[i].*name = n;
    }
  };
  mergeAndName(&Apparatus::hasTarget, &Apparatus::target, &Apparatus::targetName, "tgt_", false);
  mergeAndName(&Apparatus::hasBeam, &Apparatus::beam, &Apparatus::beamName, "beam_", false);
  mergeAndName(&Apparatus::hasDetector, &Apparatus::detector, &Apparatus::detectorName, "det_", false);

  // Parameter names of the effects ("@target" placeholders) now that the targets have names.
  for (auto &p : res_.parameterNames) {
    if (p.second.compare(0, 8, "@target.") == 0) {
      const std::string n = p.first.substr(p.first.find('_') + 1);
      for (const Apparatus &a : apparatus_)
        if (std::to_string(a.line->activeIndex) == n) p.second = a.targetName + p.second.substr(7);
    }
    if (p.second.compare(0, 7, "@sample") == 0) {
      const std::size_t dot = p.second.find('.');
      const int ai = std::atoi(p.second.substr(7, dot - 7).c_str());
      p.second = apparatus_[ai].targetName + p.second.substr(dot);
    }
  }
  // The records, in <targetInt> order, each once.
  std::set<std::string> written;
  for (const Apparatus &a : apparatus_) {
    if (a.dataKeys.empty() && a.testKeys.empty()) continue;
    const std::string where = "# <targetInt> line " + std::to_string(a.line->index);
    if (a.hasTarget && !written.count(a.targetName)) {
      io::Record r = a.target;
      r.name = a.targetName;
      r.comments.push_back(where);
      out_.targets.records.push_back(r);
      written.insert(a.targetName);
    }
    if (a.hasBeam && !written.count(a.beamName)) {
      io::Record r = a.beam;
      r.name = a.beamName;
      r.comments.push_back(where);
      out_.targets.records.push_back(r);
      written.insert(a.beamName);
    }
    if (a.hasDetector && !written.count(a.detectorName)) {
      io::Record r = a.detector;
      r.name = a.detectorName;
      r.comments.push_back(where);
      out_.targets.records.push_back(r);
      written.insert(a.detectorName);
    }
  }
  for (const io::Record &t : sampleTargets_) out_.targets.records.push_back(t);
}

}  // namespace convert
}  // namespace azure2
