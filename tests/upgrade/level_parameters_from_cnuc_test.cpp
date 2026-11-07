// physics/ level_parameters against the engine's own transformation (CNuc::TransformIn, CalcBoundaryConditions,
// TransformOut), on the upgrade test projects T1, T2, T3 (docs/upgrade/level_parameters.md section 7).
//
// Since those three CNuc members call level_parameters themselves (docs/upgrade/cnuc_split.md, step 1), the Brune,
// Park and B_c comparisons of (a) check the wiring of the engine (the System and channel map, the channel-energy
// order, the external-capture order), not the arithmetic: that was checked against the engine's own code before the
// move, and the move by comparing every output file of the engine before and after (cnuc_split.md section 6).  The
// formal comparisons still compare two eigensolvers: the engine passes GSL's, this test uses SymmetricEigen.h's.
//
// For each project the <levels> values are read from a CNuc after CNuc::Fill (the file's widths and ANCs), the channel
// functions and external-capture amplitudes from the engine (src/CNucLevelParameters.cpp), and level_parameters'
// results are compared with a CNuc after CNuc::Initialize in the same basis:
//
//  (a) Brune: amplitudes and B_c, which should be bit-identical (same arithmetic in the same order);
//      Park: amplitudes and Park's J;
//      formal (--no-brune): E_lambda, gamma_{lambda c} up to the sign of a level (the eigenvector sign is GSL's in the
//      engine and SymmetricEigen's here), B_c;
//      Park -> formal (via Park -> Brune): particle-channel amplitudes and energies against the formal CNuc;
//      formal -> Brune (brune_from_formal) against CNuc::TransformOut of the formal CNuc.
//  (b) round trips on the same parameters: observed -> Brune -> observed against the file values; Brune -> Park ->
//      Brune; formal -> Brune -> formal with the iteration tolerance at 1e-13 MeV.
//  (c) convention invariance: every block's Brune amplitudes written in every Convention (LS heavy-first, l-first,
//      i^l, photon signs; jj with either particle and either orbital order) as widths, then level_parameters in that
//      convention, gives the internal Brune amplitudes and the same formal parameters.  Parameter level only: the
//      engine reads no jj project (nothing in it is changed by this branch).
//  (d) T3's 2+ block, which lists one of two allowed channel spins, has no jj form: reported; T3's 1+ block in jj
//      with one column removed: reported and refused.
//
// The formal basis reproduces one engine behaviour that is not a convention: CNuc::TransformIn computes a photon
// channel's external-capture amplitude after the earlier J^pi groups' levels were already replaced by their formal
// parameters, so a final state in an earlier group enters with its formal E and gamma (docs/upgrade/level_parameters.md
// section 6).  The test feeds level_parameters the amplitudes in that order to reproduce the engine, and reports how
// far the photon amplitudes move with the Brune final states throughout.
//
//   level_parameters_from_cnuc_test <T1.azr> <T2.azr> <T3.azr>      (ctest: physics_level_parameters_from_cnuc)
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CNuc.h"
#include "CNucLevelParameters.h"
#include "CNucSystem.h"
#include "Config.h"
#include "CoulFuncCache.h"
#include "ECAmplitudeCache.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace phys = azure2::physics;

namespace {

int failures = 0;

struct Stat {
  long n = 0, equal = 0;
  double worst = 0.0;  // largest |diff| / scale
  void add(double got, double want, double scale) {
    n++;
    if (got == want || (std::isnan(got) && std::isnan(want))) equal++;
    const double d = std::fabs(got - want) / (scale > 0.0 ? scale : 1.0);
    if (!(d <= worst)) worst = d;
  }
};

void report(const std::string &what, const Stat &s, double tol) {
  const bool ok = s.worst <= tol;
  std::printf("  %-4s %-66s %5ld numbers, %5ld identical, largest rel. diff %.3g (tol %.0g)\n", ok ? "ok" : "FAIL",
              what.c_str(), s.n, s.equal, s.worst, tol);
  if (!ok) failures++;
}

void check(const std::string &what, bool ok) {
  std::printf("  %-4s %s\n", ok ? "ok" : "FAIL", what.c_str());
  if (!ok) failures++;
}

enum class Basis { Brune, Park, Formal };

struct Loaded {
  std::ostringstream sink;
  Config cfg;
  CNuc compound;
  Loaded() : cfg(sink) {}
};

// As AZUREAPI::Initialize loads the compound nucleus, external capture on (system_from_cnuc_test.cpp).
std::unique_ptr<Loaded> Load(const std::string &azr, Basis basis, bool initialize) {
  std::unique_ptr<Loaded> m(new Loaded);
  m->cfg.configfile = azr;
  if (basis == Basis::Park) m->cfg.paramMask |= (Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  if (basis == Basis::Formal) m->cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  m->cfg.paramMask |= Config::USE_EXTERNAL_CAPTURE;
  g_config = &m->cfg;
  if (m->cfg.ReadConfigFile() < 0) throw std::runtime_error("cannot read " + azr);
  m->cfg.ResolveMatrixRoute();
  InitializeCoulFuncCache();
  InitializeECAmplitudeCache();
  if (m->compound.Fill(m->cfg) == -1) throw std::runtime_error("CNuc::Fill failed for " + azr);
  if (initialize) m->compound.Initialize(m->cfg);
  return m;
}

double MaxAbs(const std::vector<double> &v) {
  double m = 0.0;
  for (double x : v) m = std::fmax(m, std::fabs(x));
  return m;
}

bool IsPhoton(const phys::System &s, int c) {
  return s.channel(c).kind == phys::ChannelKind::Electric || s.channel(c).kind == phys::ChannelKind::Magnetic;
}

// Overall sign that brings row la of a closest to row la of b.
double RowSign(const std::vector<double> &a, const std::vector<double> &b, std::size_t la, std::size_t nc) {
  double dot = 0.0;
  for (std::size_t k = 0; k < nc; k++) dot += a[la * nc + k] * b[la * nc + k];
  return dot < 0.0 ? -1.0 : 1.0;
}

// Amplitudes of a block as the CNuc's levels hold them.
std::vector<double> CNucGamma(CNuc &cn, const CNucSystemMap &map, std::size_t b, bool transformOut) {
  JGroup *g = cn.GetJGroup(map.groupOfBlock[b]);
  std::vector<double> out;
  for (int la : map.levelOfBlockLevel[b])
    for (int ch : map.channelOfBlockChannel[b])
      out.push_back(transformOut ? g->GetLevel(la)->GetTransformGamma(ch) : g->GetLevel(la)->GetGamma(ch));
  return out;
}
std::vector<double> CNucE(CNuc &cn, const CNucSystemMap &map, std::size_t b, bool transformOut) {
  JGroup *g = cn.GetJGroup(map.groupOfBlock[b]);
  std::vector<double> out;
  for (int la : map.levelOfBlockLevel[b]) out.push_back(transformOut ? g->GetLevel(la)->GetTransformE() : g->GetLevel(la)->GetE());
  return out;
}

void CompareB(const phys::System &sys, std::size_t b, const std::vector<double> &mine, Stat &st) {
  const phys::JpiBlock &blk = sys.block(b);
  const double scale = std::fmax(MaxAbs(blk.B()), 1e-300);
  const bool firstIsParticle = sys.channel(blk.channel(0)).kind == phys::ChannelKind::Particle;
  for (std::size_t k = 0; k < blk.num_channels(); k++) {
    // A non-particle channel copies channel 0's value in the engine; if channel 0 is not a particle channel the
    // engine's value is never set (AChannel leaves it uninitialized): not compared.
    if (sys.channel(blk.channel(k)).kind != phys::ChannelKind::Particle && !firstIsParticle) continue;
    st.add(mine[k], blk.B(k), scale);
  }
}

std::vector<std::pair<std::string, phys::Convention> > Conventions() {
  std::vector<std::pair<std::string, phys::Convention> > out;
  phys::Convention c;
  c.spin_order = phys::SpinOrder::HeavyFirst;
  out.push_back({"LS s = j2 + j1", c});
  c = phys::Convention();
  c.orbital_order = phys::OrbitalOrder::OrbitalFirst;
  c.orbital_phase = phys::OrbitalPhase::IPowL;
  c.electric_sign = -1;
  out.push_back({"LS J = l + s, i^l, E sign -", c});
  for (phys::SpinOrder so : {phys::SpinOrder::LightFirst, phys::SpinOrder::HeavyFirst})
    for (phys::OrbitalOrder oo : {phys::OrbitalOrder::SpinFirst, phys::OrbitalOrder::OrbitalFirst}) {
      c = phys::Convention();
      c.coupling = phys::CouplingScheme::JJ;
      c.spin_order = so;
      c.orbital_order = oo;
      out.push_back({std::string("jj ") + (so == phys::SpinOrder::LightFirst ? "light" : "heavy") +
                         (oo == phys::OrbitalOrder::SpinFirst ? " (j_x l)" : " (l j_x)"),
                     c});
    }
  return out;
}

void CheckProject(const std::string &name, const std::string &azr) {
  std::printf("%s\n", name.c_str());
  std::unique_ptr<Loaded> filled = Load(azr, Basis::Brune, false);
  std::unique_ptr<Loaded> refB = Load(azr, Basis::Brune, true);
  std::unique_ptr<Loaded> refP = Load(azr, Basis::Park, true);
  std::unique_ptr<Loaded> refF = Load(azr, Basis::Formal, true);
  std::unique_ptr<Loaded> emul = Load(azr, Basis::Brune, true);  // Brune CNuc, made formal group by group (below)
  const CNucSystem built = BuildSystemFromCNuc(refB->compound, refB->cfg);
  const CNucSystem builtF = BuildSystemFromCNuc(refF->compound, refF->cfg);
  const phys::System &sys = built.system;
  const CNucSystemMap &map = built.map;
  const phys::ChannelFunctions cf = ChannelFunctionsFromCNuc(filled->compound, map, filled->cfg);
  const bool ignoreZero = !!(refB->cfg.paramMask & Config::IGNORE_ZERO_WIDTHS);

  Stat bruneG, bruneB, parkG, parkJ, formalE, formalG, formalB, parkFormalE, parkFormalG, outE, outG;
  Stat rtObs, rtPark, rtFormalE, rtFormalG, photonShift;
  long signFlips = 0, levels = 0, staticMoments = 0, nonzeroExt = 0, movedExt = 0;
  int maxIterMine = 0, maxIterCNuc = 0;
  std::vector<std::vector<double> > preEC(sys.num_blocks());  // Brune amplitudes before external capture

  // Formal -> Brune by the engine: CNuc::TransformOut on the formal CNuc, its fit values set to the loaded ones.
  {
    CNuc &cn = refF->compound;
    for (int j = 1; j <= cn.NumJGroups(); j++)
      for (int la = 1; la <= cn.GetJGroup(j)->NumLevels(); la++) {
        ALevel *lev = cn.GetJGroup(j)->GetLevel(la);
        lev->SetFitE(lev->GetE());
        for (int ch = 1; ch <= cn.GetJGroup(j)->NumChannels(); ch++) lev->SetFitGamma(ch, lev->GetGamma(ch));
      }
    cn.TransformOut(refF->cfg);
  }

  for (std::size_t b = 0; b < sys.num_blocks(); b++) {
    const phys::JpiBlock &blk = sys.block(b);
    const std::size_t nc = blk.num_channels(), nl = blk.num_levels();
    const std::vector<phys::ProjectChannel> own = phys::internal_channels(sys, b);
    const phys::ProjectLevels in = ProjectLevelsFromCNuc(filled->compound, map, b);
    for (std::size_t la = 0; la < nl; la++)
      for (std::size_t k = 0; k < nc; k++)
        if (IsPhoton(sys, blk.channel(k)) && cf(b, k, in.E[la]).static_moment) staticMoments++;

    // (a) Brune.
    const std::vector<phys::ExternalWidth> extB = ExternalWidthsFromCNuc(refB->compound, map, b, refB->cfg);
    phys::LevelParameterOptions opt;
    opt.external = &extB;
    opt.ignore_zero_widths = ignoreZero;
    opt.formal = false;
    phys::Convention conv;
    const phys::InternalParameters pb = phys::level_parameters(sys, b, conv, own, in, cf, opt);
    preEC[b] = phys::amplitudes_from_observed(sys, b, own, in, phys::LevelBasis::Brune, cf).gamma;
    const std::vector<double> gB = CNucGamma(refB->compound, map, b, false);
    for (std::size_t i = 0; i < gB.size(); i++) bruneG.add(pb.basis_gamma[i], gB[i], std::fmax(std::fabs(gB[i]), 1e-300));
    CompareB(sys, b, pb.B, bruneB);

    // (a) Park.
    const std::vector<phys::ExternalWidth> extP = ExternalWidthsFromCNuc(refP->compound, map, b, refP->cfg);
    phys::Convention convP;
    convP.basis = phys::LevelBasis::Park;
    opt.external = &extP;
    const phys::InternalParameters pp = phys::level_parameters(sys, b, convP, own, in, cf, opt);
    const std::vector<double> gP = CNucGamma(refP->compound, map, b, false);
    for (std::size_t i = 0; i < gP.size(); i++) parkG.add(pp.basis_gamma[i], gP[i], std::fmax(std::fabs(gP[i]), 1e-300));
    {
      JGroup *g = refP->compound.GetJGroup(map.groupOfBlock[b]);
      for (std::size_t la = 0; la < nl; la++)
        parkJ.add(pp.park_norm[la], g->GetLevel(map.levelOfBlockLevel[b][la])->GetParkNorm(), 1.0);
    }

    // (a) formal.  External capture as the engine's TransformIn computes it: final states in earlier groups are
    // already formal (emul has groups < b replaced by the formal CNuc's values).
    const std::vector<phys::ExternalWidth> extF = ExternalWidthsFromCNuc(emul->compound, map, b, emul->cfg);
    phys::Convention convF;
    convF.basis = phys::LevelBasis::Formal;
    opt.external = &extF;
    opt.formal = true;
    const phys::InternalParameters pf = phys::level_parameters(sys, b, convF, own, in, cf, opt);
    const std::vector<double> gF = CNucGamma(refF->compound, map, b, false), EF = CNucE(refF->compound, map, b, false);
    const double scaleF = MaxAbs(gF);
    for (std::size_t la = 0; la < nl; la++) {
      levels++;
      formalE.add(pf.E[la], EF[la], std::fabs(EF[la]));
      const double s = RowSign(pf.gamma, gF, la, nc);
      if (s < 0.0) signFlips++;
      for (std::size_t k = 0; k < nc; k++) formalG.add(s * pf.gamma[la * nc + k], gF[la * nc + k], scaleF);
    }
    CompareB(builtF.system, b, pf.B, formalB);
    {  // the same with the Brune final states throughout: how much the engine's order moves the photon amplitudes
      const std::vector<phys::ExternalWidth> extAll = ExternalWidthsFromCNuc(refB->compound, map, b, refB->cfg);
      for (std::size_t i = 0; i < extF.size(); i++) {
        if (extF[i].present && std::abs(extF[i].value) > 0.0) nonzeroExt++;
        if (extF[i].value != extAll[i].value) movedExt++;
      }
      opt.external = &extAll;
      const phys::InternalParameters alt = phys::level_parameters(sys, b, convF, own, in, cf, opt);
      for (std::size_t la = 0; la < nl; la++)
        for (std::size_t k = 0; k < nc; k++)
          if (IsPhoton(sys, blk.channel(k))) {
            const double s = RowSign(alt.gamma, pf.gamma, la, nc);
            photonShift.add(s * alt.gamma[la * nc + k], pf.gamma[la * nc + k], std::fmax(std::fabs(pf.gamma[la * nc + k]), 1e-300));
          }
    }
    // emul: group b becomes formal for the later groups.
    {
      JGroup *ge = emul->compound.GetJGroup(map.groupOfBlock[b]);
      JGroup *gf = refF->compound.GetJGroup(map.groupOfBlock[b]);
      for (int la : map.levelOfBlockLevel[b]) {
        ge->GetLevel(la)->SetE(gf->GetLevel(la)->GetE());
        for (int ch : map.channelOfBlockChannel[b]) ge->GetLevel(la)->SetGamma(ch, gf->GetLevel(la)->GetGamma(ch));
      }
    }

    // (a) Park -> formal: particle channels (the photon channels' external capture is the engine's order issue).
    if (nl > 0) {
      const std::vector<double> brune = phys::brune_from_park(sys, b, in.E, pp.basis_gamma, cf);
      const phys::LevelSet f = phys::formal_from_brune(sys, b, in.E, brune, pf.B, cf);
      for (std::size_t la = 0; la < nl; la++) {
        parkFormalE.add(f.E[la], EF[la], std::fabs(EF[la]));
        const double s = RowSign(f.gamma, gF, la, nc);
        for (std::size_t k = 0; k < nc; k++)
          if (!IsPhoton(sys, blk.channel(k))) parkFormalG.add(s * f.gamma[la * nc + k], gF[la * nc + k], scaleF);
      }

      // (a) formal -> Brune against CNuc::TransformOut.
      const phys::LevelSet back = phys::brune_from_formal(sys, b, EF, gF, builtF.system.block(b).B(), cf);
      const std::vector<double> gOut = CNucGamma(refF->compound, map, b, true), EOut = CNucE(refF->compound, map, b, true);
      const double scaleO = MaxAbs(gOut);
      JGroup *g = refF->compound.GetJGroup(map.groupOfBlock[b]);
      for (std::size_t la = 0; la < nl; la++) {
        outE.add(back.E[la], EOut[la], std::fabs(EOut[la]));
        const double s = RowSign(back.gamma, gOut, la, nc);
        for (std::size_t k = 0; k < nc; k++) outG.add(s * back.gamma[la * nc + k], gOut[la * nc + k], scaleO);
        maxIterMine = std::max(maxIterMine, back.iterations[la]);
        maxIterCNuc = std::max(maxIterCNuc, g->GetLevel(map.levelOfBlockLevel[b][la])->GetTransformIterations());
      }

      // (b) round trips.
      const phys::ProjectLevels again = phys::observed_from_amplitudes(sys, b, own, in.E, preEC[b], in.is_amplitude,
                                                                        phys::LevelBasis::Brune, cf);
      for (std::size_t i = 0; i < in.value.size(); i++)
        rtObs.add(again.value[i], in.value[i], std::fmax(std::fabs(in.value[i]), 1e-300));
      const std::vector<double> park = phys::park_from_brune(sys, b, in.E, pb.basis_gamma, cf);
      const std::vector<double> back2 = phys::brune_from_park(sys, b, in.E, park, cf);
      for (std::size_t i = 0; i < back2.size(); i++)
        rtPark.add(back2[i], pb.basis_gamma[i], std::fmax(std::fabs(pb.basis_gamma[i]), 1e-300));
      phys::FormalToBruneOptions tight;
      tight.energy_tolerance_MeV = 1e-13;
      const phys::LevelSet b3 = phys::brune_from_formal(sys, b, pf.E, pf.gamma, pf.B, cf, tight);
      const phys::LevelSet f3 = phys::formal_from_brune(sys, b, b3.E, b3.gamma, pf.B, cf);
      const double sc = MaxAbs(pf.gamma);
      for (std::size_t la = 0; la < nl; la++) {
        rtFormalE.add(f3.E[la], pf.E[la], std::fabs(pf.E[la]));
        const double s = RowSign(f3.gamma, pf.gamma, la, nc);
        for (std::size_t k = 0; k < nc; k++) rtFormalG.add(s * f3.gamma[la * nc + k], pf.gamma[la * nc + k], sc);
      }
    }
  }
  report("(a) Brune amplitudes vs CNuc (Brune)", bruneG, 0.0);
  report("(a) B_c vs CNuc::CalcBoundaryConditions", bruneB, 1e-15);
  report("(a) Park amplitudes vs CNuc (--use-park)", parkG, 0.0);
  report("(a) Park's J vs CNuc", parkJ, 0.0);
  report("(a) formal E vs CNuc (--no-brune)", formalE, 1e-12);
  report("(a) formal gamma vs CNuc, up to level sign (rel. to max |gamma|)", formalG, 1e-12);
  report("(a) formal B_c vs CNuc", formalB, 1e-15);
  report("(a) Park -> formal: E vs CNuc (--no-brune)", parkFormalE, 1e-12);
  report("(a) Park -> formal: particle gamma vs CNuc, up to level sign", parkFormalG, 1e-12);
  report("(a) formal -> Brune: E vs CNuc::TransformOut", outE, 1e-12);
  report("(a) formal -> Brune: gamma vs CNuc::TransformOut, up to level sign", outG, 1e-10);
  std::printf("       formal levels whose overall sign differs from GSL's: %ld of %ld; static moments: %ld; "
              "iterations of formal -> Brune: at most %d here, %d in CNuc\n",
              signFlips, levels, staticMoments, maxIterMine, maxIterCNuc);
  report("(b) observed -> Brune -> observed (file values)", rtObs, 1e-13);
  report("(b) Brune -> Park -> Brune", rtPark, 1e-15);
  report("(b) formal -> Brune -> formal: E", rtFormalE, 1e-12);
  report("(b) formal -> Brune -> formal: gamma, up to level sign", rtFormalG, 1e-10);
  std::printf("  info external amplitudes in the formal basis: %ld non-zero, %ld differ between the engine's order and "
              "Brune final states throughout; internal photon amplitudes: %ld numbers, largest rel. change %.3g\n",
              nonzeroExt, movedExt, photonShift.n, photonShift.worst);

  // (c) convention invariance.
  for (const auto &nc_ : Conventions()) {
    Stat amp, fE, fG;
    long gapsBlocks = 0;
    std::string firstGap;
    for (std::size_t b = 0; b < sys.num_blocks(); b++) {
      const phys::JpiBlock &blk = sys.block(b);
      const std::size_t nc = blk.num_channels(), nl = blk.num_levels();
      if (nl == 0) continue;
      phys::Convention conv = nc_.second;
      conv.basis = phys::LevelBasis::Brune;
      const std::vector<phys::ProjectChannel> cols = phys::project_channels(sys, b, conv);
      const std::vector<std::string> gaps = phys::recoupling_gaps(sys, b, conv, cols);
      if (!gaps.empty()) {
        gapsBlocks++;
        if (firstGap.empty()) firstGap = gaps[0];
        continue;
      }
      const phys::ProjectLevels in = ProjectLevelsFromCNuc(filled->compound, map, b);
      const phys::ConventionTransform T = phys::convention_transform(sys, b, conv, cols);
      // The file's values written in the project's columns; amplitude-input (gammaIsRWA) entries stay amplitudes
      // only where the column maps one to one, which holds for T1-T3 (their amplitude inputs are on single channels).
      std::vector<unsigned char> rwa(nl * cols.size(), 0);
      for (std::size_t la = 0; la < nl; la++)
        for (std::size_t k = 0; k < nc; k++)
          if (in.is_amplitude[la * nc + k])
            for (std::size_t m = 0; m < cols.size(); m++)
              if (std::fabs(T(k, m)) > 0.0) rwa[la * cols.size() + m] = 1;
      const std::vector<double> gp = T.to_project(preEC[b], nl);
      const phys::ProjectLevels pin = phys::observed_from_amplitudes(sys, b, cols, in.E, gp, rwa, phys::LevelBasis::Brune, cf);
      const std::vector<phys::ExternalWidth> extB = ExternalWidthsFromCNuc(refB->compound, map, b, refB->cfg);
      phys::LevelParameterOptions opt;
      opt.external = &extB;
      opt.ignore_zero_widths = ignoreZero;
      const phys::InternalParameters ip = phys::level_parameters(sys, b, conv, cols, pin, cf, opt);
      const std::vector<double> gB = CNucGamma(refB->compound, map, b, false);
      const double scale = MaxAbs(gB);
      for (std::size_t i = 0; i < gB.size(); i++) amp.add(ip.basis_gamma[i], gB[i], scale);
      phys::Convention internal;
      opt.external = &extB;
      const phys::InternalParameters ref = phys::level_parameters(sys, b, internal, phys::internal_channels(sys, b), in, cf, opt);
      const double scaleF = MaxAbs(ref.gamma);
      for (std::size_t la = 0; la < nl; la++) {
        fE.add(ip.E[la], ref.E[la], std::fabs(ref.E[la]));
        const double s = RowSign(ip.gamma, ref.gamma, la, nc);
        for (std::size_t k = 0; k < nc; k++) fG.add(s * ip.gamma[la * nc + k], ref.gamma[la * nc + k], scaleF);
      }
    }
    report("(c) " + nc_.first + ": Brune amplitudes vs CNuc", amp, 1e-14);
    report("(c) " + nc_.first + ": formal E vs internal columns", fE, 1e-13);
    report("(c) " + nc_.first + ": formal gamma vs internal, level sign", fG, 1e-12);
    if (gapsBlocks) std::printf("       %ld block(s) not convertible, e.g. %s\n", gapsBlocks, firstGap.c_str());
  }
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: %s T1.azr T2.azr T3.azr\n", argv[0]);
    return 2;
  }
  const char *names[3] = {"T1", "T2", "T3"};
  for (int t = 0; t < 3; t++) {
    try {
      CheckProject(names[t], argv[t + 1]);
    } catch (const std::exception &e) {
      std::printf("  FAIL %s threw: %s\n", names[t], e.what());
      failures++;
    }
  }

  // (d) T3 in jj.  Its 2+ block lists only s = 2 of the two channel spins J = 2 allows for p + 7Li and n + 7Be at
  // l = 1 (j1 = 1/2, j2 = 3/2), so its jj form needs the missing s = 1: reported.  Its 1+ block lists both (p + 7Li,
  // l = 1, s = 1 and 2): complete, and reported once one of its jj columns is removed.
  std::printf("(d) missing channel\n");
  try {
    std::unique_ptr<Loaded> m = Load(argv[3], Basis::Brune, true);
    const CNucSystem built = BuildSystemFromCNuc(m->compound, m->cfg);
    const phys::System &sys = built.system;
    std::size_t b2 = sys.num_blocks(), b1 = sys.num_blocks();
    for (std::size_t i = 0; i < sys.num_blocks(); i++) {
      if (sys.block(i).two_J() == 4 && sys.block(i).parity() == 1) b2 = i;
      if (sys.block(i).two_J() == 2 && sys.block(i).parity() == 1) b1 = i;
    }
    check("T3 has a 2+ and a 1+ block", b2 < sys.num_blocks() && b1 < sys.num_blocks());
    phys::Convention jj;
    jj.coupling = phys::CouplingScheme::JJ;
    const std::vector<std::string> g2 = phys::recoupling_gaps(sys, b2, jj, phys::project_channels(sys, b2, jj));
    for (const std::string &g : g2) std::printf("       reported: %s\n", g.c_str());
    check("T3 2+ (one channel spin of two listed) in jj: reported", !g2.empty());
    std::vector<phys::ProjectChannel> cols = phys::project_channels(sys, b1, jj);
    check("T3 1+ in jj: complete, no gap", phys::recoupling_gaps(sys, b1, jj, cols).empty());
    cols.erase(cols.begin());
    const std::vector<std::string> g1 = phys::recoupling_gaps(sys, b1, jj, cols);
    for (const std::string &g : g1) std::printf("       reported: %s\n", g.c_str());
    check("T3 1+ in jj without its first column: reported", !g1.empty());
    bool threw = false;
    try {
      phys::convention_transform(sys, b1, jj, cols);
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    check("... and refused by convention_transform", threw);
  } catch (const std::exception &e) {
    std::printf("  FAIL (d) threw: %s\n", e.what());
    failures++;
  }

  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
