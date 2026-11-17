// Unit test of the lab -> c.m. frames (azure2/experiment/Frames.h), using experiment/ and physics/ alone
// (docs/upgrade/data_model.md step S6).
//
// 1. The engine's expressions as they stood in src/Unit.cpp and include/Kinematics.h before S6, written out below
//    on plain numbers, against the library on 20000 random kinematic settings each, bit for bit: an operation moved
//    (pi / 180. * x as x * pi / 180., a factor taken out of a sum) changes the last bit somewhere and fails here.
// 2. Closed forms, independent of those expressions:
//    - elastic: tan theta_lab = sin theta_cm / (cos theta_cm + m1/m2); theta_cm = 2 theta_lab for equal masses;
//    - photon aberration inverted, cos theta_lab = (cos theta_cm + beta) / (1 + beta cos theta_cm), beta = p/E of the
//      c.m. from the four-momenta;
//    - Jacobians against a numerical dOmega_lab/dOmega_cm = d cos theta_lab / d cos theta_cm (Richardson central
//      differences of the library's own angle map), to 1e-9 relative: particle (reaction and elastic, non-
//      relativistic), photon, relativistic particle; the 0 deg branch against the limit of the general formula;
//    - the relativistic branch against Lorentz boosts built here from four-vectors: E_cm from the invariant mass,
//      the decay energy from momentum balance, and theta_lab of a c.m. emission angle boosted to the lab and
//      converted back.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>

#include "azure2/experiment/Frames.h"
#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Constants.h"

using namespace azure2::experiment;

namespace {

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), detail.empty() ? "" : "  ", detail.c_str());
  if (!ok) failures++;
}

bool same(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) return true;
  return std::memcmp(&a, &b, sizeof a) == 0;
}

std::string rel(double worst) {
  char b[64];
  std::snprintf(b, sizeof b, "worst relative difference %.1e", worst);
  return b;
}

std::mt19937_64 rng(6162026);
double uni(double lo, double hi) { return std::uniform_real_distribution<double>(lo, hi)(rng); }

const double pi = 3.14159265358979323846;  // include/Constants.h
const double uconv = azure2::physics::kAtomicMassUnit_MeV;
const double hbarc = azure2::physics::kHbarC_MeV_fm;

/// A random partition: masses (u) of light-ion reactions, threshold S + E*.
Partition RandomPartition() {
  Partition p;
  p.m1 = uni(1.0, 5.0);
  p.m2 = uni(1.0, 30.0);
  p.sepE = uni(-5.0, 15.0);
  p.exE = (rng() % 3 == 0) ? uni(0.0, 4.0) : 0.0;
  return p;
}

// ---------------------------------------------------------------- the engine's expressions before S6, written out

struct OldPoint {
  double labEnergy, labAngle, cmAngle, cmEnergy;
  double labAngleMax = -1.0;
  std::string why;
};

double OldCMFromLab(double m1, double m2, double eLab) {
  const double m = m1 + m2;
  return 2.0 * m2 * eLab / (std::sqrt(m * m + 2.0 * m2 * eLab) + m);
}
double OldLabFromCM(double m1, double m2, double T) { return T * (2.0 * (m1 + m2) + T) / (2.0 * m2); }
double OldCMFromDecay(double m1, double m2, double e1) {
  const double p2 = e1 * (e1 + 2.0 * m1);
  return e1 + p2 / (std::sqrt(m2 * m2 + p2) + m2);
}

// Unit::ConvertLabEnergy, both branches: cm, Ex, bin edges
void OldConvertLabEnergy(const Partition &p, bool rel, double eLab, double bl, double bh, double out[4]) {
  if (rel) {
    const double m1 = p.m1 * uconv, m2 = p.m2 * uconv;
    out[0] = OldCMFromLab(m1, m2, eLab);
    out[1] = out[0] + p.sepE + p.exE;
    out[2] = OldCMFromLab(m1, m2, bl);
    out[3] = OldCMFromLab(m1, m2, bh);
    return;
  }
  double factor = (p.m2) / (p.m1 + p.m2);
  out[0] = eLab * factor;
  out[1] = out[0] + p.sepE + p.exE;
  out[2] = bl * factor;
  out[3] = bh * factor;
}

// Unit::ConvertDecayEnergy
double OldConvertDecayEnergy(const Partition &p, bool rel, double eLab) {
  if (rel) return OldCMFromDecay(p.m1 * uconv, p.m2 * uconv, eLab);
  return eLab / (p.m2) * (p.m1 + p.m2);
}

// Unit::ConvertLabAngle(PPair *), non-relativistic part
bool OldElastic(const Partition &p, OldPoint &pt) {
  const double ratio = p.m1 / p.m2;
  if (ratio > 1.0) {
    pt.labAngleMax = 180. / pi * asin(1.0 / ratio);
    if (ratio * sin(pi / 180. * pt.labAngle) > 1.0) return false;
  }
  pt.cmAngle = pt.labAngle + 180. / pi * asin(ratio * sin(pi / 180. * pt.labAngle));
  return true;
}

// Unit::ConvertLabAngleNR
bool OldReactionNR(const Partition &in, const Partition &out, OldPoint &pt) {
  double qValue = in.sepE + in.exE - out.sepE - out.exE;
  double a13 = (in.m1 * out.m1) * pt.labEnergy / (pt.labEnergy + qValue) / (in.m1 + in.m2) / (out.m1 + out.m2);
  double a24 = (in.m2 * out.m2) * (1 + in.m1 / in.m2 * qValue / (pt.labEnergy + qValue)) / (in.m1 + in.m2) / (out.m1 + out.m2);
  if (!(pt.labEnergy + qValue > 0.0) || !(a24 > 0.0)) {
    pt.labAngleMax = 0.0;
    return false;
  }
  if (a13 > a24) {
    pt.labAngleMax = asin(sqrt(a24 / a13)) * 180. / pi;
    if (pt.labAngleMax < pt.labAngle) return false;
  }
  double E3PerEt = a13 * pow(cos(pt.labAngle * pi / 180.) + sqrt(a24 / a13 - pow(sin(pt.labAngle * pi / 180.), 2.0)), 2.0);
  double tempE3PerEt = a13 * pow(cos((pt.labAngle + 0.001) * pi / 180.) + sqrt(a24 / a13 - pow(sin((pt.labAngle + 0.001) * pi / 180.), 2.0)), 2.0);
  double slope = sqrt(tempE3PerEt / a24) * sin((pt.labAngle + 0.001) * pi / 180.) - sqrt(E3PerEt / a24) * sin(pt.labAngle * pi / 180.);
  bool switchDomain = false;
  if (slope < 0.) switchDomain = true;
  pt.cmAngle = 180. / pi * asin(sqrt(E3PerEt / a24) * sin(pt.labAngle * pi / 180.));
  if (switchDomain) pt.cmAngle = 180. - pt.cmAngle;
  return true;
}

// Unit::RelativisticLabToCM
bool OldRelativistic(const Partition &in, const Partition &out, const OldPoint &pt, double nrCMAngle,
                     double &cmAngle, double &factor, double &maxOut) {
  const double m1 = in.m1 * uconv, m2 = in.m2 * uconv;
  const double m3 = out.m1 * uconv, m4 = out.m2 * uconv;
  const double ti = pt.cmEnergy;
  const double tf = ti + in.sepE + in.exE - out.sepE - out.exE;
  if (!(tf > 0.0) || !(ti > 0.0)) {
    maxOut = 0.0;
    return false;
  }
  const double e1 = OldLabFromCM(m1, m2, ti);
  const double sqrtS = m1 + m2 + ti, eTot = e1 + m1 + m2;
  const double pLab = sqrt(e1 * (e1 + 2.0 * m1));
  const double beta = pLab / eTot, gam = eTot / sqrtS;
  const double pf = azure2::physics::relativistic_k(m3, m4, tf) * hbarc, e3 = sqrt(m3 * m3 + pf * pf);
  const double g = beta / (pf / e3);
  const double thl = pt.labAngle * pi / 180.;
  auto beyond = [&]() {
    maxOut = (g > 1.0 && m3 > 0.0) ? asin(std::min(1.0, e3 / (m3 * gam * g))) * 180. / pi : 180.;
    return false;
  };
  double x;
  if (fabs(sin(thl)) < 1e-12) {
    x = cos(thl) > 0 ? 1.0 : -1.0;
    if (g > 1.0) x = (fabs(nrCMAngle) < 90.0) ? 1.0 : -1.0;
  } else {
    const double t2 = pow(tan(thl), 2.0), gt2 = gam * gam * t2;
    const double disc = 1.0 + gt2 * (1.0 - g * g);
    if (disc < 0.0) return beyond();
    const double r1 = (-gt2 * g + sqrt(disc)) / (1.0 + gt2), r2 = (-gt2 * g - sqrt(disc)) / (1.0 + gt2);
    const bool v1 = (r1 + g) * cos(thl) > 0.0 && fabs(r1) <= 1.0, v2 = (r2 + g) * cos(thl) > 0.0 && fabs(r2) <= 1.0;
    const double nrX = cos(nrCMAngle * pi / 180.);
    if (v1 && v2) x = (fabs(r1 - nrX) <= fabs(r2 - nrX)) ? r1 : r2;
    else if (v1) x = r1;
    else if (v2) x = r2;
    else return beyond();
  }
  x = std::max(-1.0, std::min(1.0, x));
  cmAngle = acos(x) * 180. / pi;
  factor = gam * (1.0 + g * x) / pow(gam * gam * (g + x) * (g + x) + 1.0 - x * x, 1.5);
  return true;
}

// static CMVelocity, ConvertLabAngleGammas, CalculateCrossSectionGammaConversionFactor
double OldCMVelocity(double labEnergy, const Partition &in) {
  return sqrt(labEnergy * (labEnergy + 2.0 * in.m1 * uconv)) / ((in.m1 + in.m2) * uconv + labEnergy);
}
double OldGammaAngle(const Partition &in, const OldPoint &pt) {
  double beta = OldCMVelocity(pt.labEnergy, in);
  return 180. / pi * acos((cos(pt.labAngle * pi / 180.) - beta) / (1 - beta * cos(pt.labAngle * pi / 180.)));
}
double OldGammaJacobian(const Partition &in, const OldPoint &pt) {
  double beta = OldCMVelocity(pt.labEnergy, in);
  return pow(1.0 - beta * cos(pt.labAngle * pi / 180.0), 2.0) / (1.0 - beta * beta);
}

// ConvertCrossSection / CalculateCrossSectionConversionFactor, non-relativistic (the two copies were the same text)
double OldParticleJacobian(const Partition &in, const Partition &out, const OldPoint &pt) {
  double conversionFactor;
  if (pt.labAngle == 0.0 || pt.labAngle == 180.0) {
    double m1 = in.m1;
    double m2 = in.m2;
    double m3 = out.m1;
    double m4 = out.m2;
    double e1 = pt.labEnergy;
    double qValue = in.sepE + in.exE - out.sepE - out.exE;
    double et = e1 + qValue;
    double a = m1 * m4 * e1 / (m1 + m2) / (m3 + m4) / et;
    double b = m1 * m3 * e1 / (m1 + m2) / (m3 + m4) / et;
    double c = m2 * m3 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double d = m2 * m4 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double e3et = b + d + 2 * pow(a * c, 0.5) * cos(pi / 180. * pt.cmAngle);
    conversionFactor = pow(a * c, 0.5) * pow(d / b - pow(sin(pt.labAngle * pi / 180.), 2.0), 0.5) / e3et;
  } else {
    conversionFactor = pow(sin(pi / 180. * pt.labAngle) / sin(pi / 180. * pt.cmAngle), 2.0) * cos(pi / 180. * (pt.cmAngle - pt.labAngle));
  }
  return conversionFactor;
}

// ---------------------------------------------------------------- 1. bit for bit

void TestBitForBit() {
  std::printf("the engine's expressions before S6, bit for bit (20000 settings each)\n");
  bool energy = true, decay = true, elastic = true, reaction = true, relat = true, photon = true, jac = true,
       energyMaps = true;
  int reached = 0, beyond = 0;
  for (int k = 0; k < 20000; k++) {
    const Partition in = RandomPartition(), out = RandomPartition();
    const double eLab = uni(1e-3, 30.0), bl = eLab * uni(0.9, 1.0), bh = eLab * uni(1.0, 1.1);
    double angle = uni(0.0, 180.0);
    if (k % 50 == 0) angle = 0.0;
    if (k % 50 == 1) angle = 180.0;
    const bool rel = k % 2;
    double o[4];
    OldConvertLabEnergy(in, rel, eLab, bl, bh, o);
    const BeamEnergiesCM e = BeamToCM(in, eLab, bl, bh, rel);
    energy = energy && same(e.eCM, o[0]) && same(e.eX, o[1]) && same(e.binLowCM, o[2]) && same(e.binHighCM, o[3]);
    decay = decay && same(DecayToCM(out, eLab, rel), OldConvertDecayEnergy(out, rel, eLab));
    energyMaps = energyMaps && same(CMFromLab(in.m1 * uconv, in.m2 * uconv, eLab), OldCMFromLab(in.m1 * uconv, in.m2 * uconv, eLab)) &&
                 same(LabFromCM(in.m1 * uconv, in.m2 * uconv, eLab), OldLabFromCM(in.m1 * uconv, in.m2 * uconv, eLab)) &&
                 same(CMFromDecay(in.m1 * uconv, in.m2 * uconv, eLab), OldCMFromDecay(in.m1 * uconv, in.m2 * uconv, eLab));

    OldPoint pt;
    pt.labEnergy = eLab;
    pt.labAngle = pt.cmAngle = angle;
    pt.cmEnergy = o[0];
    {
      OldPoint q = pt;
      const bool okOld = OldElastic(in, q);
      double cm = angle;
      AngleLimit lim;
      const bool okNew = ElasticLabToCM(in, angle, cm, lim);
      // The old map is the new one except where the new one refuses what the old one let pass (2026-11-17, the
      // PI's decision): a projectile no lighter than the target at a lab angle above 90 deg, which no elastic
      // scattering reaches (the old map gave theta_lab + asin((m1/m2) sin theta_lab) there); and for equal masses
      // the old map set no maximum, the new one 90 deg.
      const double ratio = in.m1 / in.m2;
      if (ratio >= 1.0 && angle > 90.0) {
        elastic = elastic && !okNew && lim.maxSet && !lim.why.empty();
      } else {
        elastic = elastic && okOld == okNew && (!okOld || same(cm, q.cmAngle)) &&
                  (lim.maxSet ? (ratio == 1.0 ? lim.max == 90.0 : same(lim.max, q.labAngleMax)) : q.labAngleMax == -1.0);
      }
    }
    {
      OldPoint q = pt;
      const bool okOld = OldReactionNR(in, out, q);
      double cm = angle;
      AngleLimit lim;
      const bool okNew = ReactionLabToCM(in, out, eLab, angle, cm, lim);
      reaction = reaction && okOld == okNew && (!okOld || same(cm, q.cmAngle)) &&
                 (lim.maxSet ? same(lim.max, q.labAngleMax) : q.labAngleMax == -1.0) && (okNew || !lim.why.empty());
      okOld ? reached++ : beyond++;
      if (okOld) {
        q.cmEnergy = o[0];
        jac = jac && same(ParticleJacobian(in, out, eLab, angle, q.cmAngle), OldParticleJacobian(in, out, q));
        // the relativistic map with the non-relativistic angle as branch
        double c1 = 0, f1 = 0, c2 = 0, f2 = 0, m1 = -1;
        OldPoint r = q;
        r.cmEnergy = BeamToCM(in, eLab, bl, bh, true).eCM;
        const bool a = OldRelativistic(in, out, r, q.cmAngle, c1, f1, m1);
        AngleLimit l2;
        const bool b = RelativisticLabToCM(in, out, r.cmEnergy, angle, q.cmAngle, c2, f2, l2);
        relat = relat && a == b && (!a || (same(c1, c2) && same(f1, f2))) && (a || same(l2.max, m1));
      }
    }
    photon = photon && same(PhotonLabToCM(CMVelocity(in, eLab), angle), OldGammaAngle(in, pt)) &&
             same(PhotonJacobian(in, eLab, angle), OldGammaJacobian(in, pt)) &&
             same(CMVelocity(in, eLab), OldCMVelocity(eLab, in));
  }
  check("beam energies E_cm, E_x, bin edges (both branches)", energy);
  check("decay energy (both branches)", decay);
  check("CMFromLab, LabFromCM, CMFromDecay (Kinematics.h)", energyMaps);
  check("elastic angle and maximum", elastic);
  char d[96];
  std::snprintf(d, sizeof d, "%d reached, %d below threshold or beyond the maximum", reached, beyond);
  check("reaction angle (Marion-Young), maximum, refusal", reaction, d);
  check("relativistic angle, Jacobian and maximum", relat);
  check("particle Jacobian (incl. 0 and 180 deg)", jac);
  check("photon velocity, aberration, Jacobian", photon);
}

// ---------------------------------------------------------------- 2. closed forms

/// d cos(theta_lab) / d cos(theta_cm) from a map theta_cm(theta_lab) (deg), Richardson-extrapolated central
/// differences in theta_lab.
template <class Map>
double NumericalJacobian(Map map, double thetaLab) {
  auto D = [&](double h) {
    const double cp = map(thetaLab + h), cm = map(thetaLab - h);
    return (cos((thetaLab + h) * pi / 180.) - cos((thetaLab - h) * pi / 180.)) / (cos(cp * pi / 180.) - cos(cm * pi / 180.));
  };
  const double h = 1e-2;  // deg; O(h^4) after the extrapolation, rounding ~1e-16 / (h sin theta)
  return (4.0 * D(h / 2) - D(h)) / 3.0;
}

void TestClosedForms() {
  std::printf("closed forms\n");
  // elastic
  {
    double worst = 0.0, worstEq = 0.0;
    for (int k = 0; k < 2000; k++) {
      Partition p = RandomPartition();
      if (p.m1 > p.m2) std::swap(p.m1, p.m2);
      const double th = uni(1.0, 179.0);
      double cm;
      AngleLimit lim;
      ElasticLabToCM(p, th, cm, lim);
      const double t = sin(cm * pi / 180.) / (cos(cm * pi / 180.) + p.m1 / p.m2);
      worst = std::max(worst, std::fabs(t - tan(th * pi / 180.)) / std::fabs(t));
      Partition q = p;
      q.m2 = q.m1;
      const double th2 = uni(1.0, 89.0);
      ElasticLabToCM(q, th2, cm, lim);
      worstEq = std::max(worstEq, std::fabs(cm - 2 * th2) / (2 * th2));
    }
    check("elastic: tan theta_lab = sin theta_cm / (cos theta_cm + m1/m2)", worst < 1e-12, rel(worst));
    check("elastic, equal masses: theta_cm = 2 theta_lab", worstEq < 1e-13, rel(worstEq));
    Partition heavy{12.0, 4.0, 0, 0};
    double cm = -1;
    AngleLimit lim;
    const bool ok = ElasticLabToCM(heavy, 30.0, cm, lim);
    check("elastic, m1 > m2: refused beyond asin(m2/m1)", !ok && lim.maxSet && std::fabs(lim.max - asin(1.0 / 3.0) * 180 / pi) < 1e-12 && !lim.why.empty());
    // A backward lab angle is beyond the maximum too, although (m1/m2) sin theta_lab < 1 there (170 deg: 0.52).
    AngleLimit back;
    check("elastic, m1 > m2: a lab angle above 90 deg is refused", !ElasticLabToCM(heavy, 170.0, cm, back) && back.maxSet && !back.why.empty());
    // Equal masses: theta_cm = 2 theta_lab reaches 180 deg at 90 deg lab; nothing is scattered beyond.
    Partition equal{4.0, 4.0, 0, 0};
    AngleLimit eq, at;
    double cmEq = -1, cmAt = -1;
    check("elastic, equal masses: a lab angle above 90 deg is refused, the maximum is 90 deg",
          !ElasticLabToCM(equal, 100.0, cmEq, eq) && eq.maxSet && std::fabs(eq.max - 90.0) < 1e-12 && !eq.why.empty());
    check("elastic, equal masses: 80 deg lab is 160 deg c.m.", ElasticLabToCM(equal, 80.0, cmAt, at) && std::fabs(cmAt - 160.0) < 1e-12);
  }
  // photon aberration, beta from four-momenta
  {
    double worst = 0.0, worstBeta = 0.0;
    for (int k = 0; k < 2000; k++) {
      const Partition p = RandomPartition();
      const double eLab = uni(0.01, 20.0), th = uni(0.0, 180.0);
      const double m1 = p.m1 * uconv, m2 = p.m2 * uconv;
      const double pz = std::sqrt((eLab + m1) * (eLab + m1) - m1 * m1), E = eLab + m1 + m2;
      const double beta = CMVelocity(p, eLab);
      worstBeta = std::max(worstBeta, std::fabs(beta - pz / E) / (pz / E));
      const double cm = PhotonLabToCM(beta, th) * pi / 180.;
      const double back = (cos(cm) + beta) / (1 + beta * cos(cm));
      worst = std::max(worst, std::fabs(back - cos(th * pi / 180.)));
    }
    check("photon: beta = p/E of the c.m. (four-momenta)", worstBeta < 1e-10, rel(worstBeta));
    check("photon aberration inverted: cos theta_lab = (cos theta_cm + beta)/(1 + beta cos theta_cm)", worst < 1e-13,
          rel(worst));
  }
  // Jacobians against the numerical dOmega_lab/dOmega_cm
  {
    double wR = 0.0, wE = 0.0, wP = 0.0, wRel = 0.0, w0 = 0.0;
    int n = 0;
    for (int k = 0; k < 400; k++) {
      const Partition in = RandomPartition(), out = RandomPartition();
      const double eLab = uni(0.5, 30.0), th = uni(5.0, 175.0);
      double cm;
      AngleLimit lim;
      auto react = [&](double t) {
        double c;
        AngleLimit l;
        return ReactionLabToCM(in, out, eLab, t, c, l) ? c : std::nan("");
      };
      if (ReactionLabToCM(in, out, eLab, th, cm, lim) && std::isfinite(react(th + 1e-3)) && std::isfinite(react(th - 1e-3)) &&
          (!lim.maxSet || th < lim.max - 0.5) && std::fabs(cm - 90.0) > 5.0) {
        const double J = ParticleJacobian(in, out, eLab, th, cm), N = NumericalJacobian(react, th);
        wR = std::max(wR, std::fabs(J - N) / std::fabs(J));
        n++;
        // the 0 deg branch: the limit of the general formula
        double c0, cs;
        AngleLimit l0;
        if (ReactionLabToCM(in, out, eLab, 0.0, c0, l0) && ReactionLabToCM(in, out, eLab, 1e-4, cs, l0) && c0 < 90.0) {
          const double J0 = ParticleJacobian(in, out, eLab, 0.0, c0), Js = ParticleJacobian(in, out, eLab, 1e-4, cs);
          w0 = std::max(w0, std::fabs(J0 - Js) / std::fabs(J0));
        }
      }
      Partition el = in;
      if (el.m1 > el.m2) std::swap(el.m1, el.m2);
      auto elastic = [&](double t) {
        double c;
        AngleLimit l;
        ElasticLabToCM(el, t, c, l);
        return c;
      };
      ElasticLabToCM(el, th, cm, lim);
      if (std::fabs(cm - th - 90.0) > 5.0) {  // asin(r sin theta) near 1: the map's own precision is sqrt(eps)
        const double JE = ParticleJacobian(el, el, eLab, th, cm), NE = NumericalJacobian(elastic, th);
        wE = std::max(wE, std::fabs(JE - NE) / std::fabs(JE));
      }
      const double beta = CMVelocity(in, eLab);
      auto phot = [&](double t) { return PhotonLabToCM(beta, t); };
      const double JP = PhotonJacobian(in, eLab, th), NP = NumericalJacobian(phot, th);
      wP = std::max(wP, std::fabs(JP - NP) / std::fabs(JP));
      // relativistic particle, at the c.m. energy of the beam
      const double eCM = BeamToCM(in, eLab, eLab, eLab, true).eCM;
      double cr, fr;
      AngleLimit lr;
      auto relmap = [&](double t) {
        double c, f;
        AngleLimit l;
        return RelativisticLabToCM(in, out, eCM, t, cm, c, f, l) ? c : std::nan("");
      };
      // away from the maximum lab angle (sqrt(disc) -> 0) and from theta* = 0, 180 (acos near +-1)
      if (RelativisticLabToCM(in, out, eCM, th, th, cr, fr, lr) && std::isfinite(relmap(th + 2.0)) &&
          std::isfinite(relmap(th - 2.0)) && cr > 5.0 && cr < 175.0) {
        // the branch nearest cr on both sides: use cr itself as the selector
        auto relmap2 = [&](double t) {
          double c, f;
          AngleLimit l;
          return RelativisticLabToCM(in, out, eCM, t, cr, c, f, l) ? c : std::nan("");
        };
        const double N = NumericalJacobian(relmap2, th);
        if (std::isfinite(N)) wRel = std::max(wRel, std::fabs(fr - N) / std::fabs(fr));
      }
    }
    char d[64];
    std::snprintf(d, sizeof d, "%d reactions, %s", n, rel(wR).c_str());
    check("particle Jacobian (reaction) = numerical dOmega_lab/dOmega_cm", wR < 1e-9, d);
    check("particle Jacobian (elastic) = numerical", wE < 1e-9, rel(wE));
    check("particle Jacobian at 0 deg = limit of the general formula", w0 < 1e-7, rel(w0));
    check("photon Jacobian = numerical", wP < 1e-9, rel(wP));
    check("relativistic particle Jacobian = numerical", wRel < 1e-9, rel(wRel));
  }
}

/// d theta_lab / d theta* (both in deg / rad as given by toLab) by a central difference.
template <class F>
double D0(F toLab, double ts) {
  return (toLab(ts + 1e-4) - toLab(ts - 1e-4)) / 2e-4 * (3.14159265358979323846 / 180.);
}

/// Relativistic kinematics from four-vectors, built here (c = 1, MeV): beam particle 1 (mass m1, kinetic energy T1)
/// on particle 2 at rest; exit particles 3, 4.
void TestRelativistic() {
  std::printf("relativistic branch against Lorentz boosts\n");
  double wE = 0.0, wD = 0.0, wA = 0.0, wJ = 0.0, wNR = 0.0;
  int n = 0;
  for (int k = 0; k < 2000; k++) {
    const Partition in = RandomPartition(), out = RandomPartition();
    const double m1 = in.m1 * uconv, m2 = in.m2 * uconv, m3 = out.m1 * uconv, m4 = out.m2 * uconv;
    const double T1 = uni(0.1, 50.0);
    // invariant mass: s = (E1 + m2)^2 - p1^2
    const double E1 = T1 + m1, p1 = std::sqrt(E1 * E1 - m1 * m1);
    const double sqrtS = std::sqrt((E1 + m2) * (E1 + m2) - p1 * p1);
    const double eCM = BeamToCM(in, T1, T1, T1, true).eCM;
    wE = std::max(wE, std::fabs(eCM - (sqrtS - m1 - m2)) / (sqrtS - m1 - m2));
    // decay at rest: particle 1 of kinetic energy e1, |p| equal and opposite
    const double e1 = uni(0.01, 10.0), pd = std::sqrt(e1 * (e1 + 2 * m1)), T = e1 + std::sqrt(m2 * m2 + pd * pd) - m2;
    wD = std::max(wD, std::fabs(DecayToCM(in, e1, true) - T) / T);
    // the exit: total c.m. kinetic energy tf = eCM + Q
    const double tf = eCM + in.sepE + in.exE - out.sepE - out.exE;
    if (!(tf > 0.1)) continue;
    const double W = m3 + m4 + tf;  // sqrt(s) of the exit, = sqrt(s) of the entrance plus Q
    const double pStar = std::sqrt((W * W - (m3 + m4) * (m3 + m4)) * (W * W - (m3 - m4) * (m3 - m4))) / (2 * W);
    const double E3s = std::sqrt(m3 * m3 + pStar * pStar);
    // boost of the c.m. in the lab, from the entrance four-momentum
    const double beta = p1 / (E1 + m2), gam = (E1 + m2) / sqrtS;
    const double thStar = uni(5.0, 175.0) * pi / 180.;
    auto toLab = [&](double ts) {
      const double pz = gam * (pStar * cos(ts) + beta * E3s), px = pStar * sin(ts);
      return atan2(px, pz) * 180. / pi;
    };
    const double thLab = toLab(thStar);
    double cm, J;
    AngleLimit lim;
    if (!RelativisticLabToCM(in, out, eCM, thLab, thStar * 180. / pi, cm, J, lim)) continue;
    // away from the maximum lab angle, where theta_lab(theta*) is stationary and J -> infinity
    if (std::fabs(D0(toLab, thStar)) < 0.05) continue;
    n++;
    wA = std::max(wA, std::fabs(cm - thStar * 180. / pi) / (thStar * 180. / pi));
    // dOmega_lab/dOmega_cm = d cos theta_lab / d cos theta*, by central differences of the boost
    auto D = [&](double h) {
      return (cos(toLab(thStar + h) * pi / 180.) - cos(toLab(thStar - h) * pi / 180.)) / (cos(thStar + h) - cos(thStar - h));
    };
    const double N = (4.0 * D(5e-6) - D(1e-5)) / 3.0;
    wJ = std::max(wJ, std::fabs(J - N) / std::fabs(J));
  }
  // non-relativistic limit: elastic scattering at 1 keV/u, the relativistic angle and Jacobian equal the
  // non-relativistic ones to O(T/m) ~ 1e-6
  for (int k = 0; k < 500; k++) {
    Partition p = RandomPartition();
    if (p.m1 > p.m2) std::swap(p.m1, p.m2);
    const double eLab = 1e-3 * p.m1, th = uni(5.0, 175.0);
    double nr, r, f;
    AngleLimit l;
    ElasticLabToCM(p, th, nr, l);
    if (!RelativisticLabToCM(p, p, BeamToCM(p, eLab, eLab, eLab, true).eCM, th, nr, r, f, l)) continue;
    wNR = std::max(wNR, std::fabs(r - nr) / nr);
    wNR = std::max(wNR, std::fabs(f - ParticleJacobian(p, p, eLab, th, nr)) / f);
  }
  // (the reference subtracts m1 + m2 ~ 10^4 MeV from sqrt(s): its own rounding is ~1e-16 m / T)
  check("E_cm = sqrt(s) - m1 - m2 (invariant mass)", wE < 1e-9, rel(wE));
  check("decay at rest: T = e1 + sqrt(m2^2 + p^2) - m2", wD < 1e-9, rel(wD));
  char d[64];
  std::snprintf(d, sizeof d, "%d settings, %s", n, rel(wA).c_str());
  check("theta* boosted to the lab and converted back", wA < 1e-9, d);
  check("relativistic Jacobian = d cos theta_lab / d cos theta* of the boost", wJ < 1e-9, rel(wJ));
  check("relativistic -> non-relativistic angle and Jacobian at 1 keV/u (elastic)", wNR < 1e-5, rel(wNR));
}

}  // namespace

// ---------------------------------------------------------------- 4. the kinematic limit

/// At the maximum lab angle dOmega_lab/dOmega_cm goes to 0: elastic 12C on p (m1 > m2), the angle at the limit as the
/// map gives it, a hair inside it, and well inside; the refusal of a datum there (E_LAB_JACOBIAN_AT_LIMIT).
void TestJacobianAtLimit() {
  Partition p;
  p.m1 = 12.0;
  p.m2 = 1.00727646688;
  FrameRule r;
  r.entrance = r.exit = p;
  r.exitKind = ExitKind::Particle;
  r.elastic = true;
  r.labDifferential = true;
  r.convertValue = true;
  double cm = 0.0;
  AngleLimit lim;
  ElasticLabToCM(p, 1.0, cm, lim);
  const double limit = lim.max;
  const auto jacobianAt = [&r](double thetaLab, bool &reached) {
    FrameRow row;
    AngleLimit l;
    reached = ConvertRow(r, 1.0, thetaLab, 0.0, 0.0, true, row, l);
    return row.jacobian;
  };
  bool atReached = false, nearReached = false, insideReached = false;
  const double at = jacobianAt(limit, atReached);
  const double near = jacobianAt(limit * (1.0 - 1e-14), nearReached);
  const double inside = jacobianAt(0.5 * limit, insideReached);
  char d[160];
  std::snprintf(d, sizeof d, "theta_max %.17g deg: J %.3g at it, %.3g a hair inside, %.3g at half", limit, at, near, inside);
  check("at the maximum lab angle the angle converts and J is below kMinLabJacobian", atReached && nearReached &&
        !(at > kMinLabJacobian) && !(near > kMinLabJacobian), d);
  check("half way to it J is far above kMinLabJacobian", insideReached && inside > 1e3 * kMinLabJacobian, d);
  const Diagnostic diag = LabJacobianAtLimit("segment #3 (data/x.dat)", 1.0, limit, at);
  check("E_LAB_JACOBIAN_AT_LIMIT: the code, the point and the way out",
        diag.code == "E_LAB_JACOBIAN_AT_LIMIT" && diag.severity == Severity::Error &&
            diag.message.find("segment #3 (data/x.dat)") != std::string::npos &&
            diag.message.find("c.m. frame") != std::string::npos);
}

int main() {
  TestBitForBit();
  TestClosedForms();
  TestRelativistic();
  TestJacobianAtLimit();
  if (failures) {
    std::printf("%d FAILURE(S)\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
