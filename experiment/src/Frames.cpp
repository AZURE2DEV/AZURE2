#include "azure2/experiment/Frames.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Constants.h"

// The expressions below are those of src/Unit.cpp (ConvertLabEnergy, ConvertDecayEnergy, ConvertLabAngle,
// ConvertLabAngleNR, RelativisticLabToCM, ConvertLabAngleGammas, ConvertCrossSection,
// CalculateCrossSectionConversionFactor, CalculateCrossSectionGammaConversionFactor) and include/Kinematics.h
// (CMFromLab, LabFromCM, CMFromDecay) as they stood, written with the same operations in the same order: the engine
// calls these, and its results must not move by a bit.  Do not "simplify" them (pi / 180. * x is not x * pi / 180.).

namespace azure2 {
namespace experiment {

namespace {
const double pi = kPi;
const double uconv = physics::kAtomicMassUnit_MeV;
const double hbarc = physics::kHbarC_MeV_fm;
}  // namespace

// ---------------------------------------------------------------- energies

double CMFromLab(double m1, double m2, double eLab) {
  const double m = m1 + m2;
  return 2.0 * m2 * eLab / (std::sqrt(m * m + 2.0 * m2 * eLab) + m);
}

double LabFromCM(double m1, double m2, double T) { return T * (2.0 * (m1 + m2) + T) / (2.0 * m2); }

double CMFromDecay(double m1, double m2, double e1) {
  const double p2 = e1 * (e1 + 2.0 * m1);
  return e1 + p2 / (std::sqrt(m2 * m2 + p2) + m2);
}

BeamEnergiesCM BeamToCM(const Partition &p, double eLab, double binLowLab, double binHighLab, bool relativistic) {
  BeamEnergiesCM r;
  if (relativistic) {
    // T = sqrt(s) - m1 - m2 with s = (m1 + m2)^2 + 2 m2 E_lab
    const double m1 = p.m1 * uconv, m2 = p.m2 * uconv;
    r.eCM = CMFromLab(m1, m2, eLab);
    r.eX = ExcitationEnergy(r.eCM, p);
    r.binLowCM = CMFromLab(m1, m2, binLowLab);
    r.binHighCM = CMFromLab(m1, m2, binHighLab);
    return r;
  }
  const double factor = BeamCMFactor(p);
  r.eCM = eLab * factor;
  r.eX = ExcitationEnergy(r.eCM, p);
  r.binLowCM = binLowLab * factor;
  r.binHighCM = binHighLab * factor;
  return r;
}

double BeamLabFromCM(const Partition &p, double eCM, bool relativistic) {
  if (relativistic) return LabFromCM(p.m1 * uconv, p.m2 * uconv, eCM);
  return eCM / BeamCMFactor(p);
}

double DecayToCM(const Partition &p, double eLab, bool relativistic) {
  // two-body decay at rest: equal momenta, T = E1 + sqrt(m2^2 + p^2) - m2
  if (relativistic) return CMFromDecay(p.m1 * uconv, p.m2 * uconv, eLab);
  return eLab / (p.m2) * (p.m1 + p.m2);
}

// ---------------------------------------------------------------- angles

bool ElasticLabToCM(const Partition &p, double thetaLab, double &thetaCM, AngleLimit &lim) {
  const double ratio = p.m1 / p.m2;
  if (ratio > 1.0) {
    lim.maxSet = true;
    lim.max = 180. / pi * asin(1.0 / ratio);
    if (ratio * sin(pi / 180. * thetaLab) > 1.0) {
      std::ostringstream why;
      why << "beyond the kinematic maximum " << lim.max << " deg of elastic scattering with a projectile heavier than "
          << "the target (only the forward branch is modelled; give the point as c.m., code 4 / 5, or drop it)";
      lim.why = why.str();
      return false;
    }
  }
  thetaCM = thetaLab + 180. / pi * asin(ratio * sin(pi / 180. * thetaLab));
  return true;
}

bool ReactionLabToCM(const Partition &in, const Partition &out, double eLab, double thetaLab, double &thetaCM,
                     AngleLimit &lim) {
  double qValue = in.sepE + in.exE - out.sepE - out.exE;
  double a13 = (in.m1 * out.m1) * eLab / (eLab + qValue) / (in.m1 + in.m2) / (out.m1 + out.m2);
  double a24 = (in.m2 * out.m2) * (1 + in.m1 / in.m2 * qValue / (eLab + qValue)) / (in.m1 + in.m2) / (out.m1 + out.m2);

  // Below the threshold of the exit pair (classically E_lab + Q <= 0 in the c.m. share, i.e. a24 <= 0) there is
  // no ejectile; the Marion-Young map then has no meaning.
  if (!(eLab + qValue > 0.0) || !(a24 > 0.0)) {
    std::ostringstream why;
    why << "below the threshold of the exit pair (Q = " << qValue << " MeV): the lab angle cannot be converted; give "
        << "the point as c.m. (code 4 / 5) or drop it";
    lim.why = why.str();
    lim.maxSet = true;
    lim.max = 0.0;
    return false;
  }
  if (a13 > a24) {
    // Two c.m. angles reach every lab angle below thetaMax and none above it; only the forward branch is modelled.
    lim.maxSet = true;
    lim.max = asin(sqrt(a24 / a13)) * 180. / pi;
    if (lim.max < thetaLab) {
      std::ostringstream why;
      why << "beyond the kinematic maximum " << lim.max << " deg at this energy (two c.m. angles reach every lab "
          << "angle below it and only the forward branch is modelled; give the point as c.m., code 4 / 5, or drop it)";
      lim.why = why.str();
      return false;
    }
  }

  double E3PerEt = a13 * pow(cos(thetaLab * pi / 180.) + sqrt(a24 / a13 - pow(sin(thetaLab * pi / 180.), 2.0)), 2.0);
  double tempE3PerEt = a13 * pow(cos((thetaLab + 0.001) * pi / 180.) + sqrt(a24 / a13 - pow(sin((thetaLab + 0.001) * pi / 180.), 2.0)), 2.0);
  double slope = sqrt(tempE3PerEt / a24) * sin((thetaLab + 0.001) * pi / 180.) - sqrt(E3PerEt / a24) * sin(thetaLab * pi / 180.);
  bool switchDomain = false;
  if (slope < 0.) switchDomain = true;

  thetaCM = 180. / pi * asin(sqrt(E3PerEt / a24) * sin(thetaLab * pi / 180.));
  if (switchDomain) thetaCM = 180. - thetaCM;
  return true;
}

bool RelativisticLabToCM(const Partition &in, const Partition &out, double ti, double thetaLab, double nrCMAngle,
                         double &cmAngle, double &factor, AngleLimit &lim) {
  const double m1 = in.m1 * uconv, m2 = in.m2 * uconv;
  const double m3 = out.m1 * uconv, m4 = out.m2 * uconv;
  const double tf = ti + in.sepE + in.exE - out.sepE - out.exE;
  if (!(tf > 0.0) || !(ti > 0.0)) {  // below the exit threshold (relativistic): no ejectile (K would take |k^2|)
    std::ostringstream why;
    why << "below the threshold of the exit pair (relativistic kinematics, c.m. energy in the exit channel " << tf
        << " MeV): the lab angle cannot be converted; give the point as c.m. (code 4 / 5) or drop it";
    lim.why = why.str();
    lim.maxSet = true;
    lim.max = 0.0;
    return false;
  }
  const double e1 = LabFromCM(m1, m2, ti);
  const double sqrtS = m1 + m2 + ti, eTot = e1 + m1 + m2;
  const double pLab = sqrt(e1 * (e1 + 2.0 * m1));
  const double beta = pLab / eTot, gam = eTot / sqrtS;
  const double pf = physics::relativistic_k(m3, m4, tf) * hbarc, e3 = sqrt(m3 * m3 + pf * pf);
  const double g = beta / (pf / e3);
  const double thl = thetaLab * pi / 180.;
  // Where two c.m. angles reach one lab angle (g > 1) the lab angle has a maximum, sin theta_max = gamma_3* /
  // (gamma g) with gamma_3* = E_3*/m_3.  It is smaller than the non-relativistic Marion-Young maximum.
  auto beyond = [&]() {
    lim.maxSet = true;
    lim.max = (g > 1.0 && m3 > 0.0) ? asin(std::min(1.0, e3 / (m3 * gam * g))) * 180. / pi : 180.;
    std::ostringstream why;
    why << "beyond the kinematic maximum " << lim.max << " deg at this energy (relativistic kinematics; two c.m. "
        << "angles reach every lab angle below it and only the forward branch is modelled; give the point as c.m., "
        << "code 4 / 5, or drop it)";
    lim.why = why.str();
    return false;
  };
  double x;  // cos theta*
  if (fabs(sin(thl)) < 1e-12) {
    x = cos(thl) > 0 ? 1.0 : -1.0;
    if (g > 1.0) x = (fabs(nrCMAngle) < 90.0) ? 1.0 : -1.0;
  } else {
    const double t2 = pow(tan(thl), 2.0), gt2 = gam * gam * t2;
    const double disc = 1.0 + gt2 * (1.0 - g * g);
    if (disc < 0.0) return beyond();
    const double r1 = (-gt2 * g + sqrt(disc)) / (1.0 + gt2), r2 = (-gt2 * g - sqrt(disc)) / (1.0 + gt2);
    // a root is valid if cos theta* + g has the sign of cos theta_l
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

double CMVelocity(const Partition &in, double labEnergy) {
  return sqrt(labEnergy * (labEnergy + 2.0 * in.m1 * uconv)) / ((in.m1 + in.m2) * uconv + labEnergy);
}

double PhotonLabToCM(double beta, double thetaLab) {
  return 180. / pi * acos((cos(thetaLab * pi / 180.) - beta) / (1 - beta * cos(thetaLab * pi / 180.)));
}

// ---------------------------------------------------------------- Jacobians

double ParticleJacobian(const Partition &in, const Partition &out, double eLab, double thetaLab, double thetaCM) {
  double conversionFactor;
  if (thetaLab == 0.0 || thetaLab == 180.0) {
    double m1 = in.m1;
    double m2 = in.m2;
    double m3 = out.m1;
    double m4 = out.m2;
    double e1 = eLab;
    double qValue = in.sepE + in.exE - out.sepE - out.exE;
    double et = e1 + qValue;
    double a = m1 * m4 * e1 / (m1 + m2) / (m3 + m4) / et;
    double b = m1 * m3 * e1 / (m1 + m2) / (m3 + m4) / et;
    double c = m2 * m3 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double d = m2 * m4 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double e3et = b + d + 2 * pow(a * c, 0.5) * cos(pi / 180. * thetaCM);
    conversionFactor = pow(a * c, 0.5) * pow(d / b - pow(sin(thetaLab * pi / 180.), 2.0), 0.5) / e3et;
  } else {
    conversionFactor = pow(sin(pi / 180. * thetaLab) / sin(pi / 180. * thetaCM), 2.0) * cos(pi / 180. * (thetaCM - thetaLab));
  }
  return conversionFactor;
}

double PhotonJacobian(const Partition &in, double eLab, double thetaLab) {
  // dOmega_lab/dOmega_cm = d cos(theta_lab)/d cos(theta_cm) = (1 - beta^2)/(1 + beta cos theta_cm)^2
  //                      = (1 - beta cos theta_lab)^2 / (1 - beta^2).
  double beta = CMVelocity(in, eLab);
  return pow(1.0 - beta * cos(thetaLab * pi / 180.0), 2.0) / (1.0 - beta * beta);
}

Diagnostic LabJacobianAtLimit(const std::string &where, double eLab, double thetaLab, double jacobian) {
  Diagnostic d;
  d.code = codes::kLabJacobianAtLimit;
  d.section = "data";
  d.recordKind = "dataset";
  std::ostringstream s;
  s.precision(12);
  s << where << ": the point at E_lab = " << eLab << " MeV, theta_lab = " << thetaLab
    << " deg lies at the kinematic maximum lab angle: dOmega_lab/dOmega_cm = " << jacobian
    << ", so the model in the lab frame diverges there; drop the point or give it in the c.m. frame (code 4 / 5)";
  std::ostringstream f;
  f << jacobian;
  d.hasFound = true;
  d.found = f.str();
  d.message = s.str();
  return d;
}

FrameRow UnconvertedRow(double eLab, double thetaLab, double binLowLab, double binHighLab) {
  FrameRow r;
  r.eLab = eLab;
  r.eCM = eLab;
  r.eX = eLab;
  r.binLowCM = binLowLab;
  r.binHighCM = binHighLab;
  r.thetaCM = thetaLab;
  r.jacobian = 1.0;
  return r;
}

void Merge(AngleLimit &into, const AngleLimit &from) {
  if (from.maxSet) {
    into.maxSet = true;
    into.max = from.max;
  }
  if (!from.why.empty()) into.why = from.why;
}

namespace {

// The relativistic map of the detected particle at the row's c.m. energy (Unit::RelativisticLabToCM before).
bool RelativisticAngle(const Partition &in, const Partition &out, const FrameRow &row, double thetaLab,
                       double nrCMAngle, double &cmAngle, double &factor, AngleLimit &lim) {
  AngleLimit l;
  const bool ok = RelativisticLabToCM(in, out, row.eCM, thetaLab, nrCMAngle, cmAngle, factor, l);
  Merge(lim, l);
  return ok;
}

// The non-relativistic reaction angle (Marion-Young), written into the row where it is reached.
bool ReactionAngleNR(const FrameRule &r, const FrameRow &row, double thetaLab, double &thetaCM, AngleLimit &lim) {
  AngleLimit l;
  double cmAngle;
  const bool ok = ReactionLabToCM(r.entrance, r.exit, row.eLab, thetaLab, cmAngle, l);
  Merge(lim, l);
  if (!ok) return false;
  thetaCM = cmAngle;
  return true;
}

// The c.m. angle of a detected particle (Unit::ConvertLabAngle, both forms, before).
bool ParticleAngle(const FrameRule &r, FrameRow &row, double thetaLab, AngleLimit &lim) {
  if (r.elastic) {
    AngleLimit l;
    double nrAngle;
    const bool ok = ElasticLabToCM(r.entrance, thetaLab, nrAngle, l);
    Merge(lim, l);
    if (!ok) return false;
    row.thetaCM = nrAngle;
    double cmAngle, factor;
    if (r.relativistic) {
      if (!RelativisticAngle(r.entrance, r.entrance, row, thetaLab, row.thetaCM, cmAngle, factor, lim)) return false;
      row.thetaCM = cmAngle;
    }
    return true;
  }
  // relativistic=1: the relativistic map decides whether the lab angle is reached and gives the angle; the
  // Marion-Young angle only picks the branch.
  if (r.relativistic) {
    const bool nrOk = ReactionAngleNR(r, row, thetaLab, row.thetaCM, lim);
    double relAngle, relFactor;
    if (!RelativisticAngle(r.entrance, r.exit, row, thetaLab, nrOk ? row.thetaCM : 0.0, relAngle, relFactor, lim))
      return false;
    row.thetaCM = relAngle;
    return true;
  }
  return ReactionAngleNR(r, row, thetaLab, row.thetaCM, lim);
}

// dOmega_lab/dOmega_cm of a detected particle at the row's angles.  Relativistic: where the map does not reach the
// lab angle the Jacobian the row started with is kept.
void ParticleRowJacobian(const FrameRule &r, FrameRow &row, double thetaLab, AngleLimit &lim) {
  if (r.relativistic) {
    double relAngle, factor;
    if (RelativisticAngle(r.entrance, r.exit, row, thetaLab, row.thetaCM, relAngle, factor, lim)) row.jacobian = factor;
    return;
  }
  row.jacobian = ParticleJacobian(r.entrance, r.exit, row.eLab, thetaLab, row.thetaCM);
}

}  // namespace

void ConvertRowEnergy(const FrameRule &r, double eLab, double binLowLab, double binHighLab, FrameRow &row) {
  row.eLab = eLab;
  if (r.betaDelayed) {
    row.eCM = DecayToCM(r.exit, eLab, r.relativistic);
    row.eX = ExcitationEnergy(row.eCM, r.exit);
  } else {
    const BeamEnergiesCM e = BeamToCM(r.entrance, eLab, binLowLab, binHighLab, r.relativistic);
    row.eCM = e.eCM;
    row.eX = e.eX;
    row.binLowCM = e.binLowCM;
    row.binHighCM = e.binHighCM;
  }
}

bool ConvertRow(const FrameRule &r, double eLab, double thetaLab, double binLowLab, double binHighLab, bool atLoad,
                FrameRow &row, AngleLimit &lim) {
  ConvertRowEnergy(r, eLab, binLowLab, binHighLab, row);
  if (r.exitKind == ExitKind::None || !r.labDifferential) return true;
  bool inRange = true;
  if (r.exitKind == ExitKind::Particle) {
    inRange = ParticleAngle(r, row, thetaLab, lim);
    if (!inRange && atLoad) return false;
    if (r.convertValue) ParticleRowJacobian(r, row, thetaLab, lim);
  }
  if (r.exitKind == ExitKind::Photon) {
    row.thetaCM = PhotonLabToCM(CMVelocity(r.entrance, eLab), thetaLab);
    if (r.convertValue) row.jacobian = PhotonJacobian(r.entrance, eLab, thetaLab);
  }
  return inRange;
}

}  // namespace experiment
}  // namespace azure2
