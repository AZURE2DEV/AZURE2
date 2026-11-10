#ifndef AZURE2_EXPERIMENT_FRAMES_H
#define AZURE2_EXPERIMENT_FRAMES_H

#include <string>

namespace azure2 {
namespace experiment {

/*!
 * Lab -> c.m. conversion of a measured point (docs/upgrade/data_model.md section 1.4, step S6).  Input data are
 * lab-frame, forward kinematics: particle 1 of the entrance partition (the beam) of lab kinetic energy E_lab (MeV)
 * on particle 2 at rest; for lab-frame differential data the detected particle is particle 1 of the exit partition
 * (a photon for a capture exit), at lab angle theta_lab (deg).  The functions give the c.m. kinetic energy E_c.m.
 * of the entrance partition (MeV), the excitation energy E_x = E_c.m. + S + E* of the compound nucleus (MeV), the
 * c.m. angle theta_c.m. (deg) and the Jacobian dOmega_lab / dOmega_c.m. (dimensionless) that takes a lab
 * differential cross section to the c.m. one, dsigma/dOmega_c.m. = J dsigma/dOmega_lab.
 *
 * Non-relativistic forms (masses in u): E_c.m. = E_lab m2 / (m1 + m2); the elastic angle theta_c.m. = theta_lab +
 * asin((m1/m2) sin theta_lab); the reaction angle by Marion and Young's a13, a24 (Nuclear Reaction Analysis, 1968);
 * J = sin^2 theta_lab / sin^2 theta_c.m. cos(theta_c.m. - theta_lab), and its limit at 0 and 180 deg.  Relativistic
 * forms (masses in MeV/c^2; "relativistic=1"): s = (m1 + m2)^2 + 2 m2 E_lab, E_c.m. = sqrt(s) - m1 - m2, and the
 * boost of the detected particle (RelativisticLabToCM).  Photons: the aberration cos theta_c.m. = (cos theta_lab -
 * beta) / (1 - beta cos theta_lab), beta the c.m. velocity, and J = (1 - beta cos theta_lab)^2 / (1 - beta^2).
 *
 * Every expression is the engine's (src/EPoint.cpp before this file, include/Kinematics.h), operation for
 * operation, so the engine's results do not move by a bit.  Only the forward branch of a double-valued map is
 * modelled: a lab angle beyond the kinematic maximum is reported (AngleLimit), not converted.
 */

/// A partition as the conversions read it (the engine's PPair): masses (u), and the separation energy S and
/// excitation energy E* of the pair (MeV); its threshold in the compound nucleus is S + E*.
struct Partition {
  double m1 = 0.0, m2 = 0.0;
  double sepE = 0.0, exE = 0.0;
};

/// The pi of the engine (include/Constants.h), so that degree conversions are the engine's bits.
constexpr double kPi = 3.14159265358979323846;

// ---------------------------------------------------------------- energies

/// Relativistic c.m. kinetic energy from the lab kinetic energy of particle 1 on particle 2 at rest (masses MeV/c^2).
double CMFromLab(double m1_MeV, double m2_MeV, double eLab);
/// Lab kinetic energy of particle 1 (particle 2 at rest) for c.m. kinetic energy T.
double LabFromCM(double m1_MeV, double m2_MeV, double T);
/// c.m. kinetic energy from the kinetic energy of particle 1 of a two-body decay at rest.
double CMFromDecay(double m1_MeV, double m2_MeV, double e1);

/// m2 / (m1 + m2): E_c.m. / E_lab of a beam (non-relativistic).
inline double BeamCMFactor(const Partition &p) { return p.m2 / (p.m1 + p.m2); }
/// (m1 + m2) / m2: E_c.m. / E_1 of a decay at rest, E_1 the light particle's energy (non-relativistic).
inline double DecayCMFactor(const Partition &p) { return (p.m1 + p.m2) / p.m2; }
/// E_x = E_c.m. + S + E* of the partition.
inline double ExcitationEnergy(double eCM, const Partition &p) { return eCM + p.sepE + p.exE; }

struct BeamEnergiesCM {
  double eCM, eX, binLowCM, binHighCM;
};
/// E_c.m., E_x and the c.m. bin edges of a beam point of lab energy eLab and bin edges binLowLab, binHighLab.
BeamEnergiesCM BeamToCM(const Partition &entrance, double eLab, double binLowLab, double binHighLab, bool relativistic);
/// The lab energy of a beam point whose E_c.m. is eCM: the inverse of BeamToCM's map (non-relativistic eCM (m1 + m2)/m2,
/// relativistic LabFromCM), for an energy grid given in the c.m. frame.  BeamToCM gives eCM back to rounding.
double BeamLabFromCM(const Partition &entrance, double eCM, bool relativistic);
/// E_c.m. of a beta-delayed point: eLab is the light particle's lab energy, the partition the decay's (exit) pair.
double DecayToCM(const Partition &exit, double eLab, bool relativistic);

// ---------------------------------------------------------------- angles

/// What an angle map reports besides the angle: the maximum lab angle (deg), if it set one, and why a lab angle was
/// not reached (empty if it was).
struct AngleLimit {
  bool maxSet = false;
  double max = 0.0;
  std::string why;
};

/// Elastic scattering, non-relativistic: theta_c.m. from theta_lab (deg).  With m1 > m2 the maximum lab angle
/// asin(m2/m1) is set, and false is returned beyond it.
bool ElasticLabToCM(const Partition &p, double thetaLab, double &thetaCM, AngleLimit &lim);
/// A reaction (entrance -> exit), non-relativistic (Marion-Young): theta_c.m. of the exit's particle 1 at lab angle
/// theta_lab and beam energy eLab.  False below the exit threshold and beyond the maximum lab angle.
bool ReactionLabToCM(const Partition &entrance, const Partition &exit, double eLab, double thetaLab, double &thetaCM,
                     AngleLimit &lim);
/// Relativistic theta_c.m. and dOmega_lab/dOmega_c.m. of the exit's particle 1 at lab angle theta_lab, for c.m.
/// energy eCM of the entrance partition; where two c.m. angles reach the lab angle, the one nearest nrCMAngle.  False
/// below the exit threshold and beyond the maximum lab angle (smaller than the non-relativistic one).
bool RelativisticLabToCM(const Partition &entrance, const Partition &exit, double eCM, double thetaLab,
                         double nrCMAngle, double &thetaCM, double &jacobian, AngleLimit &lim);

/// The c.m. velocity beta (units of c) of a beam of lab energy eLab on particle 2 at rest.
double CMVelocity(const Partition &entrance, double eLab);
/// Photon aberration: theta_c.m. (deg) of a photon at theta_lab emitted from a source moving with beta.
double PhotonLabToCM(double beta, double thetaLab);

// ---------------------------------------------------------------- Jacobians dOmega_lab / dOmega_c.m.

/// Particle, non-relativistic, at the lab and c.m. angles of the point (0 and 180 deg by the limit of a13 ... d).
double ParticleJacobian(const Partition &entrance, const Partition &exit, double eLab, double thetaLab, double thetaCM);
/// Photon: (1 - beta cos theta_lab)^2 / (1 - beta^2), beta = CMVelocity.
double PhotonJacobian(const Partition &entrance, double eLab, double thetaLab);

}  // namespace experiment
}  // namespace azure2

#endif
