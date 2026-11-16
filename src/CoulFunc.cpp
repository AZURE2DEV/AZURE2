#include "CoulFunc.h"
#include "PPair.h"
#include <iostream>
#include "azure2/physics/Derivatives.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <complex>
#include <unordered_map>
#include "Config.h"
#include "CoulombTable.h"
#include "Kinematics.h"

/*!
 * The CoulFunc object is created with reference to a PPair object.
 */

CoulFunc::CoulFunc(const PPair *pPair) :
  pair_(pPair),
  potential_(pPair->GetHybridPotential()) {
  z1_ = pPair->GetZ(1);
  z2_ = pPair->GetZ(2);
  redmass_ = (double)pPair->GetRedMass();
  kin_ = Kinematics::Channel(pPair);
  lLast_ = 0;
  radiusLast_ = 0.0;
  energyLast_ = 0.0;
  coulLast_.F = 0.0;
  coulLast_.dF = 0.0;
  coulLast_.G = 0.0;
  coulLast_.dG = 0.0;

  // The table holds plain Coulomb functions of charged channels: not neutral
  // pairs (eta = 0, already cheap) and not pairs with a nuclear potential.
  const PairFunctionSettings &settings = pPair->FunctionSettings();
  if (settings.useCoulombTable && z1_ * z2_ != 0 && !potential_)
    table_ = CoulombTable::Shared(settings.coulombTablePath, settings.coulombTableDefault);
  shiftDerivative1eV_ = settings.shiftDerivative1eV;
  tableUse_ = settings.tableUse;
}

bool CoulFunc::TableLookup(int l, double radius, double energy) {
  if (l == tLastL_ && radius == tLastR_ && energy == tLastE_) return tOk_;
  tLastL_ = l;
  tLastR_ = radius;
  tLastE_ = energy;
  const double eta = kin_.Eta(energy);
  const double rho = kin_.Rho(radius, energy);
  CoulombTable::Open o;
  tOk_ = energy > 0.0 && table_->LookupOpen(l, eta, rho, o, tableUse_);
  if (tOk_) {
    tP_ = o.P;
    tS_ = o.S;
    tPhi_ = o.phi;
  }
  return tOk_;
}

std::complex<double> CoulFunc::HardSpherePhase(int l, int two_j, double radius, double energy) {
  // The table holds pairs without a potential only, whose functions do not depend on two_j.
  if (table_ && TableLookup(l, radius, energy)) return std::complex<double>(std::cos(tPhi_), -std::sin(tPhi_));
  return CoulombHardSpherePhase(this->operator()(l, two_j, radius, energy));
}

/*!
 * Returns the atomic number of the first particle in the pair.
 */

int CoulFunc::z1() const {
  return z1_;
}

/*!
 * Returns the atomic number of the second particle in the pair.
 */

int CoulFunc::z2() const {
  return z2_;
}

/*!
 * Returns the reduced mass of the particle pair.
 */

double CoulFunc::redmass() const {
  return redmass_;
}

/*!
 * Returns the last orbital angular momentum value at which the Coulomb
 * functions were calculated.
 */

int CoulFunc::lLast() const {
  return lLast_;
}

/*!
 * Returns the last radius value at which the Coulomb
 * functions were calculated.
 */

double CoulFunc::radiusLast() const {
  return radiusLast_;
}

/*!
 * Returns the last energy value at which the Coulomb
 * functions were calculated.
 */

double CoulFunc::energyLast() const {
  return energyLast_;
}

/*!
 * Returns the last Coulomb functions which were calculated.
 */

struct CoulWaves CoulFunc::coulLast() const {
  return coulLast_;
}

/*!
 * Sets the last calculated Coulomb waves and the values for which they were calculated.
 */

void CoulFunc::setLast(int lLast, int twoJLast, double rLast, double eLast, CoulWaves coulLast) {
  lLast_ = lLast;
  twoJLast_ = twoJLast;
  radiusLast_ = rLast;
  energyLast_ = eLast;
  coulLast_.F = coulLast.F;
  coulLast_.dF = coulLast.dF;
  coulLast_.G = coulLast.G;
  coulLast_.dG = coulLast.dG;
}


/*!
 * The parenthesis operator is defined to make the class instance callable as a function.  The orbital
 * angular momentum, radius, and energy in the center of mass system are the dependent variables.
 * The function returns the Coulomb waves.
 */

CoulWaves CoulFunc::operator()(int l, int two_j, double radius, double energy) {
  struct CoulWaves result = {0.0, 0.0, 0.0, 0.0};

  // fix: 12C(a,g) has problems under 1 keV, so we put a limit
  // if(energy<=0.001) {
  //  energy = 0.001;
  //}

  // First check the local cache (last calculated value)
  if (l == lLast() && two_j == twoJLast() && radius == radiusLast() && energy == energyLast()) {
    result = coulLast();
    return result;
  }

  std::unordered_map<double, CoulWaves> *block = nullptr;
  if (radialMemo_) {
    if (radialMemo_->size >= CoulRadialMemo::kMaxRadii) {
      radialMemo_->blocks.clear();
      radialMemo_->size = 0;
    }
    block = &radialMemo_->blocks[CoulRadialMemo::Key{z1_, z2_, l, two_j, redmass_, kin_.rel ? kin_.m1 + kin_.m2 : 0.0,
                                                     energy, potential_}];
    auto w = block->find(radius);
    if (w != block->end()) {
      setLast(l, two_j, radius, energy, w->second);
      return w->second;
    }
  }

  // Otherwise computed at every call: there is no process-wide memo (docs/upgrade/coulomb.md section 3).
  struct CoulWaves newResult;

  if (potential_) {
    newResult = computeHybrid(l, two_j, radius, energy);
  } else {
    newResult = PlainWaves(l, radius, energy);
  }
  setLast(l, two_j, radius, energy, newResult);
  result = newResult;
  if (block) {
    block->emplace(radius, newResult);
    radialMemo_->size++;
  }

  return result;
}

/*!
 * Returns the penetrability as a function of orbital
 * angular momentum, radius, and energy in the center of mass system.
 */

double CoulFunc::Penetrability(int l, int two_j, double radius, double energy) {
  if (table_ && TableLookup(l, radius, energy)) return tP_;
  struct CoulWaves coul = this->operator()(l, two_j, radius, energy);
  return azure2::physics::penetrability(PhysicsWaves(coul), kin_.Rho(radius, energy));
}

/*!
 * Returns the positive energy shift function a function of orbital
 * angular momentum, radius, and energy in the center of mass system.
 */

double CoulFunc::PEShift(int l, int two_j, double radius, double energy) {
  if (table_ && TableLookup(l, radius, energy)) return tS_;
  return PEShiftExact(l, two_j, radius, energy);
}

double CoulFunc::PEShiftExact(int l, int two_j, double radius, double energy) {
  // physics/ChannelFunctions: S with the overflow guard just above a charged threshold (eta >~ 110; 11B+alpha
  // within ~1 keV, where a band of bare cross sections of 1e26-1e37 b once came from the unguarded formula).
  struct CoulWaves coul = this->operator()(l, two_j, radius, energy);
  return azure2::physics::shift(PhysicsWaves(coul), kin_.Rho(radius, energy));
}

/*!
 * Returns the energy derivative of the shift function a function of orbital
 * angular momentum, radius, and energy in the center of mass system.
 */

double CoulFunc::PEShift_dE(int l, int two_j, double radius, double energy, double step) {
  // A Richardson-extrapolated central difference with a step of 1e-3 E, accurate
  // to ~1e-10.  AZURE_SHIFT_DE_1EV (the pair's settings) restores the route used before 2026-10-19, a
  // 1 eV central difference of the exact S, which carries its rounding noise: up
  // to 3e-5 relative at a few MeV for alpha+12C.
  if (!shiftDerivative1eV_ && energy > 0.0) {
    const double h = 1e-3 * energy;
    // One route for the whole stencil: when a point of it falls outside the Coulomb table, every point is
    // taken from the exact routine.  Mixing the table (~1e-7 interpolation error) with the exact value across
    // the table's edge put a jump of order 1e-7 / (2 h) into the derivative (audit 2026-10-25, item 32).
    const bool viaTable = table_ && TableLookup(l, radius, energy + h) && TableLookup(l, radius, energy - h) &&
                          TableLookup(l, radius, energy + 0.5 * h) && TableLookup(l, radius, energy - 0.5 * h);
    // The difference itself is physics/'s (shift_derivative), with S from the route chosen here.
    return azure2::physics::shift_derivative([&](double e) { return viaTable ? PEShift(l, two_j, radius, e) : PEShiftExact(l, two_j, radius, e); }, energy);
  }

  return azure2::physics::central_derivative([&](double e) { return PEShift(l, two_j, radius, e); }, energy, step);
}

/*!
 * Returns the second energy derivative of the shift function, as the central
 * second difference of PEShift with a step of 1 keV (smaller close to the
 * threshold, where the function is only defined for positive energies).  The
 * Park level matrix needs it for the energy gradient of the overlap J.
 */

double CoulFunc::PEShift_d2E(int l, int two_j, double radius, double energy) {
  // Richardson-extrapolated central second difference, O(h^4), with a 10 keV
  // step: the Coulomb functions are good to ~1e-10, so a 1 keV step would
  // leave ~1e-4 of rounding noise in S''.  Within 4h of the threshold, where
  // the function is only defined for positive energies, a forward stencil.
  const double h = 1.0e-2;
  return azure2::physics::second_derivative([&](double e) { return PEShift(l, two_j, radius, e); }, energy, h, 1.0e-3,
                                            energy < 4.0 * h, azure2::physics::Side::Above);
}

CoulWaves CoulFunc::PlainWaves(int l, double radius, double energy) const {
  return CoulombWavesAt(l, kin_.Eta(energy), kin_.Rho(radius, energy));
}

CoulWaves CoulFunc::computeHybrid(int l, int two_j, double radius, double energy) {
  struct CoulWaves result = {0.0, 0.0, 0.0, 0.0};
  const azure2::physics::HybridPotential &potential = *potential_;

  // The engine's constants (Constants.h), as the plain Coulomb path uses them.
  const double UCONV = uconv;    // MeV/u
  const double FSTRUC = fstruc;  // the fine-structure constant
  const double HBARC = hbarc;    // hbar c in MeV fm

  try {
    // Beyond the matching radius the waves are the point-Coulomb ones.
    if (potential.matching_radius_fm <= 0.0 || radius > potential.matching_radius_fm || radius < 0.0 ||
        potential.step_fm <= 0.0)
      return PlainWaves(l, radius, energy);

    double redmass_MeV = redmass() * UCONV;
    double eta_real = z1() * z2() * FSTRUC * std::sqrt(redmass_MeV / 2.0) / std::sqrt(energy);
    double drho_dr = std::sqrt(2.0 * UCONV) / HBARC * std::sqrt(redmass() * energy);

    // H+ = G + iF of the point-Coulomb problem at r (fm), the solution beyond the matching radius.
    auto outgoing = [&](double r) -> std::complex<double> {
      const double rho_val = std::sqrt(2.0 * redmass_MeV) / HBARC * r * std::sqrt(energy);
      // physics/coulomb (coulsm) for the boundary conditions
      const CoulWaves w = CoulombWavesAt(l, eta_real, rho_val);
      return std::complex<double>(w.G, 0.0) + std::complex<double>(0.0, 1.0) * std::complex<double>(w.F, 0.0);
    };

    // H+ at the requested radius (a node of the grid) and its r-derivative; F, G and their rho-derivatives from
    // H+ = G + iF.
    const azure2::physics::RadialValue h =
        azure2::physics::IntegrateInward(potential, z1(), z2(), redmass_MeV, l, pair_->SpinOrbitFactor(l, two_j), radius, energy, outgoing);
    result.G = h.u.real();
    result.F = h.u.imag();
    result.dG = h.du_dr.real() / drho_dr;
    result.dF = h.du_dr.imag() / drho_dr;

  } catch (const std::exception &e) {
    std::cerr << "Error in hybrid method: " << e.what() << std::endl;
    std::cerr << "Falling back to standard Coulomb calculation." << std::endl;
    result = PlainWaves(l, radius, energy);
  }

  return result;
}
