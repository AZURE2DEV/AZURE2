#include "CoulFunc.h"
#include "CoulFuncCache.h"
#include "PPair.h"
#include <iostream>
#include "cwfcomp_accurate.H"
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_deriv.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include "NuclearPotential.h"
#include "NuclearPotentialManager.h"
#include "Config.h"
 
/*! 
 * The CoulFunc object is created with reference to a PPair object.
 */

CoulFunc::CoulFunc(PPair *pPair, bool useGSLFunctions) :
  useGSLFunctions_(useGSLFunctions),
  useHybridMethod_(g_config ? g_config->useHybridMethod : false), // Read from global Config
  rMatch_(15.0),
  drNumerov_(0.0275) {
  z1_=pPair->GetZ(1);
  z2_=pPair->GetZ(2);
  redmass_=(double)pPair->GetRedMass();
  lLast_=0;
  radiusLast_=0.0;
  energyLast_=0.0;
  coulLast_.F=0.0;
  coulLast_.dF=0.0;
  coulLast_.G=0.0;
  coulLast_.dG=0.0;
  dEShiftParams_.coulFunc=this;

  // Get nuclear potential from the global manager (always synchronized with GUI)
  nuclearPotential_ = NuclearPotentialManager::instance().getPotential();
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
 
void CoulFunc::setLast(int lLast, double rLast, double eLast, CoulWaves coulLast) {
  lLast_=lLast;
  radiusLast_=rLast;
  energyLast_=eLast;
  coulLast_.F=coulLast.F;
  coulLast_.dF=coulLast.dF;
  coulLast_.G=coulLast.G;
  coulLast_.dG=coulLast.dG;
}

/*!
 * The parenthesis operator is defined to make the class instance callable as a function.  The orbital
 * angular momentum, radius, and energy in the center of mass system are the dependent variables.
 * The function returns the Coulomb waves.
 */

CoulWaves CoulFunc::operator()(int l,double radius,double energy) {
  struct CoulWaves result={0.0,0.0,0.0,0.0};

  // FIXME: Cross section might be discontinuous at 0 energy, so set a minimum
  // 12C(a,g) has problem with this
  if(energy<=0.001) {
    energy = 0.001;
  }

  // First check the local cache (last calculated value)
  if(l==lLast()&&radius==radiusLast()&&energy==energyLast()) {
    result=coulLast();
    return result;
  }

  // Try to use the global cache if available
  if(g_coulFuncCache) {
    CoulFuncCache::CoulFuncKey key;
    key.z1 = z1();
    key.z2 = z2();
    key.redmass = redmass();
    key.l = l;
    key.radius = radius;

    // Check if we have cached data for this parameter set
    if(g_coulFuncCache->IsInRange(key, energy)) {
      result = g_coulFuncCache->GetInterpolatedCoulWaves(key, energy);
      setLast(l, radius, energy, result);
      return result;
    }
  }

  // If cache miss, compute directly
  struct CoulWaves newResult;

  // Check if hybrid method should be used
  if(useHybridMethod_ && nuclearPotential_) {
    newResult = computeHybrid(l, radius, energy);
  } else if(!useGSLFunctions_) {
    std::complex<double> eta(sqrt(uconv/2.)*fstruc*z1()*z2()*
			     sqrt(redmass()/energy),0.);
    std::complex<double> rho(sqrt(2.*uconv)/hbarc*radius*
			     sqrt(redmass()*energy),0.);
    std::complex<double> lValue( (double) l, 0.);
    Coulomb_wave_functions_accurate coul(true,lValue,eta);
    std::complex<double> c_F, c_dF, c_G, c_dG;
    //coul.F_dF(rho,c_F,c_dF);
    //coul.G_dG(rho,c_G,c_dG);
    coul.compute(rho,c_F,c_dF,c_G,c_dG);
    newResult.F=real(c_F);
    newResult.dF=real(c_dF);
    newResult.G=real(c_G);
    newResult.dG=real(c_dG);
  } else {
    double eta=sqrt(uconv/2.)*fstruc*z1()*z2()*
      sqrt(redmass()/energy);
    double rho=sqrt(2.*uconv)/hbarc*radius*
      sqrt(redmass()*energy);
    double lValue=double(l);
    double eF,eG;
    gsl_sf_result F,Fp,G,Gp;
    gsl_sf_coulomb_wave_FG_e(eta,rho,lValue,0,&F,&Fp,&G,&Gp,&eF,&eG);
    newResult.F=F.val*exp(eF);
    newResult.dF=Fp.val*exp(eF);
    newResult.G=G.val*exp(eG);
    newResult.dG=Gp.val*exp(eG);
  }
  setLast(l,radius,energy,newResult);
  result=newResult;

  // Add computed result to global cache for future use
  if(g_coulFuncCache) {
    CoulFuncCache::CoulFuncKey key;
    key.z1 = z1();
    key.z2 = z2();
    key.redmass = redmass();
    key.l = l;
    key.radius = radius;
    g_coulFuncCache->AddCoulWaves(key, energy, newResult);
  }

  return result;
}

/*!
 * Returns the penetrability as a function of orbital
 * angular momentum, radius, and energy in the center of mass system.
 */

double CoulFunc::Penetrability(int l,double radius,double energy)  {
  struct CoulWaves coul=this->operator()(l,radius,energy);
  double rho=sqrt(2.*uconv)/hbarc*radius*sqrt(redmass()*energy);
  return rho/(pow(coul.F,2.0)+pow(coul.G,2.0));
}

/*!
 * Returns the positive energy shift function a function of orbital
 * angular momentum, radius, and energy in the center of mass system.
 */

double CoulFunc::PEShift(int l,double radius,double energy)  {
  struct CoulWaves coul=this->operator()(l,radius,energy);
  double rho=sqrt(2.*uconv)/hbarc*radius*sqrt(redmass()*energy);
  if(pow(coul.F,2.0)==0.&&coul.F*coul.dF==0.) return rho*(coul.dG/coul.G);
  else return rho/(pow(coul.F,2.0)+pow(coul.G,2.0))*
    (coul.F*coul.dF+coul.G*coul.dG);
}

double CoulFunc::thisPEShift(double x, void *p) {
  DEShiftParams *params = (DEShiftParams*)p;
  CoulFunc *coulFunc=(params->coulFunc);
  int lValue=(params->lValue);
  double radius=(params->radius);

  return coulFunc->PEShift(lValue,radius,x);
}

/*!
 * Returns the energy derivative of the shift function a function of orbital
 * angular momentum, radius, and energy in the center of mass system.
 */

double CoulFunc::PEShift_dE(int l,double radius,double energy)  {
  double result;
  double error;

  dEShiftParams_.radius=radius;
  dEShiftParams_.lValue=l;

  gsl_function F;
  F.function=&thisPEShift;
  F.params=&dEShiftParams_;

  gsl_deriv_central (&F, energy, 1e-6, &result, &error);

  return result;
}

/*!
 * Hybrid method: Numerov integration with nuclear potential using GSL or COUL boundary conditions.
 * This method integrates the Schrödinger equation inward from a matching radius where pure
 * Coulomb boundary conditions are applied, through the nuclear potential region.
 */

CoulWaves CoulFunc::computeHybrid(int l, double radius, double energy) {
  struct CoulWaves result={0.0,0.0,0.0,0.0};

  // Physical constants (consistent with pycoulomb implementation)
  const double UCONV = 931.494;      // MeV/amu
  const double FSTRUC = 0.00729735;  // fine structure constant * hc (MeV*fm)
  const double HBARC = 197.327;      // hbar*c in MeV*fm

  try {
    // Validate inputs
    const double r_min = 0.0;
    if(rMatch_ <= r_min || radius > rMatch_ || radius < r_min || drNumerov_ <= 0.0) {
      // Fall back to standard Coulomb calculation
      if(!useGSLFunctions_) {
        std::complex<double> eta(sqrt(uconv/2.)*fstruc*z1()*z2()*
                                 sqrt(redmass()/energy),0.);
        std::complex<double> rho(sqrt(2.*uconv)/hbarc*radius*
                                 sqrt(redmass()*energy),0.);
        std::complex<double> lValue((double)l, 0.);
        Coulomb_wave_functions_accurate coul(true,lValue,eta);
        std::complex<double> c_F, c_dF, c_G, c_dG;
        coul.compute(rho,c_F,c_dF,c_G,c_dG);
        result.F=real(c_F);
        result.dF=real(c_dF);
        result.G=real(c_G);
        result.dG=real(c_dG);
      } else {
        double eta=sqrt(uconv/2.)*fstruc*z1()*z2()*sqrt(redmass()/energy);
        double rho=sqrt(2.*uconv)/hbarc*radius*sqrt(redmass()*energy);
        double lValue=double(l);
        double eF,eG;
        gsl_sf_result F,Fp,G,Gp;
        gsl_sf_coulomb_wave_FG_e(eta,rho,lValue,0,&F,&Fp,&G,&Gp,&eF,&eG);
        result.F=F.val*exp(eF);
        result.dF=Fp.val*exp(eF);
        result.G=G.val*exp(eG);
        result.dG=Gp.val*exp(eG);
      }
      return result;
    }

    // Setup: Constants and unit conversions
    std::complex<double> lValue(static_cast<double>(l), 0.0);
    double redmass_MeV = redmass() * UCONV;
    double eta_real = z1() * z2() * FSTRUC * std::sqrt(redmass_MeV / 2.0) / std::sqrt(energy);
    std::complex<double> eta(eta_real, 0.0);
    double drho_dr = std::sqrt(2.0 * UCONV) / HBARC * std::sqrt(redmass() * energy);

    // Lambda to compute H+ = G + iF from rho using Coulomb boundary conditions
    auto compute_H_plus = [&](double rho_val) -> std::complex<double> {
      if(!useGSLFunctions_) {
        // Use COUL (coulsm) for boundary conditions
        Coulomb_wave_functions_accurate coul(true, lValue, eta);
        std::complex<double> F, dF, G, dG;
        coul.compute(std::complex<double>(rho_val, 0.0), F, dF, G, dG);
        return G + std::complex<double>(0.0, 1.0) * F;
      } else {
        // Use GSL for boundary conditions
        double eta_val = eta.real();
        double lValue_val = static_cast<double>(l);
        double eF, eG;
        gsl_sf_result F, Fp, G, Gp;
        gsl_sf_coulomb_wave_FG_e(eta_val, rho_val, lValue_val, 0, &F, &Fp, &G, &Gp, &eF, &eG);
        double F_val = F.val * std::exp(eF);
        double G_val = G.val * std::exp(eG);
        return std::complex<double>(G_val, 0.0) + std::complex<double>(0.0, 1.0) * std::complex<double>(F_val, 0.0);
      }
    };

    // Lambda to compute rho from radius
    auto r_to_rho = [&](double r) -> double {
      return std::sqrt(2.0 * redmass_MeV) / HBARC * r * std::sqrt(energy);
    };

    // Create grid and set boundary conditions
    int nmax = static_cast<int>((rMatch_ - r_min) / drNumerov_) + 1;
    if(nmax % 4 != 0) nmax = nmax - (nmax % 4) + 4;
    int npts = nmax + 1;

    std::vector<double> r_grid(npts);
    std::vector<std::complex<double>> u(npts);

    for(int i = 0; i < npts; ++i) {
      r_grid[i] = r_min + i * drNumerov_;
      if(i >= npts - 2) {
        u[i] = compute_H_plus(r_to_rho(r_grid[i]));
      }
    }

    // Create combined potential (Coulomb + Nuclear + Centrifugal)
    if(!nuclearPotential_) {
      std::cerr << "Error: Nuclear potential not set for Numerov integration." << std::endl;
      return result;
    }

    auto combined_potential = std::make_shared<CoulombNuclearPotential>(
        static_cast<int>(z1()), static_cast<int>(z2()),
        *nuclearPotential_,
        l,                  // orbital angular momentum
        redmass_MeV);       // reduced mass in MeV/c^2

    // Compute effective g(r) = 2*mu/hbar^2 * (E - V(r))
    std::vector<double> g(npts);
    for(int i = 0; i < npts; ++i) {
      g[i] = 2.0 * redmass_MeV / (HBARC * HBARC) * (energy - (*combined_potential)(r_grid[i]));
    }

    // Numerov backward integration from matching radius inward
    double dr2 = drNumerov_ * drNumerov_;
    for(int i = npts - 3; i >= 0; --i) {
      double denom = 1.0 + (dr2 / 12.0) * g[i];
      if(std::abs(denom) > 1e-20) {
        u[i] = ((2.0 - (5.0 * dr2 / 6.0) * g[i+1]) * u[i+1] -
                (1.0 + (dr2 / 12.0) * g[i+2]) * u[i+2]) / denom;
      } else {
        u[i] = u[i+1];
      }
    }

    // Linear interpolation at requested radius
    int idx_lower = 0;
    for(int i = 0; i < npts - 1; ++i) {
      if(r_grid[i] <= radius && radius <= r_grid[i+1]) {
        idx_lower = i;
        break;
      }
    }

    double frac = (radius - r_grid[idx_lower]) / (r_grid[idx_lower+1] - r_grid[idx_lower]);
    std::complex<double> H_plus = u[idx_lower] + frac * (u[idx_lower+1] - u[idx_lower]);
    std::complex<double> dH_dr = (u[idx_lower+1] - u[idx_lower]) / drNumerov_;

    // Extract F, G and their derivatives from H+ = G + iF
    result.G = H_plus.real();
    result.F = H_plus.imag();
    result.dG = dH_dr.real() / drho_dr;
    result.dF = dH_dr.imag() / drho_dr;

  } catch(const std::exception& e) {
    std::cerr << "Error in hybrid method: " << e.what() << std::endl;
    std::cerr << "Falling back to standard Coulomb calculation." << std::endl;
    // Fall back to standard calculation
    if(!useGSLFunctions_) {
      std::complex<double> eta(sqrt(uconv/2.)*fstruc*z1()*z2()*
                               sqrt(redmass()/energy),0.);
      std::complex<double> rho(sqrt(2.*uconv)/hbarc*radius*
                               sqrt(redmass()*energy),0.);
      std::complex<double> lValue((double)l, 0.);
      Coulomb_wave_functions_accurate coul(true,lValue,eta);
      std::complex<double> c_F, c_dF, c_G, c_dG;
      coul.compute(rho,c_F,c_dF,c_G,c_dG);
      result.F=real(c_F);
      result.dF=real(c_dF);
      result.G=real(c_G);
      result.dG=real(c_dG);
    } else {
      double eta=sqrt(uconv/2.)*fstruc*z1()*z2()*sqrt(redmass()/energy);
      double rho=sqrt(2.*uconv)/hbarc*radius*sqrt(redmass()*energy);
      double lValue=double(l);
      double eF,eG;
      gsl_sf_result F,Fp,G,Gp;
      gsl_sf_coulomb_wave_FG_e(eta,rho,lValue,0,&F,&Fp,&G,&Gp,&eF,&eG);
      result.F=F.val*exp(eF);
      result.dF=Fp.val*exp(eF);
      result.G=G.val*exp(eG);
      result.dG=Gp.val*exp(eG);
    }
  }

  return result;
}

// Hybrid method configuration methods

/*!
 * Set the nuclear potential for hybrid method calculations.
 */
void CoulFunc::setNuclearPotential(std::shared_ptr<NuclearPotential> potential) {
  nuclearPotential_ = potential;
}

/*!
 * Get the current nuclear potential.
 */
std::shared_ptr<NuclearPotential> CoulFunc::getNuclearPotential() const {
  return nuclearPotential_;
}

/*!
 * Set the matching radius where Coulomb boundary conditions are applied.
 */
void CoulFunc::setMatchingRadius(double r_match) {
  if(r_match > 0.0) {
    rMatch_ = r_match;
  }
}

/*!
 * Get the current matching radius.
 */
double CoulFunc::getMatchingRadius() const {
  return rMatch_;
}

/*!
 * Set the Numerov grid step size.
 */
void CoulFunc::setNumerovGridStep(double dr) {
  if(dr > 0.0) {
    drNumerov_ = dr;
  }
}

/*!
 * Get the current Numerov grid step.
 */
double CoulFunc::getNumerovGridStep() const {
  return drNumerov_;
}

/*!
 * Enable/disable hybrid method for Coulomb wave calculations.
 */
void CoulFunc::setUseHybridMethod(bool useHybrid) {
  useHybridMethod_ = useHybrid;
}

/*!
 * Check if hybrid method is enabled.
 */
bool CoulFunc::getUseHybridMethod() const {
  return useHybridMethod_;
}

