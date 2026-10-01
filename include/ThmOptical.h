#ifndef THM_OPTICAL_H
#define THM_OPTICAL_H

#include <string>
#include <vector>

/*!
 * Global optical-model parametrisations for the distorted waves of a THM
 * experiment: `opticalAA=<name>` / `opticalSF=<name>` (optionally
 * `<name>:extrapolate`), docs/source/theory/thm_implementation.rst, "Global
 * optical potentials".
 *
 * Each model gives the ten numbers of the AZURE2 Woods-Saxon form
 *
 *   U(r) = V_C(r) - V f(r; R, a) - i W f(r; R_W, a_W)
 *          - 4 i W_D e^x/(1 + e^x)^2,   x = (r - R_D)/a_D,
 *
 * (V,R,a, W,RW,aW, WD,RD,aD, RC: MeV and fm, radii R = r0 A_t^(1/3) with the
 * target mass number A_t; RC a uniform sphere) for a light projectile
 * (Zp, Ap) on a target (Zt, At) at the projectile's LAB energy.  The surface
 * term is the derivative form normalized as -4 a_D W_D df/dr, which is the
 * convention of all the papers below, so W_D maps one to one.  Spin-orbit
 * terms are dropped (the distorted waves have no spin).  A negative
 * imaginary depth (Liang's volume term below 23 MeV) is set to 0.
 *
 *   ancai06      d       An & Cai, PRC 73 (2006) 054605          A 12-238,  E 0-183 MeV
 *   daehnick80   d       Daehnick, Childs & Vrcelj, PRC 21 (1980) 2253,
 *                        the global set of FRONT/TWOFNR           A 27-238,  E 11.8-90
 *   kd03         n, p    Koning & Delaroche, NPA 713 (2003) 231,
 *                        global set                               A 24-209,  E 0.001-200
 *   bg71         t, 3He  Becchetti & Greenlees (1971)             A 40-208,  E 1-40
 *   liang09      3He     Liang, Li & Cai, J. Phys. G 36 (2009) 085104
 *                                                                 A 9-208,   E 0-270
 *   mcfadden66   4He     McFadden & Satchler, NPA 84 (1966) 177   A 16-208,  E 1-25
 *   avrigeanu94  4He     Avrigeanu, Hodgson & Avrigeanu, PRC 49 (1994) 2136
 *                                                                 A 16-250,  E 1-73
 */
struct ThmGlobalOptical {
  const char *name;
  const char *projectiles;  ///< "d", "n p", ...: the light ions the model is for
  const char *reference;
  double aMin, aMax;        ///< target mass number
  double eMin, eMax;        ///< projectile lab energy (MeV)
};

/// The built-in models, in the order the documentation lists them.
const std::vector<ThmGlobalOptical> &ThmGlobalOpticals();
/// Index of a model by name, or -1.
int ThmGlobalOpticalIndex(const std::string &name);
/// The names, separated by ", ".
std::string ThmGlobalOpticalNames();
/// Whether model `index` is a potential for the projectile (Zp, Ap).
bool ThmGlobalOpticalFor(int index, int Zp, int Ap);
/// V,R,a,W,RW,aW,WD,RD,aD,RC (MeV, fm) of model `index` for (Zp, Ap) on
/// (Zt, At) at the projectile lab energy elab (MeV).  The formulas are
/// evaluated wherever asked; the validity range is the caller's business.
void ThmGlobalOpticalEvaluate(int index, int Zp, int Ap, int Zt, int At, double elab, double p[10]);

#endif
