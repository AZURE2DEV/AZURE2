#include "ThmOptical.h"

#include <algorithm>
#include <cmath>

/*
 * The formulas, transcribed from the papers as compiled in RIPL-3
 * (om-parameter-u.dat entries 6200, 2405/5405, 7100/8100, 9100, 9600; Capote
 * et al., Nucl. Data Sheets 110 (2009) 3107) and the FRONT front end of
 * TWOFNR (J.A. Tostevin, Surrey: Daehnick, Koning-Delaroche, Liang); the
 * independent check is tests/reference/thm_optical_reference.py, which reads
 * the RIPL-3 file itself.  E is the projectile lab energy, A, Z, N those of
 * the target, eta = (N - Z)/A.
 */

namespace {

const std::vector<ThmGlobalOptical> kModels = {
    {"ancai06", "d", "An & Cai, PRC 73 (2006) 054605", 12.0, 238.0, 0.0, 183.0},
    {"daehnick80", "d", "Daehnick, Childs & Vrcelj, PRC 21 (1980) 2253 (global set)", 27.0, 238.0, 11.8, 90.0},
    {"kd03", "n p", "Koning & Delaroche, NPA 713 (2003) 231 (global)", 24.0, 209.0, 0.001, 200.0},
    {"bg71", "t 3He", "Becchetti & Greenlees (1971)", 40.0, 208.0, 1.0, 40.0},
    {"liang09", "3He", "Liang, Li & Cai, J. Phys. G 36 (2009) 085104", 9.0, 208.0, 0.0, 270.0},
    {"mcfadden66", "4He", "McFadden & Satchler, NPA 84 (1966) 177", 16.0, 208.0, 1.0, 25.0},
    {"avrigeanu94", "4He", "Avrigeanu, Hodgson & Avrigeanu, PRC 49 (1994) 2136", 16.0, 250.0, 1.0, 73.0},
};

enum Species { NEUTRON, PROTON, DEUTERON, TRITON, HELION, ALPHA, OTHER };

Species Kind(int Z, int A) {
  if (Z == 0 && A == 1) return NEUTRON;
  if (Z == 1 && A == 1) return PROTON;
  if (Z == 1 && A == 2) return DEUTERON;
  if (Z == 1 && A == 3) return TRITON;
  if (Z == 2 && A == 3) return HELION;
  if (Z == 2 && A == 4) return ALPHA;
  return OTHER;
}

}  // namespace

const std::vector<ThmGlobalOptical> &ThmGlobalOpticals() { return kModels; }

int ThmGlobalOpticalIndex(const std::string &name) {
  for (size_t i = 0; i < kModels.size(); i++)
    if (name == kModels[i].name) return (int)i;
  return -1;
}

std::string ThmGlobalOpticalNames() {
  std::string out;
  for (const ThmGlobalOptical &m : kModels) out += (out.empty() ? "" : ", ") + std::string(m.name);
  return out;
}

bool ThmGlobalOpticalFor(int index, int Zp, int Ap) {
  const Species s = Kind(Zp, Ap);
  switch (index) {
    case 0:
    case 1:
      return s == DEUTERON;
    case 2:
      return s == NEUTRON || s == PROTON;
    case 3:
      return s == TRITON || s == HELION;
    case 4:
      return s == HELION;
    case 5:
    case 6:
      return s == ALPHA;
    default:
      return false;
  }
}

// The projectile enters by its charge Zp; its mass number does not.
void ThmGlobalOpticalEvaluate(int index, int Zp, int /*Ap*/, int Zt, int At, double elab, double p[10]) {
  for (int k = 0; k < 10; k++) p[k] = 0.0;
  const double A = At, Z = Zt, N = At - Zt, E = elab;
  const double a13 = std::cbrt(A), eta = (N - Z) / A;
  // Depth, reduced radius, diffuseness of the three terms; reduced Coulomb radius.
  double v = 0, rv = 0, av = 0, w = 0, rw = 0, aw = 0, wd = 0, rd = 0, ad = 0, rc = 0;
  switch (index) {
    case 0: {  // An & Cai 2006 (spin-orbit 3.557, 0.972, 1.011 dropped)
      v = 91.85 - 0.249 * E + 0.000116 * E * E + 0.642 * Z / a13;
      rv = 1.152 - 0.00776 / a13;
      av = 0.719 + 0.0126 * a13;
      w = 1.104 + 0.0622 * E;
      rw = 1.305 + 0.0997 / a13;
      aw = 0.855 - 0.100 * a13;
      wd = 10.83 - 0.0306 * E;
      rd = 1.334 + 0.152 / a13;
      ad = 0.531 + 0.062 * a13;
      rc = 1.303;
      break;
    }
    case 1: {  // Daehnick et al. 1980, the global set as FRONT21 codes it (spin-orbit dropped)
      const double g = std::exp(-(E / 100.0) * (E / 100.0));
      v = 88.5 - 0.26 * E + 0.88 * Z / a13;
      rv = 1.17;
      av = 0.709 + 0.0017 * E;
      w = (12.2 + 0.026 * E) * (1.0 - g);
      wd = (12.2 + 0.026 * E) * g;
      rw = rd = 1.325;
      double ai = 0.53 + 0.07 * a13;
      for (double magic : {8.0, 20.0, 28.0, 50.0, 82.0, 126.0}) ai -= 0.04 * std::exp(-std::pow((magic - N) / 2.0, 2));
      aw = ad = ai;
      rc = 1.30;
      break;
    }
    case 2: {  // Koning & Delaroche 2003, global (spin-orbit dropped)
      const bool proton = Zp == 1;
      const double ef = proton ? -8.4075 + 0.01378 * A : -11.2814 + 0.02646 * A;
      const double v1 = proton ? 59.30 + 21.0 * eta - 0.024 * A : 59.30 - 21.0 * eta - 0.024 * A;
      const double v2 = proton ? 0.007067 + 4.23e-6 * A : 0.007228 - 1.48e-6 * A;
      const double v3 = proton ? 1.729e-5 + 1.136e-8 * A : 1.994e-5 - 2.0e-8 * A;
      const double v4 = 7.0e-9;
      const double w1 = proton ? 14.667 + 0.009629 * A : 12.195 + 0.0167 * A;
      const double w2 = 73.55 + 0.0795 * A;
      const double d1 = proton ? 16.0 + 16.0 * eta : 16.0 - 16.0 * eta;
      const double d2 = 0.0180 + 0.003802 / (1.0 + std::exp((A - 156.0) / 8.0));
      const double d3 = 11.5;
      const double f = E - ef;
      v = v1 * (1.0 - v2 * f + v3 * f * f - v4 * f * f * f);
      if (proton) {
        rc = 1.198 + 0.697 / (a13 * a13) + 12.994 / std::pow(a13, 5);
        const double vc = 1.73 / rc * Z / a13;
        v += vc * v1 * (v2 - 2.0 * v3 * f + 3.0 * v4 * f * f);
      }
      rv = rw = 1.3039 - 0.4054 / a13;
      av = aw = 0.6778 - 1.487e-4 * A;
      w = w1 * f * f / (f * f + w2 * w2);
      wd = d1 * f * f * std::exp(-d2 * f) / (f * f + d3 * d3);
      rd = 1.3424 - 0.01585 * a13;
      ad = proton ? 0.5187 + 5.205e-4 * A : 0.5446 - 1.656e-4 * A;
      break;
    }
    case 3: {  // Becchetti & Greenlees 1971 (spin-orbit 2.5 dropped)
      const bool helion = Zp == 2;
      v = helion ? 151.9 - 0.17 * E + 50.0 * eta : 165.0 - 0.17 * E - 6.4 * eta;
      rv = 1.20;
      av = 0.72;
      w = helion ? 41.7 - 0.33 * E + 44.0 * eta : 46.0 - 0.33 * E - 110.0 * eta;
      rw = 1.40;
      aw = helion ? 0.88 : 0.84;
      rc = 1.30;
      break;
    }
    case 4: {  // Liang, Li & Cai 2009 (spin-orbit dropped)
      v = 118.36 - 0.2071 * E + 6.3961e-5 * E * E + 26.001 * eta + 0.5668 * Z / a13;
      rv = 1.1657 + 0.0401 / a13;
      av = 0.6641 + 0.0305 * a13;
      w = -6.8871 + 0.3115 * E - 6.8096e-4 * E * E;
      rw = 1.4022 + 0.0418 / a13;
      aw = 0.7732 + 0.0219 * a13;
      wd = 20.119 - 0.1626 * E - 5.4067 * eta + 1.2087 * a13;
      rd = 1.1802 + 0.0587 / a13;
      ad = 0.6292 + 0.0657 * a13;
      rc = 1.289;
      break;
    }
    case 5: {  // McFadden & Satchler 1966, energy independent
      v = 185.0;
      w = 25.0;
      rv = rw = 1.4;
      av = aw = 0.52;
      rc = 1.4;
      break;
    }
    case 6: {  // Avrigeanu, Hodgson & Avrigeanu 1994
      v = 101.1 - 0.248 * E + 6.052 * Z / a13;
      rv = 1.245;
      av = 0.817 - 0.0085 * a13;
      w = 12.64 + 0.2 * E - 1.706 * a13;
      rw = 1.57;
      aw = 0.692 - 0.02 * a13;
      rc = 1.245;
      break;
    }
    default:
      return;
  }
  p[0] = v;
  p[1] = rv * a13;
  p[2] = av;
  p[3] = std::max(0.0, w);
  p[4] = rw * a13;
  p[5] = aw;
  p[6] = std::max(0.0, wd);
  p[7] = rd * a13;
  p[8] = ad;
  p[9] = rc * a13;
  // A term that is off carries no geometry (as the ten numbers would be written).
  for (int t = 1; t < 3; t++)
    if (p[3 * t] == 0.0) p[3 * t + 1] = p[3 * t + 2] = 0.0;
}
