#include "ThmExperiment.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace {

/*!
 * Light nuclei: AME2020 atomic masses (Wang et al., Chin. Phys. C 45 (2021)
 * 030003, from the mass excesses) minus Z electron masses plus the electron
 * binding energy (Lunney, Pearson & Thibault, RMP 75 (2003) 1021, eq. A4),
 * i.e. nuclear masses in u.
 */
const ThmNuclide kNuclides[] = {
    {"n", 0, 1, 1.0086649159},    {"p", 1, 1, 1.0072764675},    {"d", 1, 2, 2.0135532134},
    {"t", 1, 3, 3.0155007169},    {"3He", 2, 3, 3.0149322434},  {"4He", 2, 4, 4.0015061756},
    {"6Li", 3, 6, 6.0134773618},  {"7Li", 3, 7, 7.0143579087},  {"9Be", 4, 9, 9.0099891714},
    {"10B", 5, 10, 10.0101946866}, {"11B", 5, 11, 11.0065629921}, {"12C", 6, 12, 11.9967096429},
    {"13C", 6, 13, 13.0000644777}, {"14N", 7, 14, 13.9992355671}, {"15N", 7, 15, 14.9962704619},
    {"16O", 8, 16, 15.9905282131}, {"17O", 8, 17, 16.9947453497}, {"18O", 8, 18, 17.9947732059},
    {"19F", 9, 19, 18.9934689019}, {"20Ne", 10, 20, 19.9869581831}, {"23Na", 11, 23, 22.9837396828},
    {"24Mg", 12, 24, 23.9784646239},
};

// A whole token as a number (operator>> that must consume everything).
bool ReadWholeDouble(const std::string &text, double &x) {
  std::istringstream s(text);
  std::string rest;
  return !!(s >> x) && !(s >> rest) && std::isfinite(x);
}
bool ReadWholeInt(const std::string &text, int &x) {
  if (text.empty()) return false;
  for (size_t i = 0; i < text.size(); i++)
    if (!std::isdigit((unsigned char)text[i]) && !(i == 0 && (text[i] == '+' || text[i] == '-'))) return false;
  std::istringstream s(text);
  return !!(s >> x);
}

std::vector<std::string> Split(const std::string &text, char sep) {
  std::vector<std::string> out;
  std::string item;
  std::istringstream s(text);
  while (std::getline(s, item, sep)) out.push_back(item);
  if (!text.empty() && text.back() == sep) out.push_back("");
  return out;
}

// "1,2,5-7" -> {1,2,5,6,7}.  "" or what is wrong.
std::string ParseSegmentList(const std::string &text, std::vector<int> &out) {
  out.clear();
  std::set<int> seen;
  for (const std::string &item : Split(text, ',')) {
    int lo, hi;
    size_t dash = item.find('-', 1);
    bool ok;
    if (dash == std::string::npos) {
      ok = ReadWholeInt(item, lo);
      hi = lo;
    } else {
      ok = ReadWholeInt(item.substr(0, dash), lo) && ReadWholeInt(item.substr(dash + 1), hi) && hi >= lo;
    }
    if (!ok || lo < 1) return "segments='" + text + "': expected segment numbers >= 1 like 1,2,4-6";
    for (int k = lo; k <= hi; k++) {
      if (!seen.insert(k).second) {
        std::ostringstream why;
        why << "segments='" << text << "': segment " << k << " is listed twice";
        return why.str();
      }
      out.push_back(k);
    }
  }
  std::sort(out.begin(), out.end());
  return "";
}

// "pmin-pmax" (MeV/c, 0 <= pmin <= pmax): the '-' whose both sides are numbers.
bool ReadWindow(const std::string &text, double &lo, double &hi) {
  for (size_t k = 1; k + 1 < text.size(); k++) {
    if (text[k] != '-') continue;
    if (ReadWholeDouble(text.substr(0, k), lo) && ReadWholeDouble(text.substr(k + 1), hi))
      return lo >= 0.0 && hi >= lo;
  }
  return false;
}

// The value of ps=: "" or what is wrong.
std::string ParsePs(const std::string &value, ThmExperiment &x) {
  const std::string usage =
      "ps='" + value +
      "': expected delta, hulthen:pmin-pmax, hulthen:a,b:pmin-pmax (a, b in fm^-1), "
      "gauss:FWHM:pmin-pmax or table:<file> (momenta in MeV/c, 0 <= pmin <= pmax)";
  x.psA = 0.2317;
  x.psB = 1.202;
  x.psFwhm = 0.0;
  x.psMin = x.psMax = 0.0;
  x.psTable.clear();
  if (value == "delta") {
    x.psKind = ThmExperiment::PS_DELTA;
    return "";
  }
  std::vector<std::string> f = Split(value, ':');
  if (f.size() >= 2 && f[0] == "table") {
    x.psKind = ThmExperiment::PS_TABLE;
    x.psTable = value.substr(6);
    return x.psTable.empty() ? usage : "";
  }
  if (f.size() == 2 && f[0] == "hulthen") {
    x.psKind = ThmExperiment::PS_HULTHEN;
    return ReadWindow(f[1], x.psMin, x.psMax) ? "" : usage;
  }
  if (f.size() == 3 && f[0] == "hulthen") {
    x.psKind = ThmExperiment::PS_HULTHEN;
    std::vector<std::string> ab = Split(f[1], ',');
    if (ab.size() != 2 || !ReadWholeDouble(ab[0], x.psA) || !ReadWholeDouble(ab[1], x.psB) ||
        !(x.psA > 0.0) || !(x.psB > x.psA))
      return "ps='" + value + "': Hulthen a,b in fm^-1 with 0 < a < b (deuteron: 0.2317,1.202)";
    return ReadWindow(f[2], x.psMin, x.psMax) ? "" : usage;
  }
  if (f.size() == 3 && f[0] == "gauss") {
    x.psKind = ThmExperiment::PS_GAUSS;
    if (!ReadWholeDouble(f[1], x.psFwhm) || !(x.psFwhm > 0.0))
      return "ps='" + value + "': the FWHM of |phi(p_s)|^2 must be a number > 0 (MeV/c)";
    return ReadWindow(f[2], x.psMin, x.psMax) ? "" : usage;
  }
  return usage;
}

// opticalAA= / opticalSF=: plane | coulomb | V,R,a,W,RW,aW,WD,RD,aD,RC.
std::string ParseOptical(const std::string &key, const std::string &value, ThmExperiment::Optical &o) {
  o = ThmExperiment::Optical();
  if (value == "plane") {
    o.kind = 0;
    return "";
  }
  if (value == "coulomb") {
    o.kind = 1;
    return "";
  }
  std::vector<std::string> f = Split(value, ',');
  bool ok = f.size() == 10;
  for (size_t k = 0; ok && k < 10; k++) ok = ReadWholeDouble(f[k], o.p[k]);
  // Radii and diffusenesses > 0 where their depth is non-zero; RC >= 0.
  for (int t = 0; ok && t < 3; t++)
    if (o.p[3 * t] != 0.0 && !(o.p[3 * t + 1] > 0.0 && o.p[3 * t + 2] > 0.0)) ok = false;
  if (ok && !(o.p[9] >= 0.0)) ok = false;
  if (!ok)
    return key + "='" + value +
           "': expected plane, coulomb or ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC (MeV and fm: real volume, "
           "imaginary volume and imaginary surface Woods-Saxon, depths > 0 attractive/absorptive, radii and "
           "diffusenesses > 0 where the depth is not 0, the Coulomb radius RC >= 0, 0 = point charge)";
  o.kind = 2;
  return "";
}

bool ValidName(const std::string &name) {
  if (name.empty()) return false;
  for (char c : name)
    if (!std::isalnum((unsigned char)c) && c != '_' && c != '-' && c != '.' && c != '+') return false;
  return true;
}

}  // namespace

std::string ThmNuclide::Parse(const std::string &text, ThmNuclide &out) {
  for (const ThmNuclide &n : kNuclides)
    if (text == n.name) {
      out = n;
      return "";
    }
  std::vector<std::string> f = Split(text, ',');
  if (f.size() == 3) {
    int Z, A;
    double mass;
    if (ReadWholeInt(f[0], Z) && ReadWholeInt(f[1], A) && ReadWholeDouble(f[2], mass) && Z >= 0 && A >= 1 &&
        Z <= A && mass > 0.0 && std::fabs(mass - A) < 1.0) {
      out.name = text;
      out.Z = Z;
      out.A = A;
      out.mass = mass;
      return "";
    }
    return "'" + text + "': expected Z,A,mass with 0 <= Z <= A and the nuclear mass in u (within 1 u of A)";
  }
  std::string names;
  for (const ThmNuclide &n : kNuclides) names += (names.empty() ? "" : " ") + n.name;
  return "unknown nuclide '" + text + "' (known: " + names + "; or Z,A,mass in u)";
}

const ThmNuclide *ThmNuclide::Find(int Z, int A) {
  for (const ThmNuclide &n : kNuclides)
    if (n.Z == Z && n.A == A) return &n;
  return nullptr;
}

const char *ThmExperiment::BackgroundName(int terms) {
  static const char *names[] = {"none", "const", "linear", "quadratic"};
  return (terms >= 0 && terms <= 3) ? names[terms] : "?";
}

std::string ParseThmExperimentLine(const std::string &line, std::vector<ThmExperiment> &experiments) {
  const std::string head = "experiment[";
  size_t close = line.find(']');
  if (line.compare(0, head.size(), head) != 0 || close == std::string::npos)
    return "expected experiment[<name>] key=value ...";
  std::string name = line.substr(head.size(), close - head.size());
  if (!ValidName(name))
    return "experiment[" + name + "]: a name is letters, digits and _ - . + only";
  std::string where = "experiment[" + name + "]: ";
  std::string rest = line.substr(close + 1);
  if (!rest.empty() && rest[0] != ' ' && rest[0] != '\t')
    return where + "expected a space after ']'";

  ThmExperiment *x = nullptr;
  for (ThmExperiment &e : experiments)
    if (e.name == name) x = &e;
  ThmExperiment fresh;
  fresh.name = name;
  ThmExperiment work = x ? *x : fresh;  // committed only if the whole line is good

  std::istringstream tokens(rest);
  std::string token;
  int count = 0;
  while (tokens >> token) {
    count++;
    size_t eq = token.find('=');
    if (eq == std::string::npos || eq == 0 || eq + 1 == token.size())
      return where + "'" + token + "' is not key=value";
    std::string key = token.substr(0, eq), value = token.substr(eq + 1);
    if (key == "theta")
      return where + "key '" + key + "' is reserved for a later version (not implemented yet)";
    if (std::find(work.keys.begin(), work.keys.end(), key) != work.keys.end())
      return where + "key '" + key + "' is given twice";
    std::string why;
    if (key == "segments") {
      why = ParseSegmentList(value, work.segments);
    } else if (key == "background") {
      int terms = -1;
      for (int t = 0; t <= 3; t++)
        if (value == ThmExperiment::BackgroundName(t)) terms = t;
      if (terms < 0) why = "background='" + value + "': expected none, const, linear or quadratic";
      work.backgroundTerms = terms;
    } else if (key == "beam" || key == "target" || key == "spectator") {
      ThmNuclide &n = key == "beam" ? work.beam : key == "target" ? work.target : work.spectator;
      why = ThmNuclide::Parse(value, n);
      work.hasKinematics = true;
      if (!why.empty()) why = key + ": " + why;
    } else if (key == "Ebeam") {
      if (!ReadWholeDouble(value, work.beamEnergy) || !(work.beamEnergy > 0.0))
        why = "Ebeam='" + value + "': expected the lab beam energy in MeV, > 0";
    } else if (key == "lineshape") {
      if (value == "on")
        work.lineshape = true;
      else if (value == "off")
        work.lineshape = false;
      else
        why = "lineshape='" + value + "': expected on or off";
    } else if (key == "ps") {
      why = ParsePs(value, work);
    } else if (key == "psNodes") {
      if (!ReadWholeInt(value, work.psNodes) || work.psNodes < 1 || work.psNodes > 64)
        why = "psNodes='" + value + "': expected a whole number of Gauss-Legendre nodes, 1 to 64";
    } else if (key == "distortion") {
      work.distortionTable.clear();
      if (value == "none")
        work.distortion = ThmExperiment::DIST_NONE;
      else if (value == "coulomb")
        work.distortion = ThmExperiment::DIST_COULOMB;
      else if (value == "optical")
        work.distortion = ThmExperiment::DIST_OPTICAL;
      else if (value.compare(0, 6, "table:") == 0 && value.size() > 6) {
        work.distortion = ThmExperiment::DIST_TABLE;
        work.distortionTable = value.substr(6);
      } else
        why = "distortion='" + value + "': expected none, coulomb, optical or table:<file>";
    } else if (key == "opticalAA" || key == "opticalSF") {
      why = ParseOptical(key, value, key == "opticalAA" ? work.opticalAA : work.opticalSF);
    } else if (key == "spectatorAngle") {
      double a = 0.0;
      if (value == "qf") {
        work.angleKind = 0;
      } else if (value.compare(0, 3, "cm:") == 0 && ReadWholeDouble(value.substr(3), a) && a >= 0.0 && a <= 180.0) {
        work.angleKind = 2;
        work.angle = a;
      } else if (ReadWholeDouble(value, a) && a >= 0.0 && a <= 180.0) {
        work.angleKind = 1;
        work.angle = a;
      } else
        why = "spectatorAngle='" + value +
              "': expected qf, a lab angle in degrees (0-180) or cm:<degrees> (0-180)";
    } else if (key == "distortionRef") {
      if (!ReadWholeDouble(value, work.distortionRef))
        why = "distortionRef='" + value + "': expected the reference energy E_ref in MeV (c.m. of x + A)";
      work.hasDistortionRef = true;
    } else if (key == "distortionRatio") {
      if (value == "dwpw")
        work.distortionRatioPW = true;
      else if (value == "dw")
        work.distortionRatioPW = false;
      else
        why = "distortionRatio='" + value + "': expected dwpw or dw";
    } else if (key == "boundState") {
      std::vector<std::string> f = Split(value, ':');
      double rmin = 0.0;
      bool ok = (f.size() == 1 || f.size() == 2) && (f[0] == "whittaker" || f[0] == "yukawa");
      if (ok && f.size() == 2) ok = ReadWholeDouble(f[1], rmin) && rmin >= 0.0 && rmin <= 50.0;
      if (ok) {
        work.boundYukawa = f[0] == "yukawa";
        work.boundRmin = rmin;
      } else
        why = "boundState='" + value + "': expected whittaker or yukawa, optionally :rmin in fm (0-50)";
    } else {
      why = "unknown key '" + key +
            "' (keys: segments, background, beam, target, spectator, Ebeam, lineshape, ps, psNodes, "
            "distortion, opticalAA, opticalSF, spectatorAngle, distortionRef, distortionRatio, boundState)";
    }
    if (!why.empty()) return where + why;
    work.keys.push_back(key);
  }
  if (count == 0) return where + "no key=value given";
  if (x)
    *x = work;
  else
    experiments.push_back(work);
  return "";
}

std::string CheckThmExperiments(const std::vector<ThmExperiment> &experiments) {
  std::map<int, std::string> owner;
  for (const ThmExperiment &x : experiments) {
    std::string where = "experiment[" + x.name + "]: ";
    auto has = [&](const char *key) { return std::find(x.keys.begin(), x.keys.end(), key) != x.keys.end(); };
    if (!has("segments")) return where + "segments= is required";
    int kin = has("beam") + has("target") + has("spectator") + has("Ebeam");
    if (kin != 0 && kin != 4) return where + "beam, target, spectator and Ebeam go together (all four or none)";
    if (x.lineshape && kin != 4)
      return where + "lineshape=on needs the kinematics of the reaction: beam, target, spectator and Ebeam";
    if (x.psKind != ThmExperiment::PS_DELTA && kin != 4)
      return where + "a ps window (ps=hulthen|gauss|table) needs the kinematics of the reaction: beam, target, "
                     "spectator and Ebeam (mu_sx from the spectator and x masses)";
    if (has("psNodes") && x.psKind == ThmExperiment::PS_DELTA)
      return where + "psNodes= needs a ps window (ps=hulthen|gauss|table)";
    const bool computed = x.distortion == ThmExperiment::DIST_COULOMB || x.distortion == ThmExperiment::DIST_OPTICAL;
    if (computed && kin != 4)
      return where + "distortion=" + (x.distortion == ThmExperiment::DIST_COULOMB ? "coulomb" : "optical") +
             " needs the kinematics of the reaction: beam, target, spectator and Ebeam";
    if ((has("opticalAA") || has("opticalSF")) && x.distortion != ThmExperiment::DIST_OPTICAL)
      return where + "opticalAA= and opticalSF= need distortion=optical";
    for (const char *key : {"spectatorAngle", "distortionRef", "distortionRatio", "boundState"})
      if (has(key) && !computed)
        return where + key + "= needs distortion=coulomb or distortion=optical";
    for (int k : x.segments) {
      auto it = owner.find(k);
      if (it != owner.end()) {
        std::ostringstream why;
        why << where << "segment " << k << " is already in experiment[" << it->second << "]";
        return why.str();
      }
      owner[k] = x.name;
    }
  }
  return "";
}

std::string ReadThmPsTable(const std::string &path, std::vector<double> &p, std::vector<double> &w) {
  p.clear();
  w.clear();
  std::ifstream in(path.c_str());
  if (!in) return "cannot read the ps table '" + path + "'";
  std::string line;
  int lineNumber = 0;
  double total = 0.0;
  while (std::getline(in, line)) {
    lineNumber++;
    size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
    std::istringstream ls(line);
    double pv, wv;
    std::string extra;
    std::ostringstream where;
    where << "'" << path << "' line " << lineNumber << ": ";
    if (!(ls >> pv >> wv) || (ls >> extra)) return where.str() + "expected two numbers, p_s (MeV/c) and w";
    if (!std::isfinite(pv) || !std::isfinite(wv) || pv < 0.0 || wv < 0.0)
      return where.str() + "p_s and the weight must be finite and >= 0";
    if (!p.empty() && !(pv > p.back())) return where.str() + "p_s must be strictly increasing";
    p.push_back(pv);
    w.push_back(wv);
    total += wv;
  }
  if (p.size() < 2) return "'" + path + "' needs at least two rows (p_s w)";
  if (!(total > 0.0)) return "'" + path + "': every weight is zero";
  return "";
}

// ---------------------------------------------------------------------------
// The linear profile

namespace {

// Inverse of a symmetric positive-definite q x q matrix (q <= 4), by Cholesky
// on the Jacobi-scaled matrix (unit diagonal), so that the columns' units
// (the model, 1, E, E^2) do not matter.  False if it is not positive definite
// to working precision.
bool InvertSPD(const std::vector<double> &Gin, int q, std::vector<double> &inv) {
  std::vector<double> dscale(q), G(q * q), L(q * q, 0.0);
  for (int i = 0; i < q; i++) {
    if (!(Gin[i * q + i] > 0.0) || !std::isfinite(Gin[i * q + i])) return false;
    dscale[i] = 1.0 / std::sqrt(Gin[i * q + i]);
  }
  for (int i = 0; i < q; i++)
    for (int j = 0; j < q; j++) G[i * q + j] = Gin[i * q + j] * dscale[i] * dscale[j];
  for (int j = 0; j < q; j++) {
    double sum = G[j * q + j];
    for (int k = 0; k < j; k++) sum -= L[j * q + k] * L[j * q + k];
    if (!(sum > 1.0e-12)) return false;  // singular (collinear columns) to working precision
    L[j * q + j] = std::sqrt(sum);
    for (int i = j + 1; i < q; i++) {
      double s = G[i * q + j];
      for (int k = 0; k < j; k++) s -= L[i * q + k] * L[j * q + k];
      L[i * q + j] = s / L[j * q + j];
    }
  }
  // inv(G) = L^-T L^-1, column by column.
  inv.assign(q * q, 0.0);
  for (int c = 0; c < q; c++) {
    std::vector<double> y(q, 0.0), x(q, 0.0);
    for (int i = 0; i < q; i++) {
      double s = (i == c) ? 1.0 : 0.0;
      for (int k = 0; k < i; k++) s -= L[i * q + k] * y[k];
      y[i] = s / L[i * q + i];
    }
    for (int i = q - 1; i >= 0; i--) {
      double s = y[i];
      for (int k = i + 1; k < q; k++) s -= L[k * q + i] * x[k];
      x[i] = s / L[i * q + i];
    }
    for (int i = 0; i < q; i++) inv[i * q + c] = x[i] * dscale[i] * dscale[c];
  }
  return true;
}

// Column j of the design matrix at a point: the model, then 1, E, E^2.
inline double Column(int j, bool withModel, double m, double energy) {
  if (withModel) {
    if (j == 0) return m;
    j--;
  }
  return j == 0 ? 1.0 : j == 1 ? energy : energy * energy;
}

}  // namespace

double ThmProfile::Background(double energy) const {
  if (terms == 0) return 0.0;
  double b = a[0] + (terms > 1 ? a[1] * energy : 0.0) + (terms > 2 ? a[2] * energy * energy : 0.0);
  return b / s;
}

double ThmProfile::Residual(double m, double d, double e, double energy) const {
  if (e == 0.0) return 0.0;
  return (s * m + s * Background(energy) - d) / e;
}

void ThmProfile::Reported(double value[4], double covariance[16]) const {
  value[0] = 1.0 / s;
  for (int k = 0; k < 3; k++) value[k + 1] = k < terms ? a[k] / s : 0.0;
  for (int i = 0; i < 16; i++) covariance[i] = 0.0;
  // Rows of d(n, b0, b1, b2)/d(c), c the profiled coefficients.
  const int q = (scaleProfiled ? 1 : 0) + (backgroundProfiled ? terms : 0);
  if (q == 0 || (int)cov.size() != q * q) return;
  std::vector<double> T(4 * q, 0.0);
  int col = 0;
  if (scaleProfiled) {
    T[0 * q + col] = -1.0 / (s * s);
    for (int k = 0; k < terms; k++) T[(k + 1) * q + col] = -a[k] / (s * s);
    col++;
  }
  if (backgroundProfiled)
    for (int k = 0; k < terms; k++) T[(k + 1) * q + col + k] = 1.0 / s;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      double sum = 0.0;
      for (int u = 0; u < q; u++)
        for (int v = 0; v < q; v++) sum += T[i * q + u] * cov[u * q + v] * T[j * q + v];
      covariance[i * 4 + j] = sum;
    }
}

ThmProfile SolveThmProfile(const std::vector<double> &m, const std::vector<double> &d,
                           const std::vector<double> &e, const std::vector<double> &energy, int terms) {
  ThmProfile p;
  p.terms = terms;
  const size_t n = m.size();
  for (size_t i = 0; i < n; i++)
    if (e[i] != 0.0) p.points++;

  // Normal equations over the columns in use; `withModel` includes the model.
  auto solve = [&](bool withModel, std::vector<double> &c, std::vector<double> &cov) -> bool {
    const int q = (withModel ? 1 : 0) + terms;
    if (q == 0) return false;
    std::vector<double> G(q * q, 0.0), h(q, 0.0), col(q);
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      const double w = 1.0 / (e[i] * e[i]);
      // Without the model column the target is d - m (s fixed at 1).
      const double y = withModel ? d[i] : d[i] - m[i];
      for (int j = 0; j < q; j++) col[j] = Column(j, withModel, m[i], energy[i]);
      for (int j = 0; j < q; j++) {
        h[j] += col[j] * y * w;
        for (int k = 0; k < q; k++) G[j * q + k] += col[j] * col[k] * w;
      }
    }
    if (!InvertSPD(G, q, cov)) return false;
    c.assign(q, 0.0);
    for (int j = 0; j < q; j++)
      for (int k = 0; k < q; k++) c[j] += cov[j * q + k] * h[k];
    return true;
  };

  std::vector<double> c, cov;
  if (terms == 0) {
    // The model column alone: s = S_md/S_mm, as ESegment::ProfileNormChiSquared.
    double Smm = 0.0, Smd = 0.0;
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      const double w = 1.0 / (e[i] * e[i]);
      Smm += m[i] * m[i] * w;
      Smd += m[i] * d[i] * w;
    }
    if (Smd > 0.0 && Smm > 0.0) {
      p.s = Smd / Smm;
      p.cov.assign(1, 1.0 / Smm);
    } else {
      p.scaleProfiled = false;
      p.s = 1.0;
      p.status = "degenerate: no positive model-data overlap, norm left at 1";
    }
    p.backgroundProfiled = false;
  } else if (solve(true, c, cov) && c[0] > 0.0) {
    p.s = c[0];
    for (int k = 0; k < terms; k++) p.a[k] = c[k + 1];
    p.cov = cov;
  } else {
    p.scaleProfiled = false;
    p.s = 1.0;
    if (solve(false, c, cov)) {
      for (int k = 0; k < terms; k++) p.a[k] = c[k];
      p.cov = cov;
      p.status = "degenerate: model and background cannot both be profiled (singular, or n <= 0); "
                 "norm left at 1, background profiled";
    } else {
      p.backgroundProfiled = false;
      p.status = "degenerate: too few points or collinear background terms; norm left at 1, no background";
    }
  }
  p.chi2 = 0.0;
  for (size_t i = 0; i < n; i++) {
    double r = p.Residual(m[i], d[i], e[i], energy[i]);
    p.chi2 += r * r;
  }
  return p;
}

void ThmProfileDerivative(const ThmProfile &p, const std::vector<double> &m, const std::vector<double> &d,
                          const std::vector<double> &e, const std::vector<double> &energy,
                          const std::vector<double> &Jm, int nCols, std::vector<double> &J,
                          std::vector<double> &G) {
  const size_t n = m.size();
  J.assign(n * nCols, 0.0);
  G.assign(n * nCols, 0.0);
  const bool withModel = p.scaleProfiled;
  const int q = (withModel ? 1 : 0) + (p.backgroundProfiled ? p.terms : 0);
  // Coefficients in the order of the profiled columns, and the fixed scale.
  std::vector<double> c;
  if (withModel) c.push_back(p.s);
  if (p.backgroundProfiled)
    for (int k = 0; k < p.terms; k++) c.push_back(p.a[k]);
  const double s = p.s;
  std::vector<double> rho(n, 0.0);  // y~ - A~ c = -r
  for (size_t i = 0; i < n; i++) rho[i] = -p.Residual(m[i], d[i], e[i], energy[i]);

  std::vector<double> v(q), dc(q);
  for (int col = 0; col < nCols; col++) {
    // v = (dA~)^T rho - A~^T (dA~) c, dA~ = [Jm/e, 0, ...]; (dA~) c = s Jm/e.
    std::fill(v.begin(), v.end(), 0.0);
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      const double jt = Jm[i * nCols + col] / e[i];
      if (withModel) v[0] += jt * rho[i];
      for (int j = 0; j < q; j++) v[j] -= Column(j, withModel, m[i], energy[i]) / e[i] * s * jt;
    }
    for (int j = 0; j < q; j++) {
      dc[j] = 0.0;
      for (int k = 0; k < q; k++) dc[j] += p.cov[j * q + k] * v[k];
    }
    for (size_t i = 0; i < n; i++) {
      // d f_i/dp = s Jm_i + sum_j A_ij dc_j (unwhitened), f the curve on the data scale.
      double df = s * Jm[i * nCols + col];
      for (int j = 0; j < q; j++) df += Column(j, withModel, m[i], energy[i]) * dc[j];
      if (e[i] != 0.0) J[i * nCols + col] = df / e[i];
      G[i * nCols + col] = df / s;
    }
  }
}
