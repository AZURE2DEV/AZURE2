// The keys of the <thm> block as two tables (ThmKeys.h), and the readers of a line that use them.
#include "ThmKeys.h"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "ThmText.h"
#include "azure2/physics/transfer/GlobalOptical.h"

using thm_text::ReadWholeDouble;
using thm_text::ReadWholeInt;
using thm_text::ReadWindow;
using thm_text::Split;

namespace {

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

// opticalAA= / opticalSF=: plane | coulomb | <global name>[:extrapolate] | V,R,a,W,RW,aW,WD,RD,aD,RC.
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
  if (!value.empty() && std::isalpha((unsigned char)value[0])) {
    const size_t colon = value.find(':');
    const std::string name = value.substr(0, colon), option = colon == std::string::npos ? "" : value.substr(colon + 1);
    namespace transfer = azure2::physics::transfer;
    const int index = transfer::GlobalOpticalIndex(name);
    if (index < 0 || (colon != std::string::npos && option != "extrapolate"))
      return key + "='" + value + "': expected plane, coulomb, a global optical potential (" +
             transfer::GlobalOpticalNames() + ", optionally :extrapolate) or ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC";
    o.kind = 3;
    o.global = index;
    o.extrapolate = colon != std::string::npos;
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

// A nucleus of the three-body reaction (beam, target, spectator): ThmNuclide::Parse, and the experiment has kinematics.
std::string ParseNuclide(const char *key, const std::string &value, ThmNuclide &n, ThmExperiment &x) {
  std::string why = ThmNuclide::Parse(value, n);
  x.hasKinematics = true;
  return why.empty() ? why : std::string(key) + ": " + why;
}

bool ParseWord(const std::string &value, std::initializer_list<const char *> words, ThmOptionLine &out) {
  out.word = value;
  for (const char *w : words)
    if (value == w) return true;
  return false;
}

std::string Present(const char *key) { return key; }

}  // namespace

const ThmKeyKind &ThmOptionKind() {
  static const ThmKeyKind kind{"<thm>", "<thm> options", "The model options of every THM experiment of the project."};
  return kind;
}

const ThmKeyKind &ThmExperimentKind() {
  static const ThmKeyKind kind{"experiment", "experiment",
                               "One THM experiment: the data sets that name it share one profiled norm and its "
                               "background; its kinematics, momentum window and distortion."};
  return kind;
}

const std::vector<ThmOptionKey> &ThmOptionKeys() {
  static const std::vector<ThmOptionKey> keys = {
      {{"vertex", "enum", {"perlevel", "constant", "onshell"}, "", "", "constant", "", "",
        "boundary in the entrance vertex (B - 1) j_l - rho j_l': perlevel S_c(E_lambda) of each level under Brune "
        "(alias real); constant B_c = S_c(E_1) of the lowest level of the J group; onshell L_c(E) = S_c(E) + i P_c(E)"},
       false, "onshell|constant|perlevel",
       [](const std::string &value, ThmOptionLine &out) {
         out.kind = ThmOptionLine::VERTEX;
         return ParseWord(value == "real" ? std::string("perlevel") : value, {"onshell", "constant", "perlevel"}, out);
       }},
      {{"kinematics", "enum", {"lacognata", "triple", "kf3body", "lambda32"}, "", "", "lacognata", "", "",
        "the kinematic factor K(E) the model carries, by how the HOES data were extracted (the experiments divide "
        "the yield by the Monte Carlo KF |phi|^2, so no further KF is applied): lacognata exit k_f/mu_f (the default; "
        "Tumino et al. 2021 Eq. 50); triple K = 1; kf3body exit Gamma_f/(mu_f k_f) and lambda32 entrance 1/k_i on "
        "shell are kept for old projects"},
       false, "lacognata|triple|kf3body|lambda32",
       [](const std::string &value, ThmOptionLine &out) {
         out.kind = ThmOptionLine::KINEMATICS;
         return ParseWord(value, {"lacognata", "triple", "kf3body", "lambda32"}, out);
       }},
      {{"entranceL", "enum", {"incoherent", "coherent"}, "", "", "incoherent", "", "",
        "entrance partial waves of one channel spin summed coherently or incoherently (for an observable integrated "
        "over the exit direction the l cross terms vanish)"},
       false, "coherent|incoherent",
       [](const std::string &value, ThmOptionLine &out) {
         out.kind = ThmOptionLine::ENTRANCE_L;
         return ParseWord(value, {"coherent", "incoherent"}, out);
       }},
      {{"coulombIntegral", "enum", {"0", "1"}, "", "", "0", "", "",
        "add the external Coulomb term 2 eta k int_a^inf O_l(kr)/O_l(ka) j_l(pr) dr to the vertex (Tribble eq. "
        "2.79)"},
       false, "0|1",
       [](const std::string &value, ThmOptionLine &out) {
         out.kind = ThmOptionLine::COULOMB_INTEGRAL;
         if (value == "1" || value == "true" || value == "on") out.flag = true;
         else if (value == "0" || value == "false" || value == "off") out.flag = false;
         else return false;
         return true;
       }},
      {{"spectatorEnergy", "number", {}, "spectatorEnergy (MeV)", "MeV", "0", "", "",
        "mean spectator kinetic energy <p_sx^2>/2mu_sx added to E + B in the half-off-shell momentum, for every THM "
        "pair (one pair: the key spectatorEnergy[<pair key>])"},
       true, "<MeV>",
       [](const std::string &value, ThmOptionLine &out) {
         out.kind = ThmOptionLine::SPECTATOR_ENERGY;
         std::istringstream vs(value);
         std::string rest;
         return !!(vs >> out.energy) && out.energy >= 0.0 && !(vs >> rest);  // "0.4junk" is not a number
       }},
  };
  return keys;
}

std::string ThmOptionUsage() {
  std::string out;
  for (const ThmOptionKey &k : ThmOptionKeys()) {
    out += (out.empty() ? "" : ", ") + k.row.key + "=" + k.usage;
    if (k.perPair) out += ", " + k.row.key + "[<pair key>]=" + k.usage;
  }
  return out;
}

const std::vector<ThmExperimentKey> &ThmExperimentKeys() {
  static const std::vector<ThmExperimentKey> keys = {
      {{"background", "enum", {"none", "const", "linear", "quadratic"}, "", "", "none", "", "",
        "the incoherent background b(E) added to the folded model, profiled with the norm"},
       [](const std::string &value, ThmExperiment &x) {
         int terms = -1;
         for (int t = 0; t <= 3; t++)
           if (value == ThmExperiment::BackgroundName(t)) terms = t;
         x.backgroundTerms = terms;
         return terms < 0 ? "background='" + value + "': expected none, const, linear or quadratic" : std::string();
       }},
      {{"beam", "text", {}, "", "", "", "<nucleus>", "", "the beam nucleus (a built-in name or Z,A,mass)"},
       [](const std::string &value, ThmExperiment &x) { return ParseNuclide("beam", value, x.beam, x); }},
      {{"target", "text", {}, "", "", "", "<nucleus>", "",
        "the target nucleus, the Trojan horse a = x + s (a built-in name or Z,A,mass)"},
       [](const std::string &value, ThmExperiment &x) { return ParseNuclide("target", value, x.target, x); }},
      {{"spectator", "text", {}, "", "", "", "<nucleus>", "", "the spectator s (a built-in name or Z,A,mass)"},
       [](const std::string &value, ThmExperiment &x) { return ParseNuclide("spectator", value, x.spectator, x); }},
      {{"Ebeam", "number", {}, "Ebeam (MeV)", "MeV", "", "", "", "the lab beam energy, > 0"},
       [](const std::string &value, ThmExperiment &x) {
         if (!ReadWholeDouble(value, x.beamEnergy) || !(x.beamEnergy > 0.0))
           return "Ebeam='" + value + "': expected the lab beam energy in MeV, > 0";
         return std::string();
       }},
      {{"lineshape", "enum", {"off", "on"}, "", "", "off", "", "",
        "the Coulomb line shape N_C of the spectator (needs beam, target, spectator, Ebeam)"},
       [](const std::string &value, ThmExperiment &x) {
         if (value != "on" && value != "off") return "lineshape='" + value + "': expected on or off";
         x.lineshape = value == "on";
         return std::string();
       }},
      {{"ps", "text", {}, "", "", "", "delta | hulthen:pmin-pmax | hulthen:a,b:pmin-pmax | gauss:FWHM:pmin-pmax | "
        "table:<file>", "", "the spectator-momentum window"},
       [](const std::string &value, ThmExperiment &x) { return ParsePs(value, x); }},
      {{"psNodes", "integer", {}, "", "", "", "", Present("ps"), "Gauss-Legendre nodes of the ps window, 1 to 64"},
       [](const std::string &value, ThmExperiment &x) {
         if (!ReadWholeInt(value, x.psNodes) || x.psNodes < 1 || x.psNodes > 64)
           return "psNodes='" + value + "': expected a whole number of Gauss-Legendre nodes, 1 to 64";
         return std::string();
       }},
      {{"distortion", "text", {}, "", "", "", "none | coulomb | optical | table:<file>", "",
        "the distortion factor R(E) multiplying the model"},
       [](const std::string &value, ThmExperiment &x) {
         x.distortionTable.clear();
         if (value == "none")
           x.distortion = ThmExperiment::DIST_NONE;
         else if (value == "coulomb")
           x.distortion = ThmExperiment::DIST_COULOMB;
         else if (value == "optical")
           x.distortion = ThmExperiment::DIST_OPTICAL;
         else if (value.compare(0, 6, "table:") == 0 && value.size() > 6) {
           x.distortion = ThmExperiment::DIST_TABLE;
           x.distortionTable = value.substr(6);
         } else
           return "distortion='" + value + "': expected none, coulomb, optical or table:<file>";
         return std::string();
       }},
      {{"opticalAA", "text", {}, "", "", "", "plane | <global set>[:extrapolate] | <parameters>", "",
        "the a + A optical potential (distortion=optical or vertexModel=dw)"},
       [](const std::string &value, ThmExperiment &x) { return ParseOptical("opticalAA", value, x.opticalAA); }},
      {{"opticalSF", "text", {}, "", "", "", "plane | <global set>[:extrapolate] | <parameters>", "",
        "the s + F optical potential (distortion=optical or vertexModel=dw)"},
       [](const std::string &value, ThmExperiment &x) { return ParseOptical("opticalSF", value, x.opticalSF); }},
      {{"spectatorAngle", "text", {}, "", "", "", "qf | <lab degrees 0-180> | cm:<degrees 0-180>", "",
        "the spectator direction of the distortion factor and the DW vertex"},
       [](const std::string &value, ThmExperiment &x) {
         double a = 0.0;
         if (value == "qf") {
           x.angleKind = 0;
         } else if (value.compare(0, 3, "cm:") == 0 && ReadWholeDouble(value.substr(3), a) && a >= 0.0 && a <= 180.0) {
           x.angleKind = 2;
           x.angle = a;
         } else if (ReadWholeDouble(value, a) && a >= 0.0 && a <= 180.0) {
           x.angleKind = 1;
           x.angle = a;
         } else
           return "spectatorAngle='" + value + "': expected qf, a lab angle in degrees (0-180) or cm:<degrees> (0-180)";
         return std::string();
       }},
      {{"distortionRef", "number", {}, "distortionRef (MeV)", "MeV", "", "", "",
        "the reference energy E_ref of R(E) = rho(E)/rho(E_ref) (c.m. of x + A)"},
       [](const std::string &value, ThmExperiment &x) {
         x.hasDistortionRef = true;
         if (!ReadWholeDouble(value, x.distortionRef))
           return "distortionRef='" + value + "': expected the reference energy E_ref in MeV (c.m. of x + A)";
         return std::string();
       }},
      {{"distortionRatio", "enum", {"dwpw", "dw"}, "", "", "dwpw", "", "",
        "rho = |M|^2/|M_PW|^2 (dwpw) or |M|^2 (dw, the published convention)"},
       [](const std::string &value, ThmExperiment &x) {
         if (value != "dwpw" && value != "dw") return "distortionRatio='" + value + "': expected dwpw or dw";
         x.distortionRatioPW = value == "dwpw";
         return std::string();
       }},
      {{"boundState", "text", {}, "", "", "", "whittaker | yukawa [:rmin in fm, 0-50]", "",
        "the s-x bound state of the distortion factor and the DW vertex"},
       [](const std::string &value, ThmExperiment &x) {
         std::vector<std::string> f = Split(value, ':');
         double rmin = 0.0;
         bool ok = (f.size() == 1 || f.size() == 2) && (f[0] == "whittaker" || f[0] == "yukawa");
         if (ok && f.size() == 2) ok = ReadWholeDouble(f[1], rmin) && rmin >= 0.0 && rmin <= 50.0;
         if (!ok) return "boundState='" + value + "': expected whittaker or yukawa, optionally :rmin in fm (0-50)";
         x.boundYukawa = f[0] == "yukawa";
         x.boundRmin = rmin;
         return std::string();
       }},
      {{"theta", "text", {}, "", "", "", "all | thmin-thmax", "",
        "the c.m. angles of the exit pair relative to p_xA (degrees)"},
       [](const std::string &value, ThmExperiment &x) {
         double lo = 0.0, hi = 0.0;
         if (value == "all") {
           x.hasTheta = false;
           x.thetaMin = 0.0;
           x.thetaMax = 180.0;
         } else if (ReadWindow(value, lo, hi) && hi <= 180.0) {
           x.hasTheta = true;
           x.thetaMin = lo;
           x.thetaMax = hi;
         } else
           return "theta='" + value +
                  "': expected all or thmin-thmax, the c.m. angles of the exit pair relative to p_xA in "
                  "degrees, 0 <= thmin <= thmax <= 180";
         return std::string();
       }},
      {{"vertexModel", "enum", {"pw", "dw", "ff"}, "", "", "pw", "", "",
        "the transfer vertex: pw the plane-wave surface term M_l; dw the surface term of the prior-form DWBA; ff the "
        "full form factor <phi_lambda|V|j_l> of each level (La Cognata et al. ApJ 708 (2010) Eqs. 20-21: M_l plus the "
        "volume term of a Woods-Saxon level function with its node count; vertex=perlevel, fixed level energies)"},
       [](const std::string &value, ThmExperiment &x) {
         if (value != "pw" && value != "dw" && value != "ff") return "vertexModel='" + value + "': expected pw, dw or ff";
         x.vertexDW = value == "dw";
         x.vertexFF = value == "ff";
         return std::string();
       }},
      {{"ffWell", "text", {}, "", "fm", "1.25,0.65", "<r0>,<a0>", "vertexModel",
        "vertexModel=ff: the Woods-Saxon radius parameter r0 (R = r0 A^(1/3), A the core's mass number) and "
        "diffuseness a0 of the level functions, fm"},
       [](const std::string &value, ThmExperiment &x) {
         std::vector<std::string> f = Split(value, ',');
         double r0 = 0.0, a0 = 0.0;
         if (f.size() != 2 || !ReadWholeDouble(f[0], r0) || !ReadWholeDouble(f[1], a0) || !(r0 > 0.0) || !(a0 > 0.0))
           return "ffWell='" + value + "': expected <r0>,<a0> (fm, > 0): the Woods-Saxon radius parameter and "
                  "diffuseness of the level functions";
         x.ffR0 = r0;
         x.ffA0 = a0;
         return std::string();
       }},
      {{"ffVolume", "enum", {"on", "off"}, "", "", "on", "", "vertexModel",
        "vertexModel=ff: off drops the volume term, so the vertex is exactly that of vertexModel=pw (a check of the "
        "limit)"},
       [](const std::string &value, ThmExperiment &x) {
         if (value != "on" && value != "off") return "ffVolume='" + value + "': expected on or off";
         x.ffVolume = value == "on";
         return std::string();
       }},
      {{"ffNodes", "text", {}, "", "", "", "<J><+|->[,k<k>][,l<l>]:<n>[;...]", "vertexModel",
        "vertexModel=ff: the radial node count n of the level functions of a J^pi group (its k-th level, its entrance "
        "channels of orbital momentum l); without a term, the Wildermuth rule of a nucleon on the core"},
       [](const std::string &value, ThmExperiment &x) { return ParseThmFfNodes(value, x.ffNodes); }},
      {{"spectatorAngles", "text", {}, "", "", "", "[cm:]thmin-thmax | [cm:]table:<file>", "",
        "a window of spectator directions (lab polar angles to the beam, or c.m.)"},
       [](const std::string &value, ThmExperiment &x) {
         std::string v = value;
         x.angleCm = v.compare(0, 3, "cm:") == 0;
         if (x.angleCm) v = v.substr(3);
         x.angleTable.clear();
         x.angleTableT.clear();
         x.angleTableW.clear();
         double lo = 0.0, hi = 0.0;
         if (v.compare(0, 6, "table:") == 0 && v.size() > 6) {
           x.angleWindow = 2;
           x.angleTable = v.substr(6);
         } else if (ReadWindow(v, lo, hi) && hi <= 180.0) {
           x.angleWindow = 1;
           x.angleMin = lo;
           x.angleMax = hi;
         } else
           return "spectatorAngles='" + value +
                  "': expected thmin-thmax (lab polar angles of the spectator to the beam, degrees, 0 <= thmin "
                  "<= thmax <= 180), cm:thmin-thmax (c.m.), table:<file> or cm:table:<file> (angle, acceptance)";
         return std::string();
       }},
      {{"spectatorAngleNodes", "integer", {}, "", "", "", "", Present("spectatorAngles"),
        "Gauss-Legendre nodes of the spectatorAngles window, 1 to 64"},
       [](const std::string &value, ThmExperiment &x) {
         if (!ReadWholeInt(value, x.angleNodes) || x.angleNodes < 1 || x.angleNodes > 64)
           return "spectatorAngleNodes='" + value + "': expected a whole number of Gauss-Legendre nodes, 1 to 64";
         return std::string();
       }},
      {{"cbackground", "text", {}, "", "", "", "<J><+|->:<exit pair>[:<s>,<l>,<s'>,<l'>][:const|:linear]"
        "[=<v>,<v>[,<v>,<v>]]", "",
        "the coherent background: a THM-only amplitude c(E) = c0 + c1 E times the entrance vertex, per J^pi and exit "
        "pair; terms separated by ';', a value followed by f is fixed"},
       [](const std::string &value, ThmExperiment &x) { return ParseThmCoherentBackground(value, x.cbackground); }},
  };
  return keys;
}

bool ParseThmOptionLine(const std::string &key, const std::string &value, ThmOptionLine &out) {
  out = ThmOptionLine();
  for (const ThmOptionKey &k : ThmOptionKeys()) {
    if (key == k.row.key) return k.parse(value, out);
    // key[<pair key>]: the option for one pair
    const std::size_t n = k.row.key.size();
    if (!k.perPair || key.compare(0, n, k.row.key) != 0) continue;
    if (!k.parse(value, out) || !(key.size() > n + 2 && key[n] == '[' && key.back() == ']')) return false;
    std::istringstream ks(key.substr(n + 1, key.size() - n - 2));
    out.perPair = true;
    return !!(ks >> out.pairKey);
  }
  return false;
}

std::string ParseThmExperimentLine(const std::string &line, std::vector<ThmExperiment> &experiments) {
  const std::string head = "experiment[";
  size_t close = line.find(']');
  if (line.compare(0, head.size(), head) != 0 || close == std::string::npos)
    return "expected experiment[<name>] key=value ...";
  std::string name = line.substr(head.size(), close - head.size());
  if (!ValidName(name))
    return "experiment[" + name + "]: a name is letters, digits and _ - . + only; rename the experiment";
  std::string where = "experiment[" + name + "]: ";
  std::string rest = line.substr(close + 1);
  if (!rest.empty() && rest[0] != ' ' && rest[0] != '\t')
    return where + "expected a space after ']'; write experiment[<name>] key=value ...";

  ThmExperiment *x = nullptr;
  for (ThmExperiment &e : experiments)
    if (e.name == name) x = &e;
  ThmExperiment fresh;
  fresh.name = name;
  ThmExperiment work = x ? *x : fresh;  // committed only if the whole line is good

  std::istringstream tokens(rest);
  std::string token;
  while (tokens >> token) {
    size_t eq = token.find('=');
    if (eq == std::string::npos || eq == 0 || eq + 1 == token.size())
      return where + "'" + token + "' is not key=value; write each key as key=value";
    std::string key = token.substr(0, eq), value = token.substr(eq + 1);
    if (std::find(work.keys.begin(), work.keys.end(), key) != work.keys.end())
      return where + "key '" + key + "' is given twice; give it once";
    std::string why;
    const ThmExperimentKey *row = nullptr;
    for (const ThmExperimentKey &k : ThmExperimentKeys())
      if (k.row.key == key) row = &k;
    if (row) {
      why = row->parse(value, work);
    } else {
      std::string keys;
      for (const ThmExperimentKey &k : ThmExperimentKeys()) keys += (keys.empty() ? "" : ", ") + k.row.key;
      why = "unknown key '" + key + "' (keys: " + keys + ")";
    }
    // How to go on, the same for every key: the table's rows say what each key takes.
    if (!why.empty()) return where + why + "; the keys and their values: pyazr.kinds()[\"thm\"]";
    work.keys.push_back(key);
    work.values.push_back(value);
  }
  if (x)
    *x = work;
  else
    experiments.push_back(work);
  return "";
}
