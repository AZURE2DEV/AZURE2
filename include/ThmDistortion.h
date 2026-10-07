#ifndef THM_DISTORTION_H
#define THM_DISTORTION_H

#include "Constants.h"

#include <atomic>
#include <memory>
#include <string>
#include <vector>

struct ThmExperiment;
struct ThmWeightTable;

/*!
 * The energy dependence of the THM transfer amplitude from the distortions of
 * the initial (a + A) and final (s + F) relative motion: `distortion=` on a
 * THM experiment line (docs/source/theory/thm_implementation.rst,
 * "Distortion factor R(E)").
 *
 * Zero-range prior-form DWBA (Mukhamedzhanov, Pang & Kadyrov, PRC 99 (2019)
 * 064618, eqs. 20-24; Mukhamedzhanov, arXiv:2609.04498, eqs. 22-30): with the
 * x-A vertex of zero range, r_sF = r_sx = r and r_aA = beta r, beta = m_s/m_a,
 * and
 *
 *   M(E) = Int d^3r chi^(-)*_{k_sF}(r) phi_sx(r) chi^(+)_{k_aA}(beta r)
 *        = 4 pi/(k_sF beta k_aA) sum_l (2l+1) e^{i(sigma_l^aA + sigma_l^sF)} P_l(x)
 *          Int_0^inf dr phi_sx(r) u_l^sF(k_sF r) u_l^aA(k_aA beta r),
 *
 * x = k_sF.k_aA/(k_sF k_aA), u_l the regular radial waves normalized to
 * u -> F_l + T_l H_l^+ (u = F_l for point Coulomb), phi_sx the s-x bound
 * state (l_sx = 0): the Whittaker tail W_{-eta_b,1/2}(2 kappa r)/r or the
 * Yukawa e^{-kappa r}/r, zero below rmin.  The plane-wave limit is the
 * Fourier transform of phi at q = k_sF - beta k_aA, the momentum distribution
 * the PWA data reduction divides by:
 *
 *   M_PW(E) = 4 pi Int dr r^2 j_0(q r) phi_sx(r).
 *
 * The weight multiplying the HOES model (before folding, as weight[k]) is
 *
 *   R(E) = rho(E)/rho(E_ref),  rho = |M|^2/|M_PW|^2  (distortionRatio=dwpw)
 *                               or   |M|^2          (dw, the papers' R),
 *
 * i.e. the published S*_PWA has to be divided by R; R grows toward low E for
 * 12C+12C.  Distorted waves: point Coulomb (distortion=coulomb) or, per
 * channel, plane / point Coulomb / Woods-Saxon real + imaginary volume +
 * imaginary surface with a uniform-sphere Coulomb term (distortion=optical,
 * opticalAA=, opticalSF=), all by the same complex Numerov integration
 * outward from the origin, matched to AZURE2's Coulomb functions (COUL)
 * beyond the turning point.  The radial integrals converge absolutely (the
 * bound state decays as e^{-kappa r}); the partial-wave sum runs until the
 * bound on the remaining terms is below 1e-13 |M|.
 */
class ThmDistortion {
 public:
  enum Kind { COULOMB, OPTICAL, TABLE };
  /// One channel's distortion: none (plane wave), point Coulomb, or a
  /// Woods-Saxon potential p[] = V,R,a, W,RW,aW, WD,RD,aD, RC (MeV, fm;
  /// RC = 0 a point charge) plus Coulomb.  A global optical potential
  /// (ThmOptical.h, global >= 0) is a Woods-Saxon channel whose p[] follow
  /// the energy: SetEnergy evaluates the model at the projectile's lab energy.
  struct Channel {
    enum Kind { PLANE, POINT_COULOMB, WOODS_SAXON };
    Kind kind = POINT_COULOMB;
    double p[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int Z1 = 0, Z2 = 0;
    double mu = 0.0;  ///< reduced mass (MeV)
    double k = 0.0;   ///< wave number (fm^-1)
    double eta = 0.0; ///< Sommerfeld parameter (0 for PLANE)
    // Global optical potential.
    int global = -1;           ///< index in ThmGlobalOpticals(), or -1
    bool extrapolate = false;  ///< allowed outside the validity range
    int Zp = 0, Ap = 0, Zt = 0, At = 0;  ///< projectile (the model's light ion) and target
    double mp = 0.0, mt = 0.0;           ///< their masses (u)
    double depthMax = 0.0;     ///< largest summed depth over the energies the channel takes (MeV)
    /// The projectile lab energy at the channel's c.m. energy ecm (MeV).
    double LabEnergy(double ecm) const { return ecm * (mp + mt) / mt; }
    /// p[] of the global model at the channel's c.m. energy ecm (no-op otherwise).
    void SetEnergy(double ecm);
    std::string Describe() const;
  };
  /// The three-body kinematics (BuildThmGroups fills it from the experiment line).
  struct Kinematics {
    int Za = 0, ZA = 0, Zs = 0, Zx = 0;  ///< Trojan horse a, the other nucleus A, spectator s, x = a - s
    int Aa = 0, AA = 0, As = 0, Ax = 0;  ///< mass numbers (global optical potentials only)
    double ma = 0.0, mA = 0.0, ms = 0.0, mx = 0.0;  ///< nuclear masses (u)
    bool horseIsBeam = true;
    double mBeam = 0.0, mTarget = 0.0;  ///< u
    double beamEnergy = 0.0;            ///< lab (MeV)
    double bind = 0.0;                  ///< B_xs (MeV)
  };
  /// Everything at one energy (Evaluate).
  struct Point {
    double energy = 0.0;  ///< E, c.m. of x + A (MeV)
    double esf = 0.0, ksf = 0.0, etasf = 0.0;
    double thetaCm = 0.0;  ///< spectator c.m. angle to the beam (deg)
    double x = 1.0;        ///< cos of the angle between k_sF and k_aA
    double q = 0.0;        ///< |k_sF - beta k_aA| (fm^-1)
    complex m = complex(0.0, 0.0);  ///< M (fm^3 up to 4 pi/(k k') conventions as above)
    double mpw = 0.0;      ///< M_PW
    int lmax = 0;          ///< highest l summed
    bool ok = false;
    std::string why;       ///< if !ok
    double tail = 0.0;  ///< |integrand(r_end)|/(kappa |M|), l = 0: the radial cutoff's size
    /// Spectator-direction window (spectatorAngles=): the nodes with a
    /// weight, and the acceptance averages <|M|^2>, <|M_PW|^2> (thetaCm, x, q
    /// are then the acceptance-weighted means, m and mpw those of the first
    /// node).  0 without a window.
    int nodes = 0;
    double m2 = 0.0, mpw2 = 0.0;
  };

  Kind kind = COULOMB;
  std::string experiment;
  std::string description;  ///< settings, for the output
  // Settings (ThmExperiment).
  enum AngleKind { QF, LAB, CM };
  AngleKind angleKind = QF;
  double angle = 0.0;  ///< deg
  /*!
   * Spectator-direction window (spectatorAngles= and/or the |p_s| cut of a
   * ps window, docs "Experimental acceptance"): instead of one direction,
   * the accepted directions at each energy.  At fixed E, |k_sF| is fixed, so the direction (azimuthal symmetry
   * about the beam: its polar angle) is the only variable, and it fixes
   * q = |k_sF - beta k_aA|; a ps window (qCut) restricts it to q in
   * [qCutLo, qCutHi].  The measure is d cos(theta_cm) (the three-body phase
   * space at fixed E) times the acceptance A(theta) (1 in a uniform window, a
   * table otherwise, in the lab or c.m. angle).  A lab window maps to one
   * c.m. interval, or two when the spectator is slower in the c.m. than the
   * c.m. itself (gamma = V_cm/v_s > 1: forward and backward branch); each
   * interval gets angNodes Gauss-Legendre nodes in cos(theta_cm).
   */
  bool angWindow = false;
  bool angAll = false;  ///< a ps window alone: every direction (c.m. 0-180) whose q is in the cut
  bool angCm = false;
  double angLo = 0.0, angHi = 0.0;  ///< deg
  std::vector<double> angT, angW;   ///< acceptance table (empty: uniform)
  int angNodes = 8;                 ///< per interval
  int angSlots = 1;                 ///< intervals (SetAngleSlots)
  bool qCut = false;
  double qCutLo = 0.0, qCutHi = 0.0;  ///< MeV/c
  std::vector<double> angGx, angGw;   ///< Gauss-Legendre rule on [-1, 1]
  struct AngleNode {
    double theta = 0.0;  ///< c.m. angle of the spectator to the beam (deg)
    double x = 1.0;      ///< cos of the angle between k_sF and k_aA
    double q = 0.0;      ///< |k_sF - beta k_aA| (fm^-1)
    double w = 0.0;      ///< Gauss-Legendre weight in cos(theta_cm) times the acceptance
  };
  bool ratioPW = true; ///< dwpw (else dw)
  bool yukawa = false;
  double rmin = 0.0;   ///< fm (as used: snapped to the grid)
  double eRef = 0.0;   ///< MeV
  /// The data range (c.m. of x + A, MeV) over which a global optical potential
  /// is checked against its validity range (Setup); dataLo > dataHi: the grid ends.
  double dataLo = 1.0, dataHi = 0.0;
  /// WARNINGs of Setup (a global potential used outside its range with
  /// :extrapolate), without the "WARNING: " prefix.
  std::vector<std::string> warnings;
  // Derived.
  Kinematics kin;
  Channel aa, sf;
  double eAA = 0.0, muSx = 0.0, kappa = 0.0, etaB = 0.0, beta = 0.0, vcm = 0.0;
  double h = 0.0;      ///< radial step for r = r_sx (fm); the a + A wave uses beta h
  int i0 = 0, n = 0;   ///< integration from node i0 to node n-1 (Simpson)
  std::vector<double> phi;  ///< phi_sx on the grid
  std::vector<double> simpson;  ///< Simpson weights (times h), 0 below i0
  std::vector<complex> sigmaAA;  ///< e^{i sigma_l^aA} ... stored as the phase factor
  std::vector<std::vector<complex>> uAA;  ///< u_l^aA(beta r_i), l = 0..
  std::vector<double> tailAA;  ///< sum_{l' >= l} (2l'+1) Int |phi u_l'^aA| (bound on what is left)
  double refRatio = 0.0;  ///< rho(E_ref)
  Point ref;
  /// ln R on a uniform energy grid (cubic Lagrange interpolation).
  double gridLo = 0.0, gridStep = 0.0;
  std::vector<double> lnR;
  double tailWorst = 0.0;  ///< largest |integrand at r_end|/(kappa |M|) seen on the grid
  /// M_PW changes sign on the grid (a node of the momentum distribution;
  /// possible with rmin > 0 at large q): dwpw is then singular there.  Only
  /// without a window of directions (one direction).
  bool pwSignChange = false;
  /// Table form (distortion=table:<file>).
  std::shared_ptr<const ThmWeightTable> table;
  mutable std::atomic<bool> warned{false};

  /// Sets up everything and fills the grid on [eLo, eHi] (MeV).  "" or what is wrong.
  std::string Build(const ThmExperiment &x, const Kinematics &k, double eLo, double eHi, double eRefDefault);
  /// The settings, kinematics, channels, bound state and radial step h (as Build
  /// uses them) without the grid of R: what the DW vertex (ThmDwVertex.h)
  /// shares.  eLo is the lowest energy (the largest k_sF).  "" or what is wrong.
  std::string Setup(const ThmExperiment &x, const Kinematics &k, double eLo);
  /// The s-x bound state phi(r) as Build tabulates it (0 below rmin).
  double Phi(double r) const;
  /// The plane-wave limit of the amplitude, M_PW(q) = 4 pi Int r^2 j_0(q r)
  /// phi(r) dr on the radial grid (Simpson), q in fm^-1.
  double PlaneWaveSource(double q) const;
  /// Sum of the Woods-Saxon depths of a channel (0 unless WOODS_SAXON; for a
  /// global potential the largest over its energies), MeV.
  static double Depth(const Channel &c);
  /// Channel 0 (a + A) or 1 (s + F) for the output; a global s + F potential
  /// at both ends of the data.
  std::string ChannelText(int which) const;
  /// The s + F channel at E (k, eta and, for a global potential, p[] at E_sF).
  Channel SfAt(double energy) const;
  /// For a global potential in channel c (0 a + A, 1 s + F): the projectile
  /// lab energies and the ten numbers at the data ends E = lo and hi (for
  /// a + A both are the same).  False if the channel is not global.
  bool GlobalEnds(int c, double lo, double hi, double elab[2], double p[2][10]) const;
  /// Direct evaluation at E (not the grid); with a spectator-direction
  /// window the acceptance averages (EvaluateWindow).
  Point Evaluate(double energy) const;
  /// The window at E: M and M_PW at every node from one set of radial
  /// integrals (only P_l(x) depends on the direction).
  Point EvaluateWindow(double energy) const;
  /// rho(E) as R uses it, from a Point.
  double Ratio(const Point &p) const {
    if (p.nodes) return ratioPW ? p.m2 / p.mpw2 : p.m2;
    return ratioPW ? std::norm(p.m) / (p.mpw * p.mpw) : std::norm(p.m);
  }
  /// R(E) = rho(E)/rho(E_ref) directly.
  double R(const Point &p) const { return Ratio(p) / refRatio; }
  /// The weight applied to the model: interpolated R (or the table).  *outside
  /// is set when E lies beyond the grid (the end value is used).
  double Weight(double energy, bool *outside = nullptr) const;
  /// Kinematics at E only (E_sF > 0, eta_sF, a reachable lab angle): "" or what is wrong.
  std::string CheckEnergy(double energy) const;
  /// cos of the angle between k_sF and k_aA at the spectator direction of
  /// spectatorAngle= (Evaluate), for the s + F wave number ksf (fm^-1).  Sets
  /// the spectator c.m. angle to the beam (deg) and whether a lab angle was
  /// beyond reach (the largest reachable one is used) when asked.
  double SpectatorCos(double ksf, double *thetaCm = nullptr, bool *clamped = nullptr) const;
  /// E_sF(E) = E_aA - B - E.
  double EsF(double energy) const { return eAA - kin.bind - energy; }
  /// The nodes of the spectator-direction window at the s + F wave number ksf
  /// (fm^-1): angSlots x angNodes entries (w = 0 in an empty interval).
  /// False if no direction is accepted.  A window of zero width is its point
  /// (one interval: equal weights; a lab angle on two branches: the limit of
  /// a shrinking window, weight |d cos(theta_cm)/d theta_lab| per branch).
  bool AngleNodes(double ksf, std::vector<AngleNode> &out) const;
  /// "" or why no direction of the window is accepted at E.
  std::string CheckWindow(double energy) const;
  /// angSlots = 2 for a lab window if gamma = V_cm/v_s exceeds 1 at E = eHi
  /// (the slowest spectator), else 1.
  void SetAngleSlots(double eHi);
  /// The window for the output ("" without one).
  std::string AngleText() const;
  /// <|M|^2> and <|M_PW|^2> of a Point (the single direction without a window).
  static double M2(const Point &p) { return p.nodes ? p.m2 : std::norm(p.m); }
  static double MPW2(const Point &p) { return p.nodes ? p.mpw2 : p.mpw * p.mpw; }

  /// u_l of channel c on the grid r_j = j step, j = 0..nStore-1 (Numerov from
  /// the origin, normalized to F_l + T_l H_l^+, T_l = (S_l - 1)/2i the nuclear
  /// part, returned in *T if asked); false if it fails.
  bool Wave(const Channel &c, int l, double step, int nStore, std::vector<complex> &u,
            complex *T = nullptr) const;
};

/// What the output file and pyazr report for one experiment's distortion.
struct ThmDistortionReport {
  std::string experiment, kind, description;
  double eRef = 0.0, eAA = 0.0, bind = 0.0, kAA = 0.0, etaAA = 0.0, kappa = 0.0, etaB = 0.0, beta = 0.0;
  std::vector<double> energy, esf, ksf, etasf, thetaCm, x, q, m2, mpw2, r, rModel;
  std::vector<int> lmax;
};

#endif
