#ifndef ISOSPINLEVELS_H
#define ISOSPINLEVELS_H

#include <Minuit2/MnUserParameters.h>

#include <set>
#include <string>
#include <vector>

#include "Constants.h"

class CNuc;
class Config;

/*!
 * Isospin-labelled levels across compound nuclei (Hale's charge-independent
 * R-matrix, as in his A=4 analysis).
 *
 * A level whose channel lines span several compound nuclei (one level ID) is
 * one isospin level.  The <isospin> block gives its isospin T.  In every
 * nucleus it lives in, its formal reduced-width amplitudes are the
 * Clebsch-Gordan projection of one amplitude per isospin channel
 * (arrangement, l, s), and its formal eigenvalue is shifted by that nucleus's
 * internal Coulomb energy shift.  The member in the reference nucleus carries
 * the free parameters; the other members' parameters are derived from them.
 * Isospin-forbidden channels (Clebsch-Gordan factor 0) are zero unless a mix
 * line gives them a free amplitude.  Every group that holds an isospin level
 * uses B_c = -l_c, so partners share their boundary values.
 *
 * Block (between </levels> or </abundances> and <segmentsData>):
 *   reference <Z> <A> <pair key>                 nucleus whose energies are the reference, and its threshold
 *   shift <Z> <A> <pair key> <value> <vary> <prior centre> <prior sigma>
 *   level <level ID> <T>
 *   mix <level ID> <pair key> <l> <s> <value> <vary> <prior sigma>
 *   factor <level ID> <pair key> <value>          overrides a computed Clebsch-Gordan factor
 *   pairT <pair key> <t1> <t2>                     particle isospins of a pair (default |N-Z|/2)
 *   soft <level ID> <sigma E> <sigma gamma>        soft tie: the level's derived slots are free,
 *                                                  with a penalty ((value - tie)/sigma)^2 each
 *
 * Everything works on the flat Minuit parameter vector: Apply() overwrites the
 * derived slots from the free ones, so every path that fills the compound
 * nucleus from a parameter vector sees consistent members.
 */

struct IsospinShift {
  int z = 0, a = 0;
  int nucleus = 0;  ///< 1-based CNuc nucleus index
  int pairKey = 0;  ///< the threshold energies are measured from
  double threshold = 0.0;  ///< Ex of that threshold (sepE + ExE)
  double value = 0.0, nominal = 0.0, sigma = 0.0;
  bool vary = false;
  bool reference = false;  ///< the reference nucleus: no shift, no parameter
  int param = -1;          ///< flat index of its parameter (shifts only)
};

struct IsospinMix {
  int levelID = 0, pairKey = 0, l = 0;
  double s = 0.0;
  double value = 0.0, sigma = 0.0;
  bool vary = false;
  int param = -1;  ///< flat index of its parameter
  int slot = -1;   ///< flat index of the width it sets
};

struct IsospinFactor {
  int levelID = 0, pairKey = 0;
  double value = 0.0;
};

struct IsospinPairT {
  int pairKey = 0;
  double t1 = 0.0, t2 = 0.0;
};

/// One width of one member: either the owner of its isospin channel (free),
/// a multiple of the owner, or isospin-forbidden (zero or a mix amplitude).
struct IsospinWidth {
  int channel = 0;      ///< 1-based channel of the member's J-group
  int slot = -1;        ///< flat index of this width
  int ownerSlot = -1;   ///< owner's flat index; -1 for the owner itself, -2 when forbidden
  double cg = 0.0;      ///< factor used
  double cgComputed = 0.0;
  bool overridden = false;
  double ratio = 1.0;   ///< cg / cg(owner)
  int mix = -1;         ///< index into mixes_, or -1
  std::string key;      ///< isospin channel (arrangement, l, s), for the report
};

struct IsospinMember {
  int nucleus = 0, jGroup = 0, level = 0;
  double t3 = 0.0;
  double threshold = 0.0;  ///< Ex of the nucleus's reference threshold
  int energySlot = -1;
  std::vector<IsospinWidth> widths;
};

/// A derived slot of the flat parameter vector as a linear function of free
/// slots: d(slot)/d(source) = coefficient.  Used to carry gradients back.
struct IsospinTie {
  int slot = -1;
  std::vector<std::pair<int, double>> sources;
};

struct IsospinLevel {
  int id = 0;
  double t = 0.0;
  int ref = 0;  ///< index of the reference member
  bool soft = false;            ///< soft tie: derived slots free, penalised
  double sigmaE = 0.0, sigmaGamma = 0.0;
  std::vector<IsospinMember> members;
};

/// One soft-tie penalty term: ((q[slot] - target)/sigma)^2 with
/// target = constant + sum coefficient * q[source].
struct IsospinSoftTerm {
  int slot = -1;
  std::vector<std::pair<int, double>> sources;
  double constant = 0.0;
  double sigma = 1.0;
};

class IsospinLevels {
 public:
  /// Reads the <isospin> block of the project file, syntax only.  0 when it
  /// is absent or read, -1 on an error (reported).
  int Read(const Config &configure);
  bool Active() const { return active_; }
  /// Whether the block declares this level ID (its lines may span nuclei).
  bool IsDeclared(int levelID) const;
  /// Resolves the declared levels against the filled compound nucleus: members,
  /// factors, owners, parameter slots, and B = -l on their groups.  Reports
  /// what it built.  -1 on an error (reported).
  int Resolve(const Config &configure, CNuc &compound);

  /// Adds the isospin parameters (shifts, then mixes) to the Minuit
  /// parameters; called right after the levels and abundances.
  void AddParams(ROOT::Minuit2::MnUserParameters &p);
  /// Number of isospin parameters.
  int NumParams() const;
  /// Flat index of isospin parameter k (0-based: the shifts, then the mixes).
  int ParamIndex(int k) const { return firstParam_ + k; }
  /// Every derived slot as a linear function of the free ones (the shifts and
  /// mixes included), for the chain rule of the gradient.
  std::vector<IsospinTie> Ties() const;
  /// Prior centre and width (0 = none) of parameter k; mixes are centred on 0.
  void ParamPrior(int k, double &centre, double &sigma) const;
  /// d(Penalty)/d(parameter k), for k = 0..NumParams()-1.
  std::vector<double> PenaltyGradient() const;
  /// f(flat index, prior centre, prior sigma) for every varied shift and mix with a prior.
  template <class F>
  void ForEachPrior(F f) const {
    if (!active_) return;
    for (const IsospinShift &e : shifts_)
      if (!e.reference && e.vary && e.sigma > 0.0) f(e.param, e.nominal, e.sigma);
    for (const IsospinMix &m : mixes_)
      if (m.vary && m.sigma > 0.0) f(m.param, 0.0, m.sigma);
  }
  /// Overwrites the derived slots of a flat parameter vector from the free ones.
  void Apply(vector_r &q) const;
  /// Stores the shift and mix values of a flat parameter vector (for the penalty and the report).
  void ReadParams(const vector_r &q);
  /// Writes the derived values into the Minuit parameters and keeps those slots fixed.
  void Sync(ROOT::Minuit2::MnUserParameters &p) const;
  /// Prior penalty on the shifts and mixes that have a prior width.
  double Penalty() const;
  bool HasPenalty() const;
  /// Soft-tie penalty at a flat parameter vector (0 without soft levels); the
  /// argument-free form uses the last ReadParams vector.
  double SoftPenalty(const vector_r &q) const;
  double SoftPenalty() const { return SoftPenalty(lastQ_); }
  bool HasSoftTies() const { return !softTerms_.empty(); }
  /// d(SoftPenalty)/d(flat index) at q, as (index, value) pairs.
  std::vector<std::pair<int, double>> SoftPenaltyGradient(const vector_r &q) const;
  std::vector<std::pair<int, double>> SoftPenaltyGradient() const { return SoftPenaltyGradient(lastQ_); }
  /// Adds the soft penalty's diagonal curvature 2 coef^2 / sigma^2 to g2 (flat indices).
  void AddSoftCurvature(std::vector<double> &g2) const;
  /// Soft levels' derived slots: set to their tie targets and freed (FillMnParams only).
  void InitSoftParams(ROOT::Minuit2::MnUserParameters &p) const;
  /// Whether a flat index is a derived (fixed, synced) slot.
  bool IsDerivedSlot(int slot) const { return derived_.count(slot) > 0; }
  /// Label of parameter flat index i if it is an isospin parameter, else "".
  std::string ParamLabel(int i) const;
  /// "ISOSPIN LEVEL <id> T=<t> (reference|derived)" for a level of a J-group, or "".
  std::string MemberLabel(int jGroup, int level) const;
  /// Writes isospin.out.
  void Write(const Config &configure, CNuc &compound) const;

 private:
  bool active_ = false;
  std::vector<IsospinShift> shifts_;
  std::vector<int> levelIDs_;
  std::vector<double> levelT_;
  std::vector<IsospinMix> mixes_;
  std::vector<IsospinFactor> factors_;
  std::vector<IsospinPairT> pairTs_;
  std::vector<IsospinLevel> levels_;
  std::set<int> derived_;
  std::vector<IsospinSoftTerm> softTerms_;
  std::vector<std::pair<int, std::pair<double, double>>> softSpec_;  // level ID -> (sigma E, sigma gamma)
  vector_r lastQ_;
  int firstParam_ = -1;

  const IsospinShift *ShiftOf(int nucleus) const;
  double ShiftValue(int nucleus, const vector_r &q) const;
};

#endif
