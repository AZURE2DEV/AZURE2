#ifndef AZURE2_EXPERIMENT_VALUES_H
#define AZURE2_EXPERIMENT_VALUES_H

#include <string>
#include <vector>

#include "azure2/experiment/Units.h"

namespace azure2 {
namespace experiment {

/*!
 * The typed values of format-3 fields (docs/upgrade/project_format3.md sections 3.3-3.10, 7.3): what the one parser
 * of each value type (Parsers.h) makes of a field's tokens, so that Schema::Build (step 3) consumes values, never text.
 * Numbers keep the token text as written next to its double (section 3.3: a writer never re-prints a number it did
 * not change; strtod of that text is what the engine uses).
 */
enum class ValueType {
  Number,          ///< a real number
  Integer,         ///< an integer (with a minimum)
  Unit,            ///< one unit token
  NumberWithUnit,  ///< <number> <unit>
  Range,           ///< <lo> <hi> <unit>
  Nuisance,        ///< <value> [<unit>] fixed|fit [prior <c> <w>[%]] [grid_max <v>] | profiled  (or a reference)
  Enum,            ///< one word of a closed list
  Reference,       ///< the name of another record (target names may carry @<epoch>)
  Name,            ///< a name that is not a record of this project (a <thm> experiment, a formula parameter)
  Path,            ///< a file path (quoted when it holds blanks)
  Text,            ///< free text, never read by the physics
  Formula,         ///< "<expression in x>" [<name>=<number> ...]
  NumberList,      ///< <number> <number> ...
  Spin,            ///< 0, 1, 3/2, ...
  SpinParity,      ///< 0+, 3-, 1/2+, ...
  Nucleus,         ///< <A><Symbol>: 13C
  Pair,            ///< <particle>+<nucleus>[<E_x>]: p+15N, a+12C[4.439]
  Composition,     ///< <Symbol> <number> ...: C 1 O 2
  Reaction,        ///< a reaction expression (section 3.5)
  Columns,         ///< the column roles of a data file (section 3.4.1)
  Grid,            ///< an energy or angle grid of an extrapolation (section 3.8)
  Windows,         ///< lab-energy windows of a fold: <lo>-<hi> ... <unit> | auto[:tol[:margin]]
  Form,            ///< a keyword-led structure of the above (stopping, udr, quadrature, ...): Alternative + parts
};

const char *ValueTypeName(ValueType t);

struct Number {
  std::string text;  ///< as written
  double value = 0;
  int col = 0;
};

enum class NuisanceMode { Fixed, Fit, Profiled };
const char *NuisanceModeName(NuisanceMode m);

struct Nuisance {
  bool hasValue = false;  ///< false only for 'profiled'
  Number value;
  const Unit *unit = nullptr;  ///< nullptr for a dimensionless nuisance (no unit token)
  NuisanceMode mode = NuisanceMode::Fixed;
  bool hasPrior = false;
  Number priorCentre;
  Number priorWidth;  ///< in the nuisance's unit, or in per cent of the centre when priorPercent
  bool priorPercent = false;
  bool hasGridMax = false;
  Number gridMax;  ///< in the nuisance's unit
  /// A free parameter without a prior (mode fit, no prior).
  bool FreeWithoutPrior() const { return mode == NuisanceMode::Fit && !hasPrior; }
};

struct Formula {
  std::string expression;  ///< the text between the quotes, in x
  struct Param {
    std::string name;
    Number value;
  };
  std::vector<Param> params;
};

struct Spin {
  int twice = 0;     ///< 2J
  int parity = 0;    ///< +1, -1, or 0 for a spin without parity
  std::string text;  ///< as written
};

/// A nucleus or light particle of a reaction or pair label, or a nuclide (13C).
struct Particle {
  std::string text;   ///< as written
  int Z = 0;
  int A = 0;
  std::string light;  ///< "p", "n", "d", "t", "h", "a", "g", "e-", "b-", "b+" when written as a letter; else ""
  bool excited = false;
  Number excitation;  ///< MeV, when excited ("a[3.0]")
  bool Same(const Particle &o) const { return Z == o.Z && A == o.A && light == o.light; }
};

/// The residual state of a reaction: ground ([0], absent, or the alias a0/p0/n0/g0), [<E_x in MeV>], or [total].
struct ResidualState {
  enum Kind { Ground, Excited, Total } kind = Ground;
  Number energy;  ///< MeV, for Excited: matched to the digits written (step 3)
};

/// One reaction <target>(<projectile>,<ejectile>)<residual>[<state>].
struct SingleReaction {
  Particle target, projectile, ejectile, residual;
  ResidualState targetState;  ///< Ground or Excited: the state of the target nucleus (an excited entrance pair)
  ResidualState state;
  std::string text;
  int col = 0;
  bool Elastic() const;           ///< same particles in and out, the residual in the target's state, no particle excitation
  bool ChargedEntrance() const;   ///< both entrance particles charged
  bool IdenticalEntrance() const; ///< projectile and target the same nucleus
  bool GammaExit() const;         ///< the ejectile is a photon
  std::string Canonical() const;  ///< the label with aliases normalised (a0 -> a, [0] dropped)
};

struct ReactionTerm {
  bool hasScale = false;
  Number scale;
  SingleReaction reaction;
  bool hasAngle = false;  ///< "@ <angle> <unit> <frame>": the term at a fixed angle
  Number angle;
  const Unit *angleUnit = nullptr;
  std::string angleFrame;  ///< "lab" or "cm"
};

/// The syntax tree of section 3.5: one reaction, a Sum of terms (each possibly Scaled), or a Ratio of two terms.
struct ReactionExpr {
  enum Op { Single, Sum, Ratio } op = Single;
  std::vector<ReactionTerm> terms;  ///< the base term first; Ratio: numerator, denominator
  /// The tree as text: "Sum(15N(p,a)12C, Scale(0.5, 20Ne(n,a)17O))", "Ratio(..., ... @ 90 deg cm)".
  std::string Tree() const;
};

struct ParticlePair {
  Particle light, heavy;
  ResidualState state;  ///< of the heavy particle
};

struct ColumnRole {
  std::string role;            ///< E_lab, E_cm, E_decay, theta_lab, theta_cm, order, value, error,
                               ///< error_percent, window_lo, window_hi, _
  const Unit *unit = nullptr;  ///< for energy and angle roles
  int col = 0;
};

struct Columns {
  bool fromFile = false;  ///< 'from_file': the roles are in the data file's header line (read in step 3)
  std::vector<ColumnRole> roles;
  const ColumnRole *Find(const std::string &role) const;
  const ColumnRole *EnergyRole() const;  ///< the E_* role, or nullptr
  const ColumnRole *AngleRole() const;   ///< the theta_* role, or nullptr
};

struct Grid {
  std::string kind;  ///< range, list, file, adaptive, none
  Number lo, hi;
  const Unit *unit = nullptr;
  std::string frame;  ///< lab or cm
  bool hasStep = false;
  Number step;  ///< range: E_k = lo + k step (decision F10)
  bool hasPoints = false;
  long points = 0;
  bool log = false;
  std::vector<Number> values;  ///< list
  std::string path;            ///< file
  Number tolerance;            ///< adaptive
  long maxPoints = 0;
  Number startStep;
};

struct FoldWindows {
  bool automatic = false;
  bool hasTolerance = false, hasMargin = false;
  Number tolerance, margin;  ///< auto:<tolerance>:<margin>
  struct Window {
    Number lo, hi;
  };
  std::vector<Window> windows;
  const Unit *unit = nullptr;
};

struct ValuePart;

/// A parsed value: `type` says which members hold it.
struct Value {
  ValueType type = ValueType::Text;
  int line = 0;
  int col = 0;
  std::string word;   ///< Enum (canonical spelling), Reference/Name, Path, Text, Form (its alternative's keyword)
  std::string epoch;  ///< Reference: the epoch of name@epoch
  long integer = 0;
  Number number;               ///< Number, Integer (its text), NumberWithUnit
  const Unit *unit = nullptr;  ///< Unit, NumberWithUnit, Range
  Number lo, hi;               ///< Range
  Nuisance nuisance;
  std::vector<Number> numbers;  ///< NumberList
  Formula formula;
  Spin spin;
  Particle nucleus;
  ParticlePair pair;
  std::vector<std::pair<std::string, Number>> composition;  ///< (element symbol, atoms per formula unit)
  ReactionExpr reaction;
  Columns columns;
  Grid grid;
  FoldWindows windows;
  std::vector<ValuePart> parts;  ///< Form: positional parts by their name, keyword parts by their keyword

  /// A part of a Form by name, or nullptr when it was not given.
  const Value *Get(const std::string &name) const;
};

struct ValuePart {
  std::string name;
  Value value;
};

}  // namespace experiment
}  // namespace azure2

#endif
