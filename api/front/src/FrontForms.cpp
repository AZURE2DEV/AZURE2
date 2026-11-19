// The forms of values (Engine.h, docs/upgrade/qt_gui.md section 11): a field's value split into the pieces a form
// edits by the engine's own parsers (experiment/Parsers.h, the ones --validate reads it with), and written back from
// them with the record grammar's token writer (io::WriteTokens) and the one writer of a nuisance (NuisanceText).  The
// pairs of <levels> and the labels of their reactions (experiment/LevelTables.h, the converter's labels).
#include "azure2/front/Engine.h"

#include <cstdlib>
#include <string>
#include <vector>

#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/Parsers.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/Quantities.h"
#include "azure2/experiment/Schema.h"
#include "azure2/experiment/Units.h"
#include "azure2/experiment/ValueText.h"
#include "azure2/io/Document.h"
#include "azure2/io/LevelRows.h"
#include "azure2/io/RecordReader.h"

namespace azure2 {
namespace front {

namespace {

namespace ex = azure2::experiment;
using Tokens = std::vector<io::Token>;

const ex::FieldDef *FindField(const std::string &kind, const std::string &key) {
  const ex::RecordKind *k = ex::Schema::FindKind(kind);
  return k ? k->Find(key) : nullptr;
}

/// The value type of a ValueInfo's name (ValueTypeName read backwards, so that no name is spelled here).
ex::ValueType TypeOf(const std::string &name) {
  for (int t = 0; t <= static_cast<int>(ex::ValueType::Form); t++)
    if (name == ex::ValueTypeName(static_cast<ex::ValueType>(t))) return static_cast<ex::ValueType>(t);
  return ex::ValueType::Text;
}

/// The word that stands alone for the columns read from the data file's header, as the columns' parser reads it.
const char *const kFromFile = "from_file";

Tokens Lex(const std::string &text) {
  std::string comment;
  std::vector<io::RecordDiagnostic> problems;
  return io::Tokenize(text, 0, comment, problems);
}

io::Token Tok(const std::string &text, bool quoted = false) {
  io::Token t;
  t.text = text;
  t.quoted = quoted;
  return t;
}

std::string Slice(const Tokens &tokens, std::size_t from, std::size_t to) {
  return io::WriteTokens(Tokens(tokens.begin() + (long)from, tokens.begin() + (long)to));
}

bool QuotedToken(const Tokens &tokens, std::size_t from, std::size_t to, const std::string &text) {
  for (std::size_t i = from; i < to; i++)
    if (tokens[i].text == text && tokens[i].quoted) return true;
  return false;
}

/// Split the tokens the cursor holds from its position, by the value spec, into `out` (the parser decides; a form's
/// alternatives and parts are walked as ParseForm walks them, each part read by ParseValue).
bool Decompose(const ex::ValueSpec &spec, ex::TokenCursor &cur, ex::ParseContext &ctx, ValueParts &out) {
  const Tokens &tokens = *cur.tokens;
  out = ValueParts();
  out.type = ex::ValueTypeName(spec.type);
  if (spec.type == ex::ValueType::Form) {
    const ex::FormSpec &form = *spec.form;
    const ex::Alternative *alt = nullptr;
    for (const ex::Alternative &a : form.alternatives)
      if (!alt && !a.keyword.empty() && !cur.AtEnd() && !cur.Peek().quoted && cur.Peek().text == a.keyword) alt = &a;
    if (alt)
      cur.Next();
    else
      for (const ex::Alternative &a : form.alternatives)
        if (a.keyword.empty()) alt = &a;
    if (!alt) return false;
    out.alternative = alt->keyword;
    for (const ex::PartSpec &p : alt->parts) {
      if (p.keyword) continue;
      if (cur.AtEnd()) {
        if (p.optional) continue;
        return false;
      }
      NamedValue part;
      part.name = p.name;
      if (!Decompose(p.spec, cur, ctx, part.value)) return false;
      out.parts.push_back(part);
    }
    while (!cur.AtEnd()) {
      const ex::PartSpec *found = nullptr;
      for (const ex::PartSpec &p : alt->parts)
        if (p.keyword && p.name == cur.Peek().text && !cur.Peek().quoted) found = &p;
      if (!found || out.Part(found->name)) return false;
      cur.Next();
      NamedValue part;
      part.name = found->name;
      if (cur.AtEnd() || !Decompose(found->spec, cur, ctx, part.value)) return false;
      out.parts.push_back(part);
    }
    return true;
  }
  const std::size_t from = cur.pos;
  ex::Value v;
  v.line = cur.line;
  if (!ex::ParseValue(spec, cur, ctx, v)) return false;
  const std::size_t to = cur.pos;
  switch (spec.type) {
    case ex::ValueType::Path:
      out.text = v.word;
      out.quoted = QuotedToken(tokens, from, to, v.word);
      break;
    case ex::ValueType::Reference:
      out.text = v.word;
      out.epoch = v.epoch;
      break;
    case ex::ValueType::NumberWithUnit:
      out.number = tokens[from].text;
      out.unit = to > from + 1 ? tokens[from + 1].text : std::string();
      break;
    case ex::ValueType::Range:
      out.lo = tokens[from].text;
      out.hi = tokens[from + 1].text;
      out.unit = to > from + 2 ? tokens[from + 2].text : std::string();
      break;
    case ex::ValueType::Nuisance: {
      if (v.type == ex::ValueType::Reference) {
        out.reference = true;
        out.text = v.word;
        break;
      }
      const ex::Nuisance &n = v.nuisance;
      out.profiled = n.mode == ex::NuisanceMode::Profiled;
      out.fit = n.mode == ex::NuisanceMode::Fit;
      if (n.hasValue) out.number = n.value.text;
      if (n.unit) out.unit = tokens[from + 1].text;
      out.hasPrior = n.hasPrior;
      if (n.hasPrior) {
        out.priorCentre = n.priorCentre.text;
        out.priorWidth = n.priorWidth.text;
        out.priorPercent = n.priorPercent;
      }
      if (n.hasGridMax) out.gridMax = n.gridMax.text;
      break;
    }
    case ex::ValueType::Formula:
      out.expression = v.formula.expression;
      for (const ex::Formula::Param &p : v.formula.params) out.params.push_back({p.name, p.value.text});
      break;
    case ex::ValueType::NumberList:
      for (const ex::Number &n : v.numbers) out.numbers.push_back(n.text);
      break;
    case ex::ValueType::Composition:
      for (const auto &c : v.composition) out.composition.push_back({c.first, c.second.text});
      break;
    case ex::ValueType::Columns:
      out.fromFile = v.columns.fromFile;
      if (!out.fromFile)
        for (std::size_t i = from; i < to; i++) {
          const std::string &t = tokens[i].text;
          const std::size_t br = t.find('[');
          if (br == std::string::npos)
            out.columns.push_back({t, ""});
          else
            out.columns.push_back({t.substr(0, br), t.substr(br + 1, t.size() - br - 2)});
        }
      break;
    case ex::ValueType::Grid: {
      const ex::Grid &g = v.grid;
      out.gridKind = g.kind;
      out.lo = g.lo.text;
      out.hi = g.hi.text;
      if (g.unit) out.unit = g.unit->token;
      out.frame = g.frame;
      out.byPoints = g.hasPoints;
      if (g.hasStep) out.step = g.step.text;
      if (g.hasPoints) out.points = std::to_string(g.points);
      out.log = g.log;
      for (const ex::Number &n : g.values) out.numbers.push_back(n.text);
      out.text = g.path;
      out.quoted = !g.path.empty() && QuotedToken(tokens, from, to, g.path);
      out.tolerance = g.tolerance.text;
      if (g.maxPoints) out.maxPoints = std::to_string(g.maxPoints);
      out.startStep = g.startStep.text;
      break;
    }
    default:
      out.text = Slice(tokens, from, to);
      break;
  }
  return true;
}

void Append(Tokens &to, const std::string &text) {
  for (const io::Token &t : Lex(text)) to.push_back(t);
}

void AppendWord(Tokens &to, const std::string &text) {
  if (!text.empty()) to.push_back(Tok(text));
}

void AppendTokens(const ValueInfo &info, const ValueParts &p, Tokens &out) {
  switch (TypeOf(info.type)) {
    case ex::ValueType::Path:
      if (!p.text.empty()) out.push_back(Tok(p.text, p.quoted));
      return;
    case ex::ValueType::Reference:
      AppendWord(out, p.text + (p.epoch.empty() ? "" : "@" + p.epoch));
      return;
    case ex::ValueType::NumberWithUnit:
      AppendWord(out, p.number);
      AppendWord(out, p.unit);
      return;
    case ex::ValueType::Range:
      AppendWord(out, p.lo);
      AppendWord(out, p.hi);
      AppendWord(out, p.unit);
      return;
    case ex::ValueType::Nuisance: {
      if (p.reference) {
        AppendWord(out, p.text);
        return;
      }
      ex::NuisanceWords w;
      w.value = p.profiled ? std::string() : p.number;
      w.unit = p.profiled ? std::string() : p.unit;
      w.mode = ex::NuisanceModeName(p.profiled ? ex::NuisanceMode::Profiled : p.fit ? ex::NuisanceMode::Fit : ex::NuisanceMode::Fixed);
      if (p.hasPrior && !p.profiled) {
        w.priorCentre = p.priorCentre;
        w.priorWidth = p.priorWidth + (p.priorPercent ? "%" : "");
      }
      if (!p.profiled) w.gridMax = p.gridMax;
      Append(out, ex::NuisanceText(w));
      return;
    }
    case ex::ValueType::Formula:
      out.push_back(Tok(p.expression, true));
      for (const auto &param : p.params) AppendWord(out, param.first + "=" + param.second);
      return;
    case ex::ValueType::NumberList:
      for (const std::string &n : p.numbers) AppendWord(out, n);
      return;
    case ex::ValueType::Composition:
      for (const auto &c : p.composition) {
        AppendWord(out, c.first);
        AppendWord(out, c.second);
      }
      return;
    case ex::ValueType::Columns:
      if (p.fromFile) {
        AppendWord(out, kFromFile);
        return;
      }
      for (const auto &c : p.columns) AppendWord(out, c.first + (c.second.empty() ? "" : "[" + c.second + "]"));
      return;
    case ex::ValueType::Grid: {
      AppendWord(out, p.gridKind);
      const std::vector<GridKind> kinds = GridKinds(true);
      const GridKind *k = nullptr;
      for (const GridKind &g : kinds)
        if (g.word == p.gridKind) k = &g;
      if (!k || !k->unit) return;
      if (k->span) {
        AppendWord(out, p.lo);
        AppendWord(out, p.hi);
      }
      if (k->list)
        for (const std::string &n : p.numbers) AppendWord(out, n);
      if (k->file && !p.text.empty()) out.push_back(Tok(p.text, p.quoted));
      AppendWord(out, p.unit);
      AppendWord(out, p.frame);
      if (k->spacing) {
        if (p.byPoints) {
          out.push_back(Tok("points"));
          AppendWord(out, p.points);
          if (p.log) out.push_back(Tok("log"));
        } else {
          out.push_back(Tok("step"));
          AppendWord(out, p.step);
        }
      }
      if (k->adaptive) {
        out.push_back(Tok("tolerance"));
        AppendWord(out, p.tolerance);
        out.push_back(Tok("max_points"));
        AppendWord(out, p.maxPoints);
        out.push_back(Tok("start_step"));
        AppendWord(out, p.startStep);
      }
      return;
    }
    case ex::ValueType::Form: {
      const AlternativeInfo *alt = nullptr;
      for (const AlternativeInfo &a : info.alternatives)
        if (a.keyword == p.alternative) alt = &a;
      if (!alt) return;
      AppendWord(out, alt->keyword);
      for (const PartInfo &part : alt->parts)
        if (!part.keyword)
          if (const ValueParts *given = p.Part(part.name)) AppendTokens(part.value, *given, out);
      for (const PartInfo &part : alt->parts)
        if (part.keyword)
          if (const ValueParts *given = p.Part(part.name)) {
            Tokens value;
            AppendTokens(part.value, *given, value);
            if (value.empty()) continue;  // an optional part left empty is not written
            out.push_back(Tok(part.name));
            out.insert(out.end(), value.begin(), value.end());
          }
      return;
    }
    default:
      Append(out, p.text);
      return;
  }
}

const FieldInfo *FieldInfoOf(const std::vector<KindInfo> &kinds, const std::string &kind, const std::string &key) {
  for (const KindInfo &k : kinds)
    if (k.kind == kind)
      for (const FieldInfo &f : k.fields)
        if (f.key == key) return &f;
  return nullptr;
}

std::string Short(const ValueParts &p) {
  switch (TypeOf(p.type)) {
    case ex::ValueType::Nuisance:
      if (p.reference) return p.text;
      if (p.profiled) return ex::NuisanceModeName(ex::NuisanceMode::Profiled);
      return p.number + (p.unit.empty() ? "" : " " + p.unit) + (p.fit ? " (" + std::string(ex::NuisanceModeName(ex::NuisanceMode::Fit)) + ")" : "");
    case ex::ValueType::Formula:
      return "\"" + p.expression + "\"";
    case ex::ValueType::Form: {
      // a form of three parts or more names them ("xi 0.34 MeV, omega 0.035 MeV, ...")
      const bool named = p.parts.size() >= 3;
      std::string s = p.alternative;
      for (std::size_t i = 0; i < p.parts.size(); i++)
      {
        // a number is named ("omega 0.035 MeV"); an expression, a unit or a word speaks for itself
        const ex::ValueType t = TypeOf(p.parts[i].value.type);
        const bool number = t == ex::ValueType::Number || t == ex::ValueType::Integer || t == ex::ValueType::NumberWithUnit ||
                            t == ex::ValueType::Nuisance;
        s += std::string(s.empty() ? "" : i && named ? ", " : " ") + (named && number ? p.parts[i].name + " " : "") +
             Short(p.parts[i].value);
      }
      return s;
    }
    case ex::ValueType::Reference:
      return p.text + (p.epoch.empty() ? "" : "@" + p.epoch);
    case ex::ValueType::NumberWithUnit:
      return p.number + " " + p.unit;
    case ex::ValueType::Range:
      return p.lo + " - " + p.hi + " " + p.unit;
    case ex::ValueType::Grid:
      if (!p.lo.empty() && !p.hi.empty()) return p.gridKind + " " + p.lo + " - " + p.hi + " " + p.unit + " " + p.frame;
      if (!p.numbers.empty()) return p.gridKind + " " + p.numbers.front() + " ... " + p.numbers.back() + " " + p.unit + " " + p.frame;
      return p.gridKind + (p.text.empty() ? "" : " " + p.text);
    case ex::ValueType::Path:
      return p.text;
    case ex::ValueType::NumberList: {
      std::string s;
      for (const std::string &n : p.numbers) s += (s.empty() ? "" : " ") + n;
      return s;
    }
    case ex::ValueType::Composition: {
      std::string s;
      for (const auto &c : p.composition) s += (s.empty() ? "" : " ") + c.first + " " + c.second;
      return s;
    }
    default:
      return p.text;
  }
}

std::string WithoutWhere(const std::string &message) {
  // the parser's messages start "dataset '', field 'norm': "; a form shows the message beside its field
  const std::size_t at = message.find("': ");
  const std::size_t field = message.find("field '");
  if (field != std::string::npos && at != std::string::npos && at > field) return message.substr(at + 3);
  return message;
}

ex::LevelTables TablesOf(const Project &project) {
  const io::Document doc = io::Document::Parse(project.Text());
  bool ok = false;
  return ex::LevelTablesOfRows(ex::LevelRowsOf(doc, ok));
}

}  // namespace

const ValueParts *ValueParts::Part(const std::string &name) const {
  for (const NamedValue &p : parts)
    if (p.name == name) return &p.value;
  return nullptr;
}

bool ParseFieldValue(const std::string &kind, const std::string &key, const std::string &text,
                     const std::string &quantity, ValueParts &out, std::string &why) {
  out = ValueParts();
  const ex::FieldDef *f = FindField(kind, key);
  if (!f) {
    why = "'" + key + "' is no field of a " + kind;
    return false;
  }
  std::string comment;
  std::vector<io::RecordDiagnostic> problems;
  const Tokens tokens = io::Tokenize(text, 1, comment, problems);
  if (!problems.empty()) {
    why = problems.front().message;
    return false;
  }
  if (!comment.empty()) {
    why = "a value holds no '#' (it starts a comment)";
    return false;
  }
  if (tokens.empty()) {
    why = "the value is missing: " + ex::Grammar(f->spec);
    return false;
  }
  std::vector<ex::Diagnostic> diagnostics;
  ex::ParseContext ctx;
  ctx.out = &diagnostics;
  ctx.base.recordKind = kind;
  ctx.base.field = key;
  ctx.quantity = quantity.empty() ? nullptr : ex::FindQuantity(quantity);
  // the engine's verdict on the whole value
  ex::TokenCursor cur;
  cur.tokens = &tokens;
  cur.line = 1;
  cur.endCol = (int)text.size() + 1;
  ex::Value whole;
  if (!ex::ParseValue(f->spec, cur, ctx, whole)) {
    why = diagnostics.empty() ? "the value cannot be read: " + ex::Grammar(f->spec) : WithoutWhere(diagnostics.front().message);
    return false;
  }
  if (!cur.AtEnd()) {
    why = "'" + cur.Peek().text + "' is one token too many (" + ex::Grammar(f->spec) + ")";
    return false;
  }
  // its pieces
  ex::TokenCursor walk;
  walk.tokens = &tokens;
  walk.line = 1;
  walk.endCol = cur.endCol;
  diagnostics.clear();
  ctx.failed = false;
  if (!Decompose(f->spec, walk, ctx, out) || !walk.AtEnd()) {
    why = "the value cannot be split into its parts";
    return false;
  }
  return true;
}

std::string ValueText(const ValueInfo &info, const ValueParts &parts) {
  Tokens tokens;
  AppendTokens(info, parts, tokens);
  return io::WriteTokens(tokens);
}

std::string FieldValueText(const std::string &kind, const std::string &key, const ValueParts &parts) {
  const std::vector<KindInfo> kinds = RecordKinds();
  const FieldInfo *f = FieldInfoOf(kinds, kind, key);
  return f ? ValueText(f->value, parts) : std::string();
}

std::string ShortValueText(const std::string &kind, const std::string &key, const std::string &text) {
  ValueParts parts;
  std::string why;
  if (!ParseFieldValue(kind, key, text, "", parts, why)) return text;
  return Short(parts);
}

std::vector<GridKind> GridKinds(bool withNone) {
  // the kinds of section 3.8, as the grid's parser (ParseGrid) reads them
  std::vector<GridKind> out;
  GridKind range{"range", "from lo to hi, by a step or by a number of points (log: evenly in log E)", true, true, false, false, false, true};
  GridKind list{"list", "the values listed", false, false, true, false, false, true};
  GridKind file{"file", "the values of a file, one per line", false, false, false, true, false, true};
  GridKind adaptive{"adaptive", "points placed where the curve bends, from lo to hi", true, false, false, false, true, true};
  out = {range, list, file, adaptive};
  if (withNone) out.push_back({"none", "no grid", false, false, false, false, false, false});
  return out;
}

std::vector<std::string> Frames() { return {"lab", "cm"}; }

std::vector<std::string> DimensionUnits(const std::string &dimension) { return ex::UnitsOfDimension(dimension); }

std::vector<std::pair<std::string, std::string>> ColumnRoles() {
  std::vector<std::pair<std::string, std::string>> out = ex::AxisRoles();
  for (const std::string &r : ex::PlainRoles()) out.push_back({r, ""});
  return out;
}

std::vector<std::pair<std::string, std::string>> DefaultColumns(const std::string &quantity) {
  std::vector<std::pair<std::string, std::string>> out;
  const ex::Quantity *q = ex::FindQuantity(quantity);
  if (!q) return out;
  auto unitOf = [](const std::string &role) {
    for (const auto &r : ex::AxisRoles())
      if (r.first == role) {
        const std::vector<std::string> units = ex::UnitsOfDimension(r.second);
        return units.empty() ? std::string() : units.front();
      }
    return std::string();
  };
  if (!q->energyRoles.empty()) out.push_back({q->energyRoles.front(), unitOf(q->energyRoles.front())});
  if (!q->angleRoles.empty()) out.push_back({q->angleRoles.front(), unitOf(q->angleRoles.front())});
  // the value and its error: the plain roles of the columns' parser
  for (const std::string &r : ex::PlainRoles())
    if (r == "value" || r == "error") out.push_back({r, ""});
  return out;
}

bool PreviewShape(const Record &record, int points, ShapeCurve &out, std::string &why) {
  out = ShapeCurve();
  // a beam's profile components and truncation, read as BuildEffects reads them
  std::vector<ex::BeamProfileComponent> components;
  double truncation = 0.0;
  auto energy = [](const ValueParts *p, double &v) {
    if (!p) return false;
    const ex::Unit *u = ex::FindUnit(p->unit);
    v = std::strtod(p->number.c_str(), nullptr);
    if (u) v = ex::ToBase(v, u);
    return u != nullptr;
  };
  for (const Field &f : record.fields) {
    ValueParts parts;
    if (record.kind == "beam" && f.key == "truncation") truncation = std::strtod(f.value.c_str(), nullptr);
    if (record.kind != "beam" || f.key != "profile") continue;
    if (!ParseFieldValue(record.kind, f.key, f.value, "", parts, why)) return false;
    ex::BeamProfileComponent c{};
    const ValueParts *alpha = parts.Part("alpha"), *weight = parts.Part("weight");
    if (!energy(parts.Part("xi"), c.xi) || !energy(parts.Part("omega"), c.omega) || !alpha || !weight) {
      why = "a profile component is incomplete";
      return false;
    }
    c.alpha = std::strtod(alpha->text.c_str(), nullptr);
    c.weight = std::strtod(weight->text.c_str(), nullptr);
    components.push_back(c);
  }
  if (components.empty()) {
    why = "no shape the engine can draw from this record alone";
    return false;
  }
  ex::KernelParameters p;
  p.profile = &components;
  p.truncation = truncation;
  double low = 0, high = 0;
  ex::BeamProfileSupport(p, low, high);
  if (!(high > low) || points < 2) {
    why = "the profile has no width";
    return false;
  }
  out.title = "beam energy profile";
  out.xLabel = "E_lab (MeV)";
  for (int i = 0; i < points; i++) {
    const double e = low + (high - low) * i / (points - 1);
    out.x.push_back(e);
    out.y.push_back(ex::BeamProfileWeight(p, e));
  }
  return true;
}

std::vector<PairChoice> Pairs(const Project &project) {
  std::vector<PairChoice> out;
  for (const ex::PairInfo &p : TablesOf(project).pairs) {
    PairChoice c;
    c.key = p.key;
    c.label = ex::PairLabel(p);
    c.kind = p.kind;
    // the compound nucleus: a photon pair's heavy particle; a beta pair's daughter (the beta's charge -1 or +1 taken
    // from the parent); else the sum of the two
    c.nucleus = p.kind == 10   ? ex::NuclideLabel(p.heavyZ, p.heavyA)
                : p.kind == 20 ? ex::NuclideLabel(p.heavyZ - p.lightZ, p.heavyA)
                               : ex::NuclideLabel(p.lightZ + p.heavyZ, p.lightA + p.heavyA);
    c.threshold = p.threshold;
    out.push_back(c);
  }
  return out;
}

std::string ReactionLabel(const Project &project, int entrance, int exit) {
  const ex::LevelTables tables = TablesOf(project);
  const ex::PairInfo *in = ex::PairByKey(entrance, tables.pairs);
  const ex::PairInfo *out = exit > 0 ? ex::PairByKey(exit, tables.pairs) : nullptr;
  if (!in || (exit > 0 && !out)) return std::string();
  return ex::ReactionLabel(*in, out);
}

bool ReactionPairs(const Project &project, const std::string &reaction, int &entrance, int &exit) {
  entrance = exit = 0;
  std::string comment;
  std::vector<io::RecordDiagnostic> problems;
  const Tokens tokens = io::Tokenize(reaction, 1, comment, problems);
  if (!problems.empty() || tokens.empty()) return false;
  std::vector<ex::Diagnostic> diagnostics;
  ex::ParseContext ctx;
  ctx.out = &diagnostics;
  ex::TokenCursor cur;
  cur.tokens = &tokens;
  cur.line = 1;
  cur.endCol = (int)reaction.size() + 1;
  ex::ReactionExpr expr;
  if (!ex::ParseReaction(cur, ctx, expr) || !cur.AtEnd() || expr.op != ex::ReactionExpr::Single || expr.terms.size() != 1 ||
      expr.terms[0].hasAngle || expr.terms[0].hasScale)
    return false;
  const ex::SingleReaction &r = expr.terms[0].reaction;
  const ex::LevelTables tables = TablesOf(project);
  const std::vector<int> in = ex::PairsNamed(r.projectile, r.target, r.targetState, tables.pairs);
  entrance = in.size() == 1 ? in[0] : 0;
  if (r.state.kind == ex::ResidualState::Total) {
    exit = -1;
  } else {
    const std::vector<int> out = ex::PairsNamed(r.ejectile, r.residual, r.state, tables.pairs);
    exit = out.size() == 1 ? out[0] : 0;
  }
  return true;
}

}  // namespace front
}  // namespace azure2
