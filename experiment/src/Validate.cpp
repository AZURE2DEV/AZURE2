// Schema::Validate: the one table-driven validator of format-3 records (Schema.h; docs/upgrade/project_format3.md
// sections 6.2, 7.3).  Pass 1 parses every field of every record with its value type; pass 2 decides with the whole
// project in hand: which fields apply and are required (kind rows, quantity rows, conditions on other fields and on
// referenced records), units, reactions, column roles and norm modes of the quantity, names and references.
#include <algorithm>
#include <map>
#include <set>

#include "azure2/experiment/Parsers.h"
#include "azure2/experiment/Schema.h"

namespace azure2 {
namespace experiment {

namespace {

std::string Joined(const std::vector<std::string> &v, const char *sep = ", ") {
  std::string s;
  for (std::size_t i = 0; i < v.size(); i++) s += (i ? sep : "") + v[i];
  return s;
}

std::size_t Distance(const std::string &a, const std::string &b) {
  std::vector<std::size_t> row(b.size() + 1);
  for (std::size_t j = 0; j <= b.size(); j++) row[j] = j;
  for (std::size_t i = 1; i <= a.size(); i++) {
    std::size_t diag = row[0];
    row[0] = i;
    for (std::size_t j = 1; j <= b.size(); j++) {
      const std::size_t up = row[j];
      row[j] = std::min({row[j] + 1, row[j - 1] + 1, diag + (a[i - 1] == b[j - 1] ? 0 : 1)});
      diag = up;
    }
  }
  return row[b.size()];
}

// The candidate nearest to `word`, if it is close enough to be a misspelling.
std::string Nearest(const std::string &word, const std::vector<std::string> &candidates) {
  std::string best;
  std::size_t bestD = std::max<std::size_t>(2, word.size() / 3) + 1;
  for (const std::string &c : candidates) {
    const std::size_t d = Distance(word, c);
    if (d < bestD) {
      bestD = d;
      best = c;
    }
  }
  return best;
}

// One record during validation.
struct Work {
  ValidatedRecord vr;
  const io::Record *src = nullptr;
  const RecordKind *def = nullptr;
  std::set<std::string> given;   // keys written (parsed or not)
  std::set<std::string> failed;  // keys written whose value did not parse
};

class Validator {
 public:
  explicit Validator(std::vector<Diagnostic> &out) : out_(out) {}

  void Parse(const std::string &section, const io::Record &rec);
  void CheckNames();
  void CheckRecord(Work &w);
  std::vector<Work> works;
  std::vector<RecordReference> references;  ///< collected by CheckReferences

 private:
  Diagnostic Base(const Work &w, const std::string &field, int line, int col) const {
    Diagnostic d;
    d.line = line;
    d.column = col;
    d.section = w.vr.section;
    d.recordKind = w.vr.kind;
    d.record = w.vr.name;
    d.field = field;
    return d;
  }
  void Emit(Diagnostic d, const char *code, const std::string &message) {
    d.code = code;
    d.message = message;
    out_.push_back(d);
  }
  std::string Where(const Work &w, const std::string &field = "") const {
    return w.vr.kind + " '" + w.vr.name + "'" + (field.empty() ? "" : ", field '" + field + "'");
  }
  const io::Field *Source(const Work &w, const std::string &key) const {
    for (const io::Field &f : w.src->fields)
      if (f.key == key) return &f;
    return nullptr;
  }
  const Work *Resolve(const std::string &name) const {
    const auto it = byName_.find(name);
    return it == byName_.end() ? nullptr : &works[it->second];
  }
  // A referenced record named by `field`: 1 found, 0 not given, -1 given but unknown or unparsed, or a required
  // field that is missing (reported already; what depends on it is not judged).
  int Referenced(const Work &w, const std::string &field, const Work *&out) const;
  // `field` is not given although it is required unconditionally (by the kind or the quantity row).
  bool Missing(const Work &w, const std::string &field) const;
  // 1 holds, 0 does not, -1 cannot be decided (a field it depends on did not parse or names nothing).
  int Holds(const Condition &c, const Work &w) const;
  void CheckFields(Work &w);
  void CheckQuantity(Work &w);
  void CheckReaction(Work &w, const Value &v, const std::string &field);
  void CheckColumns(Work &w);
  void CheckNorm(Work &w);
  void CheckReferences(Work &w);
  void CollectRefs(const ValueSpec &spec, const Value &v, std::vector<std::pair<const Value *, std::vector<std::string>>> &out) const;

  std::vector<Diagnostic> &out_;
  std::map<std::string, std::size_t> byName_;
};

void Validator::Parse(const std::string &section, const io::Record &rec) {
  Work w;
  w.src = &rec;
  w.vr.kind = rec.kind;
  w.vr.name = rec.name;
  w.vr.section = section;
  w.vr.line = w.vr.lastLine = rec.line;
  for (const io::Field &f : rec.fields) w.vr.lastLine = std::max(w.vr.lastLine, f.line);
  w.def = Schema::FindKind(rec.kind);
  if (!w.def) {
    std::vector<std::string> kinds;
    for (const RecordKind &k : Schema::Kinds()) kinds.push_back(k.kind);
    Diagnostic d = Base(w, "", rec.line, rec.col);
    d.hasFound = true;
    d.found = rec.kind;
    d.allowed = kinds;
    Emit(d, io::record_codes::kUnknownKind, "'" + rec.kind + "' is not a record kind (" + Joined(kinds) + ")");
    return;
  }
  const std::vector<std::string> &sections = w.def->sections;
  if (std::find(sections.begin(), sections.end(), section) == sections.end()) {
    Diagnostic d = Base(w, "", rec.line, rec.col);
    d.hasFound = true;
    d.found = section;
    d.allowed = sections;
    Emit(d, codes::kWrongSection, Where(w) + " is in <" + section + ">; a " + rec.kind + " record belongs in <" +
                                      Joined(sections, "> or <") + ">");
  }
  if (!rec.name.empty() && !IsName(rec.name)) {
    Diagnostic d = Base(w, "", rec.line, rec.nameCol);
    d.hasFound = true;
    d.found = rec.name;
    d.allowed = {"[A-Za-z][A-Za-z0-9_-]*"};
    Emit(d, codes::kBadName, rec.kind + " '" + rec.name + "': a name is a letter followed by letters, digits, '_' and "
                                                          "'-' ('.' and '@' are reserved)");
  }
  std::vector<std::string> keys, sectionKinds = Schema::KindsOfSection(section);
  for (const FieldDef &f : w.def->fields) keys.push_back(f.key);
  for (const io::Field &f : rec.fields) {
    const FieldDef *def = w.def->Find(f.key);
    Diagnostic base = Base(w, f.key, f.line, f.col);
    if (!def) {
      std::string refused;
      for (const auto &r : w.def->refused)
        if (r.first == f.key) refused = r.second;
      base.hasFound = true;
      base.found = f.key;
      if (!refused.empty()) {
        Emit(base, codes::kFieldNotApplicable, Where(w, f.key) + ": " + refused);
        continue;
      }
      std::vector<std::string> candidates = keys;
      candidates.insert(candidates.end(), sectionKinds.begin(), sectionKinds.end());
      const std::string near = Nearest(f.key, candidates);
      const bool kind = !near.empty() && std::find(sectionKinds.begin(), sectionKinds.end(), near) != sectionKinds.end() &&
                        f.tokens.size() == 1;
      base.allowed = kind ? sectionKinds : keys;
      if (!near.empty()) {
        base.hasFix = true;
        base.fix = {"replace_key", f.line, f.col, 0, near};
      }
      if (kind)
        Emit(base, io::record_codes::kUnknownKind, "'" + f.key + " " + f.tokens[0].text + "' (line " +
                                                       std::to_string(f.line) + ") looks like a header: '" + f.key +
                                                       "' is not a record kind; did you mean '" + near + "'?");
      else
        Emit(base, codes::kUnknownField, Where(w) + ": unknown field '" + f.key + "'" +
                                             (near.empty() ? "" : "; did you mean '" + near + "'?") + " (fields of " +
                                             w.vr.kind + ": " + Joined(keys) + ")");
      continue;
    }
    if (w.given.count(f.key) && !def->repeatable) {
      base.hasFound = true;
      base.found = f.key;
      const io::Field *first = Source(w, f.key);
      Emit(base, codes::kFieldTwice, Where(w, f.key) + ": given twice (first on line " + std::to_string(first->line) +
                                         "); a field is written once");
      continue;
    }
    w.given.insert(f.key);
    if (f.tokens.empty()) {
      w.failed.insert(f.key);
      base.allowed = {Grammar(def->spec)};
      Emit(base, codes::kValueMissing, Where(w, f.key) + ": the value is missing: " + Grammar(def->spec));
      continue;
    }
    {
      bool placeholder = false;
      for (const io::Token &t : f.tokens)
        if (t.text == "?" && !t.quoted) {
          placeholder = true;
          base.column = t.col;
        }
      if (placeholder) {
        w.failed.insert(f.key);
        base.hasFound = true;
        base.found = "?";
        base.allowed = {Grammar(def->spec)};
        Emit(base, codes::kPlaceholder, Where(w, f.key) + " is '?': a value is needed (" + Grammar(def->spec) + ")");
        continue;
      }
    }
    {
      ParseContext ctx;
      ctx.out = &out_;
      ctx.base = base;
      TokenCursor cur;
      cur.tokens = &f.tokens;
      cur.line = f.line;
      cur.endCol = f.tokens.back().col + static_cast<int>(f.tokens.back().text.size()) + (f.tokens.back().quoted ? 2 : 0);
      ValidatedField vf;
      vf.key = f.key;
      vf.line = f.line;
      vf.col = f.col;
      bool ok = ParseValue(def->spec, cur, ctx, vf.value);
      if (ok && !cur.AtEnd()) {
        std::string rest;
        for (std::size_t i = cur.pos; i < f.tokens.size(); i++) rest += (rest.empty() ? "" : " ") + f.tokens[i].text;
        ctx.Error(codes::kTooManyTokens, cur.Peek().col, rest, {Grammar(def->spec)},
                  Where(w, f.key) + ": '" + rest + "' after the value; the value is " + Grammar(def->spec));
        ok = false;
      }
      if (ok)
        w.vr.fields.push_back(vf);
      else
        w.failed.insert(f.key);
    }
  }
  if (const Value *q = w.vr.Get("quantity")) w.vr.quantity = FindQuantity(q->word);
  works.push_back(std::move(w));
}

void Validator::CheckNames() {
  for (std::size_t i = 0; i < works.size(); i++) {
    const Work &w = works[i];
    if (w.vr.name.empty() || !w.def) continue;
    const auto it = byName_.find(w.vr.name);
    if (it == byName_.end()) {
      byName_[w.vr.name] = i;
      continue;
    }
    const Work &first = works[it->second];
    Diagnostic d = Base(w, "", w.vr.line, w.src->nameCol);
    d.hasFound = true;
    d.found = w.vr.name;
    Emit(d, codes::kDuplicateName, w.vr.kind + " '" + w.vr.name + "': the name is already used by the " + first.vr.kind +
                                       " on line " + std::to_string(first.vr.line) +
                                       "; names are unique across the project (F3)");
  }
}

bool Validator::Missing(const Work &w, const std::string &field) const {
  if (w.given.count(field)) return false;
  const FieldDef *f = w.def->Find(field);
  if (!f) return false;
  if (!f->byQuantity) return f->applies.kind == Condition::Always && !f->applies.negate && f->required.kind == Condition::Always && !f->required.negate;
  if (!w.vr.quantity) return false;
  for (const QuantityField &qf : w.vr.quantity->fields)
    if (qf.field == field && qf.required && qf.when.kind == Condition::Always && !qf.when.negate) return true;
  return false;
}

int Validator::Referenced(const Work &w, const std::string &field, const Work *&out) const {
  out = nullptr;
  if (!w.given.count(field)) return Missing(w, field) ? -1 : 0;
  const Value *v = w.vr.Get(field);
  if (!v || v->type != ValueType::Reference) return -1;
  out = Resolve(v->word);
  return out ? 1 : -1;
}

int Validator::Holds(const Condition &c, const Work &w) const {
  auto result = [&](int r) { return r < 0 ? -1 : ((r == 1) != c.negate ? 1 : 0); };
  switch (c.kind) {
    case Condition::Always:
      return result(1);
    case Condition::FieldIs: {
      if (!w.given.count(c.field)) return Missing(w, c.field) ? -1 : result(0);
      const Value *v = w.vr.Get(c.field);
      if (!v) return -1;
      return result(std::find(c.values.begin(), c.values.end(), v->word) != c.values.end() ? 1 : 0);
    }
    case Condition::FieldPresent:
      return result(w.given.count(c.field) ? 1 : 0);
    case Condition::UnitDimension: {
      if (!w.given.count("unit")) return result(0);
      const Value *v = w.vr.Get("unit");
      if (!v || !v->unit) return -1;
      return result(std::find(c.values.begin(), c.values.end(), v->unit->dimension) != c.values.end() ? 1 : 0);
    }
    case Condition::RatioReaction:
    case Condition::SingleReaction:
    case Condition::GammaExit: {
      const Value *v = w.vr.Get("reaction");
      if (!v) return w.given.count("reaction") ? -1 : result(0);
      if (c.kind == Condition::RatioReaction) return result(v->reaction.op == ReactionExpr::Ratio ? 1 : 0);
      if (c.kind == Condition::SingleReaction) return result(v->reaction.terms.size() == 1 ? 1 : 0);
      bool all = true;
      for (const ReactionTerm &t : v->reaction.terms) all = all && t.reaction.GammaExit();
      return result(all ? 1 : 0);
    }
    case Condition::AngleAxis: {
      if (w.vr.kind == "extrapolation") {
        const Value *g = w.vr.Get("angles");
        if (!g) return w.given.count("angles") ? -1 : result(0);
        return result(g->grid.kind != "none" ? 1 : 0);
      }
      const Value *v = w.vr.Get("columns");
      if (!v) return w.given.count("columns") ? -1 : result(0);
      if (v->columns.fromFile) return -1;  // the roles are in the data file (step 3)
      return result(v->columns.AngleRole() ? 1 : 0);
    }
    case Condition::ThmExperiment: {
      if (!w.given.count("thm")) return result(0);
      const Value *v = w.vr.Get("thm");
      if (!v) return -1;
      return result(v->word != "global" ? 1 : 0);
    }
    case Condition::Fold: {
      const Work *t = nullptr, *b = nullptr, *d = nullptr;
      const int rt = Referenced(w, "target", t), rb = Referenced(w, "beam", b), rd = Referenced(w, "detector", d);
      if (rt < 0 || rb < 0 || rd < 0) return -1;
      bool fold = rb == 1;
      if (t) {
        const Value *kind = t->vr.Get("kind");
        if (!kind) return -1;
        // A layer or an infinite target integrates over the energy loss; a neutron sample enters through the
        // transmission or yield transform (no sub-points) unless it Doppler-broadens.
        fold = fold || kind->word == "layer" || kind->word == "infinite" || t->given.count("doppler");
      }
      if (d) fold = fold || d->given.count("energy_resolution") || d->given.count("udr");
      return result(fold ? 1 : 0);
    }
  }
  return -1;
}

void Validator::CheckFields(Work &w) {
  const Quantity *q = w.vr.quantity;
  for (const FieldDef &f : w.def->fields) {
    int applies = 1, required = 0;
    std::string why;
    if (f.byQuantity) {
      if (!q) {
        continue;  // the quantity is missing or unknown, which is reported: what depends on it is not judged
      } else {
        const QuantityField *rule = nullptr;
        for (const QuantityField &qf : q->fields)
          if (qf.field == f.key) rule = &qf;
        if (!rule) {
          applies = 0;
          std::vector<std::string> allowed;
          for (const QuantityField &qf : q->fields)
            if (w.def->Find(qf.field)) allowed.push_back(qf.field);
          why = "quantity '" + q->name + "' does not use it" + (allowed.empty() ? "" : " (it uses " + Joined(allowed) + ")");
        } else {
          applies = Holds(rule->when, w);
          required = rule->required ? applies : 0;
          why = "for quantity '" + q->name + "' it applies only when " + Describe(rule->when);
        }
      }
    } else {
      applies = Holds(f.applies, w);
      required = applies == 1 ? Holds(f.required, w) : 0;
      why = "it applies only when " + Describe(f.applies);
    }
    const bool given = w.given.count(f.key) > 0;
    if (given && applies == 0) {
      const io::Field *src = Source(w, f.key);
      Diagnostic d = Base(w, f.key, src->line, src->col);
      d.hasFound = true;
      d.found = f.key;
      Emit(d, codes::kFieldNotApplicable, Where(w, f.key) + " does not belong here: " + why +
                                              " (a value the engine would not use is refused, not ignored)");
    }
    if (!given && applies == 1 && required == 1) {
      Diagnostic d = Base(w, f.key, w.vr.line, 1);
      d.allowed = {Grammar(f.spec)};
      d.hasFix = true;
      std::string key = f.key;
      d.fix = {"insert_field", 0, 0, w.vr.lastLine, "    " + key + std::string(key.size() < 15 ? 15 - key.size() : 1, ' ') + "?"};
      std::string reason = f.byQuantity ? "quantity '" + q->name + "' requires it" : "it is required";
      if (f.byQuantity) {
        for (const QuantityField &qf : q->fields)
          if (qf.field == f.key && qf.when.kind != Condition::Always) reason += " when " + Describe(qf.when);
      } else if (f.applies.kind != Condition::Always) {
        reason += " when " + Describe(f.applies);
      } else if (f.required.kind != Condition::Always) {
        reason += " when " + Describe(f.required);
      }
      for (const Refusal &r : q ? q->refusals : std::vector<Refusal>{})
        if (r.code == codes::kFieldRequired) reason += ": " + r.message;
      Emit(d, codes::kFieldRequired, Where(w) + ": field '" + f.key + "' is missing; " + reason + " (" + Grammar(f.spec) + ")");
    }
  }
  for (const std::vector<std::string> &group : w.def->oneOf) {
    bool any = false;
    for (const std::string &k : group) any = any || w.given.count(k);
    if (any) continue;
    Diagnostic d = Base(w, group[0], w.vr.line, 1);
    d.allowed = group;
    Emit(d, codes::kFieldRequired, Where(w) + " has none of " + Joined(group) + "; it needs at least one");
  }
}

void Validator::CheckReaction(Work &w, const Value &v, const std::string &field) {
  const Quantity &q = *w.vr.quantity;
  const ReactionRule &rule = q.reaction;
  auto refuse = [&](const std::string &ruleName, const SingleReaction *r, const std::string &fallback) {
    std::string message = fallback;
    for (const Refusal &x : q.refusals)
      if (x.rule == ruleName) message = x.message;
    Diagnostic d = Base(w, field, v.line, r ? r->col : v.col);
    d.hasFound = true;
    d.found = r ? r->text : v.reaction.Tree();
    d.allowed = {rule.text.empty() ? "one reaction" : rule.text};
    Emit(d, codes::kReactionNotAllowed, Where(w, field) + ": quantity '" + q.name + "' with " + d.found + ": " + message);
  };
  if ((v.reaction.op == ReactionExpr::Sum && !rule.Has(ReactionRule::kSum)) ||
      (v.reaction.op == ReactionExpr::Ratio && !rule.Has(ReactionRule::kRatio)))
    return refuse("composite", nullptr, std::string("takes one reaction, not a ") +
                                             (v.reaction.op == ReactionExpr::Sum ? "sum" : "ratio"));
  for (const ReactionTerm &t : v.reaction.terms) {
    const SingleReaction &r = t.reaction;
    if (rule.Has(ReactionRule::kElastic) && !r.Elastic()) return refuse("inelastic", &r, "is defined for elastic scattering only");
    if (rule.Has(ReactionRule::kChargedEntrance) && !r.ChargedEntrance())
      return refuse("neutral", &r, "needs two charged entrance particles");
    if (rule.Has(ReactionRule::kNeutralEntrance) && r.ChargedEntrance())
      return refuse("charged_entrance", &r, "needs an uncharged entrance particle");
    if (rule.Has(ReactionRule::kIdentical) && !r.IdenticalEntrance())
      return refuse("non_identical", &r, "needs identical particles");
    if (rule.Has(ReactionRule::kNotIdentical) && r.IdenticalEntrance())
      return refuse("identical", &r, "is not defined for identical particles");
    if (rule.Has(ReactionRule::kNoChargedElastic) && r.Elastic() && r.ChargedEntrance())
      return refuse("charged_elastic", &r, "is not defined for charged-particle elastic scattering");
    if (rule.Has(ReactionRule::kNoTotal) && r.state.kind == ResidualState::Total)
      return refuse("total", &r, "does not take the [total] exit");
    if (rule.Has(ReactionRule::kSpinHalfProjectile) &&
        (r.projectile.light == "d" || r.projectile.light == "a" || r.projectile.light == "g"))
      return refuse("spin", &r, "needs a spin-1/2 projectile");
  }
}

void Validator::CheckColumns(Work &w) {
  const Quantity &q = *w.vr.quantity;
  const Value *v = w.vr.Get("columns");
  if (!v || v->columns.fromFile) return;
  const Columns &c = v->columns;
  const io::Field *src = Source(w, "columns");
  auto problem = [&](const char *code, int col, const std::string &found, std::vector<std::string> allowed,
                     const std::string &message) {
    Diagnostic d = Base(w, "columns", src->line, col ? col : src->col);
    d.hasFound = !found.empty();
    d.found = found;
    d.allowed = std::move(allowed);
    Emit(d, code, Where(w, "columns") + ": " + message);
  };
  for (const ColumnRole &r : c.roles) {
    const bool energy = r.role.compare(0, 2, "E_") == 0, angle = r.role.compare(0, 6, "theta_") == 0;
    if (energy && std::find(q.energyRoles.begin(), q.energyRoles.end(), r.role) == q.energyRoles.end())
      problem(codes::kRoleNotAllowed, r.col, r.role, q.energyRoles,
              "quantity '" + q.name + "' takes the energy column " + Joined(q.energyRoles, " or ") + ", not " + r.role);
    if (angle && (q.angle == Quantity::NoAngle ||
                  std::find(q.angleRoles.begin(), q.angleRoles.end(), r.role) == q.angleRoles.end()))
      problem(codes::kRoleNotAllowed, r.col, r.role, q.angleRoles,
              q.angle == Quantity::NoAngle ? "quantity '" + q.name + "' is angle-integrated; write '_' for a column that is not read"
                                           : "quantity '" + q.name + "' takes the angle column " + Joined(q.angleRoles, " or "));
    if (r.role == "order" && !q.orderRole)
      problem(codes::kRoleNotAllowed, r.col, r.role, {}, "only legendre data have an order column");
    if ((r.role == "window_lo" || r.role == "window_hi") && !q.windowRoles)
      problem(codes::kRoleNotAllowed, r.col, r.role, {}, "quantity '" + q.name + "' has no beam-profile window columns");
  }
  if (!c.EnergyRole()) problem(codes::kRoleRequired, 0, "", q.energyRoles, "an energy column (" + Joined(q.energyRoles, " or ") + "[<unit>]) is required");
  if (q.angle == Quantity::RequiredAngle && !c.AngleRole())
    problem(codes::kRoleRequired, 0, "", q.angleRoles, "quantity '" + q.name + "' needs an angle column (" + Joined(q.angleRoles, " or ") + "[<unit>])");
  if (q.orderRole && !c.Find("order")) problem(codes::kRoleRequired, 0, "", {"order"}, "legendre data need the column 'order'");
  if (!c.Find("value")) problem(codes::kRoleRequired, 0, "", {"value"}, "the column 'value' is required");
  if (!c.Find("error") && !c.Find("error_percent"))
    problem(codes::kRoleRequired, 0, "", {"error", "error_percent"}, "an uncertainty column (error or error_percent) is required");
  if (c.Find("error") && c.Find("error_percent"))
    problem(codes::kRoleNotAllowed, c.Find("error_percent")->col, "error_percent", {"error"}, "give 'error' or 'error_percent', not both");
  if (q.windowRoles && (c.Find("window_lo") != nullptr) != (c.Find("window_hi") != nullptr))
    problem(codes::kRoleRequired, 0, "", {"window_lo", "window_hi"}, "window_lo and window_hi come together");
  if (c.Find("window_lo") && q.windowRoles) {
    const Work *beam = nullptr;
    if (Referenced(w, "beam", beam) >= 0 && (!beam || !beam->given.count("profile")))
      problem(codes::kRoleNotAllowed, c.Find("window_lo")->col, "window_lo", {},
              "window columns are read by a beam-profile kernel: the data set needs a beam with 'profile'");
  }
}

void Validator::CheckNorm(Work &w) {
  const Value *v = w.vr.Get("norm");
  if (!v) return;
  const Nuisance *n = &v->nuisance;
  std::string via;
  if (v->type == ValueType::Reference) {
    const Work *r = Resolve(v->word);
    const Value *rv = r && r->vr.kind == "normalization" ? r->vr.Get("value") : nullptr;
    if (!rv) return;
    n = &rv->nuisance;
    via = " (normalization record '" + v->word + "')";
  }
  const io::Field *src = Source(w, "norm");
  Diagnostic d = Base(w, "norm", src->line, src->col);
  const std::string mode = NuisanceModeName(n->mode);
  // Every THM data set with a fitted norm has it profiled (ESegment::IsProfiledNorm), in an experiment or not.
  if (n->mode == NuisanceMode::Profiled && Holds(Condition::Present("thm"), w) == 0) {
    d.hasFound = true;
    d.found = mode;
    d.allowed = {"fixed", "fit"};
    Emit(d, codes::kModeNotAllowed, Where(w, "norm") + ": 'profiled' is the norm mode of THM data sets "
                                                       "(thm <experiment> or thm global, section 3.9)");
    return;
  }
  std::vector<const NormRule *> rules;
  for (const NormRule &r : w.vr.quantity->norm) rules.push_back(&r);
  for (const NormRule &r : GeneralNormRules()) rules.push_back(&r);
  for (const NormRule *r : rules) {
    if (Holds(r->when, w) != 1) continue;
    const bool modeOk = std::find(r->modes.begin(), r->modes.end(), mode) != r->modes.end();
    if (!modeOk || (n->hasPrior && !r->prior)) {
      d.hasFound = true;
      d.found = mode + (n->hasPrior ? " prior" : "");
      d.allowed = r->modes;
      Emit(d, codes::kModeNotAllowed, Where(w, "norm") + via + ": " + mode + (n->hasPrior ? " with a prior" : "") +
                                          " is not allowed: " + r->message);
    }
    return;
  }
}

void Validator::CollectRefs(const ValueSpec &spec, const Value &v,
                            std::vector<std::pair<const Value *, std::vector<std::string>>> &out) const {
  if (v.type == ValueType::Reference) out.push_back({&v, spec.type == ValueType::Nuisance ? spec.refKinds : spec.values});
  if (spec.type != ValueType::Form) return;
  for (const Alternative &a : spec.form->alternatives) {
    if (a.keyword != v.word) continue;
    for (const PartSpec &p : a.parts)
      if (const Value *pv = v.Get(p.name)) CollectRefs(p.spec, *pv, out);
  }
}

void Validator::CheckReferences(Work &w) {
  for (const ValidatedField &f : w.vr.fields) {
    const FieldDef *def = w.def->Find(f.key);
    std::vector<std::pair<const Value *, std::vector<std::string>>> refs;
    CollectRefs(def->spec, f.value, refs);
    for (const auto &ref : refs) {
      const Value &v = *ref.first;
      references.push_back({w.vr.name, f.key, v.line, v.col, v.word, v.epoch});
      Diagnostic d = Base(w, f.key, v.line, v.col);
      d.hasFound = true;
      d.found = v.word;
      const Work *r = Resolve(v.word);
      std::vector<std::string> names;
      for (const Work &x : works)
        if (std::find(ref.second.begin(), ref.second.end(), x.vr.kind) != ref.second.end() && !x.vr.name.empty())
          names.push_back(x.vr.name);
      d.allowed = names;
      if (!r) {
        const std::string near = Nearest(v.word, names);
        if (!near.empty()) {
          d.hasFix = true;
          d.fix = {"replace_token", v.line, v.col, 0, near};
        }
        Emit(d, codes::kUnknownReference, Where(w, f.key) + ": no " + Joined(ref.second, " or ") + " record is named '" +
                                              v.word + "'" + (names.empty() ? "" : " (there are: " + Joined(names) + ")"));
        continue;
      }
      if (std::find(ref.second.begin(), ref.second.end(), r->vr.kind) == ref.second.end()) {
        Emit(d, codes::kReferenceKind, Where(w, f.key) + ": '" + v.word + "' is a " + r->vr.kind + " (line " +
                                           std::to_string(r->vr.line) + "); this field names a " + Joined(ref.second, " or "));
        continue;
      }
      if (!v.epoch.empty()) {
        std::vector<std::string> epochs;
        for (const Value *e : r->vr.GetAll("epoch"))
          if (const Value *name = e->Get("name")) epochs.push_back(name->word);
        if (std::find(epochs.begin(), epochs.end(), v.epoch) == epochs.end()) {
          d.allowed = epochs;
          d.found = v.word + "@" + v.epoch;
          Emit(d, codes::kUnknownEpoch, Where(w, f.key) + ": target '" + v.word + "' has no epoch '" + v.epoch + "'" +
                                            (epochs.empty() ? " (it has none)" : " (it has " + Joined(epochs) + ")"));
        }
      }
    }
  }
}

void Validator::CheckQuantity(Work &w) {
  const Quantity &q = *w.vr.quantity;
  // Units of the value: the unit field (data sets, extrapolations) and the value fields of a datum.
  std::vector<std::string> allowed, conditional;
  for (const UnitRule &u : q.units) {
    const std::vector<std::string> tokens = u.dimension.empty() ? std::vector<std::string>{u.unit} : UnitsOfDimension(u.dimension);
    const int holds = Holds(u.when, w);
    if (holds != 0) allowed.insert(allowed.end(), tokens.begin(), tokens.end());
    if (holds == 0) conditional.push_back(Joined(tokens) + " when " + Describe(u.when));
  }
  for (const ValidatedField &f : w.vr.fields) {
    const FieldDef *def = w.def->Find(f.key);
    if (def->spec.dimension != "@quantity" || !f.value.unit) continue;
    if (std::find(allowed.begin(), allowed.end(), f.value.unit->token) != allowed.end()) continue;
    const io::Field *src = Source(w, f.key);
    const int col = f.value.type == ValueType::Unit ? src->tokens[0].col : src->tokens.size() > 1 ? src->tokens[1].col : src->col;
    Diagnostic d = Base(w, f.key, f.line, col);
    d.hasFound = true;
    d.found = f.value.unit->token;
    d.allowed = allowed;
    Emit(d, codes::kUnitNotAllowed, Where(w, f.key) + ": unit '" + d.found + "' is not a unit of quantity '" + q.name +
                                        "' (" + Joined(allowed) + (conditional.empty() ? "" : "; " + Joined(conditional, "; ")) + ")");
  }
  if (const Value *v = w.vr.Get("reaction")) CheckReaction(w, *v, "reaction");
  if (q.reaction.Has(ReactionRule::kNeutralPair))
    if (const Value *p = w.vr.Get("pair"))
      if (p->pair.light.Z > 0 && p->pair.heavy.Z > 0) {
        Diagnostic d = Base(w, "pair", p->line, p->col);
        d.hasFound = true;
        d.found = Source(w, "pair")->tokens[0].text;
        std::string message = "needs a neutral pair";
        for (const Refusal &r : q.refusals)
          if (r.rule == "charged_pair") message = r.message;
        Emit(d, codes::kReactionNotAllowed, Where(w, "pair") + ": " + d.found + " is a charged pair; " + message);
      }
  if (w.vr.kind == "dataset") {
    CheckColumns(w);
    CheckNorm(w);
    if (!q.datasetWarning.empty()) {  // a quantity computed for tests and comparison, not measured (W_NOT_MEASURABLE)
      const Value *qv = w.vr.Get("quantity");
      Diagnostic d = Base(w, "quantity", qv->line, qv->col);
      d.severity = Severity::Warning;
      d.hasFound = true;
      d.found = q.name;
      Emit(d, codes::kNotMeasurable, Where(w, "quantity") + ": " + q.datasetWarning);
    }
  }
  if (w.vr.kind == "extrapolation") {
    const Value *g = w.vr.Get("angles");
    const bool grid = g && g->grid.kind != "none";
    if (q.angle == Quantity::RequiredAngle && !grid && !(w.given.count("angles") && !g)) {
      Diagnostic d = Base(w, "angles", g ? g->line : w.vr.line, g ? g->col : 1);
      d.allowed = {"range|list|file|adaptive"};
      if (!g) {
        d.hasFix = true;
        d.fix = {"insert_field", 0, 0, w.vr.lastLine, "    angles         ?"};
      }
      Emit(d, codes::kFieldRequired, Where(w) + ": quantity '" + q.name + "' needs an angle grid in 'angles'");
    }
    if (q.angle == Quantity::NoAngle && grid) {
      Diagnostic d = Base(w, "angles", g->line, g->col);
      d.hasFound = true;
      d.found = g->grid.kind;
      d.allowed = {"none"};
      Emit(d, codes::kFieldNotApplicable, Where(w, "angles") + ": quantity '" + q.name + "' is angle-integrated; write "
                                                                "'angles none' or leave the field out");
    }
  }
  // The target's kind.
  const Work *t = nullptr;
  if (Referenced(w, "target", t) == 1 && t->vr.kind == "target" && !q.targetKinds.empty()) {
    const Value *kind = t->vr.Get("kind");
    if (kind && std::find(q.targetKinds.begin(), q.targetKinds.end(), kind->word) == q.targetKinds.end()) {
      const io::Field *src = Source(w, "target");
      Diagnostic d = Base(w, "target", src->line, src->tokens[0].col);
      d.hasFound = true;
      d.found = kind->word;
      d.allowed = q.targetKinds;
      std::string message = "quantity '" + q.name + "' takes a target of kind " + Joined(q.targetKinds, " or ");
      for (const Refusal &r : q.refusals)
        if (r.code == codes::kTargetKind && r.rule.find(kind->word) != std::string::npos) message += ": " + r.message;
      Emit(d, codes::kTargetKind, Where(w, "target") + ": target '" + t->vr.name + "' is of kind " + kind->word + "; " + message);
    }
  }
  // F9: efficiency and norm multiply; both free without a prior is degenerate.
  const Value *eff = w.vr.Get("efficiency"), *norm = w.vr.Get("norm");
  if (eff && norm && eff->type == ValueType::Nuisance && eff->nuisance.FreeWithoutPrior()) {
    const Nuisance *n = norm->type == ValueType::Nuisance ? &norm->nuisance : nullptr;
    if (!n)
      if (const Work *r = Resolve(norm->word))
        if (const Value *rv = r->vr.Get("value")) n = &rv->nuisance;
    if (n && n->FreeWithoutPrior()) {
      const io::Field *src = Source(w, "efficiency");
      Diagnostic d = Base(w, "efficiency", src->line, src->col);
      Emit(d, codes::kDegenerateNuisance, Where(w) + ": 'efficiency' and 'norm' multiply and are both free without a "
                                                     "prior, so only their product is determined; give one of them a prior "
                                                     "or fix it (F9)");
    }
  }
}

void Validator::CheckRecord(Work &w) {
  if (!w.def) return;
  CheckFields(w);
  if (w.vr.quantity) CheckQuantity(w);
  CheckReferences(w);
}

}  // namespace

ValidationResult Schema::Validate(const std::vector<SectionRecords> &sections) {
  ValidationResult result;
  Validator v(result.diagnostics);
  for (const SectionRecords &s : sections)
    for (const io::Record &r : s.records) v.Parse(s.section, r);
  v.CheckNames();
  for (Work &w : v.works) v.CheckRecord(w);
  std::stable_sort(result.diagnostics.begin(), result.diagnostics.end(),
                   [](const Diagnostic &a, const Diagnostic &b) { return a.line < b.line; });
  for (Work &w : v.works)
    if (w.def) result.records.push_back(std::move(w.vr));
  result.references = std::move(v.references);
  return result;
}

}  // namespace experiment
}  // namespace azure2
