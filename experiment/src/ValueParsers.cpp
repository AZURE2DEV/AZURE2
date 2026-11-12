// The parsers of the scalar and structured value types of format 3 (Parsers.h; docs/upgrade/project_format3.md
// section 3.3).  Reactions are in ReactionParser.cpp, column roles and grids in ColumnsGrid.cpp.
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>

#include "azure2/experiment/Parsers.h"
#include "azure2/physics/Elements.h"

namespace azure2 {
namespace experiment {

namespace {

std::string Joined(const std::vector<std::string> &v, const char *sep) {
  std::string s;
  for (std::size_t i = 0; i < v.size(); i++) s += (i ? sep : "") + v[i];
  return s;
}

// Light particles of reaction and pair labels (section 3.5): letter, Z, A.
struct Light {
  const char *letter;
  int Z, A;
};
const Light kLight[] = {{"n", 0, 1}, {"p", 1, 1}, {"d", 1, 2}, {"t", 1, 3}, {"h", 2, 3},
                        {"a", 2, 4}, {"g", 0, 0}, {"e-", -1, 0}, {"b-", -1, 0}, {"b+", 1, 0}};

bool ParseUnitToken(const std::string &dimension, TokenCursor &cur, ParseContext &ctx, const Unit *&out,
                    const std::string &what) {
  if (cur.AtEnd()) {
    ctx.Error(codes::kUnitMissing, cur.endCol, "", UnitsOfDimension(dimension),
              ctx.Where() + ": " + what + " needs a unit" +
                  (dimension.empty() || dimension[0] == '@' ? std::string() : " of " + dimension) +
                  "; there is no default unit");
    return false;
  }
  const io::Token &t = cur.Peek();
  const Unit *u = FindUnit(t.text);
  if (!u) {
    const Unit *fixed = FindUnitIgnoringCase(t.text);
    const bool known = fixed && (dimension.empty() || dimension[0] == '@' || dimension == fixed->dimension);
    Fix fix{"replace_token", cur.line, t.col, 0, known ? fixed->token : ""};
    std::vector<std::string> allowed =
        (dimension.empty() || dimension[0] == '@') ? std::vector<std::string>{} : UnitsOfDimension(dimension);
    // A dimensioned value followed by a mode word or nothing unit-like: the unit is missing, not unknown.
    if (!known && (t.text == "fixed" || t.text == "fit" || t.text == "profiled" || t.text == "lab" || t.text == "cm")) {
      ctx.Error(codes::kUnitMissing, t.col, t.text, allowed,
                ctx.Where() + ": " + what + " needs a unit" + (allowed.empty() ? "" : " (" + Joined(allowed, ", ") + ")") +
                    " before '" + t.text + "'; there is no default unit");
      return false;
    }
    ctx.Error(codes::kUnknownUnit, t.col, t.text, allowed,
              ctx.Where() + ": '" + t.text + "' is not a unit" +
                  (allowed.empty() ? "" : "; " + what + " takes " + Joined(allowed, ", ")) +
                  (known ? " (did you mean '" + std::string(fixed->token) + "'?)" : ""),
              known ? &fix : nullptr);
    return false;
  }
  if (!dimension.empty() && dimension[0] != '@' && dimension != u->dimension) {
    ctx.Error(codes::kUnitNotAllowed, t.col, t.text, UnitsOfDimension(dimension),
              ctx.Where() + ": unit '" + t.text + "' is a " + u->dimension + " unit; " + what + " is a " + dimension +
                  " (" + Joined(UnitsOfDimension(dimension), ", ") + ")");
    return false;
  }
  cur.Next();
  out = u;
  return true;
}

bool NeedToken(TokenCursor &cur, ParseContext &ctx, const std::string &what) {
  if (!cur.AtEnd()) return true;
  ctx.Error(codes::kTooFewTokens, cur.endCol, "", {}, ctx.Where() + ": " + what + " is missing");
  return false;
}

bool ParseNumberToken(TokenCursor &cur, ParseContext &ctx, Number &out, const std::string &what, bool positive) {
  if (!NeedToken(cur, ctx, what)) return false;
  const io::Token &t = cur.Peek();
  if (!ParseNumberText(t.text, out)) {
    ctx.Error(codes::kNotANumber, t.col, t.text, {}, ctx.Where() + ": " + what + " must be a number, not '" + t.text + "'");
    return false;
  }
  if (positive && !(out.value > 0)) {
    ctx.Error(codes::kOutOfRange, t.col, t.text, {"> 0"}, ctx.Where() + ": " + what + " must be positive, not " + t.text);
    return false;
  }
  out.col = t.col;
  cur.Next();
  return true;
}

bool ParseNuisance(const ValueSpec &spec, TokenCursor &cur, ParseContext &ctx, Value &out) {
  if (!NeedToken(cur, ctx, "the value")) return false;
  Nuisance &n = out.nuisance;
  const io::Token &first = cur.Peek();
  Number probe;
  const bool modeWord = first.text == "fixed" || first.text == "fit" || first.text == "profiled";
  if (!spec.refKinds.empty() && !modeWord && !ParseNumberText(first.text, probe)) {
    if (!IsName(first.text)) {
      ctx.Error(codes::kBadName, first.col, first.text, {},
                ctx.Where() + ": '" + first.text + "' is neither a number nor the name of a " +
                    Joined(spec.refKinds, " or ") + " record");
      return false;
    }
    out.type = ValueType::Reference;
    out.word = first.text;
    cur.Next();
    return true;
  }
  if (first.text == "profiled") {
    if (!spec.profiled) {
      ctx.Error(codes::kModeNotAllowed, first.col, first.text, {"fixed", "fit"},
                ctx.Where() + ": 'profiled' is the norm mode of THM data sets of an experiment; here the mode is "
                              "fixed or fit");
      return false;
    }
    cur.Next();
    n.mode = NuisanceMode::Profiled;
    return true;
  }
  if (!ParseNumberToken(cur, ctx, n.value, "the value", false)) return false;
  n.hasValue = true;
  if (!spec.dimension.empty() && !ParseUnitToken(spec.dimension, cur, ctx, n.unit, "the value")) return false;
  if (cur.AtEnd()) {
    ctx.Error(codes::kNuisanceMode, cur.endCol, "", {"fixed", "fit"},
              ctx.Where() + ": the mode is missing after the value: write 'fixed' or 'fit' (there is no default mode)");
    return false;
  }
  const io::Token &m = cur.Peek();
  if (m.text == "fixed") {
    n.mode = NuisanceMode::Fixed;
  } else if (m.text == "fit") {
    n.mode = NuisanceMode::Fit;
  } else {
    const bool unitHere = spec.dimension.empty() && FindUnit(m.text);
    ctx.Error(codes::kNuisanceMode, m.col, m.text, {"fixed", "fit"},
              ctx.Where() + ": '" + m.text + "' is not a mode: write 'fixed' or 'fit' after the value" +
                  (m.text == "profiled" ? " ('profiled' takes no value)" : "") +
                  (unitHere ? " (this value is dimensionless and takes no unit)" : ""));
    return false;
  }
  cur.Next();
  while (!cur.AtEnd()) {
    const io::Token &k = cur.Peek();
    if (k.text == "prior" && !n.hasPrior) {
      cur.Next();
      if (cur.pos + 2 > cur.tokens->size()) {
        ctx.Error(codes::kPriorSyntax, k.col, k.text, {"prior <centre> <width>[%]"},
                  ctx.Where() + ": 'prior' needs its centre and width: prior <centre> <width>[%]");
        return false;
      }
      const io::Token &c = cur.Next();
      io::Token w = cur.Next();
      if (!w.text.empty() && w.text.back() == '%') {
        n.priorPercent = true;
        w.text.pop_back();
      }
      if (!ParseNumberText(c.text, n.priorCentre) || !ParseNumberText(w.text, n.priorWidth) || !(n.priorWidth.value > 0)) {
        ctx.Error(codes::kPriorSyntax, c.col, c.text + " " + w.text + (n.priorPercent ? "%" : ""),
                  {"prior <centre> <width>[%]"},
                  ctx.Where() + ": the prior is 'prior <centre> <width>[%]' with numbers and a width > 0");
        return false;
      }
      n.priorCentre.col = c.col;
      n.priorWidth.col = w.col;
      n.hasPrior = true;
    } else if (k.text == "grid_max" && !n.hasGridMax) {
      cur.Next();
      if (!ParseNumberToken(cur, ctx, n.gridMax, "grid_max", false)) return false;
      n.hasGridMax = true;
    } else {
      break;
    }
  }
  return true;
}

bool ParseForm(const ValueSpec &spec, TokenCursor &cur, ParseContext &ctx, Value &out) {
  const FormSpec &form = *spec.form;
  const Alternative *alt = nullptr;
  std::vector<std::string> keywords;
  for (const Alternative &a : form.alternatives) {
    if (a.keyword.empty()) continue;
    keywords.push_back(a.keyword);
    if (!alt && !cur.AtEnd() && !cur.Peek().quoted && cur.Peek().text == a.keyword) alt = &a;
  }
  if (alt) {
    cur.Next();
  } else {
    for (const Alternative &a : form.alternatives)
      if (a.keyword.empty()) alt = &a;
  }
  if (!alt) {
    if (cur.AtEnd()) {
      ctx.Error(codes::kTooFewTokens, cur.endCol, "", keywords, ctx.Where() + ": the value is missing: " + Grammar(spec));
    } else {
      const io::Token &t = cur.Peek();
      ctx.Error(codes::kNotInList, t.col, t.text, keywords,
                ctx.Where() + ": '" + t.text + "' is not one of " + Joined(keywords, ", "));
    }
    return false;
  }
  out.word = alt->keyword;
  const std::string head = alt->keyword.empty() ? std::string() : "'" + alt->keyword + "' ";
  std::vector<std::string> partKeywords;
  for (const PartSpec &p : alt->parts) {
    if (p.keyword) {
      partKeywords.push_back(p.name);
      continue;
    }
    if (cur.AtEnd()) {
      if (p.optional) continue;
      ctx.Error(codes::kTooFewTokens, cur.endCol, "", {}, ctx.Where() + ": " + head + "needs <" + p.name + "> (" + Grammar(spec) + ")");
      return false;
    }
    Value v;
    v.line = cur.line;
    v.col = cur.Peek().col;
    if (!ParseValue(p.spec, cur, ctx, v)) return false;
    out.parts.push_back({p.name, v});
  }
  while (!cur.AtEnd() && !partKeywords.empty()) {
    const io::Token &t = cur.Peek();
    const PartSpec *part = nullptr;
    for (const PartSpec &p : alt->parts)
      if (p.keyword && p.name == t.text && !t.quoted) part = &p;
    if (!part) {
      ctx.Error(codes::kNotInList, t.col, t.text, partKeywords,
                ctx.Where() + ": '" + t.text + "' is not a keyword of " + head + "(" + Joined(partKeywords, ", ") + ")");
      return false;
    }
    if (out.Get(part->name)) {
      ctx.Error(codes::kPartTwice, t.col, t.text, {}, ctx.Where() + ": '" + t.text + "' is given twice");
      return false;
    }
    cur.Next();
    Value v;
    v.line = cur.line;
    v.col = t.col;
    if (cur.AtEnd()) {
      ctx.Error(codes::kTooFewTokens, cur.endCol, "", {}, ctx.Where() + ": '" + part->name + "' needs a value (" + Grammar(part->spec) + ")");
      return false;
    }
    if (!ParseValue(part->spec, cur, ctx, v)) return false;
    out.parts.push_back({part->name, v});
  }
  for (const PartSpec &p : alt->parts) {
    if (!p.keyword || p.optional || out.Get(p.name)) continue;
    ctx.Error(codes::kPartMissing, cur.endCol, "", {p.name},
              ctx.Where() + ": " + head + "needs '" + p.name + " " + Grammar(p.spec) + "'");
    return false;
  }
  return true;
}

}  // namespace

void ParseContext::Error(const char *code, int col, const std::string &found, std::vector<std::string> allowed,
                         const std::string &message, const Fix *fix) {
  if (failed) return;
  failed = true;
  Diagnostic d = base;
  d.code = code;
  d.column = col;
  d.hasFound = !found.empty();
  d.found = found;
  d.allowed = std::move(allowed);
  d.message = message;
  if (fix && !fix->text.empty()) {
    d.hasFix = true;
    d.fix = *fix;
  }
  out->push_back(d);
}

std::string ParseContext::Where() const {
  std::string s = base.recordKind + " '" + base.record + "'";
  if (!base.field.empty()) s += ", field '" + base.field + "'";
  return s;
}

bool ParseNumberText(const std::string &text, Number &out) {
  if (text.empty()) return false;
  const char c = text[0];
  if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '+' || c == '-' || c == '.')) return false;
  errno = 0;
  char *end = nullptr;
  const double v = std::strtod(text.c_str(), &end);
  if (end != text.c_str() + text.size() || !std::isfinite(v) || errno == ERANGE) return false;
  out.text = text;
  out.value = v;
  return true;
}

bool IsName(const std::string &text) {
  if (text.empty() || !std::isalpha(static_cast<unsigned char>(text[0]))) return false;
  for (char c : text)
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-') return false;
  return true;
}

bool ParseSpinText(const std::string &text, bool withParity, Spin &out) {
  std::string s = text;
  out = Spin();
  out.text = text;
  if (withParity) {
    if (s.empty() || (s.back() != '+' && s.back() != '-')) return false;
    out.parity = s.back() == '+' ? 1 : -1;
    s.pop_back();
  }
  const std::size_t slash = s.find('/');
  auto digits = [](const std::string &d) {
    if (d.empty() || d.size() > 4) return false;
    for (char c : d)
      if (!std::isdigit(static_cast<unsigned char>(c))) return false;
    return true;
  };
  if (slash == std::string::npos) {
    if (!digits(s)) return false;
    out.twice = 2 * std::atoi(s.c_str());
    return true;
  }
  const std::string num = s.substr(0, slash), den = s.substr(slash + 1);
  if (!digits(num) || den != "2" || std::atoi(num.c_str()) % 2 != 1) return false;
  out.twice = std::atoi(num.c_str());
  return true;
}

bool ParseParticleText(const std::string &text, bool light, Particle &out, std::string &error, const char *&code) {
  out = Particle();
  out.text = text;
  code = codes::kNucleusSyntax;
  if (light) {
    for (const Light &l : kLight) {
      const std::string letter = l.letter;
      const bool alias = text.size() == letter.size() + 1 && text.compare(0, letter.size(), letter) == 0 &&
                         text.back() == '0' && letter.size() == 1 && letter != "e";
      if (text == letter || alias) {
        out.light = letter;
        out.Z = l.Z;
        out.A = l.A;
        return true;
      }
    }
    if (text.size() >= 2 && std::string("npdtha").find(text[0]) != std::string::npos) {
      bool ordinal = true;
      for (std::size_t i = 1; i < text.size(); i++) ordinal = ordinal && std::isdigit(static_cast<unsigned char>(text[i]));
      if (ordinal) {
        code = codes::kReactionOrdinal;
        error = "'" + text + "' is an ordinal label, which counts the states the project happens to include; name the "
                "residual state by its excitation energy in MeV instead ('" + std::string(1, text[0]) + "' and " +
                "[<E_x>] after the residual nucleus)";
        return false;
      }
    }
  }
  std::size_t i = 0;
  while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) i++;
  const std::string symbol = text.substr(i);
  if (i == 0 || i > 3 || symbol.empty()) {
    error = "'" + text + "' is not a nuclide <A><Symbol> (e.g. 13C)" +
            std::string(light ? " or a light particle (n p d t h a g e- b- b+)" : "");
    return false;
  }
  const int z = physics::ElementZ(symbol);
  if (z < 0) {
    code = codes::kUnknownElement;
    error = "'" + symbol + "' in '" + text + "' is not an element symbol (He, Li, C, ...; case matters)";
    return false;
  }
  out.Z = z;
  out.A = std::atoi(text.substr(0, i).c_str());
  if (out.A < 1 || out.A < out.Z) {
    error = "'" + text + "': the mass number A = " + std::to_string(out.A) + " is smaller than Z = " + std::to_string(z);
    return false;
  }
  return true;
}

bool ParseValue(const ValueSpec &spec, TokenCursor &cur, ParseContext &ctx, Value &out) {
  out.type = spec.type;
  if (!cur.AtEnd()) {
    out.line = cur.line;
    out.col = cur.Peek().col;
  }
  switch (spec.type) {
    case ValueType::Number:
      return ParseNumberToken(cur, ctx, out.number, "the value", spec.positive);
    case ValueType::Integer: {
      if (!NeedToken(cur, ctx, "the value")) return false;
      const io::Token &t = cur.Peek();
      std::size_t k = (t.text[0] == '-' || t.text[0] == '+') ? 1 : 0;
      bool ok = t.text.size() > k && t.text.size() < 10;
      for (std::size_t i = k; ok && i < t.text.size(); i++) ok = std::isdigit(static_cast<unsigned char>(t.text[i]));
      if (!ok) {
        ctx.Error(codes::kNotAnInteger, t.col, t.text, {}, ctx.Where() + ": '" + t.text + "' is not an integer");
        return false;
      }
      out.integer = std::atol(t.text.c_str());
      out.number = {t.text, static_cast<double>(out.integer), t.col};
      if (out.integer < spec.min) {
        ctx.Error(codes::kOutOfRange, t.col, t.text, {">= " + std::to_string(spec.min)},
                  ctx.Where() + ": " + t.text + " is out of range; the value is an integer >= " + std::to_string(spec.min));
        return false;
      }
      cur.Next();
      return true;
    }
    case ValueType::Unit:
      return ParseUnitToken(spec.dimension, cur, ctx, out.unit, "the value");
    case ValueType::NumberWithUnit:
      return ParseNumberToken(cur, ctx, out.number, "the value", spec.positive) &&
             ParseUnitToken(spec.dimension, cur, ctx, out.unit, "the value");
    case ValueType::Range: {
      if (!ParseNumberToken(cur, ctx, out.lo, "<lo>", false) || !ParseNumberToken(cur, ctx, out.hi, "<hi>", false))
        return false;
      if (out.lo.value > out.hi.value) {
        ctx.Error(codes::kEmptyRange, out.lo.col, out.lo.text + " " + out.hi.text, {"<lo> <hi> with lo <= hi"},
                  ctx.Where() + ": the range " + out.lo.text + " " + out.hi.text + " is empty (lo > hi)");
        return false;
      }
      return ParseUnitToken(spec.dimension, cur, ctx, out.unit, "the range");
    }
    case ValueType::Nuisance:
      return ParseNuisance(spec, cur, ctx, out);
    case ValueType::Enum: {
      if (!NeedToken(cur, ctx, "the value")) return false;
      const io::Token &t = cur.Peek();
      for (const std::string &v : spec.values)
        if (t.text == v) out.word = v;
      for (const auto &a : spec.aliases)
        if (t.text == a.first) out.word = a.second;
      if (out.word.empty() || t.quoted) {
        ctx.Error(codes::kNotInList, t.col, t.text, spec.values,
                  ctx.Where() + ": '" + t.text + "' is not allowed; the value is one of " + Joined(spec.values, ", "));
        return false;
      }
      cur.Next();
      return true;
    }
    case ValueType::Reference:
    case ValueType::Name: {
      if (!NeedToken(cur, ctx, "the name")) return false;
      const io::Token &t = cur.Peek();
      std::string name = t.text;
      const std::size_t at = spec.epoch ? name.find('@') : std::string::npos;
      if (at != std::string::npos) {
        out.epoch = name.substr(at + 1);
        name = name.substr(0, at);
      }
      if (!IsName(name) || (at != std::string::npos && !IsName(out.epoch))) {
        ctx.Error(codes::kBadName, t.col, t.text, {"[A-Za-z][A-Za-z0-9_-]*" + std::string(spec.epoch ? "[@<epoch>]" : "")},
                  ctx.Where() + ": '" + t.text + "' is not a name (a letter, then letters, digits, '_' or '-')");
        return false;
      }
      out.word = name;
      cur.Next();
      return true;
    }
    case ValueType::Path:
      if (!NeedToken(cur, ctx, "the path")) return false;
      out.word = cur.Next().text;
      return true;
    case ValueType::Text:
      if (!NeedToken(cur, ctx, "the text")) return false;
      while (!cur.AtEnd()) out.word += (out.word.empty() ? "" : " ") + cur.Next().text;
      return true;
    case ValueType::Formula: {
      if (!NeedToken(cur, ctx, "the formula")) return false;
      const io::Token &t = cur.Peek();
      if (!t.quoted) {
        ctx.Error(codes::kFormulaSyntax, t.col, t.text, {"\"<expression in x>\" [<name>=<number> ...]"},
                  ctx.Where() + ": the formula is written in quotes: \"" + t.text + "\"");
        return false;
      }
      out.formula.expression = t.text;
      cur.Next();
      while (!cur.AtEnd() && cur.Peek().text.find('=') != std::string::npos) {
        const io::Token &p = cur.Peek();
        const std::size_t eq = p.text.find('=');
        Formula::Param param;
        param.name = p.text.substr(0, eq);
        if (!IsName(param.name) || !ParseNumberText(p.text.substr(eq + 1), param.value)) {
          ctx.Error(codes::kFormulaSyntax, p.col, p.text, {"<name>=<number>"},
                    ctx.Where() + ": a formula parameter is <name>=<number>, not '" + p.text + "'");
          return false;
        }
        param.value.col = p.col;
        out.formula.params.push_back(param);
        cur.Next();
      }
      return true;
    }
    case ValueType::NumberList: {
      Number n;
      while (!cur.AtEnd() && ParseNumberText(cur.Peek().text, n)) {
        n.col = cur.Next().col;
        out.numbers.push_back(n);
      }
      if (out.numbers.empty()) {
        if (!NeedToken(cur, ctx, "the numbers")) return false;
        ctx.Error(codes::kNotANumber, cur.Peek().col, cur.Peek().text, {},
                  ctx.Where() + ": '" + cur.Peek().text + "' is not a number");
        return false;
      }
      return true;
    }
    case ValueType::Spin:
    case ValueType::SpinParity: {
      if (!NeedToken(cur, ctx, "the spin")) return false;
      const io::Token &t = cur.Peek();
      const bool parity = spec.type == ValueType::SpinParity;
      if (!ParseSpinText(t.text, parity, out.spin)) {
        ctx.Error(codes::kSpinSyntax, t.col, t.text, {parity ? "<J><+|->: 0+, 3-, 1/2+" : "<J>: 0, 1, 3/2"},
                  ctx.Where() + ": '" + t.text + "' is not a spin" + (parity ? " with parity (0+, 3-, 1/2+)" : " (0, 1, 3/2)") +
                      "; half-integers are written n/2");
        return false;
      }
      cur.Next();
      return true;
    }
    case ValueType::Nucleus: {
      if (!NeedToken(cur, ctx, "the nucleus")) return false;
      const io::Token &t = cur.Peek();
      std::string error;
      const char *code = nullptr;
      if (!ParseParticleText(t.text, false, out.nucleus, error, code)) {
        ctx.Error(code, t.col, t.text, {"<A><Symbol>"}, ctx.Where() + ": " + error);
        return false;
      }
      cur.Next();
      return true;
    }
    case ValueType::Pair: {
      if (!NeedToken(cur, ctx, "the pair")) return false;
      const io::Token &t = cur.Peek();
      std::size_t plus = std::string::npos;
      for (std::size_t p = 1; p + 1 < t.text.size() && plus == std::string::npos; p++)
        if (t.text[p] == '+' && std::isdigit(static_cast<unsigned char>(t.text[p + 1]))) plus = p;
      std::string heavy = plus == std::string::npos ? "" : t.text.substr(plus + 1), error;
      std::string state;
      const std::size_t br = heavy.find('[');
      if (br != std::string::npos && heavy.back() == ']') {
        state = heavy.substr(br + 1, heavy.size() - br - 2);
        heavy = heavy.substr(0, br);
      }
      // the light particle's excitation, as in a reaction label (a[3.0], section 3.5)
      std::string light = plus == std::string::npos ? "" : t.text.substr(0, plus), lightExc;
      const std::size_t lbr = light.find('[');
      if (lbr != std::string::npos && light.back() == ']') {
        lightExc = light.substr(lbr + 1, light.size() - lbr - 2);
        light = light.substr(0, lbr);
      }
      const char *code = nullptr;
      Number e, le;
      const bool ok = plus != std::string::npos &&
                      ParseParticleText(light, true, out.pair.light, error, code) &&
                      (lightExc.empty() || (ParseNumberText(lightExc, le) && le.value >= 0)) &&
                      ParseParticleText(heavy, false, out.pair.heavy, error, code) &&
                      (state.empty() || state == "0" || (ParseNumberText(state, e) && e.value > 0));
      if (!ok) {
        ctx.Error(codes::kPairSyntax, t.col, t.text, {"<particle>+<nucleus>[<E_x in MeV>]"},
                  ctx.Where() + ": '" + t.text + "' is not a particle pair <particle>+<nucleus>[<E_x>] (p+15N, a+12C[4.439])" +
                      (error.empty() ? "" : ": " + error));
        return false;
      }
      if (!lightExc.empty() && le.value > 0) {
        out.pair.light.excited = true;
        out.pair.light.excitation = le;
      }
      if (!state.empty() && state != "0") {
        out.pair.state.kind = ResidualState::Excited;
        out.pair.state.energy = e;
      }
      cur.Next();
      return true;
    }
    case ValueType::Composition: {
      if (!NeedToken(cur, ctx, "the composition")) return false;
      while (!cur.AtEnd()) {
        const io::Token &s = cur.Peek();
        const int z = physics::ElementZ(s.text);
        if (z < 1) {
          ctx.Error(codes::kUnknownElement, s.col, s.text, {}, ctx.Where() + ": '" + s.text +
                    "' is not an element symbol; the composition is <Symbol> <atoms> ... (C 1 O 2)");
          return false;
        }
        const std::string symbol = s.text;
        cur.Next();
        Number n;
        if (!ParseNumberToken(cur, ctx, n, "the number of " + symbol + " atoms", true)) return false;
        out.composition.push_back({symbol, n});
      }
      return true;
    }
    case ValueType::Reaction:
      return ParseReaction(cur, ctx, out.reaction);
    case ValueType::Columns:
      return ParseColumns(cur, ctx, out.columns);
    case ValueType::Grid:
      return ParseGrid(spec, cur, ctx, out.grid);
    case ValueType::Form:
      return ParseForm(spec, cur, ctx, out);
  }
  return false;
}

}  // namespace experiment
}  // namespace azure2
