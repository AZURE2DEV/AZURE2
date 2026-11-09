// Column roles of a data file (docs/upgrade/project_format3.md section 3.4.1) and extrapolation grids (section 3.8).
// Which roles a quantity requires or refuses is the validator's (Quantity::energyRoles, ...); here only the grammar.
#include "azure2/experiment/Parsers.h"

namespace azure2 {
namespace experiment {

const std::vector<std::pair<std::string, std::string>> &AxisRoles() {
  static const std::vector<std::pair<std::string, std::string>> roles = {
      {"E_lab", "energy"}, {"E_cm", "energy"}, {"E_x", "energy"}, {"E_decay", "energy"},
      {"theta_lab", "angle"}, {"theta_cm", "angle"}};
  return roles;
}

const std::vector<std::string> &PlainRoles() {
  static const std::vector<std::string> roles = {"order", "value", "error", "error_percent", "window_lo", "window_hi", "_"};
  return roles;
}

bool ParseColumns(TokenCursor &cur, ParseContext &ctx, Columns &out) {
  if (cur.AtEnd()) {
    ctx.Error(codes::kTooFewTokens, cur.endCol, "", {}, ctx.Where() + ": the column roles are missing");
    return false;
  }
  if (cur.Peek().text == "from_file" && cur.tokens->size() - cur.pos == 1) {
    cur.Next();
    out.fromFile = true;
    return true;
  }
  std::vector<std::string> all;
  for (const auto &r : AxisRoles()) all.push_back(r.first + "[<unit>]");
  for (const std::string &r : PlainRoles()) all.push_back(r);
  while (!cur.AtEnd()) {
    const io::Token &t = cur.Peek();
    ColumnRole role;
    role.col = t.col;
    const std::size_t br = t.text.find('[');
    role.role = t.text.substr(0, br);
    std::string unit;
    const bool bracket = br != std::string::npos;
    if (bracket) {
      if (t.text.back() != ']') {
        ctx.Error(codes::kRoleUnit, t.col, t.text, {}, ctx.Where() + ": '" + t.text + "': the unit is written in brackets, as E_lab[MeV]");
        return false;
      }
      unit = t.text.substr(br + 1, t.text.size() - br - 2);
    }
    std::string dimension;
    bool known = false;
    for (const auto &r : AxisRoles())
      if (r.first == role.role) {
        known = true;
        dimension = r.second;
      }
    for (const std::string &r : PlainRoles()) known = known || r == role.role;
    if (!known) {
      ctx.Error(codes::kUnknownRole, t.col, t.text, all,
                ctx.Where() + ": '" + t.text + "' is not a column role" +
                    (t.text == "from_file" ? " ('from_file' stands alone: the roles are then in the data file)" : ""));
      return false;
    }
    if (!dimension.empty()) {
      role.unit = FindUnit(unit);
      if (!bracket || !role.unit || dimension != role.unit->dimension) {
        ctx.Error(codes::kRoleUnit, t.col, t.text, UnitsOfDimension(dimension),
                  ctx.Where() + ": the column " + role.role + " needs its " + dimension + " unit in brackets (" + role.role +
                      "[" + UnitsOfDimension(dimension)[0] + "])" + (bracket ? ", not '" + unit + "'" : "") +
                      "; there is no default unit");
        return false;
      }
    } else if (bracket) {
      ctx.Error(codes::kRoleUnit, t.col, t.text, {role.role},
                ctx.Where() + ": the column " + role.role + " takes no unit (its unit is the record's)");
      return false;
    }
    if (role.role != "_" && out.Find(role.role)) {
      ctx.Error(codes::kRoleTwice, t.col, t.text, {}, ctx.Where() + ": the column role " + role.role + " is given twice");
      return false;
    }
    out.roles.push_back(role);
    cur.Next();
  }
  return true;
}

bool ParseGrid(const ValueSpec &spec, TokenCursor &cur, ParseContext &ctx, Grid &out) {
  std::vector<std::string> kinds = {"range", "list", "file", "adaptive"};
  if (spec.none) kinds.push_back("none");
  if (cur.AtEnd()) {
    ctx.Error(codes::kTooFewTokens, cur.endCol, "", kinds, ctx.Where() + ": the grid is missing: " + Grammar(spec));
    return false;
  }
  const io::Token &k = cur.Peek();
  bool known = false;
  for (const std::string &x : kinds) known = known || (k.text == x && !k.quoted);
  if (!known) {
    ctx.Error(codes::kNotInList, k.col, k.text, kinds, ctx.Where() + ": '" + k.text + "' is not a grid kind (" + Grammar(spec) + ")");
    return false;
  }
  out.kind = cur.Next().text;
  if (out.kind == "none") return true;
  ValueSpec number(ValueType::Number), unit(ValueType::Unit, spec.dimension), frame(ValueType::Enum), integer(ValueType::Integer);
  frame.values = {"lab", "cm"};
  Value v;
  auto num = [&](Number &n, bool positive) {
    number.positive = positive;
    if (!ParseValue(number, cur, ctx, v)) return false;
    n = v.number;
    return true;
  };
  auto unitFrame = [&]() {
    if (!ParseValue(unit, cur, ctx, v)) return false;
    out.unit = v.unit;
    if (!ParseValue(frame, cur, ctx, v)) return false;
    out.frame = v.word;
    return true;
  };
  if (out.kind == "list") {
    Number n;
    while (!cur.AtEnd() && ParseNumberText(cur.Peek().text, n)) {
      n.col = cur.Next().col;
      out.values.push_back(n);
    }
    if (out.values.empty()) return num(n, false);  // reports the missing or non-numeric token
    return unitFrame();
  }
  if (out.kind == "file") {
    ValueSpec path(ValueType::Path);
    if (!ParseValue(path, cur, ctx, v)) return false;
    out.path = v.word;
    return unitFrame();
  }
  if (!num(out.lo, false) || !num(out.hi, false)) return false;
  if (!(out.lo.value < out.hi.value)) {
    ctx.Error(codes::kEmptyRange, out.lo.col, out.lo.text + " " + out.hi.text, {"<lo> < <hi>"},
              ctx.Where() + ": the grid " + out.lo.text + " to " + out.hi.text + " is empty (lo >= hi)");
    return false;
  }
  if (!unitFrame()) return false;
  if (out.kind == "range") {
    if (cur.AtEnd() || (cur.Peek().text != "step" && cur.Peek().text != "points")) {
      ctx.Error(codes::kPartMissing, cur.AtEnd() ? cur.endCol : cur.Peek().col, cur.AtEnd() ? "" : cur.Peek().text,
                {"step <s>", "points <n> [log]"}, ctx.Where() + ": a range grid needs 'step <s>' or 'points <n> [log]'");
      return false;
    }
    if (cur.Next().text == "step") return (out.hasStep = num(out.step, true));
    integer.min = 2;
    if (!ParseValue(integer, cur, ctx, v)) return false;
    out.hasPoints = true;
    out.points = v.integer;
    if (!cur.AtEnd() && cur.Peek().text == "log") {
      cur.Next();
      out.log = true;
    }
    return true;
  }
  // adaptive: tolerance, max_points, start_step, in any order, each once.
  bool tol = false, maxp = false, start = false;
  while (!cur.AtEnd()) {
    const io::Token &t = cur.Peek();
    bool *seen = t.text == "tolerance" ? &tol : t.text == "max_points" ? &maxp : t.text == "start_step" ? &start : nullptr;
    if (!seen) {
      ctx.Error(codes::kNotInList, t.col, t.text, {"tolerance", "max_points", "start_step"},
                ctx.Where() + ": '" + t.text + "' is not a keyword of an adaptive grid (tolerance, max_points, start_step)");
      return false;
    }
    if (*seen) {
      ctx.Error(codes::kPartTwice, t.col, t.text, {}, ctx.Where() + ": '" + t.text + "' is given twice");
      return false;
    }
    *seen = true;
    cur.Next();
    if (seen == &maxp) {
      integer.min = 2;
      if (!ParseValue(integer, cur, ctx, v)) return false;
      out.maxPoints = v.integer;
    } else if (!num(seen == &tol ? out.tolerance : out.startStep, true)) {
      return false;
    }
  }
  if (!tol || !maxp || !start) {
    ctx.Error(codes::kPartMissing, cur.endCol, "", {"tolerance", "max_points", "start_step"},
              ctx.Where() + ": an adaptive grid needs tolerance <t>, max_points <n> and start_step <s>");
    return false;
  }
  return true;
}

}  // namespace experiment
}  // namespace azure2
