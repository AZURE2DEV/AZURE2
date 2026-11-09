// Reaction expressions of format 3 (Parsers.h; docs/upgrade/project_format3.md section 3.5):
//   expression := term { + term } | term / term
//   term       := [<number> *] <reaction> [@ <angle> <unit> <frame>]
//   reaction   := <A><Symbol>(<particle>,<particle>)<A><Symbol>[<state>]
// with light particles n p d t h a g e- b- b+ (aliases a0 p0 n0 d0 t0 h0 g0 for the ground state; ordinals a1 refused),
// a particle excitation <particle>[<E>], and the residual state absent, [0], [<E_x in MeV>] or [total] (photons).
// The pair matching against <levels> is step 3's (it needs the project's pairs).
#include <cstddef>

#include "azure2/experiment/Parsers.h"

namespace azure2 {
namespace experiment {

namespace {

// <particle>[<E>] -> particle and excitation text; false on unbalanced brackets.
bool SplitExcitation(const std::string &text, std::string &particle, std::string &excitation, bool &has) {
  const std::size_t br = text.find('[');
  has = br != std::string::npos;
  if (!has) {
    particle = text;
    return text.find(']') == std::string::npos;
  }
  if (text.back() != ']' || text.find('[', br + 1) != std::string::npos) return false;
  particle = text.substr(0, br);
  excitation = text.substr(br + 1, text.size() - br - 2);
  return true;
}

bool IsGroundAlias(const Particle &p) { return !p.light.empty() && p.text.size() == 2 && p.text.back() == '0'; }

bool ParseSingle(const std::string &text, int col, ParseContext &ctx, SingleReaction &r) {
  r.text = text;
  r.col = col;
  auto fail = [&](const char *code, const std::string &why) {
    ctx.Error(code, col, text, {"<A><Symbol>(<particle>,<particle>)<A><Symbol>[<state>]"},
              ctx.Where() + ": '" + text + "' " + why);
    return false;
  };
  const std::size_t open = text.find('('), comma = text.find(',', open == std::string::npos ? 0 : open),
                    close = text.find(')', comma == std::string::npos ? 0 : comma);
  if (text.find('(', open + 1) != std::string::npos || text.find('/') != std::string::npos ||
      (close != std::string::npos && text.find('+', close) != std::string::npos))
    return fail(codes::kReactionSyntax, "holds more than one reaction: terms are joined by ' + ' or ' / ' with blanks");
  if (open == std::string::npos || comma == std::string::npos || close == std::string::npos)
    return fail(codes::kReactionSyntax, "is not a reaction <target>(<projectile>,<ejectile>)<residual>[<state>], e.g. "
                                        "15N(p,a)12C or 12C(p,g)13N[total]");
  std::string residual, state, error;
  bool hasState = false;
  if (!SplitExcitation(text.substr(close + 1), residual, state, hasState))
    return fail(codes::kReactionSyntax, "has unbalanced brackets in its residual state");
  const char *code = nullptr;
  if (!ParseParticleText(text.substr(0, open), false, r.target, error, code) ||
      !ParseParticleText(residual, false, r.residual, error, code))
    return fail(code == codes::kNucleusSyntax ? codes::kReactionSyntax : code, "-- " + error);
  struct Light {
    std::string text;
    Particle *out;
  } inner[] = {{text.substr(open + 1, comma - open - 1), &r.projectile}, {text.substr(comma + 1, close - comma - 1), &r.ejectile}};
  for (Light &l : inner) {
    std::string name, exc;
    bool hasExc = false;
    if (!SplitExcitation(l.text, name, exc, hasExc)) return fail(codes::kReactionSyntax, "has unbalanced brackets in '" + l.text + "'");
    if (!ParseParticleText(name, true, *l.out, error, code))
      return fail(code == codes::kNucleusSyntax ? codes::kReactionSyntax : code, "-- " + error);
    if (l.out->light.empty()) {
      for (const char *letter : {"n", "p", "d", "t", "h", "a"}) {
        Particle lp;
        ParseParticleText(letter, true, lp, error, code);
        if (lp.Z == l.out->Z && lp.A == l.out->A)
          return fail(codes::kReactionSyntax, "-- inside the parentheses " + name + " is written '" + letter + "'");
      }
    }
    if (hasExc) {
      if (!ParseNumberText(exc, l.out->excitation) || l.out->excitation.value < 0)
        return fail(codes::kReactionSyntax, "-- the excitation of '" + name + "' is [<E in MeV>], not [" + exc + "]");
      l.out->excited = l.out->excitation.value > 0;
    }
  }
  if (hasState) {
    if (state == "total") {
      if (r.ejectile.light != "g")
        return fail(codes::kReactionSyntax, "-- [total] is the sum over the photon exits; the ejectile must be 'g'");
      r.state.kind = ResidualState::Total;
    } else if (!ParseNumberText(state, r.state.energy) || r.state.energy.value < 0) {
      return fail(codes::kReactionSyntax, "-- the residual state is [<E_x in MeV>], [0] or [total], not [" + state + "]");
    } else if (r.state.energy.value > 0) {
      r.state.kind = ResidualState::Excited;
    }
  }
  if (IsGroundAlias(r.ejectile) && r.state.kind != ResidualState::Ground)
    return fail(codes::kReactionSyntax, "-- '" + r.ejectile.text + "' names the ground state, the brackets another state");
  bool beta = false;
  for (const Particle *p : {&r.projectile, &r.ejectile})
    beta = beta || p->light == "e-" || p->light == "b-" || p->light == "b+";
  if (!beta && (r.target.A + r.projectile.A != r.ejectile.A + r.residual.A ||
                r.target.Z + r.projectile.Z != r.ejectile.Z + r.residual.Z)) {
    ctx.Error(codes::kReactionBalance, col, text, {},
              ctx.Where() + ": '" + text + "' does not conserve charge and mass number: Z " +
                  std::to_string(r.target.Z + r.projectile.Z) + " -> " + std::to_string(r.ejectile.Z + r.residual.Z) +
                  ", A " + std::to_string(r.target.A + r.projectile.A) + " -> " + std::to_string(r.ejectile.A + r.residual.A));
    return false;
  }
  return true;
}

}  // namespace

bool ParseReaction(TokenCursor &cur, ParseContext &ctx, ReactionExpr &out) {
  char op = 0;
  while (true) {
    if (cur.AtEnd()) {
      ctx.Error(codes::kTooFewTokens, cur.endCol, "", {}, ctx.Where() + ": a reaction is missing" + std::string(op ? " after '" + std::string(1, op) + "'" : ""));
      return false;
    }
    ReactionTerm term;
    io::Token tok = cur.Next();
    std::string text = tok.text;
    int col = tok.col;
    const std::size_t star = text.find('*');
    std::string scale;
    if (star != std::string::npos) {
      scale = text.substr(0, star);
      text = text.substr(star + 1);
      col += static_cast<int>(star) + 1;
    } else if (!cur.AtEnd() && cur.Peek().text == "*") {
      scale = text;
      cur.Next();
      text.clear();
    }
    if (text.empty() && !scale.empty()) {
      if (cur.AtEnd()) {
        ctx.Error(codes::kTooFewTokens, cur.endCol, "", {}, ctx.Where() + ": a reaction is missing after '" + scale + "*'");
        return false;
      }
      text = cur.Peek().text;
      col = cur.Next().col;
    }
    if (!scale.empty()) {
      if (!ParseNumberText(scale, term.scale)) {
        ctx.Error(codes::kReactionSyntax, tok.col, tok.text, {"<number>*<reaction>"},
                  ctx.Where() + ": the scale of a term is a number, as 0.5*20Ne(n,a)17O, not '" + scale + "'");
        return false;
      }
      term.scale.col = tok.col;
      term.hasScale = true;
    }
    if (!ParseSingle(text, col, ctx, term.reaction)) return false;
    if (!cur.AtEnd() && cur.Peek().text == "@") {
      const io::Token at = cur.Next();
      if (cur.pos + 3 > cur.tokens->size()) {
        ctx.Error(codes::kTooFewTokens, cur.endCol, "", {"@ <angle> <unit> lab|cm"},
                  ctx.Where() + ": a fixed angle is '@ <angle> <unit> lab|cm', e.g. @ 90 deg cm");
        return false;
      }
      const io::Token &a = cur.Next();
      const io::Token &u = cur.Next();
      const io::Token &f = cur.Next();
      term.angleUnit = FindUnit(u.text);
      if (!ParseNumberText(a.text, term.angle) || !term.angleUnit || std::string(term.angleUnit->dimension) != "angle" ||
          (f.text != "lab" && f.text != "cm")) {
        ctx.Error(codes::kReactionSyntax, at.col, a.text + " " + u.text + " " + f.text, {"@ <angle> deg|rad lab|cm"},
                  ctx.Where() + ": a fixed angle is '@ <angle> <unit> lab|cm' with the unit deg or rad, e.g. @ 90 deg cm");
        return false;
      }
      term.angle.col = a.col;
      term.angleFrame = f.text;
      term.hasAngle = true;
    }
    out.terms.push_back(term);
    if (cur.AtEnd()) break;
    const io::Token &o = cur.Peek();
    if (o.text != "+" && o.text != "/") {
      ctx.Error(codes::kReactionSyntax, o.col, o.text, {"+", "/", "@"},
                ctx.Where() + ": '" + o.text + "' follows a reaction; terms are joined by ' + ' or ' / ' (with blanks)");
      return false;
    }
    if (op && op != o.text[0]) {
      ctx.Error(codes::kReactionSyntax, o.col, o.text, {},
                ctx.Where() + ": a reaction expression is a sum or a ratio of two terms, not both (no ratio of sums)");
      return false;
    }
    op = o.text[0];
    cur.Next();
  }
  if (op == '/' && out.terms.size() != 2) {
    ctx.Error(codes::kReactionSyntax, out.terms[2].reaction.col, out.terms[2].reaction.text, {},
              ctx.Where() + ": a ratio has one numerator and one denominator");
    return false;
  }
  out.op = op == '+' ? ReactionExpr::Sum : op == '/' ? ReactionExpr::Ratio : ReactionExpr::Single;
  return true;
}

}  // namespace experiment
}  // namespace azure2
