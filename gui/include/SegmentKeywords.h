#ifndef SEGMENTKEYWORDS_H
#define SEGMENTKEYWORDS_H

#include <algorithm>

#include <QString>
#include <QStringList>

/*!
 * The keyword blocks a segment line may end with (SegLine.h / ExtrapLine.h):
 * "sqrtshift <b> <bError> <vary>", "sample <1|name>" (alias "natural", the keyword before
 * samples existed), "shiftgroup <label>", "shiftfngroup <label>", "normgroup <label>", "covariance <file>", "total <0|1>", "transmission <n> [<vary> <sigma>]", "yield <n> [<vary> <sigma>]" (a sample transmission or self-shielded yield, n in atoms/barn; with vary and sigma n is the fit parameter segment_<key>_thickness), "thicknessgroup <label>", "absolute <0|1>" (angular-coefficient data are the absolute A_k) and "shiftfn \"<f(x)>\" N v s vary ...".  The GUI keeps the tail verbatim in
 * SegmentsDataData::extraTokens / SegmentsTestData::extraTokens and edits only the sample, covariance, spin,
 * transmission / yield and group (normgroup / shiftgroup / shiftfngroup / thicknessgroup) blocks through these
 * helpers, so a block it does not know survives a load/save untouched.
 */
namespace SegmentKeywords {

inline bool isKeyword(const QString &token) {
  // sqrtshift, sample (alias natural), shiftgroup, shiftfngroup, normgroup, covariance, shiftfn, total, transmission,
  // yield and absolute are
  // the engine's (SegLine.h);
  // mixture is the retired 2026-10-12 block, still carried through a load/save of an old file.
  return token.compare("sqrtshift", Qt::CaseInsensitive) == 0 ||
         token.compare("sample", Qt::CaseInsensitive) == 0 ||
         token.compare("natural", Qt::CaseInsensitive) == 0 ||
         token.compare("shiftgroup", Qt::CaseInsensitive) == 0 ||
         token.compare("shiftfngroup", Qt::CaseInsensitive) == 0 ||
         token.compare("normgroup", Qt::CaseInsensitive) == 0 ||
         token.compare("covariance", Qt::CaseInsensitive) == 0 ||
         token.compare("shiftfn", Qt::CaseInsensitive) == 0 ||
         token.compare("total", Qt::CaseInsensitive) == 0 ||
         token.compare("absolute", Qt::CaseInsensitive) == 0 ||
         token.compare("spin", Qt::CaseInsensitive) == 0 ||
         token.compare("transmission", Qt::CaseInsensitive) == 0 ||
         token.compare("yield", Qt::CaseInsensitive) == 0 ||
         token.compare("thicknessgroup", Qt::CaseInsensitive) == 0 ||
         token.compare("mixture", Qt::CaseInsensitive) == 0;
}

inline bool isSampleKeyword(const QString &token) {
  return token.compare("sample", Qt::CaseInsensitive) == 0 || token.compare("natural", Qt::CaseInsensitive) == 0;
}

/*! The channel spin of a phase-shift segment's tail ("spin <s>"); "" when none. */
inline QString phaseSpin(const QString &extras) {
  QStringList tokens = extras.split(' ', Qt::SkipEmptyParts);
  for (int i = 0; i + 1 < tokens.size(); i++)
    if (tokens[i].compare("spin", Qt::CaseInsensitive) == 0) return tokens[i + 1];
  return QString();
}

/*! The tail with its spin block set to \p s or removed (empty \p s); an existing block is replaced where it
 *  stands, a new one goes before a sample block (writers put that last) or at the end. */
inline QString withPhaseSpin(const QString &extras, const QString &s) {
  QStringList tokens = extras.split(' ', Qt::SkipEmptyParts);
  int at = -1;
  for (int i = 0; i < tokens.size(); i++)
    if (tokens[i].compare("spin", Qt::CaseInsensitive) == 0) {
      at = i;
      tokens.erase(tokens.begin() + i, tokens.begin() + std::min<int>(i + 2, tokens.size()));
      break;
    }
  const QString v = s.trimmed();
  if (!v.isEmpty()) {
    if (at < 0) {
      at = tokens.size();
      for (int i = 0; i < tokens.size(); i++)
        if (isSampleKeyword(tokens[i])) {
          at = i;
          break;
        }
    }
    tokens.insert(at, v);
    tokens.insert(at, "spin");
  }
  return tokens.join(" ");
}

/*! The sample token of the tail: "" when there is none (or "sample 0"), "1" for the default
 *  sample, else the sample's name.  Reads "sample <t>" and the older "natural <t>". */
inline QString naturalSample(const QString &extras) {
  QStringList tokens = extras.split(' ', Qt::SkipEmptyParts);
  for (int i = 0; i + 1 < tokens.size(); i++)
    if (isSampleKeyword(tokens[i])) return tokens[i + 1] == "0" ? QString() : tokens[i + 1];
  return QString();
}

/*! True when the tail carries a sample block ("sample 1", "sample <name>" or the older "natural ..."). */
inline bool isNatural(const QString &extras) { return !naturalSample(extras).isEmpty(); }

/*! The tail with its sample block removed (every other block kept in order). */
inline QString withoutNatural(const QString &extras) {
  QStringList tokens = extras.split(' ', Qt::SkipEmptyParts);
  QStringList out;
  for (int i = 0; i < tokens.size(); i++) {
    if (isSampleKeyword(tokens[i])) {
      i++;  // skip the flag
      continue;
    }
    out.append(tokens[i]);
  }
  return out.join(' ');
}

/*! The tail with the sample token set ("1" or a sample name) or cleared (""); the block is
 *  written last as "sample <token>", the order the engine's writers use, and a cleared token
 *  leaves no block (an older "natural" block is replaced). */
inline QString withNaturalToken(const QString &extras, const QString &token) {
  QString rest = withoutNatural(extras);
  if (token.isEmpty()) return rest;
  return rest.isEmpty() ? "sample " + token : rest + " sample " + token;
}

/*! The tail with the default sample set or cleared. */
inline QString withNatural(const QString &extras, bool natural) {
  return withNaturalToken(extras, natural ? QString("1") : QString());
}

/*! The tail split into tokens; a double-quoted path with spaces stays one token (quotes kept). */
inline QStringList tailTokens(const QString &extras) {
  QStringList out;
  QString cur;
  bool quoted = false;
  for (const QChar c : extras) {
    if (c == '"') quoted = !quoted;
    if (c.isSpace() && !quoted) {
      if (!cur.isEmpty()) out.append(cur);
      cur.clear();
    } else {
      cur.append(c);
    }
  }
  if (!cur.isEmpty()) out.append(cur);
  return out;
}

/*! The data covariance file of the tail ("covariance <file>"), without surrounding quotes; "" when none. */
inline QString covarianceFile(const QString &extras) {
  const QStringList tokens = tailTokens(extras);
  for (int i = 0; i + 1 < tokens.size(); i++)
    if (tokens[i].compare("covariance", Qt::CaseInsensitive) == 0) {
      QString f = tokens[i + 1];
      if (f.size() >= 2 && f.startsWith('"') && f.endsWith('"')) f = f.mid(1, f.size() - 2);
      return f;
    }
  return QString();
}

/*! The tail with its covariance block set to \p file or removed (empty \p file).  An existing block is
 *  replaced where it stands; a new one goes before a sample block (the engine's writers put that last) or
 *  at the end.  A path with a space is written in double quotes, which the engine reads (ReadPathToken). */
inline QString withCovarianceFile(const QString &extras, const QString &file) {
  QStringList tokens = tailTokens(extras);
  int at = -1;
  for (int i = 0; i < tokens.size(); i++)
    if (tokens[i].compare("covariance", Qt::CaseInsensitive) == 0) {
      at = i;
      tokens.erase(tokens.begin() + i, tokens.begin() + std::min<int>(i + 2, tokens.size()));
      break;
    }
  const QString f = file.trimmed();
  if (!f.isEmpty()) {
    if (at < 0) {
      at = tokens.size();
      for (int i = 0; i < tokens.size(); i++)
        if (isSampleKeyword(tokens[i])) {
          at = i;
          break;
        }
    }
    tokens.insert(at, f.contains(' ') ? '"' + f + '"' : f);
    tokens.insert(at, "covariance");
  }
  return tokens.join(' ');
}

/*! The group keywords of a data-segment line, each followed by one label (SegLine.h, EData::ResolveNormGroups /
 *  ResolveShiftGroups / ResolveShiftFunctions): every active segment with the same label shares the normalization
 *  (normgroup), the constant energy shift (shiftgroup) or the shift function (shiftfngroup) of the group's first
 *  active segment in file order. */
inline QStringList groupKeywords() {
  return QStringList() << "normgroup" << "shiftgroup" << "shiftfngroup" << "thicknessgroup";
}

/*! True when \p token reads as a number (the engine's WholeNumber). */
inline bool isNumberToken(const QString &token) {
  bool ok = false;
  token.toDouble(&ok);
  return ok;
}

/*! Tokens of the keyword block that starts at tokens[i], the keyword included: 4 for sqrtshift, 3 + 3N for
 *  "shiftfn "<f>" N ...", 3 + 2n for the retired "mixture <w_self> <n> {"file" <w>}...", 4 for the fit form
 *  "transmission n vary sigma" / "yield n vary sigma" (a number after n), 2 for every other keyword, 1 for a token
 *  that is no keyword (never past the end).  Walking the tail block by block keeps a label
 *  or a value that happens to read like a keyword from being taken for one. */
inline int blockLength(const QStringList &tokens, int i) {
  const int left = tokens.size() - i;
  if (left <= 0) return 0;
  const QString k = tokens[i].toLower();
  int n = 1;
  if (k == "sqrtshift") {
    n = 4;
  } else if (k == "shiftfn") {
    bool ok = false;
    const int params = (i + 2 < tokens.size()) ? tokens[i + 2].toInt(&ok) : 0;
    n = 3 + (ok && params > 0 ? 3 * params : 0);
  } else if (k == "mixture") {
    bool ok = false;
    const int files = (i + 2 < tokens.size()) ? tokens[i + 2].toInt(&ok) : 0;
    n = 3 + (ok && files > 0 ? 2 * files : 0);
  } else if (k == "transmission" || k == "yield") {
    n = (i + 2 < tokens.size() && isNumberToken(tokens[i + 2])) ? 4 : 2;
  } else if (isKeyword(tokens[i])) {
    n = 2;
  }
  return std::min(n, left);
}

/*! Index of the token that starts the \p keyword block among the tail's tokens, or -1.  The engine reads the
 *  lower-case keyword only, so the match is exact. */
inline int blockIndex(const QStringList &tokens, const QString &keyword) {
  for (int i = 0; i < tokens.size(); i += std::max(1, blockLength(tokens, i)))
    if (tokens[i] == keyword) return i;
  return -1;
}

/*! The label of the tail's \p keyword block ("normgroup", "shiftgroup" or "shiftfngroup"); "" when none. */
inline QString groupLabel(const QString &extras, const QString &keyword) {
  const QStringList tokens = tailTokens(extras);
  const int i = blockIndex(tokens, keyword);
  return (i >= 0 && i + 1 < tokens.size()) ? tokens[i + 1] : QString();
}

/*! True when \p label can be written as a group label: one token, no blank and no double quote. */
inline bool isValidGroupLabel(const QString &label) {
  for (const QChar c : label)
    if (c.isSpace() || c == '"') return false;
  return true;
}

/*! The tail with its \p keyword block set to \p label or removed (empty \p label).  An existing block has its label
 *  replaced where it stands; a new one goes before a sample block (the engine's writers put that last) or at the
 *  end (the engine reads the blocks in any order).  Every other token is kept as it was. */
inline QString withGroupLabel(const QString &extras, const QString &keyword, const QString &label) {
  QStringList tokens = tailTokens(extras);
  const QString v = label.trimmed();
  const int at = blockIndex(tokens, keyword);
  if (at >= 0) {
    const int n = blockLength(tokens, at);
    if (v.isEmpty())
      tokens.erase(tokens.begin() + at, tokens.begin() + at + n);
    else if (n >= 2)
      tokens[at + 1] = v;
    else
      tokens.insert(at + 1, v);  // a bare keyword at the end of the tail: give it its label
    return tokens.join(' ');
  }
  if (v.isEmpty()) return tokens.join(' ');
  int where = tokens.size();
  for (int i = 0; i < tokens.size(); i += std::max(1, blockLength(tokens, i)))
    if (isSampleKeyword(tokens[i])) {
      where = i;
      break;
    }
  tokens.insert(where, v);
  tokens.insert(where, keyword);
  return tokens.join(' ');
}

/*!
 * The transmission / self-shielded-yield block of a segment line ("transmission n [vary sigma]", "yield n [vary
 * sigma]", SegLine.h ReadThicknessBlock): kind "" (none), "transmission" or "yield", the thickness n (atoms/barn) as
 * written, and the fit form (vary 0/1 and the prior width sigma, atoms/barn, 0 = none), which makes n the fit parameter
 * segment_<key>_thickness.  Texts are kept as written, so an untouched block compares equal and is not rewritten.
 */
struct SampleBlock {
  QString kind;
  QString n;
  bool fit = false;
  bool vary = false;
  QString sigma;
  bool operator==(const SampleBlock &o) const {
    return kind == o.kind && n == o.n && fit == o.fit && (!fit || (vary == o.vary && sigma == o.sigma));
  }
  bool operator!=(const SampleBlock &o) const { return !(*this == o); }
};

/*! The tail's transmission / yield block (kind "" when there is none). */
inline SampleBlock sampleBlock(const QString &extras) {
  const QStringList tokens = tailTokens(extras);
  SampleBlock b;
  for (const QString kw : {QString("transmission"), QString("yield")}) {
    const int i = blockIndex(tokens, kw);
    if (i < 0) continue;
    b.kind = kw;
    b.n = (i + 1 < tokens.size()) ? tokens[i + 1] : QString();
    if (blockLength(tokens, i) >= 4) {
      b.fit = true;
      b.vary = tokens[i + 2].toDouble() == 1.0;
      b.sigma = tokens[i + 3];
    }
    return b;
  }
  return b;
}

/*! The tail with its transmission / yield block set to \p b, or removed (empty kind).  An existing block is replaced
 *  where it stands; a new one goes before a sample block (the engine's writers put that last) or at the end.  Every
 *  other token is kept as it was. */
inline QString withSampleBlock(const QString &extras, const SampleBlock &b) {
  QStringList tokens = tailTokens(extras);
  int at = -1;
  for (const QString kw : {QString("transmission"), QString("yield")}) {
    const int i = blockIndex(tokens, kw);
    if (i < 0) continue;
    const int n = blockLength(tokens, i);
    tokens.erase(tokens.begin() + i, tokens.begin() + i + n);
    if (at < 0) at = i;
  }
  if (b.kind.isEmpty()) return tokens.join(' ');
  QStringList block;
  block << b.kind << b.n.trimmed();
  if (b.fit) block << (b.vary ? "1" : "0") << (b.sigma.trimmed().isEmpty() ? QString("0") : b.sigma.trimmed());
  if (at < 0) {
    at = tokens.size();
    for (int i = 0; i < tokens.size(); i += std::max(1, blockLength(tokens, i)))
      if (isSampleKeyword(tokens[i])) {
        at = i;
        break;
      }
  }
  for (int k = block.size() - 1; k >= 0; k--) tokens.insert(at, block[k]);
  return tokens.join(' ');
}

}  // namespace SegmentKeywords

#endif
