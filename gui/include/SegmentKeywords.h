#ifndef SEGMENTKEYWORDS_H
#define SEGMENTKEYWORDS_H

#include <algorithm>

#include <QString>
#include <QStringList>

/*!
 * The keyword blocks a segment line may end with (SegLine.h / ExtrapLine.h):
 * "sqrtshift <b> <bError> <vary>", "sample <1|name>" (alias "natural", the keyword before
 * samples existed), "shiftgroup <label>", "shiftfngroup <label>", "normgroup <label>", "covariance <file>", "total <0|1>", "absolute <0|1>" (angular-coefficient data are the absolute A_k) and "shiftfn \"<f(x)>\" N v s vary ...".  The GUI keeps the tail verbatim in
 * SegmentsDataData::extraTokens / SegmentsTestData::extraTokens and edits only the sample
 * token through these helpers, so a block it does not know survives a load/save untouched.
 */
namespace SegmentKeywords {

inline bool isKeyword(const QString &token) {
  // sqrtshift, sample (alias natural), shiftgroup, shiftfngroup, normgroup, covariance, shiftfn, total and absolute are
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
         token.compare("mixture", Qt::CaseInsensitive) == 0;
}

inline bool isSampleKeyword(const QString &token) {
  return token.compare("sample", Qt::CaseInsensitive) == 0 || token.compare("natural", Qt::CaseInsensitive) == 0;
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

}  // namespace SegmentKeywords

#endif
