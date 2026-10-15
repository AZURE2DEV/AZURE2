#ifndef SEGMENTKEYWORDS_H
#define SEGMENTKEYWORDS_H

#include <QString>
#include <QStringList>

/*!
 * The keyword blocks a segment line may end with (SegLine.h / ExtrapLine.h):
 * "sqrtshift <b> <bError> <vary>", "sample <1|name>" (alias "natural", the keyword before
 * samples existed), "shiftgroup <label>", "normgroup <label>" and "shiftfn \"<f(x)>\" N v s vary ...".  The GUI keeps the tail verbatim in
 * SegmentsDataData::extraTokens / SegmentsTestData::extraTokens and edits only the sample
 * token through these helpers, so a block it does not know survives a load/save untouched.
 */
namespace SegmentKeywords {

inline bool isKeyword(const QString &token) {
  // sqrtshift, sample (alias natural), shiftgroup, normgroup and shiftfn are the engine's (SegLine.h);
  // mixture is the retired 2026-10-12 block, still carried through a load/save of an old file.
  return token.compare("sqrtshift", Qt::CaseInsensitive) == 0 ||
         token.compare("sample", Qt::CaseInsensitive) == 0 ||
         token.compare("natural", Qt::CaseInsensitive) == 0 ||
         token.compare("shiftgroup", Qt::CaseInsensitive) == 0 ||
         token.compare("normgroup", Qt::CaseInsensitive) == 0 ||
         token.compare("shiftfn", Qt::CaseInsensitive) == 0 ||
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

}  // namespace SegmentKeywords

#endif
