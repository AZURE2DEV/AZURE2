#ifndef SEGMENTKEYWORDS_H
#define SEGMENTKEYWORDS_H

#include <QString>
#include <QStringList>

/*!
 * The keyword blocks a segment line may end with (SegLine.h / ExtrapLine.h):
 * "sqrtshift <b> <bError> <vary>", "natural <0|1>" and "shiftgroup <label>".  The GUI keeps the
 * tail verbatim in SegmentsDataData::extraTokens / SegmentsTestData::extraTokens
 * and edits only the natural flag through these helpers, so a block it does
 * not know survives a load/save untouched.
 */
namespace SegmentKeywords {

inline bool isKeyword(const QString &token) {
  // sqrtshift, natural and shiftgroup are the engine's (SegLine.h); mixture is the
  // retired 2026-10-12 block, still carried through a load/save of an old file.
  return token.compare("sqrtshift", Qt::CaseInsensitive) == 0 ||
         token.compare("natural", Qt::CaseInsensitive) == 0 ||
         token.compare("shiftgroup", Qt::CaseInsensitive) == 0 ||
         token.compare("mixture", Qt::CaseInsensitive) == 0;
}

/*! The natural-target token of the tail: "" when there is none (or "natural 0"), "1" for the default
 *  sample, else the sample's name. */
inline QString naturalSample(const QString &extras) {
  QStringList tokens = extras.split(' ', Qt::SkipEmptyParts);
  for (int i = 0; i + 1 < tokens.size(); i++)
    if (tokens[i].compare("natural", Qt::CaseInsensitive) == 0) return tokens[i + 1] == "0" ? QString() : tokens[i + 1];
  return QString();
}

/*! True when the tail carries a natural block ("natural 1" or "natural <sample>"). */
inline bool isNatural(const QString &extras) { return !naturalSample(extras).isEmpty(); }

/*! The tail with its natural block removed (every other block kept in order). */
inline QString withoutNatural(const QString &extras) {
  QStringList tokens = extras.split(' ', Qt::SkipEmptyParts);
  QStringList out;
  for (int i = 0; i < tokens.size(); i++) {
    if (tokens[i].compare("natural", Qt::CaseInsensitive) == 0) {
      i++;  // skip the flag
      continue;
    }
    out.append(tokens[i]);
  }
  return out.join(' ');
}

/*! The tail with the natural token set ("1" or a sample name) or cleared (""); the block is
 *  written last, the order the engine's writers use, and a cleared token leaves no block. */
inline QString withNaturalToken(const QString &extras, const QString &token) {
  QString rest = withoutNatural(extras);
  if (token.isEmpty()) return rest;
  return rest.isEmpty() ? "natural " + token : rest + " natural " + token;
}

/*! The tail with the natural flag set (default sample) or cleared. */
inline QString withNatural(const QString &extras, bool natural) {
  return withNaturalToken(extras, natural ? QString("1") : QString());
}

}  // namespace SegmentKeywords

#endif
