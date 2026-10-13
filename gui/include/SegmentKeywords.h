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

/*! True when the tail carries "natural 1". */
inline bool isNatural(const QString &extras) {
  QStringList tokens = extras.split(' ', Qt::SkipEmptyParts);
  for (int i = 0; i + 1 < tokens.size(); i++)
    if (tokens[i].compare("natural", Qt::CaseInsensitive) == 0) return tokens[i + 1].toInt() != 0;
  return false;
}

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

/*! The tail with the natural flag set or cleared; "natural 1" is written last,
 *  the order the engine's writers use, and a cleared flag leaves no block. */
inline QString withNatural(const QString &extras, bool natural) {
  QString rest = withoutNatural(extras);
  if (!natural) return rest;
  return rest.isEmpty() ? QString("natural 1") : rest + " natural 1";
}

}  // namespace SegmentKeywords

#endif
