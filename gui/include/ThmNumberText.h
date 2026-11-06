#ifndef THMNUMBERTEXT_H
#define THMNUMBERTEXT_H

#include <QLocale>
#include <QString>
#include <cmath>
#include <sstream>
#include <string>

/*!
 * Numbers as the THM pages of the GUI write and read them: written with the
 * shortest text that reads back as the same double (Qt's
 * FloatingPointShortest, which pyazr's _thm_number ports), read with
 * operator>> on an istringstream, as the engine reads them, so that "0.5 " or
 * "5" are taken exactly as AZURE2 takes them.
 */
namespace ThmText {

/// Shortest text that reads back as the same double.
inline QString number(double x) { return QString::number(x, 'g', QLocale::FloatingPointShortest); }

/// A leading number (strict: and nothing after it).
inline bool readDouble(const QString &text, double &x, bool strict = false) {
  std::istringstream s(text.toStdString());
  if (!(s >> x)) return false;
  std::string rest;
  return !strict || !(s >> rest);
}

inline bool readInt(const QString &text, int &x, bool strict = false) {
  std::istringstream s(text.toStdString());
  if (!(s >> x)) return false;
  std::string rest;
  return !strict || !(s >> rest);
}

/// A whole token as a finite number, as the engine reads Ebeam or a <levels> column.
inline bool readWholeDouble(const QString &text, double &x) {
  return readDouble(text.trimmed(), x, true) && std::isfinite(x);
}

}  // namespace ThmText

#endif
