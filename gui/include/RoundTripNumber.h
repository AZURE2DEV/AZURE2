#ifndef ROUNDTRIPNUMBER_H
#define ROUNDTRIPNUMBER_H

#include <QString>

/*!
 * Text for a number the engine will read back from the .azr: the fewest
 * significant digits (six or more, in QTextStream's 'g' style) that parse back
 * to exactly the same double.
 *
 * QTextStream writes doubles with six significant digits, so every save used to
 * round the model -- a level at 27.49431 MeV came back at 27.4943, a width of
 * 55469980 eV as 5.547e+07 -- and moved the project's physics.  Starting at six
 * digits keeps every value that already round-tripped byte-identical to what
 * earlier versions wrote.
 */
inline QString roundTripNumber(double x) {
  for (int precision = 6; precision < 17; precision++) {
    QString text = QString::number(x, 'g', precision);
    if (text.toDouble() == x) return text;
  }
  return QString::number(x, 'g', 17);
}

/*!
 * roundTripNumber for a fixed-width column of a left-aligned QTextStream
 * (qSetFieldWidth(width)): a field pads only up to its width, so text that
 * fills it gets a trailing space -- 0.30000000000000004 would otherwise run
 * into the next column and the engine would read one fused token.
 */
inline QString roundTripField(double x, int width = 15) {
  QString text = roundTripNumber(x);
  return text.size() >= width ? text + " " : text;
}

#endif
