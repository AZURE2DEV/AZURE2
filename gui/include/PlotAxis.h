#ifndef PLOTAXIS_H
#define PLOTAXIS_H

#include <algorithm>
#include <cmath>

#include <QString>
#include <QVector>

/// The axis ticks shared by the plot (PlotWidget) and the map (HeatmapWidget): one copy of the tick rules.
namespace plotaxis {

/// Linear ticks: about n nice steps (1, 2, 5 x 10^k) over [a, b].
inline QVector<double> LinearTicks(double a, double b, int n) {
  QVector<double> out;
  if (!(b > a) || !std::isfinite(a) || !std::isfinite(b)) return out;
  // a range at the doubles' resolution (|U| = 1 to 1e-16) has no ticks: t += step would not move
  if (b - a < 1e-12 * std::max(std::fabs(a), std::fabs(b))) return out;
  const double raw = (b - a) / n, mag = std::pow(10.0, std::floor(std::log10(raw)));
  const double r = raw / mag, step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
  const double first = std::ceil(a / step);
  for (int k = 0; k <= 4 * n + 4; k++) {
    const double t = (first + k) * step;
    if (t > b + 1e-9 * step) break;
    out << (std::fabs(t) < 1e-12 * step ? 0 : t);
  }
  return out;
}

/// Decade ticks over [a, b] (a, b > 0); with fewer than two decades, the 2 and 5 of each decade too.
inline QVector<double> LogTicks(double a, double b) {
  QVector<double> out;
  if (!(b > a) || a <= 0) return out;
  const int lo = (int)std::floor(std::log10(a)), hi = (int)std::ceil(std::log10(b));
  const int every = std::max(1, (hi - lo) / 8);
  for (int k = lo; k <= hi; k += every)
    for (double m : (hi - lo <= 2 ? QVector<double>{1, 2, 5} : QVector<double>{1})) {
      const double t = m * std::pow(10.0, k);
      if (t >= a * (1 - 1e-12) && t <= b * (1 + 1e-12)) out << t;
    }
  return out;
}

/// A tick's text: plain up to 5 digits, else 1e-4 style.
inline QString TickText(double v) {
  const double a = std::fabs(v);
  if (a != 0 && (a < 1e-3 || a >= 1e5)) return QString::number(v, 'e', 0).replace("e-0", "e-").replace("e+0", "e");
  return QString::number(v, 'g', 6);
}

}  // namespace plotaxis

#endif
