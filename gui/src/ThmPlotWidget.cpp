#include "ThmPlotWidget.h"

#include <QAbstractTextDocumentLayout>
#include <QPainter>
#include <QPainterPath>
#include <QTextDocument>
#include <algorithm>
#include <cmath>

QColor thmPlotColor(int i) {
  static const QColor colors[] = {QColor(0x00, 0x72, 0xB2), QColor(0xD5, 0x5E, 0x00), QColor(0x00, 0x9E, 0x73),
                                  QColor(0xCC, 0x79, 0xA7), QColor(0xE6, 0x9F, 0x00), QColor(0x56, 0xB4, 0xE9)};
  return colors[((i % 6) + 6) % 6];
}

ThmPlotWidget::ThmPlotWidget(QWidget *parent) : QWidget(parent) {
  setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  setAutoFillBackground(true);
  setBackgroundRole(QPalette::Base);
}

void ThmPlotWidget::clear() {
  series_.clear();
  markers_.clear();
  message_.clear();
  update();
}

void ThmPlotWidget::setTitle(const QString &html) {
  title_ = html;
  update();
}

void ThmPlotWidget::setAxisLabels(const QString &x, const QString &y) {
  xLabel_ = x;
  yLabel_ = y;
  update();
}

void ThmPlotWidget::setLogY(bool log) {
  logY_ = log;
  update();
}

void ThmPlotWidget::addSeries(const Series &s) {
  series_ << s;
  update();
}

void ThmPlotWidget::addMarker(const Marker &m) {
  markers_ << m;
  update();
}

void ThmPlotWidget::setMessage(const QString &text) {
  message_ = text;
  update();
}

QVector<double> ThmPlotWidget::linearTicks(double lo, double hi, int target) {
  QVector<double> ticks;
  if (!(hi > lo) || !std::isfinite(lo) || !std::isfinite(hi)) return ticks;
  const double raw = (hi - lo) / std::max(target, 1);
  const double mag = std::pow(10.0, std::floor(std::log10(raw)));
  double step = mag;
  for (double f : {1.0, 2.0, 5.0, 10.0})
    if (f * mag >= raw) {
      step = f * mag;
      break;
    }
  for (double t = std::ceil(lo / step - 1e-9) * step; t <= hi + 1e-9 * step; t += step)
    ticks << (std::fabs(t) < 1e-12 * step ? 0.0 : t);
  return ticks;
}

namespace {

QString tickText(double v, double step) {
  // Enough decimals for the step, no more.
  int decimals = 0;
  while (decimals < 6 && std::fabs(step * std::pow(10.0, decimals) - std::round(step * std::pow(10.0, decimals))) >
                             1e-6 * step * std::pow(10.0, decimals))
    decimals++;
  return QString::number(v, 'f', decimals);
}

}  // namespace

void ThmPlotWidget::paintEvent(QPaintEvent *) {
  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing, true);
  const QPalette pal = palette();
  const QColor ink = pal.color(QPalette::Text);
  p.fillRect(rect(), pal.color(QPalette::Base));

  QFont tickFont = font();
  tickFont.setPointSizeF(tickFont.pointSizeF() > 0 ? tickFont.pointSizeF() * 1.05 : 10.0);
  QFont labelFont = tickFont;
  const QFontMetrics fm(tickFont);

  // Data ranges.
  double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
  for (const Series &s : series_)
    for (int i = 0; i < s.x.size() && i < s.y.size(); i++) {
      const double y = s.y[i];
      if (!std::isfinite(s.x[i]) || !std::isfinite(y) || (logY_ && !(y > 0.0))) continue;
      x0 = std::min(x0, s.x[i]);
      x1 = std::max(x1, s.x[i]);
      y0 = std::min(y0, y);
      y1 = std::max(y1, y);
    }
  const bool empty = !(x1 > x0) || !(y1 >= y0);

  // Layout: room for tick labels and axis labels, and for the legend above
  // the axes (a strip of entries that wraps), so it never covers a curve or
  // the label inside.
  const int left = fm.horizontalAdvance("-0.000") + fm.height() + 14;
  const int bottom = 2 * fm.height() + 14;
  const int lineLen = 24;
  struct Entry {
    const Series *s;
    QRect r;
  };
  QList<Entry> legend;
  {
    int x = left, y = 6;
    for (const Series &s : series_) {
      if (s.label.isEmpty()) continue;
      const int w = lineLen + 6 + fm.horizontalAdvance(s.label) + 12;
      if (x > left && x + w > width() - 8) {
        x = left;
        y += fm.height() + 2;
      }
      legend << Entry{&s, QRect(x, y, w, fm.height())};
      x += w;
    }
  }
  const int top = legend.isEmpty() ? 10 : legend.last().r.bottom() + 8;
  const QRect frame(left, top, std::max(10, width() - left - 12), std::max(10, height() - bottom - top));

  p.setPen(QPen(ink, 1));
  p.drawRect(frame);
  p.setFont(labelFont);
  if (!xLabel_.isEmpty())
    p.drawText(QRect(frame.left(), height() - fm.height() - 4, frame.width(), fm.height()), Qt::AlignCenter, xLabel_);
  if (!yLabel_.isEmpty()) {
    p.save();
    p.translate(fm.height() / 2 + 2, frame.center().y());
    p.rotate(-90);
    p.drawText(QRect(-frame.height() / 2, -fm.height() / 2, frame.height(), fm.height()), Qt::AlignCenter, yLabel_);
    p.restore();
  }

  if (empty || !message_.isEmpty()) {
    p.setPen(pal.color(QPalette::PlaceholderText));
    p.drawText(frame.adjusted(8, 8, -8, -8), Qt::AlignCenter | Qt::TextWordWrap,
               message_.isEmpty() ? tr("nothing to show") : message_);
    return;
  }

  // Axis ranges with a margin; log axes span whole decades when short.
  double ylo = y0, yhi = y1;
  if (logY_) {
    ylo = std::log10(y0);
    yhi = std::log10(y1);
    if (yhi - ylo < 1.0) {
      const double mid = 0.5 * (ylo + yhi);
      ylo = mid - 0.5;
      yhi = mid + 0.5;
    }
    const double pad = 0.06 * (yhi - ylo);
    ylo -= pad;
    yhi += pad + 0.14 * (yhi - ylo);  // room for the label
  } else {
    if (yhi - ylo < 1e-12 * std::max(1.0, std::fabs(yhi))) {
      ylo -= 0.5 * std::max(1e-3, std::fabs(ylo));
      yhi += 0.5 * std::max(1e-3, std::fabs(yhi));
    }
    const double pad = 0.06 * (yhi - ylo);
    ylo -= pad;
    yhi += pad + 0.15 * (yhi - ylo);  // room for the label
  }
  const double xlo = x0, xhi = x1;
  auto X = [&](double x) { return frame.left() + (x - xlo) / (xhi - xlo) * frame.width(); };
  auto Y = [&](double y) {
    const double v = logY_ ? std::log10(y) : y;
    return frame.bottom() - (v - ylo) / (yhi - ylo) * frame.height();
  };

  // Ticks, inward, with labels outside.
  p.setFont(tickFont);
  const int tick = 6;
  const QVector<double> xt = linearTicks(xlo, xhi, std::max(4, frame.width() / 80));
  const double xstep = xt.size() > 1 ? xt[1] - xt[0] : 1.0;
  for (double t : xt) {
    const double px = X(t);
    p.drawLine(QPointF(px, frame.bottom()), QPointF(px, frame.bottom() - tick));
    p.drawLine(QPointF(px, frame.top()), QPointF(px, frame.top() + tick));
    p.drawText(QRectF(px - 40, frame.bottom() + 3, 80, fm.height()), Qt::AlignHCenter | Qt::AlignTop,
               tickText(t, xstep));
  }
  if (logY_) {
    const int d0 = (int)std::ceil(ylo - 1e-9), d1 = (int)std::floor(yhi + 1e-9);
    const int every = std::max(1, (d1 - d0 + 1) / 6 + ((d1 - d0 + 1) % 6 ? 1 : 0));
    for (int d = d0; d <= d1; d++) {
      const double py = frame.bottom() - (d - ylo) / (yhi - ylo) * frame.height();
      p.drawLine(QPointF(frame.left(), py), QPointF(frame.left() + tick, py));
      p.drawLine(QPointF(frame.right(), py), QPointF(frame.right() - tick, py));
      if ((d - d0) % every == 0)
        p.drawText(QRectF(0, py - fm.height() / 2.0, frame.left() - 4, fm.height()), Qt::AlignRight | Qt::AlignVCenter,
                   QString("1e%1").arg(d));
    }
  } else {
    const QVector<double> yt = linearTicks(ylo, yhi, std::max(4, frame.height() / 45));
    const double ystep = yt.size() > 1 ? yt[1] - yt[0] : 1.0;
    for (double t : yt) {
      const double py = frame.bottom() - (t - ylo) / (yhi - ylo) * frame.height();
      p.drawLine(QPointF(frame.left(), py), QPointF(frame.left() + tick, py));
      p.drawLine(QPointF(frame.right(), py), QPointF(frame.right() - tick, py));
      p.drawText(QRectF(0, py - fm.height() / 2.0, frame.left() - 4, fm.height()), Qt::AlignRight | Qt::AlignVCenter,
                 std::fabs(ystep) < 1e-3 || std::fabs(t) >= 1e5 ? QString::number(t, 'g', 3) : tickText(t, ystep));
    }
  }

  // Markers and curves, clipped to the frame.
  p.save();
  p.setClipRect(frame.adjusted(1, 1, -1, -1));
  for (const Marker &m : markers_) {
    if (m.x < xlo || m.x > xhi) continue;
    p.setPen(QPen(m.color.isValid() ? m.color : ink, 1.2, Qt::DashLine));
    p.drawLine(QPointF(X(m.x), frame.top()), QPointF(X(m.x), frame.bottom()));
  }
  for (const Series &s : series_) {
    QPainterPath path;
    bool pen = false;
    for (int i = 0; i < s.x.size() && i < s.y.size(); i++) {
      if (!std::isfinite(s.x[i]) || !std::isfinite(s.y[i]) || (logY_ && !(s.y[i] > 0.0))) {
        pen = false;
        continue;
      }
      const QPointF q(X(s.x[i]), Y(s.y[i]));
      if (pen)
        path.lineTo(q);
      else
        path.moveTo(q);
      pen = true;
    }
    p.setPen(QPen(s.color.isValid() ? s.color : ink, 2.0, s.style));
    p.drawPath(path);
  }
  p.restore();
  for (const Marker &m : markers_) {
    if (m.x < xlo || m.x > xhi || m.label.isEmpty()) continue;
    p.setPen(m.color.isValid() ? m.color : ink);
    p.drawText(QRectF(X(m.x) + 3, frame.bottom() - fm.height() - tick - 2, 120, fm.height()),
               Qt::AlignLeft | Qt::AlignVCenter, m.label);
  }

  // The reaction, bold, inside the axes, in the corner the curves and marker
  // lines cross least (upper left on a tie), over a light box so it stays
  // readable where no corner is free.
  if (!title_.isEmpty()) {
    QTextDocument doc;
    QFont f = tickFont;
    f.setBold(true);
    f.setPointSizeF(f.pointSizeF() * 1.1);
    doc.setDefaultFont(f);
    doc.setHtml(QString("<b>%1</b>").arg(title_));
    const QSizeF ts = doc.size();
    const double pad = 6;
    const QRectF inner = QRectF(frame).adjusted(pad, pad, -pad, -pad);
    const QRectF corners[4] = {
        QRectF(inner.left(), inner.top(), ts.width(), ts.height()),
        QRectF(inner.right() - ts.width(), inner.top(), ts.width(), ts.height()),
        QRectF(inner.left(), inner.bottom() - ts.height(), ts.width(), ts.height()),
        QRectF(inner.right() - ts.width(), inner.bottom() - ts.height(), ts.width(), ts.height())};
    // Occupancy: sample every drawn segment finely and count samples in each
    // corner box, grown by a small margin; marker lines count too.
    int hits[4] = {0, 0, 0, 0};
    auto count = [&](const QPointF &q) {
      for (int c = 0; c < 4; c++)
        if (corners[c].adjusted(-4, -4, 4, 4).contains(q)) hits[c]++;
    };
    for (const Series &s : series_) {
      bool have = false;
      QPointF prev;
      for (int i = 0; i < s.x.size() && i < s.y.size(); i++) {
        if (!std::isfinite(s.x[i]) || !std::isfinite(s.y[i]) || (logY_ && !(s.y[i] > 0.0))) {
          have = false;
          continue;
        }
        const QPointF q(X(s.x[i]), Y(s.y[i]));
        if (have) {
          const int n = std::max(1, (int)(std::hypot(q.x() - prev.x(), q.y() - prev.y()) / 2.0));
          for (int k = 0; k <= n; k++) count(prev + (q - prev) * (double(k) / n));
        } else {
          count(q);
        }
        prev = q;
        have = true;
      }
    }
    for (const Marker &m : markers_) {
      if (m.x < xlo || m.x > xhi) continue;
      for (double yy = frame.top(); yy <= frame.bottom(); yy += 2.0) count(QPointF(X(m.x), yy));
    }
    int best = 0;
    for (int c = 1; c < 4; c++)
      if (hits[c] < hits[best]) best = c;
    const QRectF box = corners[best].adjusted(-3, -1, 3, 1);
    QColor bg = pal.color(QPalette::Base);
    bg.setAlpha(225);
    p.fillRect(box, bg);
    p.save();
    p.translate(corners[best].topLeft());
    QAbstractTextDocumentLayout::PaintContext ctx;
    ctx.palette.setColor(QPalette::Text, ink);
    doc.documentLayout()->draw(&p, ctx);
    p.restore();
  }

  // Legend, above the axes.
  for (const Entry &e : legend) {
    const double yy = e.r.center().y();
    p.setPen(QPen(e.s->color.isValid() ? e.s->color : ink, 2.0, e.s->style));
    p.drawLine(QPointF(e.r.left(), yy), QPointF(e.r.left() + lineLen, yy));
    p.setPen(ink);
    p.drawText(QRect(e.r.left() + lineLen + 6, e.r.top(), e.r.width() - lineLen - 6, e.r.height()),
               Qt::AlignLeft | Qt::AlignVCenter, e.s->label);
  }
}
