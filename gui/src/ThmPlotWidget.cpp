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

void ThmPlotWidget::changed() {
  // The peers' frames depend on this one's legend and tick labels.
  update();
  for (ThmPlotWidget *p : rowPeers_ + columnPeers_)
    if (p) p->update();
}

void ThmPlotWidget::clear() {
  series_.clear();
  legendOnly_.clear();
  markers_.clear();
  bands_.clear();
  message_.clear();
  changed();
}

void ThmPlotWidget::setTitle(const QString &html) {
  title_ = html;
  changed();
}

void ThmPlotWidget::setAxisLabels(const QString &x, const QString &y) {
  xLabel_ = x;
  yLabel_ = y;
  changed();
}

void ThmPlotWidget::setLogY(bool log) {
  logY_ = log;
  changed();
}

void ThmPlotWidget::addSeries(const Series &s) {
  series_ << s;
  changed();
}

void ThmPlotWidget::addMarker(const Marker &m) {
  markers_ << m;
  changed();
}

void ThmPlotWidget::addBand(const Band &b) {
  bands_ << b;
  changed();
}

void ThmPlotWidget::setMessage(const QString &text) {
  message_ = text;
  changed();
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

void ThmPlotWidget::addLegendEntry(const QString &label, Qt::PenStyle style, const QColor &color) {
  Series s;
  s.label = label;
  s.style = style;
  s.color = color;
  legendOnly_ << s;
  changed();
}

void ThmPlotWidget::setAlignedWith(const QList<ThmPlotWidget *> &row, const QList<ThmPlotWidget *> &column) {
  rowPeers_ = row;
  columnPeers_ = column;
  rowPeers_.removeAll(this);
  columnPeers_.removeAll(this);
  update();
}

QFont ThmPlotWidget::tickFont() const {
  QFont f = font();
  if (f.pointSizeF() <= 0) f.setPointSizeF(9.0);
  return f;
}

ThmPlotWidget::Axes ThmPlotWidget::axes() const {
  Axes a;
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
  a.empty = !(x1 > x0) || !(y1 >= y0) || !message_.isEmpty();
  if (!(x1 > x0) || !(y1 >= y0)) return a;
  // A margin; log axes span whole decades when short; room at the top for the label.
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
    yhi += pad + 0.14 * (yhi - ylo);
  } else {
    if (yhi - ylo < 1e-12 * std::max(1.0, std::fabs(yhi))) {
      ylo -= 0.5 * std::max(1e-3, std::fabs(ylo));
      yhi += 0.5 * std::max(1e-3, std::fabs(yhi));
    }
    const double pad = 0.06 * (yhi - ylo);
    ylo -= pad;
    yhi += pad + 0.15 * (yhi - ylo);
  }
  a.xlo = x0;
  a.xhi = x1;
  a.ylo = ylo;
  a.yhi = yhi;
  return a;
}

QList<ThmPlotWidget::Tick> ThmPlotWidget::yTicks(const Axes &a, int frameHeight) const {
  QList<Tick> out;
  if (a.empty) return out;
  if (logY_) {
    const int d0 = (int)std::ceil(a.ylo - 1e-9), d1 = (int)std::floor(a.yhi + 1e-9);
    const int maxLabels = std::max(2, frameHeight / 28);
    const int n = d1 - d0 + 1;
    const int every = std::max(1, n / maxLabels + (n % maxLabels ? 1 : 0));
    for (int d = d0; d <= d1; d++) out << Tick{double(d), (d - d0) % every == 0 ? QString("1e%1").arg(d) : QString()};
  } else {
    // About one label per 35 px, and at least three.
    int target = std::max(4, frameHeight / 35);
    QVector<double> yt = linearTicks(a.ylo, a.yhi, target);
    while (yt.size() < 3 && target < 12) yt = linearTicks(a.ylo, a.yhi, target += 2);
    const double ystep = yt.size() > 1 ? yt[1] - yt[0] : 1.0;
    for (double t : yt)
      out << Tick{t, std::fabs(ystep) < 1e-3 || std::fabs(t) >= 1e5 ? QString::number(t, 'g', 3) : tickText(t, ystep)};
  }
  return out;
}

QList<ThmPlotWidget::Series> ThmPlotWidget::legendSeries() const {
  QList<Series> out;
  for (const Series &s : series_)
    if (!s.label.isEmpty()) out << s;
  for (const Series &s : legendOnly_)
    if (!s.label.isEmpty()) out << s;
  return out;
}

namespace {
const int kLineLen = 16;  // legend sample
const int kGap = 10;      // between legend entries
const int kPad = 8;       // widget edge to legend
QFont legendFont(QFont f) {
  f.setPointSizeF(f.pointSizeF() * 0.9);
  return f;
}
}  // namespace

QList<QList<QPair<ThmPlotWidget::Series, int>>> ThmPlotWidget::legendRows() const {
  const QFontMetrics fm(legendFont(tickFont()));
  QList<QList<QPair<Series, int>>> rows;
  int x = 0;
  const int avail = std::max(50, width() - 2 * kPad);
  for (const Series &s : legendSeries()) {
    const int w = kLineLen + 4 + fm.horizontalAdvance(s.label) + kGap;
    if (rows.isEmpty() || (x > 0 && x + w > avail)) {
      rows << QList<QPair<Series, int>>();
      x = 0;
    }
    rows.last() << qMakePair(s, w);
    x += w;
  }
  return rows;
}

int ThmPlotWidget::ownTop() const {
  const QFontMetrics fm(legendFont(tickFont()));
  const int n = legendRows().size();
  return n == 0 ? kPad : 4 + n * (fm.height() + 2) + 6;
}

int ThmPlotWidget::top() const {
  int t = ownTop();
  for (const ThmPlotWidget *p : rowPeers_)
    if (p && !p->isHidden()) t = std::max(t, p->ownTop());
  return t;
}

int ThmPlotWidget::bottomMargin() const {
  const QFontMetrics fm(tickFont());
  return 2 * fm.height() + 10;
}

int ThmPlotWidget::ownLeft() const {
  // Room for the rotated y label and the widest y tick label.
  const QFontMetrics fm(tickFont());
  int widest = fm.horizontalAdvance("0.00");
  const Axes a = axes();
  for (const Tick &t : yTicks(a, height() - top() - bottomMargin())) widest = std::max(widest, fm.horizontalAdvance(t.text));
  return fm.height() + 6 + widest + 6;
}

int ThmPlotWidget::left() const {
  int l = ownLeft();
  for (const ThmPlotWidget *p : columnPeers_)
    if (p && !p->isHidden()) l = std::max(l, p->ownLeft());
  return l;
}

QRect ThmPlotWidget::frameRect() const {
  const int l = left(), t = top();
  return QRect(l, t, std::max(10, width() - l - 14), std::max(10, height() - bottomMargin() - t));
}

void ThmPlotWidget::paintEvent(QPaintEvent *) {
  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing, true);
  const QPalette pal = palette();
  const QColor ink = pal.color(QPalette::Text);
  p.fillRect(rect(), pal.color(QPalette::Base));

  const QFont font = tickFont();
  const QFontMetrics fm(font);
  p.setFont(font);
  const Axes a = axes();
  const QRect frame = frameRect();

  p.setPen(QPen(ink, 1));
  p.drawRect(frame);
  if (!xLabel_.isEmpty())
    p.drawText(QRect(frame.left(), height() - fm.height() - 3, frame.width(), fm.height()), Qt::AlignCenter, xLabel_);
  if (!yLabel_.isEmpty()) {
    p.save();
    p.translate(4 + fm.height() / 2, frame.center().y());
    p.rotate(-90);
    p.drawText(QRect(-frame.height() / 2, -fm.height() / 2, frame.height(), fm.height()), Qt::AlignCenter, yLabel_);
    p.restore();
  }

  // Legend, above the axes, each row centred on the frame (inside the widget).
  {
    const QList<QList<QPair<Series, int>>> rows = legendRows();
    const QFont lf = legendFont(font);
    const QFontMetrics fm(lf);
    p.setFont(lf);
    int y = 4;
    for (const QList<QPair<Series, int>> &row : rows) {
      int total = 0;
      for (const auto &e : row) total += e.second;
      total -= kGap;  // no gap after the last entry
      int x = frame.center().x() - total / 2;
      x = std::max(kPad, std::min(x, width() - kPad - total));
      for (const auto &e : row) {
        const Series &s = e.first;
        const double yy = y + fm.height() / 2.0;
        const QColor c = s.color.isValid() ? s.color : ink;
        if (s.symbols) {
          p.setPen(QPen(c, 1.0));
          p.setBrush(c);
          p.drawEllipse(QPointF(x + kLineLen / 2.0, yy), 3.5, 3.5);
          p.setBrush(Qt::NoBrush);
        } else {
          p.setPen(QPen(c, 2.0, s.style));
          p.drawLine(QPointF(x, yy), QPointF(x + kLineLen, yy));
        }
        p.setPen(ink);
        p.drawText(QRect(x + kLineLen + 4, y, e.second - kLineLen - 4, fm.height()), Qt::AlignLeft | Qt::AlignVCenter,
                   s.label);
        x += e.second;
      }
      y += fm.height() + 2;
    }
    p.setFont(font);
  }

  if (a.empty) {
    p.setPen(pal.color(QPalette::PlaceholderText));
    p.drawText(frame.adjusted(8, 8, -8, -8), Qt::AlignCenter | Qt::TextWordWrap,
               message_.isEmpty() ? tr("nothing to show") : message_);
    return;
  }

  const double xlo = a.xlo, xhi = a.xhi, ylo = a.ylo, yhi = a.yhi;
  auto X = [&](double x) { return frame.left() + (x - xlo) / (xhi - xlo) * frame.width(); };
  auto Y = [&](double y) {
    const double v = logY_ ? std::log10(y) : y;
    return frame.bottom() - (v - ylo) / (yhi - ylo) * frame.height();
  };

  // Ticks, inward, with labels outside.
  const int tick = 5;
  p.setPen(QPen(ink, 1));
  int xTarget = std::max(4, frame.width() / 60);
  QVector<double> xt = linearTicks(xlo, xhi, xTarget);
  while (xt.size() < 3 && xTarget < 12) xt = linearTicks(xlo, xhi, xTarget += 2);
  const double xstep = xt.size() > 1 ? xt[1] - xt[0] : 1.0;
  for (double t : xt) {
    const double px = X(t);
    p.drawLine(QPointF(px, frame.bottom()), QPointF(px, frame.bottom() - tick));
    p.drawLine(QPointF(px, frame.top()), QPointF(px, frame.top() + tick));
    // Centred under the tick, but kept inside the widget at the right end.
    const QString text = tickText(t, xstep);
    const double half = fm.horizontalAdvance(text) / 2.0 + 1.0;
    const double cx = std::min(px, width() - 2.0 - half);
    p.drawText(QRectF(cx - half, frame.bottom() + 3, 2.0 * half, fm.height()), Qt::AlignHCenter | Qt::AlignTop, text);
  }
  for (const Tick &t : yTicks(a, frame.height())) {
    const double py = frame.bottom() - (t.at - ylo) / (yhi - ylo) * frame.height();
    p.drawLine(QPointF(frame.left(), py), QPointF(frame.left() + tick, py));
    p.drawLine(QPointF(frame.right(), py), QPointF(frame.right() - tick, py));
    if (!t.text.isEmpty())
      p.drawText(QRectF(0, py - fm.height() / 2.0, frame.left() - 5, fm.height()), Qt::AlignRight | Qt::AlignVCenter,
                 t.text);
  }

  // Markers and curves, clipped to the frame.
  p.save();
  p.setClipRect(frame.adjusted(1, 1, -1, -1));
  for (const Band &b : bands_) {
    QColor c = b.color.isValid() ? b.color : QColor(Qt::gray);
    c.setAlpha(60);
    const double x0 = X(std::max(std::min(b.x0, b.x1), xlo)), x1 = X(std::min(std::max(b.x0, b.x1), xhi));
    if (x1 < x0) continue;
    p.fillRect(QRectF(x0, frame.top(), std::max(2.0, x1 - x0), frame.height()), c);
  }
  for (const Marker &m : markers_) {
    if (m.x < xlo || m.x > xhi) continue;
    p.setPen(QPen(m.color.isValid() ? m.color : ink, 1.0, Qt::DashLine));
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
    const QColor c = s.color.isValid() ? s.color : ink;
    if (s.symbols) {
      p.setPen(QPen(c, 1.0));
      p.setBrush(c);
      for (int i = 0; i < s.x.size() && i < s.y.size(); i++)
        if (std::isfinite(s.x[i]) && std::isfinite(s.y[i]) && (!logY_ || s.y[i] > 0.0))
          p.drawEllipse(QPointF(X(s.x[i]), Y(s.y[i])), 3.0, 3.0);
      p.setBrush(Qt::NoBrush);
      continue;
    }
    p.setPen(QPen(c, 1.8, s.style));
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
    QFont f = font;
    f.setBold(true);
    doc.setDefaultFont(f);
    doc.setDocumentMargin(1);
    doc.setHtml(QString("<b>%1</b>").arg(title_));
    const QSizeF ts = doc.size();
    const double pad = 5;
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
        if (s.symbols) {
          for (int k = 0; k < 8; k++) count(q + QPointF(4.0 * std::cos(k * M_PI / 4), 4.0 * std::sin(k * M_PI / 4)));
          continue;
        }
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
    const QRectF box = corners[best].adjusted(-2, 0, 2, 0);
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
}
