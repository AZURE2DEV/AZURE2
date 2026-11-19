// The Plot tab's drawing with QPainter (PlotWidget.h).
#include "PlotWidget.h"

#include <algorithm>
#include <cmath>

#include <QFileInfo>
#include <QImage>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QPdfWriter>
#include <QRubberBand>
#include <QSvgGenerator>

namespace {

const int kLeft = 78, kRight = 16, kTop = 14, kBottom = 44, kGap = 30;

double Log(double v) { return std::log10(v); }

/// Linear ticks: about n nice steps (1, 2, 5 x 10^k) over [a, b].
QVector<double> LinearTicks(double a, double b, int n) {
  QVector<double> out;
  if (!(b > a)) return out;
  const double raw = (b - a) / n, mag = std::pow(10.0, std::floor(std::log10(raw)));
  const double r = raw / mag, step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
  for (double t = std::ceil(a / step) * step; t <= b + 1e-9 * step; t += step) out << (std::fabs(t) < 1e-12 * step ? 0 : t);
  return out;
}

/// Decade ticks over [a, b] (a, b > 0); with fewer than two decades, the 2 and 5 of each decade too.
QVector<double> LogTicks(double a, double b) {
  QVector<double> out;
  if (!(b > a) || a <= 0) return out;
  const int lo = (int)std::floor(Log(a)), hi = (int)std::ceil(Log(b));
  const int every = std::max(1, (hi - lo) / 8);
  for (int k = lo; k <= hi; k += every)
    for (double m : (hi - lo <= 2 ? QVector<double>{1, 2, 5} : QVector<double>{1})) {
      const double t = m * std::pow(10.0, k);
      if (t >= a * (1 - 1e-12) && t <= b * (1 + 1e-12)) out << t;
    }
  return out;
}

QString TickText(double v) {
  const double a = std::fabs(v);
  if (a != 0 && (a < 1e-3 || a >= 1e5)) return QString::number(v, 'e', 0).replace("e-0", "e-").replace("e+0", "e");
  return QString::number(v, 'g', 6);
}

void Marker(QPainter &p, int shape, QPointF c, double r) {
  switch (shape % 5) {
    case 1: p.drawRect(QRectF(c.x() - r, c.y() - r, 2 * r, 2 * r)); break;
    case 2: p.drawPolygon(QPolygonF({c + QPointF(0, -r - 1), c + QPointF(r + 1, 0), c + QPointF(0, r + 1), c + QPointF(-r - 1, 0)})); break;
    case 3: p.drawPolygon(QPolygonF({c + QPointF(0, -r - 1), c + QPointF(r + 1, r), c + QPointF(-r - 1, r)})); break;
    case 4: p.drawPolygon(QPolygonF({c + QPointF(0, r + 1), c + QPointF(r + 1, -r), c + QPointF(-r - 1, -r)})); break;
    default: p.drawEllipse(c, r, r);
  }
}

}  // namespace

PlotWidget::PlotWidget(QWidget *parent) : QWidget(parent) {
  setObjectName("plotWidget");
  setMinimumSize(320, 240);
  setAutoFillBackground(true);
  setToolTip(tr("Drag a rectangle to zoom; double-click to see the whole range."));
}

void PlotWidget::setSeries(const QList<PlotSeries> &series) {
  series_ = series;
  resetZoom();
}

void PlotWidget::setAxisTitles(const QString &x, const QString &y) {
  xTitle_ = x;
  yTitle_ = y;
  update();
}

void PlotWidget::setLogX(bool on) {
  logX_ = on;
  resetZoom();
}

void PlotWidget::setLogY(bool on) {
  logY_ = on;
  resetZoom();
}

void PlotWidget::setShowBand(bool on) {
  band_ = on;
  resetZoom();
}

void PlotWidget::setShowResiduals(bool on) {
  residuals_ = on;
  update();
}

void PlotWidget::setShowLegend(bool on) {
  legend_ = on;
  update();
}

bool PlotWidget::usable(double x, double y) const {
  return std::isfinite(x) && std::isfinite(y) && (!logX_ || x > 0) && (!logY_ || y > 0);
}

QRectF PlotWidget::dataRange() const {
  double x0 = INFINITY, x1 = -INFINITY, y0 = INFINITY, y1 = -INFINITY;
  auto take = [&](double x, double y) {
    if (!usable(x, y)) return;
    x0 = std::min(x0, x);
    x1 = std::max(x1, x);
    y0 = std::min(y0, y);
    y1 = std::max(y1, y);
  };
  for (const PlotSeries &s : series_) {
    for (int i = 0; i < s.dataX.size(); i++) {
      take(s.dataX[i], s.dataY[i]);
      if (i < s.dataError.size() && std::isfinite(s.dataError[i])) {
        take(s.dataX[i], s.dataY[i] + s.dataError[i]);
        if (!logY_ || s.dataY[i] - s.dataError[i] > 0) take(s.dataX[i], s.dataY[i] - s.dataError[i]);
      }
    }
    for (int i = 0; i < s.modelX.size(); i++) {
      take(s.modelX[i], s.modelY[i]);
      if (band_ && i < s.bandLow.size()) {
        take(s.modelX[i], s.bandHigh[i]);
        take(s.modelX[i], s.bandLow[i]);
      }
    }
  }
  if (!(x1 >= x0)) return QRectF(logX_ ? 0.1 : 0, logY_ ? 0.1 : 0, 1, 1);
  auto widen = [](double &a, double &b, bool log) {
    if (log) {
      const double la = Log(a), lb = Log(b), m = lb > la ? 0.05 * (lb - la) : 0.5;
      a = std::pow(10.0, la - m);
      b = std::pow(10.0, lb + m);
    } else {
      const double m = b > a ? 0.05 * (b - a) : (a != 0 ? 0.1 * std::fabs(a) : 1.0);
      a -= m;
      b += m;
    }
  };
  widen(x0, x1, logX_);
  widen(y0, y1, logY_);
  return QRectF(QPointF(x0, y0), QPointF(x1, y1));
}

void PlotWidget::zoomTo(const QRectF &data) {
  view_ = data.normalized();
  update();
}

void PlotWidget::resetZoom() {
  view_ = dataRange();
  update();
}

double PlotWidget::mapX(double x, const QRect &r) const {
  const double a = logX_ ? Log(view_.left()) : view_.left(), b = logX_ ? Log(view_.right()) : view_.right();
  const double v = logX_ ? Log(x) : x;
  return r.left() + (v - a) / (b - a) * r.width();
}

double PlotWidget::mapY(double y, const QRect &r, const QRectF &view, bool log) const {
  const double a = log ? Log(view.top()) : view.top(), b = log ? Log(view.bottom()) : view.bottom();
  const double v = log ? Log(y) : y;
  return r.bottom() - (v - a) / (b - a) * r.height();
}

double PlotWidget::unmapX(double px, const QRect &r) const {
  const double a = logX_ ? Log(view_.left()) : view_.left(), b = logX_ ? Log(view_.right()) : view_.right();
  const double v = a + (px - r.left()) / r.width() * (b - a);
  return logX_ ? std::pow(10.0, v) : v;
}

double PlotWidget::unmapY(double py, const QRect &r) const {
  const double a = logY_ ? Log(view_.top()) : view_.top(), b = logY_ ? Log(view_.bottom()) : view_.bottom();
  const double v = a + (r.bottom() - py) / r.height() * (b - a);
  return logY_ ? std::pow(10.0, v) : v;
}

PlotWidget::Frame PlotWidget::framesOf(const QRect &area) const {
  bool anyResidual = false;
  for (const PlotSeries &s : series_) anyResidual = anyResidual || !s.residual.isEmpty();
  const QRect inner = area.adjusted(kLeft, kTop, -kRight, -kBottom);
  Frame f;
  if (residuals_ && anyResidual && inner.height() > 200) {
    const int lowerHeight = inner.height() / 4;
    f.main = QRect(inner.left(), inner.top(), inner.width(), inner.height() - lowerHeight - kGap);
    f.lower = QRect(inner.left(), f.main.bottom() + kGap, inner.width(), lowerHeight);
  } else {
    f.main = inner;
  }
  return f;
}

int PlotWidget::dataPointsDrawn() const {
  int n = 0;
  for (const PlotSeries &s : series_)
    for (int i = 0; i < s.dataX.size(); i++) n += usable(s.dataX[i], s.dataY[i]);
  return n;
}

void PlotWidget::render(QPainter &p, const QRect &area) const {
  p.save();
  p.setRenderHint(QPainter::Antialiasing);
  p.fillRect(area, Qt::white);
  const Frame f = framesOf(area);
  const QRect &r = f.main;
  const QFontMetrics fm(p.font());

  // axes and ticks of the main panel
  p.setPen(QPen(QColor(0xdd, 0xdd, 0xdd), 1, Qt::DotLine));
  const QVector<double> xt = logX_ ? LogTicks(view_.left(), view_.right()) : LinearTicks(view_.left(), view_.right(), 7);
  const QVector<double> yt = logY_ ? LogTicks(view_.top(), view_.bottom()) : LinearTicks(view_.top(), view_.bottom(), 6);
  for (double t : xt) p.drawLine(QPointF(mapX(t, r), r.top()), QPointF(mapX(t, r), r.bottom()));
  for (double t : yt) p.drawLine(QPointF(r.left(), mapY(t, r, view_, logY_)), QPointF(r.right(), mapY(t, r, view_, logY_)));
  p.setPen(Qt::black);
  p.drawRect(r);
  const QRect xLabels = f.lower.isEmpty() ? r : f.lower;
  for (double t : xt) {
    const double x = mapX(t, r);
    p.drawLine(QPointF(x, xLabels.bottom()), QPointF(x, xLabels.bottom() - 5));
    p.drawLine(QPointF(x, r.bottom()), QPointF(x, r.bottom() - 5));
    const QString s = TickText(t);
    p.drawText(QPointF(x - fm.horizontalAdvance(s) / 2.0, xLabels.bottom() + fm.height()), s);
  }
  for (double t : yt) {
    const double y = mapY(t, r, view_, logY_);
    p.drawLine(QPointF(r.left(), y), QPointF(r.left() + 5, y));
    const QString s = TickText(t);
    p.drawText(QPointF(r.left() - 6 - fm.horizontalAdvance(s), y + fm.ascent() / 2.0 - 1), s);
  }
  p.drawText(QRect(area.left(), xLabels.bottom() + fm.height() + 2, area.width(), fm.height() + 4), Qt::AlignHCenter,
             xTitle_);
  p.save();
  p.translate(area.left() + fm.height() - 2, r.center().y());
  p.rotate(-90);
  p.drawText(QRect(-r.height() / 2, -fm.height(), r.height(), fm.height() + 4), Qt::AlignHCenter, yTitle_);
  p.restore();

  // the series, clipped to the panel
  p.save();
  p.setClipRect(r);
  for (const PlotSeries &s : series_) {
    QColor shade = s.color;
    shade.setAlpha(60);
    if (band_ && s.bandLow.size() == s.modelX.size() && !s.modelX.isEmpty()) {
      if (s.modelX.size() == 1) {
        if (usable(s.modelX[0], s.bandLow[0]) && usable(s.modelX[0], s.bandHigh[0])) {
          p.setPen(QPen(shade, 8));
          p.drawLine(QPointF(mapX(s.modelX[0], r), mapY(s.bandLow[0], r, view_, logY_)),
                     QPointF(mapX(s.modelX[0], r), mapY(s.bandHigh[0], r, view_, logY_)));
        }
      } else {
        QPolygonF poly;
        for (int i = 0; i < s.modelX.size(); i++)
          if (usable(s.modelX[i], s.bandHigh[i])) poly << QPointF(mapX(s.modelX[i], r), mapY(s.bandHigh[i], r, view_, logY_));
        for (int i = s.modelX.size() - 1; i >= 0; i--) {
          const double low = logY_ && !(s.bandLow[i] > 0) ? view_.top() : s.bandLow[i];
          if (usable(s.modelX[i], low)) poly << QPointF(mapX(s.modelX[i], r), mapY(low, r, view_, logY_));
        }
        p.setPen(Qt::NoPen);
        p.setBrush(shade);
        p.drawPolygon(poly);
        p.setBrush(Qt::NoBrush);
      }
    }
    // the calculated curve, a marker for a single point
    p.setPen(QPen(s.color, 2));
    if (s.modelX.size() == 1) {
      if (usable(s.modelX[0], s.modelY[0])) {
        const QPointF c(mapX(s.modelX[0], r), mapY(s.modelY[0], r, view_, logY_));
        p.drawLine(c - QPointF(9, 0), c + QPointF(9, 0));
      }
    } else {
      QPainterPath path;
      bool open = false;
      for (int i = 0; i < s.modelX.size(); i++) {
        if (!usable(s.modelX[i], s.modelY[i])) {
          open = false;
          continue;
        }
        const QPointF c(mapX(s.modelX[i], r), mapY(s.modelY[i], r, view_, logY_));
        if (open) path.lineTo(c);
        else path.moveTo(c);
        open = true;
      }
      p.drawPath(path);
    }
    // the data with error bars
    p.setPen(QPen(s.color.darker(130), 1.2));
    p.setBrush(s.color);
    for (int i = 0; i < s.dataX.size(); i++) {
      if (!usable(s.dataX[i], s.dataY[i])) continue;
      const QPointF c(mapX(s.dataX[i], r), mapY(s.dataY[i], r, view_, logY_));
      const double e = i < s.dataError.size() ? s.dataError[i] : NAN;
      if (std::isfinite(e) && e > 0) {
        const double lo = s.dataY[i] - e;
        const double ylo = logY_ && !(lo > 0) ? r.bottom() + 2 : mapY(lo, r, view_, logY_);
        const double yhi = mapY(s.dataY[i] + e, r, view_, logY_);
        p.drawLine(QPointF(c.x(), ylo), QPointF(c.x(), yhi));
        p.drawLine(QPointF(c.x() - 3, ylo), QPointF(c.x() + 3, ylo));
        p.drawLine(QPointF(c.x() - 3, yhi), QPointF(c.x() + 3, yhi));
      }
      Marker(p, s.marker, c, 3.2);
    }
    p.setBrush(Qt::NoBrush);
  }
  p.restore();

  // the residual panel: (data - model)/error on a linear axis
  if (!f.lower.isEmpty()) {
    const QRect &lr = f.lower;
    double m = 1;
    for (const PlotSeries &s : series_)
      for (int i = 0; i < s.residual.size(); i++)
        if (std::isfinite(s.residual[i]) && s.residualX[i] >= std::min(view_.left(), view_.right()) &&
            s.residualX[i] <= std::max(view_.left(), view_.right()))
          m = std::max(m, std::fabs(s.residual[i]));
    m *= 1.1;
    const QRectF lv(QPointF(view_.left(), -m), QPointF(view_.right(), m));
    p.setPen(Qt::black);
    p.drawRect(lr);
    p.setPen(QPen(Qt::gray, 1, Qt::DashLine));
    p.drawLine(QPointF(lr.left(), mapY(0, lr, lv, false)), QPointF(lr.right(), mapY(0, lr, lv, false)));
    p.setPen(Qt::black);
    for (double t : LinearTicks(-m, m, 4)) {
      const double y = mapY(t, lr, lv, false);
      p.drawLine(QPointF(lr.left(), y), QPointF(lr.left() + 5, y));
      const QString s = TickText(t);
      p.drawText(QPointF(lr.left() - 6 - fm.horizontalAdvance(s), y + fm.ascent() / 2.0 - 1), s);
    }
    p.save();
    p.translate(area.left() + fm.height() - 2, lr.center().y());
    p.rotate(-90);
    p.drawText(QRect(-lr.height(), -fm.height(), 2 * lr.height(), fm.height() + 4), Qt::AlignHCenter,
               tr("(data - calc)/error"));
    p.restore();
    p.setClipRect(lr);
    for (const PlotSeries &s : series_) {
      p.setPen(QPen(s.color.darker(130), 1));
      p.setBrush(s.color);
      for (int i = 0; i < s.residual.size(); i++) {
        if (!std::isfinite(s.residual[i]) || !std::isfinite(s.residualX[i]) || (logX_ && s.residualX[i] <= 0)) continue;
        Marker(p, s.marker, QPointF(mapX(s.residualX[i], r), mapY(s.residual[i], lr, lv, false)), 2.5);
      }
    }
    p.setClipping(false);
    p.setBrush(Qt::NoBrush);
  }

  // the legend
  if (legend_ && !series_.isEmpty()) {
    int width = 0;
    for (const PlotSeries &s : series_) width = std::max(width, fm.horizontalAdvance(s.label));
    const int row = fm.height() + 2;
    QRect box(r.right() - width - 46, r.top() + 6, width + 40, row * (int)series_.size() + 8);
    p.setPen(QColor(0x99, 0x99, 0x99));
    p.setBrush(QColor(255, 255, 255, 225));
    p.drawRect(box);
    for (int i = 0; i < series_.size(); i++) {
      const PlotSeries &s = series_[i];
      const int y = box.top() + 4 + i * row + row / 2;
      p.setPen(QPen(s.color, 2));
      p.drawLine(box.left() + 6, y, box.left() + 26, y);
      p.setBrush(s.color);
      if (!s.dataX.isEmpty()) Marker(p, s.marker, QPointF(box.left() + 16, y), 3);
      p.setBrush(Qt::NoBrush);
      p.setPen(Qt::black);
      p.drawText(QPointF(box.left() + 32, y + fm.ascent() / 2.0 - 1), s.label);
    }
  }
  p.restore();
}

void PlotWidget::paintEvent(QPaintEvent *) {
  QPainter p(this);
  render(p, rect());
}

bool PlotWidget::exportTo(const QString &file, QString &why) const {
  const QString suffix = QFileInfo(file).suffix().toLower();
  const QSize size(1000, 720);
  if (suffix == "png") {
    QImage image(size, QImage::Format_ARGB32);
    image.fill(Qt::white);
    QPainter p(&image);
    render(p, QRect(QPoint(0, 0), size));
    p.end();
    if (!image.save(file)) {
      why = tr("cannot write %1").arg(file);
      return false;
    }
    return true;
  }
  if (suffix == "svg") {
    QSvgGenerator svg;
    svg.setFileName(file);
    svg.setSize(size);
    svg.setViewBox(QRect(QPoint(0, 0), size));
    svg.setTitle(tr("AZURE2 plot"));
    QPainter p(&svg);
    render(p, QRect(QPoint(0, 0), size));
    return p.end();
  }
  if (suffix == "pdf") {
    QPdfWriter pdf(file);
    pdf.setPageSize(QPageSize(QSizeF(size.width() * 0.75, size.height() * 0.75), QPageSize::Point));
    pdf.setPageMargins(QMarginsF(0, 0, 0, 0));
    pdf.setResolution(96);
    QPainter p(&pdf);
    render(p, QRect(QPoint(0, 0), size));
    return p.end();
  }
  why = tr("%1: the plot is written as .png, .svg or .pdf").arg(file);
  return false;
}

void PlotWidget::mousePressEvent(QMouseEvent *e) {
  if (e->button() != Qt::LeftButton) return;
  origin_ = e->position().toPoint();
  if (!rubber_) rubber_ = new QRubberBand(QRubberBand::Rectangle, this);
  rubber_->setGeometry(QRect(origin_, QSize()));
  rubber_->show();
}

void PlotWidget::mouseMoveEvent(QMouseEvent *e) {
  if (rubber_ && rubber_->isVisible()) rubber_->setGeometry(QRect(origin_, e->position().toPoint()).normalized());
}

void PlotWidget::mouseReleaseEvent(QMouseEvent *e) {
  if (!rubber_ || !rubber_->isVisible()) return;
  rubber_->hide();
  const QRect band = QRect(origin_, e->position().toPoint()).normalized();
  if (band.width() < 6 || band.height() < 6) return;
  const QRect r = framesOf(rect()).main;
  zoomTo(QRectF(QPointF(unmapX(band.left(), r), unmapY(band.bottom(), r)),
                QPointF(unmapX(band.right(), r), unmapY(band.top(), r))));
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent *) { resetZoom(); }
