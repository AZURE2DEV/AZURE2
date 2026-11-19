// The pole map's drawing with QPainter (HeatmapWidget.h, docs/upgrade/qt_gui.md section 15.5).
#include "HeatmapWidget.h"

#include <algorithm>
#include <cmath>

#include <QFile>
#include <QFileInfo>
#include <QImage>
#include <QMouseEvent>
#include <QPainter>
#include <QPdfWriter>
#include <QRubberBand>
#include <QSvgGenerator>
#include <QWheelEvent>

#include "PlotAxis.h"

namespace {

const int kRight = 12, kTop = 12, kBottom = 44, kBarGap = 14, kBarWidth = 16, kBarLabels = 64;

/// The colour map: nine anchors of a viridis-like ramp (dark blue-violet, blue, teal, green, yellow), lightness
/// rising monotonically from the first to the last, so that an order of values reads as an order of lightness in
/// colour and in grey, on light and dark paper alike.  Linear in RGB between anchors.
const int kRamp[9][3] = {{68, 1, 84},    {72, 40, 120},  {62, 74, 137},  {49, 104, 142}, {38, 130, 142},
                         {31, 158, 137}, {53, 183, 121}, {109, 205, 89}, {253, 231, 37}};

}  // namespace

QColor HeatmapWidget::colourAt(double t) {
  if (!std::isfinite(t)) return QColor();
  t = std::clamp(t, 0.0, 1.0) * 8;
  const int i = std::min(7, (int)std::floor(t));
  const double f = t - i;
  return QColor((int)std::lround(kRamp[i][0] + f * (kRamp[i + 1][0] - kRamp[i][0])),
                (int)std::lround(kRamp[i][1] + f * (kRamp[i + 1][1] - kRamp[i][1])),
                (int)std::lround(kRamp[i][2] + f * (kRamp[i + 1][2] - kRamp[i][2])));
}

HeatmapWidget::HeatmapWidget(QWidget *parent) : QWidget(parent) {
  setObjectName("heatmapWidget");
  setMinimumSize(320, 200);
  setMouseTracking(true);
  setFocusPolicy(Qt::WheelFocus);
  setToolTip(tr("Wheel: zoom in energy about the cursor (Shift: the rows).  Left drag: zoom to a rectangle.  Right or "
                "middle drag: pan.  Double-click: the whole map.  Click a row: its cut below."));
}

void HeatmapWidget::setMap(const QVector<double> &x, const QStringList &rows, const QVector<QVector<double>> &values,
                           const QString &xTitle, const QString &valueTitle) {
  x_ = x, rows_ = rows, values_ = values, xTitle_ = xTitle, valueTitle_ = valueTitle;
  rowTips_.clear();
  points_.clear();
  vertical_.clear();
  rangeOfColours();
  resetZoom();
}

void HeatmapWidget::setLogColour(bool on) {
  log_ = on;
  rangeOfColours();
  update();
}

double HeatmapWidget::colourValue(double v) const {
  if (!std::isfinite(v)) return NAN;
  if (log_) return v > 0 ? std::log10(v) : NAN;
  return v;
}

void HeatmapWidget::rangeOfColours() {
  low_ = INFINITY, high_ = -INFINITY;
  for (const QVector<double> &row : values_)
    for (double v : row) {
      const double c = colourValue(v);
      if (!std::isfinite(c)) continue;
      low_ = std::min(low_, c), high_ = std::max(high_, c);
    }
  if (!(high_ >= low_)) low_ = 0, high_ = 1;
  // a map of one value (|U| = 1 of a single open channel, to round-off): a range of +-0.5 around it
  if (high_ - low_ < 1e-9 * std::max(1.0, std::fabs(high_))) {
    const double mid = 0.5 * (low_ + high_);
    low_ = mid - 0.5, high_ = mid + 0.5;
  }
}

int HeatmapWidget::leftMargin(const QFontMetrics &fm) const {
  int w = 0;
  for (const QString &r : rows_) w = std::max(w, fm.horizontalAdvance(r));
  return std::min(w, 220) + 16;
}

HeatmapWidget::Frames HeatmapWidget::framesOf(const QRect &area) const {
  const QFontMetrics fm(font());
  Frames f;
  const int left = leftMargin(fm);
  const int barSpace = kBarGap + kBarWidth + kBarLabels + fm.height() + 4;
  f.map = QRect(area.left() + left, area.top() + kTop, std::max(10, area.width() - left - kRight - barSpace),
                std::max(10, area.height() - kTop - kBottom));
  f.bar = QRect(f.map.right() + kBarGap, f.map.top(), kBarWidth, f.map.height());
  return f;
}

double HeatmapWidget::mapX(double x, const QRect &r) const {
  return r.left() + (x - view_.left()) / (view_.right() - view_.left()) * r.width();
}
double HeatmapWidget::mapRow(double row, const QRect &r) const {
  return r.top() + (row - view_.top()) / (view_.bottom() - view_.top()) * r.height();
}
double HeatmapWidget::unmapX(double px, const QRect &r) const {
  return view_.left() + (px - r.left()) / r.width() * (view_.right() - view_.left());
}
double HeatmapWidget::unmapRow(double py, const QRect &r) const {
  return view_.top() + (py - r.top()) / r.height() * (view_.bottom() - view_.top());
}

void HeatmapWidget::zoomTo(double x0, double x1, double row0, double row1) {
  if (x1 < x0) std::swap(x0, x1);
  if (row1 < row0) std::swap(row0, row1);
  if (!(x1 > x0) || !(row1 > row0)) return;
  view_ = QRectF(QPointF(x0, row0), QPointF(x1, row1));
  update();
}

void HeatmapWidget::resetZoom() {
  if (x_.isEmpty()) view_ = QRectF(QPointF(0, -0.5), QPointF(1, 0.5));
  else {
    const double a = x_.front(), b = x_.back();
    view_ = QRectF(QPointF(std::min(a, b), -0.5), QPointF(b > a ? b : a + 1, rows_.size() - 0.5));
  }
  update();
}

bool HeatmapWidget::cellAt(const QPointF &pixel, int &column, int &row, double &value) const {
  const QRect r = framesOf(rect()).map;
  if (!r.contains(pixel.toPoint()) || x_.isEmpty() || rows_.isEmpty()) return false;
  const double x = unmapX(pixel.x(), r);
  row = (int)std::lround(unmapRow(pixel.y(), r));
  if (row < 0 || row >= rows_.size()) return false;
  // the nearest column (the cells reach half way to their neighbours)
  column = (int)(std::lower_bound(x_.begin(), x_.end(), x) - x_.begin());
  if (column >= x_.size()) column = x_.size() - 1;
  if (column > 0 && std::fabs(x_[column - 1] - x) < std::fabs(x_[column] - x)) column--;
  value = row < values_.size() && column < values_[row].size() ? values_[row][column] : NAN;
  return true;
}

QString HeatmapWidget::readoutAt(const QPointF &pixel) const {
  int c = 0, r = 0;
  double v = NAN;
  if (!cellAt(pixel, c, r, v)) return QString();
  const QString name = r < rowTips_.size() && !rowTips_[r].isEmpty() ? rowTips_[r] : rows_.value(r);
  return QString("x = %1, %2, value = %3").arg(x_[c], 0, 'g', 6).arg(name).arg(v, 0, 'g', 5);
}

void HeatmapWidget::render(QPainter &p, const QRect &area, bool onScreen) const {
  p.save();
  QColor paper = onScreen ? palette().color(QPalette::Base) : QColor(Qt::white);
  paper.setAlpha(255);
  const QColor ink = onScreen ? palette().color(QPalette::Text) : QColor(Qt::black);
  p.fillRect(area, paper);
  const Frames f = framesOf(area);
  const QRect &r = f.map;
  const QFontMetrics fm(p.font());

  // the cells: column i spans half way to its neighbours
  p.save();
  p.setClipRect(r);
  const int nx = x_.size();
  for (int row = 0; row < rows_.size() && row < values_.size(); row++) {
    const double top = mapRow(row - 0.5, r), bottom = mapRow(row + 0.5, r);
    if (bottom < r.top() || top > r.bottom()) continue;
    for (int i = 0; i < nx && i < values_[row].size(); i++) {
      const double c = colourValue(values_[row][i]);
      if (!std::isfinite(c)) continue;
      const double a = i > 0 ? 0.5 * (x_[i - 1] + x_[i]) : x_[i] - (nx > 1 ? 0.5 * (x_[1] - x_[0]) : 0.5);
      const double b = i + 1 < nx ? 0.5 * (x_[i] + x_[i + 1]) : x_[i] + (nx > 1 ? 0.5 * (x_[i] - x_[i - 1]) : 0.5);
      const double pa = mapX(a, r), pb = mapX(b, r);
      if (pb < r.left() || pa > r.right()) continue;
      p.fillRect(QRectF(QPointF(pa, top), QPointF(pb + 0.6, bottom)), colourAt((c - low_) / (high_ - low_)));
    }
  }
  // the poles: vertical lines, then points ringed so that they read on every colour
  p.setRenderHint(QPainter::Antialiasing);
  for (double x : vertical_) {
    const double px = mapX(x, r);
    p.setPen(QPen(QColor(0, 0, 0, 160), 3));
    p.drawLine(QPointF(px, r.top()), QPointF(px, r.bottom()));
    p.setPen(QPen(Qt::white, 1, Qt::DashLine));
    p.drawLine(QPointF(px, r.top()), QPointF(px, r.bottom()));
  }
  for (const QPointF &pt : points_) {
    const QPointF c(mapX(pt.x(), r), mapRow(pt.y(), r));
    p.setPen(QPen(Qt::black, 1.5));
    p.setBrush(Qt::white);
    p.drawEllipse(c, 4.5, 4.5);
  }
  p.setBrush(Qt::NoBrush);
  p.restore();

  // the frame, the x ticks, the row names
  p.setPen(ink);
  p.drawRect(r);
  for (double t : plotaxis::LinearTicks(view_.left(), view_.right(), 7)) {
    const double x = mapX(t, r);
    p.drawLine(QPointF(x, r.bottom()), QPointF(x, r.bottom() + 4));
    const QString s = plotaxis::TickText(t);
    p.drawText(QPointF(x - fm.horizontalAdvance(s) / 2.0, r.bottom() + 4 + fm.ascent()), s);
  }
  p.drawText(QRect(r.left(), r.bottom() + fm.height() + 6, r.width(), fm.height() + 4), Qt::AlignHCenter, xTitle_);
  const double rowHeight = r.height() / (view_.bottom() - view_.top());
  const int every = std::max(1, (int)std::ceil(fm.height() / std::max(1.0, rowHeight)));
  for (int row = 0; row < rows_.size(); row += every) {
    const double y = mapRow(row, r);
    if (y < r.top() - 1 || y > r.bottom() + 1) continue;
    p.drawLine(QPointF(r.left() - 4, y), QPointF(r.left(), y));
    const QString s = fm.elidedText(rows_[row], Qt::ElideMiddle, r.left() - area.left() - 10);
    p.drawText(QPointF(r.left() - 7 - fm.horizontalAdvance(s), y + fm.ascent() / 2.0 - 1), s);
  }

  // the colour bar with its numbers and title
  const QRect &bar = f.bar;
  for (int y = bar.top(); y <= bar.bottom(); y++)
    p.fillRect(QRect(bar.left(), y, bar.width(), 1), colourAt(1.0 - double(y - bar.top()) / std::max(1, bar.height())));
  p.setPen(ink);
  p.drawRect(bar);
  for (double t : plotaxis::LinearTicks(low_, high_, 5)) {
    const double y = bar.bottom() - (t - low_) / (high_ - low_) * bar.height();
    p.drawLine(QPointF(bar.right(), y), QPointF(bar.right() + 4, y));
    p.drawText(QPointF(bar.right() + 7, y + fm.ascent() / 2.0 - 1), plotaxis::TickText(t));
  }
  p.save();
  p.translate(bar.right() + kBarLabels + fm.height(), bar.center().y());
  p.rotate(-90);
  const QString title = log_ ? tr("log10 %1").arg(valueTitle_) : valueTitle_;
  p.drawText(QRect(-bar.height() / 2, -fm.height(), bar.height(), fm.height() + 4), Qt::AlignHCenter, title);
  p.restore();

  if (onScreen && !cursorText_.isEmpty()) {
    const QString shown = fm.elidedText(cursorText_, Qt::ElideMiddle, r.width() - 16);
    const QRect box(r.left() + 6, r.top() + 6, fm.horizontalAdvance(shown) + 12, fm.height() + 6);
    QColor back = paper;
    back.setAlpha(225);
    p.setPen(QColor(0x99, 0x99, 0x99));
    p.setBrush(back);
    p.drawRoundedRect(box, 3, 3);
    p.setPen(ink);
    p.drawText(box, Qt::AlignCenter, shown);
  }
  p.restore();
}

void HeatmapWidget::paintEvent(QPaintEvent *) {
  QPainter p(this);
  render(p, rect(), true);
}

bool HeatmapWidget::exportTo(const QString &file, QString &why) const {
  const QString suffix = QFileInfo(file).suffix().toLower();
  const QSize size(1000, 600);
  if (suffix == "png") {
    QImage image(size, QImage::Format_ARGB32);
    image.fill(Qt::white);
    QPainter p(&image);
    render(p, QRect(QPoint(0, 0), size));
    p.end();
    if (image.save(file)) return true;
    why = tr("cannot write %1").arg(file);
    return false;
  }
  if (suffix == "svg") {
    QSvgGenerator svg;
    svg.setFileName(file);
    svg.setSize(size);
    svg.setViewBox(QRect(QPoint(0, 0), size));
    svg.setTitle(tr("AZURE2 map"));
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
  if (suffix == "csv") {
    QFile out(file);
    if (!out.open(QIODevice::WriteOnly | QIODevice::Text)) {
      why = tr("cannot write %1").arg(file);
      return false;
    }
    out.write(csv().toUtf8());
    return true;
  }
  why = tr("%1: the map is written as .png, .svg, .pdf or .csv").arg(file);
  return false;
}

QString HeatmapWidget::csv() const {
  auto num = [](double v) { return std::isfinite(v) ? QString::number(v, 'g', 10) : QString(); };
  QString out = QString("# %1; %2; one column per row of the map\nx").arg(xTitle_, valueTitle_);
  for (int r = 0; r < rows_.size(); r++)
    out += ",\"" + QString(r < rowTips_.size() && !rowTips_[r].isEmpty() ? rowTips_[r] : rows_[r]).replace('"', '\'') + "\"";
  out += '\n';
  for (int i = 0; i < x_.size(); i++) {
    out += num(x_[i]);
    for (int r = 0; r < rows_.size(); r++) out += "," + num(r < values_.size() ? values_[r].value(i, NAN) : NAN);
    out += '\n';
  }
  return out;
}

void HeatmapWidget::wheelZoom(const QPointF &pixel, double steps, Qt::KeyboardModifiers modifiers) {
  const QRect r = framesOf(rect()).map;
  const double factor = std::pow(1.25, -steps);
  QRectF v = view_;
  if (modifiers & Qt::ShiftModifier) {
    const double c = r.contains(pixel.toPoint()) ? unmapRow(pixel.y(), r) : view_.center().y();
    v.setTop(c + (view_.top() - c) * factor);
    v.setBottom(c + (view_.bottom() - c) * factor);
  } else {
    const double c = r.contains(pixel.toPoint()) ? unmapX(pixel.x(), r) : view_.center().x();
    v.setLeft(c + (view_.left() - c) * factor);
    v.setRight(c + (view_.right() - c) * factor);
  }
  view_ = v;
  update();
}

void HeatmapWidget::panBy(const QPoint &pixels) {
  const QRect r = framesOf(rect()).map;
  const double dx = double(pixels.x()) / r.width() * view_.width(), dy = double(pixels.y()) / r.height() * view_.height();
  view_.translate(-dx, -dy);
  update();
}

void HeatmapWidget::wheelEvent(QWheelEvent *e) {
  const double steps = e->angleDelta().y() / 120.0;
  if (steps == 0) return;
  wheelZoom(e->position(), steps, e->modifiers());
  e->accept();
}

void HeatmapWidget::leaveEvent(QEvent *) {
  cursorText_.clear();
  emit cursorMoved(QString());
  update();
}

void HeatmapWidget::mousePressEvent(QMouseEvent *e) {
  moved_ = false;
  if (e->button() == Qt::RightButton || e->button() == Qt::MiddleButton) {
    panning_ = true;
    panLast_ = e->position().toPoint();
    setCursor(Qt::ClosedHandCursor);
    return;
  }
  if (e->button() != Qt::LeftButton) return;
  origin_ = e->position().toPoint();
  if (!rubber_) rubber_ = new QRubberBand(QRubberBand::Rectangle, this);
  rubber_->setGeometry(QRect(origin_, QSize()));
  rubber_->show();
}

void HeatmapWidget::mouseMoveEvent(QMouseEvent *e) {
  if (panning_) {
    const QPoint now = e->position().toPoint();
    panBy(now - panLast_);
    panLast_ = now;
  }
  if (rubber_ && rubber_->isVisible()) {
    rubber_->setGeometry(QRect(origin_, e->position().toPoint()).normalized());
    if ((e->position().toPoint() - origin_).manhattanLength() > 5) moved_ = true;
  }
  cursorText_ = readoutAt(e->position());
  emit cursorMoved(cursorText_);
  update();
}

void HeatmapWidget::mouseReleaseEvent(QMouseEvent *e) {
  if (panning_) {
    panning_ = false;
    unsetCursor();
    return;
  }
  if (!rubber_ || !rubber_->isVisible()) return;
  rubber_->hide();
  const QRect band = QRect(origin_, e->position().toPoint()).normalized();
  const QRect r = framesOf(rect()).map;
  if (!moved_ || band.width() < 6 || band.height() < 6) {
    int c = 0, row = 0;
    double v = 0;
    if (cellAt(e->position(), c, row, v)) emit rowClicked(row);
    return;
  }
  zoomTo(unmapX(band.left(), r), unmapX(band.right(), r), unmapRow(band.top(), r), unmapRow(band.bottom(), r));
}

void HeatmapWidget::mouseDoubleClickEvent(QMouseEvent *) { resetZoom(); }
