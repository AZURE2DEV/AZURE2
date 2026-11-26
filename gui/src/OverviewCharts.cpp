// The Overview page's bar charts (OverviewCharts.h).
#include "OverviewCharts.h"

#include <QPainter>
#include <QPainterPath>
#include <algorithm>
#include <cmath>
#include <limits>

#include "AzureTheme.h"

BarChart::BarChart(bool horizontal, QWidget *parent)
    : QWidget(parent), horizontal_(horizontal), reference_(std::numeric_limits<double>::quiet_NaN()) {
  setSizePolicy(QSizePolicy::Expanding, horizontal ? QSizePolicy::Fixed : QSizePolicy::Expanding);
}

QColor BarChart::warningColor() { return azureThemeIsDark() ? QColor("#f08a7e") : QColor("#c0392b"); }

double BarChart::ScaleTop(const QVector<double> &values, bool *capped) {
  QVector<double> positive;
  double top = 0;
  for (double v : values)
    if (std::isfinite(v) && v > 0) {
      positive << v;
      top = std::max(top, v);
    }
  if (capped) *capped = false;
  if (positive.size() < 3) return top;
  std::sort(positive.begin(), positive.end());
  const double median = positive[positive.size() / 2];
  if (median > 0 && top > 10 * median) {
    if (capped) *capped = true;
    return 10 * median;
  }
  return top;
}

void BarChart::setBars(const QStringList &labels, const QVector<double> &values, const QVector<bool> &highlight) {
  labels_ = labels;
  values_ = values;
  highlight_ = highlight;
  bool cut = false;
  ScaleTop(values, &cut);
  capped_ = cut && !horizontal_;
  // each bar's name and value in the tip (the bars are numbered when their names do not fit under them)
  QStringList tip;
  for (int i = 0; i < values.size() && i < 40; i++)
    tip << QString("%1: %2  %3").arg(i + 1).arg(labels.value(i)).arg(QString::number(values[i], 'g', 4));
  if (capped_)
    tip.prepend(tr("The scale is cut at ten times the median: a cut bar (marked) is longer than the axis."));
  setToolTip(tip.join("\n"));
  updateGeometry();
  update();
}

void BarChart::setReference(double value, const QString &text) {
  reference_ = value;
  referenceText_ = text;
  update();
}

void BarChart::setEmptyText(const QString &text) {
  empty_ = text;
  update();
}

QSize BarChart::sizeHint() const {
  const int line = fontMetrics().height();
  if (horizontal_) return QSize(16 * line, std::max(2, (int)values_.size()) * (line + 4) + 4);
  return QSize(16 * line, 9 * line);
}

QSize BarChart::minimumSizeHint() const {
  const int line = fontMetrics().height();
  return horizontal_ ? QSize(10 * line, sizeHint().height()) : QSize(10 * line, 6 * line);
}

void BarChart::paintEvent(QPaintEvent *) {
  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing);
  const QPalette &pal = palette();
  const QFontMetrics fm = fontMetrics();
  const int line = fm.height();
  if (values_.isEmpty()) {
    p.setPen(pal.color(QPalette::PlaceholderText));
    p.drawText(rect(), Qt::AlignCenter | Qt::TextWordWrap, empty_);
    return;
  }
  double top = horizontal_ ? 0 : ScaleTop(values_);
  if (horizontal_)
    for (double v : values_)
      if (std::isfinite(v)) top = std::max(top, v);
  if (std::isfinite(reference_)) top = std::max(top, reference_);
  if (top <= 0) top = 1;
  const QColor bar = pal.color(QPalette::Highlight);
  const QColor warn = warningColor();

  if (horizontal_) {
    int labelWidth = 0;
    for (const QString &l : labels_) labelWidth = std::max(labelWidth, fm.horizontalAdvance(l));
    const int countWidth = fm.horizontalAdvance("0000");
    const int x0 = labelWidth + line / 2;
    const int span = std::max(10, width() - x0 - countWidth - line / 2);
    for (int i = 0; i < values_.size(); i++) {
      const int y = 2 + i * (line + 4);
      p.setPen(pal.color(QPalette::WindowText));
      p.drawText(QRect(0, y, labelWidth, line), Qt::AlignRight | Qt::AlignVCenter, labels_.value(i));
      const int w = std::max(2, (int)std::lround(span * values_[i] / top));
      p.setPen(Qt::NoPen);
      p.setBrush(highlight_.value(i) ? warn : bar);
      p.drawRoundedRect(QRectF(x0, y + 2, w, line - 4), 2, 2);
      p.setPen(pal.color(QPalette::WindowText));
      p.drawText(QRect(x0 + w + 4, y, countWidth, line), Qt::AlignLeft | Qt::AlignVCenter,
                 QString::number(values_[i], 'g', 6));
    }
    return;
  }

  // vertical: an axis on the left with 0 and the top, a bar per item, its label under it when there is room
  const QString topText = QString::number(top, 'g', 3);
  const int left = fm.horizontalAdvance(topText) + line / 2;
  const QRect area(left, line / 2, std::max(10, width() - left - 2), std::max(10, height() - line / 2 - line - 4));
  p.setPen(pal.color(QPalette::Mid));
  p.drawLine(area.bottomLeft(), area.bottomRight());
  p.drawLine(area.bottomLeft(), area.topLeft());
  p.setPen(pal.color(QPalette::PlaceholderText));
  p.drawText(QRect(0, area.top() - line / 2, left - 3, line), Qt::AlignRight | Qt::AlignVCenter, topText);
  p.drawText(QRect(0, area.bottom() - line / 2, left - 3, line), Qt::AlignRight | Qt::AlignVCenter, "0");
  const int n = values_.size();
  const double slot = (double)area.width() / n;
  const double w = std::max(2.0, slot * 0.7);
  bool labels = true;
  for (const QString &l : labels_) labels = labels && fm.horizontalAdvance(l) < slot - 2;
  for (int i = 0; i < n; i++) {
    const double v = std::isfinite(values_[i]) ? values_[i] : 0;
    const bool cut = v > top;
    const double h = area.height() * std::min(v, top) / top;  // never above the axis
    const double x = area.left() + slot * i + (slot - w) / 2;
    p.setPen(Qt::NoPen);
    p.setBrush(highlight_.value(i) ? warn : bar);
    p.drawRect(QRectF(x, area.bottom() - h, w, h));
    if (cut) {  // the bar is longer than the scale: two breaks near its top
      p.setPen(QPen(pal.color(QPalette::Base), 2));
      for (int k = 0; k < 2; k++) {
        const double y = area.top() + 6 + 5 * k;
        p.drawLine(QPointF(x - 1, y + 3), QPointF(x + w + 1, y - 3));
      }
    }
    p.setPen(pal.color(QPalette::PlaceholderText));
    const QString label = labels ? labels_.value(i) : QString::number(i + 1);
    if (labels || slot > fm.horizontalAdvance(label) + 1)
      p.drawText(QRectF(area.left() + slot * i, area.bottom() + 2, slot, line), Qt::AlignHCenter | Qt::AlignTop, label);
  }
  if (std::isfinite(reference_)) {
    const double y = area.bottom() - area.height() * reference_ / top;
    QPen dash(pal.color(QPalette::WindowText));
    dash.setStyle(Qt::DashLine);
    p.setPen(dash);
    p.drawLine(QPointF(area.left(), y), QPointF(area.right(), y));
    // its text left of the axis, where no bar is
    if (area.bottom() - y > line && y - area.top() > line)  // not over the axis's own numbers
      p.drawText(QRectF(0, y - line / 2.0, left - 3, line), Qt::AlignRight | Qt::AlignVCenter, referenceText_);
  }
}
