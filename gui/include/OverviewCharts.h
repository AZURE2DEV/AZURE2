#ifndef OVERVIEWCHARTS_H
#define OVERVIEWCHARTS_H

#include <QColor>
#include <QStringList>
#include <QVector>
#include <QWidget>

/*!
 * The Overview page's small charts (docs/upgrade/qt_gui.md section 15), drawn with QPainter in the window's palette.
 *
 * BarChart: one bar per item, vertical (chi^2/N per data set, a dashed line at a reference value) or horizontal
 * (levels per J^pi group, the count at the bar's end).  The numbers are the caller's; nothing is computed here but
 * the bars' lengths.
 */
class BarChart : public QWidget {
  Q_OBJECT

 public:
  explicit BarChart(bool horizontal, QWidget *parent = nullptr);
  /// The bars: a label and a value each; `highlight` marks bars drawn in the warning colour.
  void setBars(const QStringList &labels, const QVector<double> &values, const QVector<bool> &highlight = {});
  /// A dashed line at `value` (vertical charts; NaN: none) with its text.
  void setReference(double value, const QString &text);
  void setEmptyText(const QString &text);
  int bars() const { return (int)values_.size(); }
  QSize sizeHint() const override;
  QSize minimumSizeHint() const override;

  /// The colour of a value to look at, readable on the light and the dark palette.
  static QColor warningColor();
  /// The top of the bars' scale: the largest finite value, or, when it is more than ten times the median of the
  /// positive values (one outlier would flatten the others), ten times the median; `capped` says which.  Bars above
  /// it are drawn full and cut (qt_gui.md section 19.9).
  static double ScaleTop(const QVector<double> &values, bool *capped = nullptr);
  /// Whether the last drawing cut the bars at a capped scale (the tests).
  bool capped() const { return capped_; }

 protected:
  void paintEvent(QPaintEvent *) override;

 private:
  bool horizontal_;
  QStringList labels_;
  QVector<double> values_;
  QVector<bool> highlight_;
  double reference_;
  QString referenceText_, empty_;
  bool capped_ = false;
};

#endif
