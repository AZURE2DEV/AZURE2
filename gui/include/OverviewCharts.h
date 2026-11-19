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

 protected:
  void paintEvent(QPaintEvent *) override;

 private:
  bool horizontal_;
  QStringList labels_;
  QVector<double> values_;
  QVector<bool> highlight_;
  double reference_;
  QString referenceText_, empty_;
};

#endif
