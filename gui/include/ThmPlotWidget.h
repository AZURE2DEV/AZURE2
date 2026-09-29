#ifndef THMPLOTWIDGET_H
#define THMPLOTWIDGET_H

#include <QColor>
#include <QList>
#include <QString>
#include <QVector>
#include <QWidget>

/*!
 * A small line plot drawn with QPainter, for the THM workspace's Diagnostics
 * page.  It needs no plotting library, so the page works in every build (the
 * Plot tab's Qwt is optional).  Deliberately plain: one frame, inward ticks,
 * large tick labels, a bold label inside the axes (the reaction), a legend of
 * a few entries, optional vertical markers (e.g. the nodes of a form factor)
 * and a linear or logarithmic y axis.
 */
class ThmPlotWidget : public QWidget {
  Q_OBJECT

 public:
  struct Series {
    QVector<double> x, y;
    QColor color;
    Qt::PenStyle style = Qt::SolidLine;
    QString label;  ///< legend text; "" = not in the legend
  };
  struct Marker {
    double x = 0.0;
    QColor color;
    QString label;  ///< drawn at the top of the line
  };

  explicit ThmPlotWidget(QWidget *parent = 0);

  void clear();
  /// Rich text (HTML subset of QTextDocument), drawn bold in the upper left.
  void setTitle(const QString &html);
  void setAxisLabels(const QString &x, const QString &y);
  void setLogY(bool log);
  void addSeries(const Series &s);
  void addMarker(const Marker &m);
  /// Shown in the middle instead of the curves (e.g. "not computed yet").
  void setMessage(const QString &text);

  const QList<Series> &series() const { return series_; }
  const QList<Marker> &markers() const { return markers_; }
  QString title() const { return title_; }
  bool logY() const { return logY_; }

  QSize sizeHint() const override { return QSize(420, 300); }
  QSize minimumSizeHint() const override { return QSize(260, 200); }

  /// "Nice" tick positions covering [lo, hi] (about `target` of them).
  static QVector<double> linearTicks(double lo, double hi, int target = 5);

 protected:
  void paintEvent(QPaintEvent *event) override;

 private:
  QList<Series> series_;
  QList<Marker> markers_;
  QString title_, xLabel_, yLabel_, message_;
  bool logY_ = false;
};

/// Okabe-Ito colours, for curves that must stay distinguishable.
QColor thmPlotColor(int i);

#endif
