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
 * a few entries, optional vertical markers (e.g. the nodes of a form factor),
 * series drawn as points (e.g. quadrature nodes) and a linear or logarithmic
 * y axis.  Plots side by side can be aligned (setAlignedWith): the frames of
 * a row share their top (the tallest legend of the row) and those of a column
 * their left edge (the widest tick labels of the column), so that the axes
 * line up in a grid of equally sized panels.
 */
class ThmPlotWidget : public QWidget {
  Q_OBJECT

 public:
  struct Series {
    QVector<double> x, y;
    QColor color;
    Qt::PenStyle style = Qt::SolidLine;
    QString label;  ///< legend text; "" = not in the legend
    bool symbols = false;  ///< filled circles at the points, no line (e.g. quadrature nodes)
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
  /// A legend entry that is no curve (e.g. "dashed: window average").
  void addLegendEntry(const QString &label, Qt::PenStyle style, const QColor &color = QColor());
  /// The plots whose frame top (`row`) and left edge (`column`) this one
  /// shares; hidden ones are ignored.  Empty lists: aligned with nothing.
  void setAlignedWith(const QList<ThmPlotWidget *> &row, const QList<ThmPlotWidget *> &column);
  /// Shown in the middle instead of the curves (e.g. "not computed yet").
  void setMessage(const QString &text);

  const QList<Series> &series() const { return series_; }
  const QList<Marker> &markers() const { return markers_; }
  QString title() const { return title_; }
  bool logY() const { return logY_; }

  QSize sizeHint() const override { return QSize(400, 260); }
  QSize minimumSizeHint() const override { return QSize(240, 180); }
  /// The frame (the axes box) as painted, in widget coordinates.
  QRect frameRect() const;

  /// "Nice" tick positions covering [lo, hi] (about `target` of them).
  static QVector<double> linearTicks(double lo, double hi, int target = 5);

 protected:
  void paintEvent(QPaintEvent *event) override;

 private:
  struct Axes {
    bool empty = true;
    double xlo = 0, xhi = 1, ylo = 0, yhi = 1;  ///< y in log10 for a log axis
  };
  struct Tick {
    double at;     ///< axis units (log10 for a log axis)
    QString text;  ///< "" = unlabelled
  };
  QFont tickFont() const;
  Axes axes() const;
  QList<Tick> yTicks(const Axes &a, int frameHeight) const;
  QList<Series> legendSeries() const;
  /// The legend's rows, each a list of (series, width), wrapped to the widget width.
  QList<QList<QPair<Series, int>>> legendRows() const;
  int ownTop() const;
  int top() const;
  int ownLeft() const;
  int left() const;
  int bottomMargin() const;
  void changed();

  QList<ThmPlotWidget *> rowPeers_, columnPeers_;
  QList<Series> legendOnly_;
  QList<Series> series_;
  QList<Marker> markers_;
  QString title_, xLabel_, yLabel_, message_;
  bool logY_ = false;
};

/// Okabe-Ito colours, for curves that must stay distinguishable.
QColor thmPlotColor(int i);

#endif
