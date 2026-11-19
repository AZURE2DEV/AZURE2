#ifndef PLOTWIDGET_H
#define PLOTWIDGET_H

#include <QColor>
#include <QList>
#include <QRectF>
#include <QString>
#include <QVector>
#include <QWidget>

QT_BEGIN_NAMESPACE
class QPainter;
class QRubberBand;
QT_END_NAMESPACE

/// One series of the plot: data points with error bars, a calculated curve, its uncertainty band; all in the units of
/// the axes the tab chose.  An empty vector leaves that part out.
struct PlotSeries {
  QString label;
  QColor color;
  int marker = 0;  ///< 0 circle, 1 square, 2 diamond, 3 triangle, 4 inverted triangle
  QVector<double> dataX, dataY, dataError;
  QVector<double> modelX, modelY;
  QVector<double> bandLow, bandHigh;  ///< beside modelX, the model -+ 1 sigma
  QVector<double> residualX, residual;  ///< (data - model) / error, for the lower panel
};

/*!
 * The Plot tab's drawing (docs/upgrade/qt_gui.md section 5): a plot drawn with QPainter, so that every GUI build has
 * it without a plotting library.  Linear or logarithmic axes with decade ticks, data with error bars, the calculated
 * curve (a single point drawn as a marker), the band as a shaded area, a legend, and a panel of normalised residuals
 * below.  Zoom: drag a rectangle with the left button; double-click or resetZoom() returns to the whole range.
 * render() draws into any QPainter: the widget, a PNG, an SVG (QSvgGenerator) or a PDF (QPdfWriter).
 */
class PlotWidget : public QWidget {
  Q_OBJECT

 public:
  explicit PlotWidget(QWidget *parent = nullptr);
  void setSeries(const QList<PlotSeries> &series);
  const QList<PlotSeries> &series() const { return series_; }
  void setAxisTitles(const QString &x, const QString &y);
  void setLogX(bool on);
  void setLogY(bool on);
  void setShowBand(bool on);
  void setShowResiduals(bool on);
  void setShowLegend(bool on);
  /// Zoom to a rectangle in data coordinates; resetZoom() to the data's range.
  void zoomTo(const QRectF &data);
  void resetZoom();
  QRectF view() const { return view_; }
  /// The number of data points drawn (finite and, on a log axis, positive), for the tests.
  int dataPointsDrawn() const;
  /// Draw the whole plot into `p` over `area`.
  void render(QPainter &p, const QRect &area) const;
  /// Write the plot to a file: .png, .svg or .pdf by its extension.  False with `why`.
  bool exportTo(const QString &file, QString &why) const;
  QSize sizeHint() const override { return QSize(700, 480); }

 protected:
  void paintEvent(QPaintEvent *) override;
  void mousePressEvent(QMouseEvent *) override;
  void mouseMoveEvent(QMouseEvent *) override;
  void mouseReleaseEvent(QMouseEvent *) override;
  void mouseDoubleClickEvent(QMouseEvent *) override;

 private:
  struct Frame {
    QRect main, lower;  ///< the main panel; the residual panel (empty when off)
  };
  Frame framesOf(const QRect &area) const;
  QRectF dataRange() const;
  double mapX(double x, const QRect &r) const;
  double mapY(double y, const QRect &r, const QRectF &view, bool log) const;
  double unmapX(double px, const QRect &r) const;
  double unmapY(double py, const QRect &r) const;
  bool usable(double x, double y) const;

  QList<PlotSeries> series_;
  QString xTitle_, yTitle_;
  bool logX_ = false, logY_ = false, band_ = true, residuals_ = true, legend_ = true;
  QRectF view_;
  QRubberBand *rubber_ = nullptr;
  QPoint origin_;
};

#endif
