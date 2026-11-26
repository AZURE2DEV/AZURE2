#ifndef PLOTWIDGET_H
#define PLOTWIDGET_H

#include <QColor>
#include <QList>
#include <QRectF>
#include <QString>
#include <QStringList>
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
 * below.  Zoom: drag a rectangle with the left button, or turn the wheel about the cursor (Ctrl: x only, Shift: y
 * only); pan: drag with the right or middle button, or Shift + left; double-click or resetZoom() returns to the whole
 * range.  The cursor's coordinates are drawn in the panel's corner and sent as cursorMoved.  render() draws into any
 * QPainter: the widget (in the window's palette), a PNG, an SVG (QSvgGenerator) or a PDF (QPdfWriter) (on white);
 * exportTo() also writes the points drawn as CSV.
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
  /// Vertical dashed lines at `x` (data units) with a label each, drawn over the curves (the Physics tab's poles).
  void setMarkers(const QVector<double> &x, const QStringList &labels = QStringList());
  const QVector<double> &markers() const { return markersX_; }
  /// Zoom to a rectangle in data coordinates; resetZoom() to the data's range.
  void zoomTo(const QRectF &data);
  void resetZoom();
  QRectF view() const { return view_; }
  /// The number of data points drawn (finite and, on a log axis, positive), for the tests.
  int dataPointsDrawn() const;
  /// Draw the whole plot into `p` over `area`: on white (a file) or, `onScreen`, in the widget's palette.
  void render(QPainter &p, const QRect &area, bool onScreen = false) const;
  /// Where the legend is drawn in the widget (empty without one), and whether it is below the axes (the tests).
  QRect legendRect() const { return framesOf(rect()).legend; }
  bool legendBelowAxes() const;
  /// The pole labels drawn at the last paint, as rectangles (the tests: none overlaps another).
  QList<QRect> markerLabelRects() const;
  /// Write the plot to a file: .png, .svg, .pdf or .csv (the points drawn) by its extension.  False with `why`.
  bool exportTo(const QString &file, QString &why) const;
  /// The points of the series as CSV: per series a block "# <label>: data" (x, y, error) and "# <label>: calculated"
  /// (x, y[, band low, band high]).
  QString csv() const;
  /// The data coordinates under a pixel of the main panel (false outside it).
  bool dataAt(const QPointF &pixel, QPointF &data) const;
  /// Zoom about a pixel by wheel steps (> 0 in): what a wheel turn does.  Ctrl: x only, Shift: y only.
  void wheelZoom(const QPointF &pixel, double steps, Qt::KeyboardModifiers modifiers = Qt::NoModifier);
  /// Pan by a drag of `pixels` (the content follows the mouse).
  void panBy(const QPoint &pixels);
  /// The cursor's read-out as drawn ("x = 0.52, y = 0.031"), empty outside the panel.
  QString cursorText() const { return cursorText_; }
  QSize sizeHint() const override { return QSize(700, 480); }

 signals:
  void cursorMoved(double x, double y, bool inside);

 protected:
  void paintEvent(QPaintEvent *) override;
  void mousePressEvent(QMouseEvent *) override;
  void mouseMoveEvent(QMouseEvent *) override;
  void mouseReleaseEvent(QMouseEvent *) override;
  void mouseDoubleClickEvent(QMouseEvent *) override;
  void wheelEvent(QWheelEvent *) override;
  void leaveEvent(QEvent *) override;

 private:
  struct Frame {
    QRect main, lower;  ///< the main panel; the residual panel (empty when off)
    QRect legend;       ///< the legend's box (empty without one): a corner of the main panel that covers no point, else
                        ///< below the axes
    int legendColumns = 1;
  };
  Frame framesOf(const QRect &area) const;
  Frame framesOfInner(const QRect &area) const;  ///< the panels alone, without the legend's room
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
  bool panning_ = false;
  QPoint panLast_;
  QString cursorText_;
  QVector<double> markersX_;
  QStringList markerLabels_;
};

#endif
