#ifndef HEATMAPWIDGET_H
#define HEATMAPWIDGET_H

#include <QColor>
#include <QPointF>
#include <QRectF>
#include <QString>
#include <QStringList>
#include <QVector>
#include <QWidget>

class QRubberBand;

/*!
 * A map of values over (x, row): the Physics tab's pole map (docs/upgrade/qt_gui.md section 15.5).  x is a number
 * axis (the excitation energy), the rows are named lines (J^pi groups, channel pairs, level pairs), the colour is the
 * value on a colour map of our own (colourAt: a viridis-like ramp, dark blue to yellow) with a colour bar and its
 * numbers.  The ticks are the plot's (PlotAxis.h).  Zoom: the wheel about the cursor (x; Shift: rows), a rectangle
 * dragged with the left button; pan with the right or middle button; double click: the whole map.  The cursor's
 * read-out (x, row, value) is drawn in the corner and sent as cursorMoved; a click on a row sends rowClicked.  Points
 * (x, row) and vertical lines at x are drawn over the map (the poles).  Drawn in the window's palette on screen, on
 * white in a file (PNG, SVG, PDF; CSV the values).
 */
class HeatmapWidget : public QWidget {
  Q_OBJECT

 public:
  explicit HeatmapWidget(QWidget *parent = nullptr);
  /// values[row][column] at x[column]; NaN cells are left empty.
  void setMap(const QVector<double> &x, const QStringList &rows, const QVector<QVector<double>> &values,
              const QString &xTitle, const QString &valueTitle);
  /// Longer names of the rows, for the read-out (the short ones label the axis).
  void setRowTips(const QStringList &tips) { rowTips_ = tips; }
  /// Colour by log10 of the value (cells with values <= 0 left empty).
  void setLogColour(bool on);
  void setPoints(const QVector<QPointF> &points) { points_ = points, update(); }  ///< (x, row index)
  void setVerticalMarkers(const QVector<double> &x) { vertical_ = x, update(); }
  const QVector<double> &x() const { return x_; }
  QString valueTitle() const { return valueTitle_; }
  const QStringList &rows() const { return rows_; }
  const QVector<QVector<double>> &values() const { return values_; }
  const QVector<QPointF> &points() const { return points_; }
  /// The colour range (in log10 when log colour is on).
  double colourLow() const { return low_; }
  double colourHigh() const { return high_; }
  /// The cell under a pixel: its column, row and value (false outside the map).
  bool cellAt(const QPointF &pixel, int &column, int &row, double &value) const;
  /// The read-out at a pixel as drawn ("x = 2.365, 1/2+, value = -3.21"), empty outside.
  QString readoutAt(const QPointF &pixel) const;
  QString cursorText() const { return cursorText_; }
  /// The rectangle of the map's cells inside the widget (for the tests: where a pixel lies).
  QRect mapRect() const { return framesOf(rect()).map; }
  void zoomTo(double x0, double x1, double row0, double row1);
  void resetZoom();
  QRectF view() const { return view_; }
  void wheelZoom(const QPointF &pixel, double steps, Qt::KeyboardModifiers modifiers = Qt::NoModifier);
  void panBy(const QPoint &pixels);
  void render(QPainter &p, const QRect &area, bool onScreen = false) const;
  bool exportTo(const QString &file, QString &why) const;
  QString csv() const;
  /// The colour map at t in [0, 1].
  static QColor colourAt(double t);
  QSize sizeHint() const override { return QSize(700, 360); }

 signals:
  void cursorMoved(const QString &readout);
  void rowClicked(int row);

 protected:
  void paintEvent(QPaintEvent *) override;
  void mousePressEvent(QMouseEvent *) override;
  void mouseMoveEvent(QMouseEvent *) override;
  void mouseReleaseEvent(QMouseEvent *) override;
  void mouseDoubleClickEvent(QMouseEvent *) override;
  void wheelEvent(QWheelEvent *) override;
  void leaveEvent(QEvent *) override;

 private:
  struct Frames {
    QRect map, bar;
  };
  Frames framesOf(const QRect &area) const;
  int leftMargin(const QFontMetrics &fm) const;
  double mapX(double x, const QRect &r) const;
  double mapRow(double row, const QRect &r) const;
  double unmapX(double px, const QRect &r) const;
  double unmapRow(double py, const QRect &r) const;
  double colourValue(double v) const;  ///< the value on the colour scale (log10 or itself); NaN when not drawn
  void rangeOfColours();

  QVector<double> x_;
  QStringList rows_, rowTips_;
  QVector<QVector<double>> values_;
  QString xTitle_, valueTitle_;
  bool log_ = false;
  double low_ = 0, high_ = 1;
  QVector<QPointF> points_;
  QVector<double> vertical_;
  QRectF view_;  ///< x from left to right; rows from top (the first row) to bottom, row i over [i - 0.5, i + 0.5]
  QRubberBand *rubber_ = nullptr;
  QPoint origin_;
  bool panning_ = false, moved_ = false;
  QPoint panLast_;
  QString cursorText_;
};

#endif
