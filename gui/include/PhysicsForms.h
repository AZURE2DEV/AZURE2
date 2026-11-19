#ifndef PHYSICSFORMS_H
#define PHYSICSFORMS_H

#include <vector>

#include <QColor>
#include <QList>
#include <QString>
#include <QVector>
#include <QWidget>

#include "PlotWidget.h"
#include "azure2/front/Engine.h"

class QCheckBox;
class QComboBox;
class QFormLayout;
class QLabel;
class QLineEdit;
class QSpinBox;

/// The Physics tab's pieces shared by its pages (docs/upgrade/qt_gui.md section 15).
namespace physicsqt {

/// The colour of the n-th curve (the Plot tab's ten colours).
QColor CurveColour(int n);
/// n values from a to b, evenly spaced or (log) evenly in log10 (a, b > 0).
QVector<double> Grid(double a, double b, int n, bool log);
/// A number as a field shows it.
QString Number(double v);

/*!
 * The fields of a sweep: from, to, the number of points and log spacing.  Read by the page that owns it.
 */
class SweepFields : public QWidget {
  Q_OBJECT
 public:
  explicit SweepFields(const QString &unit, QWidget *parent = nullptr);
  void set(double from, double to, int points, bool log);
  void setUnit(const QString &unit);
  /// The grid; false with `why` when a field does not read or log spacing meets a value <= 0.
  bool grid(QVector<double> &out, QString &why) const;
  QLineEdit *from, *to;
  QSpinBox *points;
  QCheckBox *log;

 signals:
  void changed();

 private:
  QLabel *unit_;
};

/*!
 * A form made from a row of the physics or stopping table: one field per argument (a real a number field, an integer
 * a spin box with the row's range, a flag a check box, a choice a drop-down, a text a line), its start value the
 * default, else the example; its unit beside it, its doc and domain as the tool tip.  The swept argument is chosen in
 * a drop-down of the arguments that may be an axis (`hasSweep`): its field is disabled and the sweep's fields (from,
 * to, points, log, started at the row's sweep range) stand below.
 */
class ArgumentForm : public QWidget {
  Q_OBJECT
 public:
  explicit ArgumentForm(QWidget *parent = nullptr);
  void setFunction(const azure2::front::PhysicsFunction &function);
  const azure2::front::PhysicsFunction &function() const { return function_; }
  int count() const { return (int)function_.arguments.size(); }
  /// Argument k's value (a choice its index, a flag 0 or 1); NaN when its field does not read.
  double value(int k) const;
  QString text(int k) const;  ///< a "text" argument's value
  void setValue(int k, double v);
  void setText(int k, const QString &text);
  QWidget *editor(int k) const { return editors_.value(k); }
  /// The swept argument (its index in the row), -1 when the row has none.
  int swept() const;
  void setSwept(int k);
  SweepFields *sweep() const { return sweep_; }
  /// The columns EvaluatePhysics takes: the swept argument's grid, every other a column of one value.  The numeric
  /// arguments only when `skipText` (EvaluateStopping).  False with `why` when a field does not read.
  bool columns(std::vector<std::vector<double>> &out, QVector<double> &grid, QString &why, bool skipText) const;
  /// The fixed arguments in words ("l = 1, a_c = 5 fm"): a curve's label.
  QString fixedText() const;
  /// The same as a list ("l = 1", ...), one entry per fixed argument in the row's order.
  QStringList fixedParts() const;

 signals:
  void changed();

 private:
  void sweptChanged();
  azure2::front::PhysicsFunction function_;
  QFormLayout *fields_;
  QList<QWidget *> editors_;
  QComboBox *sweepChoice_;
  SweepFields *sweep_;
};

/*!
 * The plot of a page with its strip: log x, log y, legend, Reset zoom, Export (PNG, SVG, PDF, CSV), the cursor's
 * read-out; a status line under the plot (what was left out and why) and, for the stopping page, a citation line.
 */
class PlotPanel : public QWidget {
  Q_OBJECT
 public:
  explicit PlotPanel(const QString &name, QWidget *parent = nullptr);
  PlotWidget *plot() const { return plot_; }
  /// Draw these curves with these axis titles; logY suggested (only when the user has not chosen).
  void show(const QList<PlotSeries> &curves, const QString &xTitle, const QString &yTitle);
  void setStatus(const QString &text);
  QString status() const;
  void setCitation(const QString &text);
  QCheckBox *logX() const { return logX_; }
  QCheckBox *logY() const { return logY_; }
  QString cursor() const;

 private:
  PlotWidget *plot_;
  QCheckBox *logX_, *logY_, *legend_;
  QLabel *cursor_, *status_, *citation_;
};

}  // namespace physicsqt

#endif
