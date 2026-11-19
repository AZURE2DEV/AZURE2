#ifndef PLOTTAB_H
#define PLOTTAB_H

#include <QWidget>

#include "ProjectHolder.h"

class PlotWidget;

QT_BEGIN_NAMESPACE
class QCheckBox;
class QComboBox;
class QLabel;
class QListWidget;
class QPushButton;
class QRadioButton;
QT_END_NAMESPACE

/*!
 * The Plot tab (docs/upgrade/qt_gui.md section 5): the data sets and extrapolations of the last run, as the facade
 * reads them from the run's output files (azure2::front::ReadCurves: the files matched to the records by the engine's
 * build of the project; nothing is parsed here).  Choose data sets and extrapolations (several overlaid, each with
 * its colour and marker in the legend); x the c.m. energy, the excitation energy or the c.m. angle; y the cross
 * section (the observable as written) or the S factor where it is one; linear or log axes; the uncertainty band where
 * the run wrote one (--covariance-band); the residuals (data - calc)/error in a panel below; zoom by a rectangle;
 * export to PNG, SVG or PDF.  The window refreshes it after every run.
 */
class PlotTab : public QWidget {
  Q_OBJECT

 public:
  explicit PlotTab(ProjectHolder *holder, QWidget *parent = nullptr);
  PlotWidget *plot() const { return plot_; }
  /// Select the curves of these records (data sets or extrapolations) and draw them.
  void selectRecords(const QStringList &names);
  /// The records the last run's files hold, in the lists' order.
  QStringList availableRecords() const;
  QString notes() const;
  bool exportPlot(const QString &file, QString &why) const;

 public slots:
  /// Read the last run's curves again (after a run, after opening a project) and redraw.
  void refresh();
  void draw();

 private slots:
  void exportDialog();

 private:
  ProjectHolder *holder_;
  azure2::front::Curves curves_;
  QListWidget *dataList_, *extrapList_;
  QComboBox *xAxis_;
  QRadioButton *ySigma_, *yS_;
  QCheckBox *logX_, *logY_, *band_, *residuals_, *legend_;
  QLabel *notes_;
  PlotWidget *plot_;
};

#endif
