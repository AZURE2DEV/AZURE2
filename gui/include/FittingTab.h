#ifndef FITTINGTAB_H
#define FITTINGTAB_H

#include <QWidget>

#include "ProjectHolder.h"

class SectionTextPanel;

QT_BEGIN_NAMESPACE
class QComboBox;
class QLineEdit;
class QSpinBox;
class QTableWidget;
class QLabel;
QT_END_NAMESPACE

/*!
 * Fitting Settings: the fitter of a fit and its controls (the engine's table of fitters, azure2::front::Fitters: a
 * fitter gone from the engine is gone from the list), the priors (the prior covariance file of <config>; a prior on
 * one parameter is a field of its record, "norm 1 fit prior 1 0.05"), the parameter settings (<parameterSettings>)
 * and the sampler's priors (<mcmc>), the last two as text checked by the engine.
 *
 * The parameter table (docs/upgrade/qt_gui.md section 6): every parameter of a fit of the saved project as the engine
 * sets it up (azure2::front::ReadParameters: name, kind, start value in fit units, fixed or free, the fit's limits),
 * with its <parameterSettings> row (limits and nuisance prior, edited in the grammar of the engine's reader through
 * azure2::front::SetParameterSetting) and, after a fit, the values and errors of the run's param.sav beside the
 * start.  "Take the fitted values as the new start" copies param.sav beside the project and hands it to the next
 * run as --parameter-file, the engine's own way of starting from a parameter file.
 */
class FittingTab : public QWidget {
  Q_OBJECT

 public:
  explicit FittingTab(ProjectHolder *holder, QWidget *parent = nullptr);
  /// The program's switches of a run of `command` (the fitter's switch, --max-calls, --refreeze-covariance).
  QStringList programOptions(const QString &command) const;
  /// Write the text sections and the prior covariance key into the project (before a save).
  bool apply(QString &why);
  QComboBox *fitterCombo() const { return fitter_; }
  QTableWidget *parameterTable() const { return table_; }
  /// Fill the table from a session of the saved project (`projectFile`); false with `why`.
  bool readParameters(const QString &projectFile, QString &why);
  /// The values and errors of a parameter file (the run's param.sav) beside the start values; false with `why`.
  bool readFitted(const QString &file, QString &why);
  /// Copy `fitted` (a param.sav) to `<project>.start.sav` beside the project and ask for it as the next start.
  bool takeFittedAsStart(const QString &fitted, QString &why);
  /// The run's param.sav (the project's output directory).
  QString fittedFile() const;
  /// Write the limits and the prior of one row of the table as its <parameterSettings> row (an empty row removed).
  bool applyRow(int row, QString &why);

 signals:
  /// The parameter file the next run starts from (RunTab: --parameter-file).
  void startFileChosen(const QString &file);
  /// The table is asked for: the window saves the project and calls readParameters.
  void parametersRequested();

 public slots:
  void reload();

 private:
  ProjectHolder *holder_;
  QComboBox *fitter_;
  QSpinBox *maxCalls_;
  QSpinBox *refreeze_;
  QLineEdit *priorCovariance_;
  SectionTextPanel *parameterSettings_;
  SectionTextPanel *mcmc_;
  QTableWidget *table_;
  QLabel *tableNote_;
  bool filling_ = false;
};

#endif
