#ifndef FITTINGTAB_H
#define FITTINGTAB_H

#include <QWidget>

#include "ProjectHolder.h"

class SectionTextPanel;

QT_BEGIN_NAMESPACE
class QComboBox;
class QLineEdit;
class QSpinBox;
QT_END_NAMESPACE

/*!
 * Fitting Settings: the fitter of a fit and its controls (the engine's table of fitters, azure2::front::Fitters: a
 * fitter gone from the engine is gone from the list), the priors (the prior covariance file of <config>; a prior on
 * one parameter is a field of its record, "norm 1 fit prior 1 0.05"), the parameter settings (<parameterSettings>)
 * and the sampler's priors (<mcmc>), the last two as text checked by the engine.
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
};

#endif
