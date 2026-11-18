#include "FittingTab.h"

#include <QComboBox>
#include <QDir>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QTabWidget>
#include <QToolButton>
#include <QVBoxLayout>

#include "SectionTextPanel.h"

FittingTab::FittingTab(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  fitter_ = new QComboBox;
  fitter_->setObjectName("fitterCombo");
  for (const azure2::front::Fitter &f : azure2::front::Fitters()) {
    fitter_->addItem(QString::fromStdString(f.name), QString::fromStdString(f.option));
    fitter_->setItemData(fitter_->count() - 1, QString::fromStdString(f.doc), Qt::ToolTipRole);
  }
  maxCalls_ = new QSpinBox;
  maxCalls_->setRange(0, 1000000000);
  maxCalls_->setSpecialValueText(tr("none"));
  maxCalls_->setToolTip(tr("The fit's call limit (--max-calls); a fit stopped by it ends at the lowest objective seen."));
  refreeze_ = new QSpinBox;
  refreeze_->setRange(0, 1000);
  refreeze_->setSpecialValueText(tr("none"));
  refreeze_->setToolTip(tr("Refit up to N times with relvector covariance terms re-frozen at the fitted parameters "
                           "(--refreeze-covariance)."));
  QFormLayout *fitForm = new QFormLayout;
  fitForm->addRow(tr("Minimizer:"), fitter_);
  fitForm->addRow(tr("Call limit:"), maxCalls_);
  fitForm->addRow(tr("Covariance refreeze passes:"), refreeze_);
  QGroupBox *fitBox = new QGroupBox(tr("Fit"));
  fitBox->setLayout(fitForm);

  priorCovariance_ = new QLineEdit;
  priorCovariance_->setObjectName("priorCovariance");
  priorCovariance_->setPlaceholderText(tr("none"));
  priorCovariance_->setToolTip(tr("The prior parameter covariance file, relative to the project (<config> "
                                  "priorCovariance): a correlated Gaussian prior on named parameters."));
  QToolButton *browse = new QToolButton;
  browse->setText("...");
  connect(browse, &QToolButton::clicked, this, [this]() {
    const QString f = QFileDialog::getOpenFileName(this, tr("Prior covariance file"), holder_->directory());
    if (!f.isEmpty()) priorCovariance_->setText(QDir(holder_->directory()).relativeFilePath(f));
  });
  QHBoxLayout *priorRow = new QHBoxLayout;
  priorRow->addWidget(priorCovariance_, 1);
  priorRow->addWidget(browse);
  QLabel *priorNote = new QLabel(tr("A Gaussian prior on one parameter is written on the parameter's own field: "
                                    "\"norm 1 fit prior 1 0.05\" on a data set, \"areal_density ... fit prior ...\" on a "
                                    "target, an energy scale's or a normalization record's value.  Nuisance parameters "
                                    "are the fields written \"fit\"."));
  priorNote->setWordWrap(true);
  QVBoxLayout *priorLayout = new QVBoxLayout;
  QFormLayout *priorForm = new QFormLayout;
  priorForm->addRow(tr("Prior covariance file:"), priorRow);
  priorLayout->addLayout(priorForm);
  priorLayout->addWidget(priorNote);
  QGroupBox *priorBox = new QGroupBox(tr("Priors"));
  priorBox->setLayout(priorLayout);

  parameterSettings_ = new SectionTextPanel(
      holder, "parameterSettings",
      tr("Parameter settings (<parameterSettings>): per parameter, by its name, the limits and whether it is fixed, "
         "overriding the levels' and records' flags."));
  mcmc_ = new SectionTextPanel(holder, "mcmc",
                               tr("The sampler's priors and settings (<mcmc>), read by the command mcmc where the build "
                                  "has the sampler."));
  QTabWidget *texts = new QTabWidget;
  texts->addTab(parameterSettings_, tr("Parameter Settings"));
  texts->addTab(mcmc_, tr("Sampler (MCMC)"));

  QHBoxLayout *top = new QHBoxLayout;
  top->addWidget(fitBox);
  top->addWidget(priorBox, 1);
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addLayout(top);
  layout->addWidget(texts, 1);
  connect(holder_, &ProjectHolder::opened, this, &FittingTab::reload);
  reload();
}

void FittingTab::reload() { priorCovariance_->setText(holder_->configValue("priorCovariance")); }

QStringList FittingTab::programOptions(const QString &command) const {
  QStringList out;
  if (command != "fit" && command != "minos") return out;
  const QString option = fitter_->currentData().toString();
  if (!option.isEmpty()) out << option;
  if (maxCalls_->value() > 0) out << QString("--max-calls=%1").arg(maxCalls_->value());
  if (refreeze_->value() > 0) out << QString("--refreeze-covariance=%1").arg(refreeze_->value());
  return out;
}

bool FittingTab::apply(QString &why) {
  parameterSettings_->apply();
  mcmc_->apply();
  const QString value = priorCovariance_->text().trimmed();
  if (value == holder_->configValue("priorCovariance")) return true;
  return holder_->setConfig("priorCovariance", value, why, value.isEmpty());
}
