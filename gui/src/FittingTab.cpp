#include "FittingTab.h"
#include "ExpertMode.h"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QFile>
#include <QFileInfo>
#include <QHeaderView>
#include <QMessageBox>
#include <QPushButton>
#include <QTableWidget>
#include <algorithm>
#include <cmath>
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

#include <map>

FittingTab::FittingTab(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  fitter_ = new QComboBox;
  fitter_->setObjectName("fitterCombo");
  for (const azure2::front::Fitter &f : azure2::front::Fitters()) {
    fitter_->addItem(QString::fromStdString(f.name), QString::fromStdString(f.option));
    fitter_->setItemData(fitter_->count() - 1, QString::fromStdString(f.doc), Qt::ToolTipRole);
  }
  // the program's own choice without a switch (the row with no option: MIGRAD), not the table's first row
  fitter_->setCurrentIndex(std::max(0, fitter_->findData(QString())));
  maxCalls_ = new QSpinBox;
  maxCalls_->setRange(0, 1000000000);
  maxCalls_->setSpecialValueText(tr("none"));
  maxCalls_->setToolTip(tr("The fit's call limit (--max-calls); a fit stopped by it ends at the lowest objective seen."));
  refreeze_ = new QSpinBox;
  refreeze_->setRange(0, 1000);
  refreeze_->setSpecialValueText(tr("none"));
  refreeze_->setToolTip(tr("Refit up to N times with relvector covariance terms re-frozen at the fitted parameters "
                           "(--refreeze-covariance)."));
  expert_ = expert::Enabled();
  QFormLayout *fitForm = new QFormLayout;
  fitForm_ = fitForm;
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
  QVBoxLayout *priorLayout = new QVBoxLayout;
  QFormLayout *priorForm = new QFormLayout;
  priorForm->addRow(tr("Prior covariance file:"), priorRow);
  priorLayout->addLayout(priorForm);
  QGroupBox *priorBox = new QGroupBox(tr("Priors"));
  priorBox->setLayout(priorLayout);
  // where a single parameter's prior is written: the box's tip, not a paragraph (qt_gui.md section 19.6)
  priorBox->setToolTip(tr("A Gaussian prior on one parameter is written on the parameter's own field: \"norm 1 fit "
                          "prior 1 0.05\" on a data set, \"areal_density ... fit prior ...\" on a target, an energy "
                          "scale's or a normalization record's value.  Nuisance parameters are the fields written "
                          "\"fit\"."));
  priorBox_ = priorBox;

  parameterSettings_ = new SectionTextPanel(
      holder, "parameterSettings",
      tr("Parameter settings (<parameterSettings>): per parameter, by its name, the limits and whether it is fixed, "
         "overriding the levels' and records' flags."));
  mcmc_ = new SectionTextPanel(holder, "mcmc",
                               tr("The sampler's priors and settings (<mcmc>), read by the command mcmc where the build "
                                  "has the sampler."));
  // the parameter table
  table_ = new QTableWidget(0, 11);
  table_->setObjectName("parameterTable");
  table_->setHorizontalHeaderLabels({tr("Name"), tr("Kind"), tr("Start"), tr("Fixed"), tr("Lower Limit"),
                                     tr("Upper Limit"), tr("Prior"), tr("Prior Centre"), tr("Prior Sigma"),
                                     tr("Fitted"), tr("Fit Error")});
  table_->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
  table_->horizontalHeader()->setStretchLastSection(true);
  table_->verticalHeader()->setVisible(false);
  table_->setToolTip(tr("Every parameter of a fit as the engine sets it up.  Start, Fitted and Fit Error are in the "
                        "fit's units (as param.sav: a width as its reduced-width amplitude, MeV^1/2).  Limits and "
                        "the prior are the parameter's <parameterSettings> row: a level width's in eV (partial "
                        "width), converted by the engine; empty limits: none.  Fixed or free is set where the "
                        "parameter is written (the Levels tab, a record's 'fixed'/'fit')."));
  connect(table_, &QTableWidget::itemChanged, this, [this](QTableWidgetItem *item) {
    if (filling_ || item->column() < 4 || item->column() > 8) return;
    QString why;
    if (!applyRow(item->row(), why)) QMessageBox::warning(this, tr("Parameter settings"), why);
  });
  tableNote_ = new QLabel(tr("Read the parameters to fill the table (the project is saved first)."));
  tableNote_->setObjectName("parameterNote");
  tableNote_->setWordWrap(true);
  QPushButton *read = new QPushButton(tr("&Read Parameters"));
  read->setObjectName("readParameters");
  read->setToolTip(tr("Save the project and set up a session of it to list its parameters."));
  connect(read, &QPushButton::clicked, this, &FittingTab::parametersRequested);
  QPushButton *fitted = new QPushButton(tr("Show &Fitted Values"));
  fitted->setToolTip(tr("The values and errors of the last fit's param.sav, beside the start values."));
  connect(fitted, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!readFitted(fittedFile(), why)) QMessageBox::warning(this, tr("Fitted values"), why);
  });
  QPushButton *take = new QPushButton(tr("&Take the Fitted Values as the New Start"));
  take->setObjectName("takeFitted");
  take->setToolTip(tr("Copy the last fit's param.sav beside the project (<project>.start.sav) and start the next "
                      "run from it (--parameter-file).  The project's own values are not changed."));
  connect(take, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!takeFittedAsStart(fittedFile(), why)) QMessageBox::warning(this, tr("New start"), why);
  });
  QHBoxLayout *tableButtons = new QHBoxLayout;
  tableButtons->addWidget(read);
  tableButtons->addWidget(fitted);
  tableButtons->addWidget(take);
  tableButtons->addStretch(1);
  QWidget *tablePage = new QWidget;
  QVBoxLayout *tableLayout = new QVBoxLayout(tablePage);
  tableLayout->addLayout(tableButtons);
  tableLayout->addWidget(table_, 1);
  tableLayout->addWidget(tableNote_);

  QTabWidget *texts = new QTabWidget;
  pages_ = texts;
  texts->setObjectName("fittingPages");
  texts->addTab(tablePage, tr("Parameters"));
  texts->addTab(parameterSettings_, tr("Parameter Settings (text)"));
  texts->addTab(mcmc_, tr("Sampler (MCMC)"));

  QHBoxLayout *top = new QHBoxLayout;
  top->addWidget(fitBox);
  top->addWidget(priorBox, 1);
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addLayout(top);
  layout->addWidget(texts, 1);
  connect(holder_, &ProjectHolder::opened, this, &FittingTab::reload);
  // a row of the Parameters table writes <parameterSettings>: its text page shown at once in standard view (a hidden
  // page must never hold lines in effect)
  connect(holder_, &ProjectHolder::sectionEdited, this, [this](const QString &s) {
    if (s == "parameterSettings" || s == "mcmc" || s == "config") updateExpertRows();
  });
  reload();
}

void FittingTab::reload() {
  priorCovariance_->setText(holder_->configValue("priorCovariance"));
  updateExpertRows();
}

void FittingTab::setExpertMode(bool on) {
  expert_ = on;
  updateExpertRows();
}

void FittingTab::updateExpertRows() {
  auto row = [this](QWidget *field, bool shown) {
    field->setVisible(shown);
    if (QWidget *label = fitForm_->labelForField(field)) label->setVisible(shown);
  };
  row(maxCalls_, expert_ || maxCalls_->value() > 0);
  row(refreeze_, expert_ || refreeze_->value() > 0);
  priorBox_->setVisible(expert_ || !priorCovariance_->text().trimmed().isEmpty());
  // a text page with lines other than comments is in effect: shown in either view
  auto inUse = [this](const QString &section) {
    for (const QString &line : holder_->sectionBody(section)) {
      const QString t = line.trimmed();
      if (!t.isEmpty() && !t.startsWith('#')) return true;
    }
    return false;
  };
  pages_->setTabVisible(pages_->indexOf(parameterSettings_), expert_ || inUse("parameterSettings"));
  pages_->setTabVisible(pages_->indexOf(mcmc_), expert_ || inUse("mcmc"));
}

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

namespace {

QTableWidgetItem *Fixed(const QString &text) {
  QTableWidgetItem *i = new QTableWidgetItem(text);
  i->setFlags(i->flags() & ~Qt::ItemIsEditable);
  return i;
}

QString Num(double v) { return std::isfinite(v) ? QString::number(v, 'g', 8) : QString(); }

}  // namespace

QString FittingTab::fittedFile() const {
  return QDir(QDir(holder_->directory()).filePath(holder_->configValue("output_dir"))).filePath("param.sav");
}

bool FittingTab::readParameters(const QString &projectFile, QString &why) {
  std::vector<azure2::front::Parameter> params;
  std::string reason;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const bool ok = azure2::front::ReadParameters(projectFile.toStdString(), params, reason);
  QApplication::restoreOverrideCursor();
  if (!ok) {
    why = QString::fromStdString(reason);
    tableNote_->setText(why);
    return false;
  }
  std::map<std::string, azure2::front::ParameterSetting> settings;
  for (const azure2::front::ParameterSetting &r : azure2::front::ParameterSettings(holder_->project())) settings[r.name] = r;
  filling_ = true;
  table_->setRowCount(0);
  int free = 0;
  for (const azure2::front::Parameter &p : params) {
    const int row = table_->rowCount();
    table_->insertRow(row);
    table_->setItem(row, 0, Fixed(QString::fromStdString(p.name)));
    table_->setItem(row, 1, Fixed(QString::fromStdString(p.kind)));
    QTableWidgetItem *start = Fixed(Num(p.value));
    start->setToolTip(QString::fromStdString(p.unit));
    table_->setItem(row, 2, start);
    table_->setItem(row, 3, Fixed(p.fixed ? tr("fixed") : tr("free")));
    free += !p.fixed;
    const auto s = settings.find(p.name);
    const bool has = s != settings.end();
    const bool limits = has && (s->second.lower != 0 || s->second.upper != 0);
    QTableWidgetItem *lower = new QTableWidgetItem(limits ? Num(s->second.lower) : QString());
    QTableWidgetItem *upper = new QTableWidgetItem(limits ? Num(s->second.upper) : QString());
    const QString now = tr("the fit's limits now (fit units): %1 .. %2")
                            .arg(std::isfinite(p.lower) ? Num(p.lower) : tr("none"),
                                 std::isfinite(p.upper) ? Num(p.upper) : tr("none"));
    lower->setToolTip(now);
    upper->setToolTip(now);
    table_->setItem(row, 4, lower);
    table_->setItem(row, 5, upper);
    QTableWidgetItem *prior = new QTableWidgetItem;
    prior->setFlags((prior->flags() | Qt::ItemIsUserCheckable) & ~Qt::ItemIsEditable);
    prior->setCheckState(has && s->second.nuisance ? Qt::Checked : Qt::Unchecked);
    table_->setItem(row, 6, prior);
    table_->setItem(row, 7, new QTableWidgetItem(has && s->second.nuisance ? Num(s->second.nominal) : QString()));
    table_->setItem(row, 8, new QTableWidgetItem(has && s->second.nuisance ? Num(s->second.error) : QString()));
    table_->setItem(row, 9, Fixed(QString()));
    table_->setItem(row, 10, Fixed(QString()));
  }
  filling_ = false;
  table_->resizeColumnsToContents();
  tableNote_->setText(tr("%1 parameters, %2 free.  Limits and prior: the parameter's <parameterSettings> row (a level "
                         "width's in eV); they apply to free parameters.")
                          .arg(params.size())
                          .arg(free));
  return true;
}

bool FittingTab::readFitted(const QString &file, QString &why) {
  std::vector<azure2::front::ParameterValue> values;
  std::string reason;
  if (!azure2::front::ReadParameterFile(file.toStdString(), values, reason)) {
    why = QString::fromStdString(reason) + tr(" (no fit has written it yet)");
    return false;
  }
  filling_ = true;
  int shown = 0;
  for (const azure2::front::ParameterValue &v : values)
    for (int row = 0; row < table_->rowCount(); row++)
      if (table_->item(row, 0)->text().toStdString() == v.name) {
        table_->item(row, 9)->setText(Num(v.value));
        table_->item(row, 10)->setText(Num(v.error));
        shown++;
      }
  filling_ = false;
  tableNote_->setText(tr("%1 fitted values of %2 beside the start values.").arg(shown).arg(file));
  return true;
}

bool FittingTab::takeFittedAsStart(const QString &fitted, QString &why) {
  if (!QFileInfo::exists(fitted)) {
    why = tr("%1: no fit has written it yet").arg(fitted);
    return false;
  }
  const QFileInfo project(holder_->path());
  const QString start = project.dir().filePath(project.completeBaseName() + ".start.sav");
  QFile::remove(start);
  if (!QFile::copy(fitted, start)) {
    why = tr("cannot write %1").arg(start);
    return false;
  }
  emit startFileChosen(start);
  tableNote_->setText(tr("The next run starts from %1 (--parameter-file).").arg(start));
  return true;
}

bool FittingTab::applyRow(int row, QString &why) {
  if (row < 0 || row >= table_->rowCount()) return false;
  auto number = [this, row](int column, bool &ok) {
    const QString t = table_->item(row, column)->text().trimmed();
    ok = true;
    if (t.isEmpty()) return 0.0;
    return t.toDouble(&ok);
  };
  bool okLow, okHigh, okCentre, okSigma;
  azure2::front::ParameterSetting s;
  s.name = table_->item(row, 0)->text().toStdString();
  s.lower = number(4, okLow);
  s.upper = number(5, okHigh);
  s.nuisance = table_->item(row, 6)->checkState() == Qt::Checked;
  s.nominal = number(7, okCentre);
  s.error = number(8, okSigma);
  if (!okLow || !okHigh || !okCentre || !okSigma) {
    why = tr("%1: a limit or a prior number does not read as a number").arg(QString::fromStdString(s.name));
    return false;
  }
  const QString kind = table_->item(row, 1)->text();
  s.category = kind == "energy" || kind == "width" ? "level" : kind.toStdString();
  const bool empty = s.lower == 0 && s.upper == 0 && !s.nuisance;
  std::string reason;
  if (!azure2::front::SetParameterSetting(holder_->project(), s, reason, empty)) {
    why = QString::fromStdString(reason);
    return false;
  }
  emit holder_->sectionEdited("parameterSettings");
  parameterSettings_->reload();
  return true;
}
