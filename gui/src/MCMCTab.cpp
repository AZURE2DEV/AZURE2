#include "MCMCTab.h"

#include <cmath>
#include <map>

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QSpinBox>
#include <QTabWidget>
#include <QTableWidget>
#include <QThread>
#include <QVBoxLayout>
#include <QSplitter>
#include <algorithm>
#include <vector>

#include "PlotWidget.h"

namespace {

QString Num(double v) { return std::isfinite(v) ? QString::number(v, 'g', 8) : QString(); }

}  // namespace

MCMCTab::MCMCTab(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  // Parameters: the priors table
  table_ = new QTableWidget(0, 5);
  table_->setObjectName("mcmcPriorTable");
  table_->setHorizontalHeaderLabels({tr("Parameter"), tr("Current Value"), tr("Prior Mean"), tr("Prior Std"),
                                     tr("Use Gaussian Prior")});
  table_->horizontalHeader()->setStretchLastSection(true);
  table_->verticalHeader()->setVisible(false);
  table_->setToolTip(tr("The sampler's parameter rows: a Gaussian prior on a free parameter, in the sampler's basis "
                        "(--basis): a width as its reduced-width amplitude (rwa) or as the project writes it "
                        "(physical).  Only checked rows are written."));
  connect(table_, &QTableWidget::itemChanged, this, [this]() {
    if (!filling_) dirty_ = true;
  });
  note_ = new QLabel(tr("Load the parameters to fill the table (the project is saved first)."));
  note_->setWordWrap(true);
  QPushButton *load = new QPushButton(tr("&Load Parameters"));
  connect(load, &QPushButton::clicked, this, &MCMCTab::parametersRequested);
  QWidget *paramsPage = new QWidget;
  QVBoxLayout *pl = new QVBoxLayout(paramsPage);
  QHBoxLayout *pb = new QHBoxLayout;
  pb->addWidget(load);
  pb->addStretch(1);
  pl->addLayout(pb);
  pl->addWidget(table_, 1);
  pl->addWidget(note_);

  // Sampling settings: the command's options, with the ranges and defaults of the engine's table
  walkers_ = new QSpinBox;
  walkers_->setObjectName("mcmcWalkers");
  steps_ = new QSpinBox;
  steps_->setObjectName("mcmcSteps");
  threads_ = new QSpinBox;
  spread_ = new QDoubleSpinBox;
  energySpread_ = new QDoubleSpinBox;
  basis_ = new QComboBox;
  basis_->setObjectName("mcmcBasis");
  fresh_ = new QCheckBox(tr("Fresh start (overwrite samples.mcmc; without it the sampler appends)"));
  walkers_->setRange(2, 1000);
  walkers_->setValue(50);
  steps_->setRange(100, 10000000);
  steps_->setValue(10000);
  threads_->setRange(1, std::max(1, QThread::idealThreadCount()));
  threads_->setValue(std::min(4, threads_->maximum()));
  spread_->setDecimals(3);
  spread_->setRange(0.001, 50);
  spread_->setValue(5);
  spread_->setSuffix(" %");
  energySpread_->setDecimals(3);
  energySpread_->setRange(0.001, 1);
  energySpread_->setValue(1);
  energySpread_->setSuffix(" keV");
  QFormLayout *form = new QFormLayout;
  for (const azure2::front::Option &o : azure2::front::Options()) {
    bool mcmc = false;
    for (const std::string &c : o.commands) mcmc = mcmc || c == "mcmc";
    if (!mcmc) continue;
    const QString option = QString::fromStdString(o.option), help = QString::fromStdString(o.help);
    QWidget *editor = nullptr;
    QString label;
    if (option == "--walkers=") editor = walkers_, label = tr("Number of Walkers:");
    else if (option == "--steps=") editor = steps_, label = tr("Number of Steps:");
    else if (option == "--parameter-spread=") editor = spread_, label = tr("Initial Chain Spread:");
    else if (option == "--energy-spread=") editor = energySpread_, label = tr("Level Energy Spread:");
    else if (option == "--sampler-threads=") editor = threads_, label = tr("Number of Threads:");
    else if (option == "--basis=") {
      for (const std::string &w : o.words) basis_->addItem(QString::fromStdString(w));
      editor = basis_, label = tr("Sampled Parameters:");
    } else if (option == "--overwrite-samples") {
      fresh_->setToolTip(help);
      continue;
    } else {
      continue;
    }
    if (o.high > o.low && qobject_cast<QSpinBox *>(editor))
      qobject_cast<QSpinBox *>(editor)->setRange((int)o.low, (int)o.high);
    editor->setToolTip(option + " " + help);
    // the row's help in one line beside the editor
    QWidget *row = new QWidget;
    QHBoxLayout *h = new QHBoxLayout(row);
    h->setContentsMargins(0, 0, 0, 0);
    h->addWidget(editor);
    QLabel *hint = new QLabel(help);
    hint->setProperty("role", "hint");
    hint->setWordWrap(true);
    h->addWidget(hint, 1);
    form->addRow(label, row);
  }
  form->addRow(QString(), fresh_);
  QGroupBox *settings = new QGroupBox(tr("Sampling Settings"));
  settings->setLayout(form);

  run_ = new QPushButton(tr("&Run MCMC Sampling"));
  run_->setObjectName("mcmcRun");
  run_->setProperty("primary", true);
  stop_ = new QPushButton(QIcon(":/actions/stop.svg"), tr("Stop Sampling"));
  stop_->setProperty("danger", true);
  stop_->setEnabled(false);
  connect(run_, &QPushButton::clicked, this, &MCMCTab::runRequested);
  connect(stop_, &QPushButton::clicked, this, &MCMCTab::stopRequested);
  bar_ = new QProgressBar;
  bar_->setObjectName("mcmcProgress");
  bar_->setRange(0, 1);
  bar_->setValue(0);
  status_ = new QLabel;
  status_->setObjectName("mcmcStatus");
  acceptance_ = new QLabel;
  acceptance_->setObjectName("mcmcAcceptance");
  log_ = new QPlainTextEdit;
  log_->setReadOnly(true);
  log_->setMaximumBlockCount(5000);
  log_->setFont(QFont("Monospace"));
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(run_);
  buttons->addWidget(stop_);
  buttons->addStretch(1);
  QWidget *runPage = new QWidget;
  QVBoxLayout *rl = new QVBoxLayout(runPage);
  rl->addWidget(settings);
  rl->addLayout(buttons);
  rl->addWidget(bar_);
  rl->addWidget(status_);
  rl->addWidget(acceptance_);
  rl->addWidget(log_, 1);

  // Results: the chain of the last run
  accepted_ = NAN;
  QPushButton *readButton = new QPushButton(tr("Read &Results"));
  readButton->setObjectName("mcmcReadResults");
  readButton->setProperty("primary", true);
  readButton->setToolTip(tr("Read samples.mcmc of the output directory (every step of the chain)."));
  connect(readButton, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!readResults(why)) resultsNote_->setText(why);
  });
  burnIn_ = new QSpinBox;
  burnIn_->setObjectName("mcmcBurnIn");
  burnIn_->setRange(0, 10000000);
  burnIn_->setToolTip(tr("The first steps of the chain left out of the statistics and the histogram."));
  connect(burnIn_, &QSpinBox::valueChanged, this, [this]() { fillResults(); });
  plotChoice_ = new QComboBox;
  plotChoice_->setObjectName("mcmcPlotChoice");
  connect(plotChoice_, &QComboBox::currentIndexChanged, this,
          [this]() { showResultPlot(plotChoice_->currentData().toString()); });
  resultsNote_ = new QLabel(tr("Run the sampler, then read its results."));
  resultsNote_->setObjectName("mcmcResultsNote");
  resultsNote_->setWordWrap(true);
  results_ = new QTableWidget(0, 6);
  results_->setObjectName("mcmcResults");
  results_->setHorizontalHeaderLabels({tr("Parameter"), tr("Mean"), tr("Std. dev."), tr("16 %"), tr("50 %"), tr("84 %")});
  results_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  results_->horizontalHeader()->setStretchLastSection(true);
  results_->verticalHeader()->setVisible(false);
  results_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  results_->setSelectionBehavior(QAbstractItemView::SelectRows);
  connect(results_, &QTableWidget::currentCellChanged, this, [this](int row) {
    if (row >= 0) {
      const int at = plotChoice_->findData(QString("trace:%1").arg(row));
      if (at >= 0) plotChoice_->setCurrentIndex(at);
    }
  });
  resultsPlot_ = new PlotWidget;
  resultsPlot_->setObjectName("mcmcResultsPlot");
  resultsPlot_->setShowLegend(false);
  resultsPlot_->setShowResiduals(false);
  QHBoxLayout *resultsTop = new QHBoxLayout;
  resultsTop->addWidget(readButton);
  resultsTop->addSpacing(12);
  resultsTop->addWidget(new QLabel(tr("Burn-in (steps):")));
  resultsTop->addWidget(burnIn_);
  resultsTop->addSpacing(12);
  resultsTop->addWidget(new QLabel(tr("Plot:")));
  resultsTop->addWidget(plotChoice_, 1);
  QSplitter *resultsSplit = new QSplitter(Qt::Horizontal);
  resultsSplit->addWidget(results_);
  resultsSplit->addWidget(resultsPlot_);
  resultsSplit->setStretchFactor(0, 2);
  resultsSplit->setStretchFactor(1, 3);
  resultsSplit->setSizes({520, 780});
  QWidget *resultsPage = new QWidget;
  QVBoxLayout *resl = new QVBoxLayout(resultsPage);
  resl->addLayout(resultsTop);
  resl->addWidget(resultsNote_);
  resl->addWidget(resultsSplit, 1);

  QTabWidget *pages = new QTabWidget;
  pages->setObjectName("mcmcPages");
  pages->addTab(paramsPage, tr("Parameters"));
  pages->addTab(runPage, tr("Sampling Settings"));
  pages->addTab(resultsPage, tr("Results"));
  (new QVBoxLayout(this))->addWidget(pages);
  connect(holder_, &ProjectHolder::opened, this, &MCMCTab::reload);
  reload();
}

void MCMCTab::reload() {
  filling_ = true;
  table_->setRowCount(0);
  for (const azure2::front::SamplerPrior &p : azure2::front::SamplerPriors(holder_->project())) {
    const int row = table_->rowCount();
    table_->insertRow(row);
    QTableWidgetItem *name = new QTableWidgetItem(QString::fromStdString(p.name));
    name->setFlags(name->flags() & ~Qt::ItemIsEditable);
    table_->setItem(row, 0, name);
    table_->setItem(row, 1, new QTableWidgetItem(Num(p.value)));
    table_->setItem(row, 2, new QTableWidgetItem(Num(p.mean)));
    table_->setItem(row, 3, new QTableWidgetItem(Num(p.sigma)));
    QTableWidgetItem *use = new QTableWidgetItem;
    use->setFlags((use->flags() | Qt::ItemIsUserCheckable) & ~Qt::ItemIsEditable);
    use->setCheckState(p.use ? Qt::Checked : Qt::Unchecked);
    table_->setItem(row, 4, use);
  }
  filling_ = false;
  dirty_ = false;
}

bool MCMCTab::readParameters(const QString &projectFile, QString &why) {
  std::vector<azure2::front::Parameter> params;
  std::string reason;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const bool ok = azure2::front::ReadParameters(projectFile.toStdString(), params, reason);
  QApplication::restoreOverrideCursor();
  if (!ok) {
    why = QString::fromStdString(reason);
    note_->setText(why);
    return false;
  }
  std::map<std::string, azure2::front::SamplerPrior> rows;
  for (const azure2::front::SamplerPrior &p : azure2::front::SamplerPriors(holder_->project())) rows[p.name] = p;
  filling_ = true;
  table_->setRowCount(0);
  for (const azure2::front::Parameter &p : params) {
    if (p.fixed) continue;
    const auto r = rows.find(p.name);
    const bool has = r != rows.end();
    const int row = table_->rowCount();
    table_->insertRow(row);
    QTableWidgetItem *name = new QTableWidgetItem(QString::fromStdString(p.name));
    name->setFlags(name->flags() & ~Qt::ItemIsEditable);
    name->setToolTip(QString::fromStdString(p.kind + ", " + p.unit));
    table_->setItem(row, 0, name);
    table_->setItem(row, 1, new QTableWidgetItem(Num(p.value)));
    table_->setItem(row, 2, new QTableWidgetItem(Num(has ? r->second.mean : p.value)));
    table_->setItem(row, 3, new QTableWidgetItem(Num(has ? r->second.sigma : 0.1 * std::fabs(p.value))));
    QTableWidgetItem *use = new QTableWidgetItem;
    use->setFlags((use->flags() | Qt::ItemIsUserCheckable) & ~Qt::ItemIsEditable);
    use->setCheckState(has && r->second.use ? Qt::Checked : Qt::Unchecked);
    table_->setItem(row, 4, use);
  }
  filling_ = false;
  note_->setText(tr("%1 free parameters (values in fit units, as param.sav).").arg(table_->rowCount()));
  return true;
}

void MCMCTab::apply() {
  if (!dirty_) return;
  std::vector<azure2::front::SamplerPrior> rows;
  for (int row = 0; row < table_->rowCount(); row++) {
    if (table_->item(row, 4)->checkState() != Qt::Checked) continue;
    azure2::front::SamplerPrior p;
    p.name = table_->item(row, 0)->text().toStdString();
    p.value = table_->item(row, 1)->text().toDouble();
    p.mean = table_->item(row, 2)->text().toDouble();
    p.sigma = table_->item(row, 3)->text().toDouble();
    p.use = true;
    rows.push_back(p);
  }
  azure2::front::SetSamplerPriors(holder_->project(), rows);
  dirty_ = false;
  emit holder_->sectionEdited("mcmc");
}

QStringList MCMCTab::options() const {
  QStringList out;
  out << QString("--walkers=%1").arg(walkers_->value()) << QString("--steps=%1").arg(steps_->value())
      << QString("--parameter-spread=%1").arg(spread_->value()) << QString("--energy-spread=%1").arg(energySpread_->value())
      << QString("--sampler-threads=%1").arg(threads_->value()) << QString("--basis=%1").arg(basis_->currentText());
  if (fresh_->isChecked()) out << "--overwrite-samples";
  return out;
}

void MCMCTab::setWalkers(int walkers) { walkers_->setValue(walkers); }
void MCMCTab::setSteps(int steps) { steps_->setValue(steps); }
void MCMCTab::setSpread(double percent) { spread_->setValue(percent); }
void MCMCTab::setThreads(int threads) { threads_->setValue(threads); }
bool MCMCTab::setBasis(const QString &basis) {
  const int at = basis_->findText(basis);
  if (at >= 0) basis_->setCurrentIndex(at);
  return at >= 0;
}
void MCMCTab::setFreshStart(bool fresh) { fresh_->setChecked(fresh); }

void MCMCTab::setRunning(bool running) {
  run_->setEnabled(!running);
  stop_->setEnabled(running);
  if (running) {
    log_->clear();
    acceptance_->clear();
    bar_->setRange(0, steps_->value());
    bar_->setValue(0);
    lastStep_ = 0;
  }
}

QString MCMCTab::acceptance() const { return acceptance_->text(); }

void MCMCTab::showProgress(const QString &stage, qlonglong step, qlonglong total, double logPosterior,
                           double logLikelihood) {
  if (stage != "mcmc") return;
  if (total > 0) bar_->setRange(0, (int)total);
  bar_->setValue((int)step);
  lastStep_ = (int)step;
  status_->setText(tr("step %1 of %2: log posterior %3, log likelihood %4")
                       .arg(step)
                       .arg(total)
                       .arg(logPosterior, 0, 'g', 10)
                       .arg(logLikelihood, 0, 'g', 10));
}

void MCMCTab::appendLine(const QString &line, bool transient) {
  if (transient) return;
  log_->appendPlainText(line);
  // the sampler's own line, shown as the program prints it
  // (it follows the last progress line, which begins with '\r' and so shares its line)
  const int at = line.indexOf("Acceptance fraction:");
  if (at >= 0) acceptance_->setText(line.mid(at).trimmed());
}

// ------------------------------------------------------------------------------------------------- the Results page

bool MCMCTab::readResults(QString &why) {
  std::string reason;
  azure2::front::SamplerChain chain;
  if (!azure2::front::ReadSamplerChain(holder_->project(), chain, reason)) {
    why = QString::fromStdString(reason);
    return false;
  }
  chain_ = std::move(chain);
  // the parameters' names: the priors table's rows when they are the run's free parameters, else by place
  chainNames_.clear();
  for (int k = 0; k < chain_.parameters; k++)
    chainNames_ << (table_->rowCount() == chain_.parameters && table_->item(k, 0) ? table_->item(k, 0)->text()
                                                                                    : tr("parameter %1").arg(k));
  plotChoice_->blockSignals(true);
  plotChoice_->clear();
  plotChoice_->addItem(tr("Log posterior (mean of the walkers) against the step"), "logprob");
  for (int k = 0; k < chain_.parameters; k++) plotChoice_->addItem(tr("Trace: %1").arg(chainNames_[k]), QString("trace:%1").arg(k));
  for (int k = 0; k < chain_.parameters; k++)
    plotChoice_->addItem(tr("Histogram: %1").arg(chainNames_[k]), QString("histogram:%1").arg(k));
  plotChoice_->blockSignals(false);
  burnIn_->blockSignals(true);
  burnIn_->setMaximum(std::max(0, chain_.steps - 1));
  if (burnIn_->value() >= chain_.steps) burnIn_->setValue(chain_.steps / 5);
  burnIn_->blockSignals(false);
  fillResults();
  return true;
}

void MCMCTab::setBurnIn(int steps) { burnIn_->setValue(steps); }

void MCMCTab::fillResults() {
  const int S = chain_.steps, W = chain_.walkers, P = chain_.parameters;
  results_->setRowCount(0);
  if (S == 0 || W == 0) return;
  const int burn = std::min(burnIn_->value(), S - 1);
  auto at = [&](int step, int walker, int k) { return chain_.chain[((size_t)step * W + walker) * P + k]; };
  // the acceptance fraction: moves after the burn-in that changed the position
  long moves = 0, changed = 0;
  for (int t = std::max(1, burn); t < S; t++)
    for (int w = 0; w < W; w++) {
      bool same = true;
      for (int k = 0; k < P && same; k++) same = at(t, w, k) == at(t - 1, w, k);
      moves++;
      changed += !same;
    }
  accepted_ = moves ? (double)changed / moves : NAN;
  for (int k = 0; k < P; k++) {
    std::vector<double> v;
    v.reserve((size_t)(S - burn) * W);
    for (int t = burn; t < S; t++)
      for (int w = 0; w < W; w++) v.push_back(at(t, w, k));
    double mean = 0;
    for (double x : v) mean += x;
    mean /= v.size();
    double var = 0;
    for (double x : v) var += (x - mean) * (x - mean);
    const double sd = v.size() > 1 ? std::sqrt(var / (v.size() - 1)) : 0;
    std::sort(v.begin(), v.end());
    auto quantile = [&v](double q) {  // linear between order statistics
      const double pos = q * (v.size() - 1);
      const size_t lo = (size_t)std::floor(pos), hi = std::min(lo + 1, v.size() - 1);
      return v[lo] + (pos - lo) * (v[hi] - v[lo]);
    };
    const int row = results_->rowCount();
    results_->insertRow(row);
    results_->setItem(row, 0, new QTableWidgetItem(chainNames_.value(k)));
    const double values[] = {mean, sd, quantile(0.16), quantile(0.5), quantile(0.84)};
    for (int c = 0; c < 5; c++) {
      QTableWidgetItem *item = new QTableWidgetItem(QString::number(values[c], 'g', 7));
      item->setData(Qt::UserRole, values[c]);
      item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
      results_->setItem(row, c + 1, item);
    }
  }
  resultsNote_->setText(tr("%1 steps of %2 walkers, %3 parameters (%4); the first %5 steps left out.  Acceptance "
                           "fraction %6.")
                            .arg(S)
                            .arg(W)
                            .arg(P)
                            .arg(QString::fromStdString(chain_.file))
                            .arg(burn)
                            .arg(accepted_, 0, 'f', 3));
  showResultPlot(plotChoice_->currentData().toString());
}

bool MCMCTab::showResultPlot(const QString &what) {
  const int S = chain_.steps, W = chain_.walkers, P = chain_.parameters;
  QList<PlotSeries> series;
  const int burn = std::min(burnIn_->value(), std::max(0, S - 1));
  const int k = what.section(':', 1).toInt();
  const QColor colours[] = {QColor(0x1f, 0x77, 0xb4), QColor(0xd6, 0x27, 0x28), QColor(0x2c, 0xa0, 0x2c),
                            QColor(0xff, 0x7f, 0x0e), QColor(0x94, 0x67, 0xbd), QColor(0x8c, 0x56, 0x4b)};
  if (S == 0) {
    resultsPlot_->setSeries(series);
    return false;
  }
  if (what == "logprob") {
    PlotSeries s;
    s.label = tr("mean log posterior");
    s.color = colours[0];
    for (int t = 0; t < S; t++) {
      double m = 0;
      for (int w = 0; w < W; w++) m += chain_.logProbability[(size_t)t * W + w];
      s.modelX << t;
      s.modelY << m / W;
    }
    series << s;
    resultsPlot_->setAxisTitles(tr("step"), tr("log posterior"));
  } else if (what.startsWith("trace:") && k >= 0 && k < P) {
    for (int w = 0; w < W; w++) {
      PlotSeries s;
      s.label = tr("walker %1").arg(w);
      QColor c = colours[w % 6];
      c.setAlpha(150);
      s.color = c;
      for (int t = 0; t < S; t++) {
        s.modelX << t;
        s.modelY << chain_.chain[((size_t)t * W + w) * P + k];
      }
      series << s;
    }
    resultsPlot_->setAxisTitles(tr("step"), chainNames_.value(k));
  } else if (what.startsWith("histogram:") && k >= 0 && k < P) {
    std::vector<double> v;
    for (int t = burn; t < S; t++)
      for (int w = 0; w < W; w++) v.push_back(chain_.chain[((size_t)t * W + w) * P + k]);
    const auto [lo, hi] = std::minmax_element(v.begin(), v.end());
    const int bins = std::max(5, std::min(60, (int)std::sqrt((double)v.size())));
    const double a = *lo, width = *hi > *lo ? (*hi - *lo) / bins : 1.0;
    std::vector<int> count(bins, 0);
    for (double x : v) count[std::min(bins - 1, (int)((x - a) / width))]++;
    PlotSeries s;
    s.label = chainNames_.value(k);
    s.color = colours[0];
    for (int b = 0; b < bins; b++) {  // a step outline
      s.modelX << a + b * width << a + (b + 1) * width;
      s.modelY << count[b] << count[b];
    }
    series << s;
    resultsPlot_->setAxisTitles(chainNames_.value(k), tr("samples per bin (after the burn-in)"));
  } else {
    return false;
  }
  resultsPlot_->setSeries(series);
  const int choice = plotChoice_->findData(what);
  if (choice >= 0 && choice != plotChoice_->currentIndex()) {
    plotChoice_->blockSignals(true);
    plotChoice_->setCurrentIndex(choice);
    plotChoice_->blockSignals(false);
  }
  return true;
}
