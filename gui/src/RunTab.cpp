#include "RunTab.h"
#include "ExpertMode.h"
#include <QListView>

#include <algorithm>
#include <cmath>

#include <QCheckBox>
#include <QComboBox>
#include <QDesktopServices>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFontDatabase>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPainter>
#include <QPainterPath>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QRegularExpression>
#include <QSpinBox>
#include <QSplitter>
#include <QSyntaxHighlighter>
#include <QTableWidget>
#include <QThread>
#include <QToolButton>
#include <QUrl>
#include <QVBoxLayout>

#include "azure2/front/Engine.h"

namespace {

/// The old program's words for the commands (its Calculation Type list), the engine's name otherwise.
QString Title(const QString &command) {
  if (command == "calculate") return QObject::tr("Calculate Segments From Data");
  if (command == "fit") return QObject::tr("Fit Segments From Data");
  if (command == "extrapolate") return QObject::tr("Calculate Segments Without Data");
  if (command == "minos") return QObject::tr("Perform MINOS Error Analysis");
  if (command == "rate") return QObject::tr("Calculate Reaction Rate");
  if (command == "mcmc") return QObject::tr("MCMC Bayesian Sampling");
  return command;
}

/// A number of the option table: the value written as the command line would ('g', no trailing zeros); the box's
/// minimum, one below the row's low end, means "not set" (the option left out: the engine's default applies).
class NumberBox : public QDoubleSpinBox {
 public:
  using QDoubleSpinBox::QDoubleSpinBox;
  QString textFromValue(double v) const override { return QString::number(v, 'g', 10); }
};

double Low(const azure2::front::Option &o) { return std::isfinite(o.low) ? o.low : -1e9; }
double High(const azure2::front::Option &o) { return std::isfinite(o.high) && o.high > o.low ? o.high : 1e9; }

/// A table item that sorts by its number (UserRole) where it has one.
class NumberItem : public QTableWidgetItem {
 public:
  NumberItem(const QString &text, double value) : QTableWidgetItem(text) {
    setData(Qt::UserRole, value);
    setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
  }
  bool operator<(const QTableWidgetItem &other) const override {
    const QVariant a = data(Qt::UserRole), b = other.data(Qt::UserRole);
    if (a.isValid() && b.isValid()) return a.toDouble() < b.toDouble();
    return QTableWidgetItem::operator<(other);
  }
};

/// The log's colours: a warning in amber, an error in red, bold (RunTab::lineSeverity).
class LogHighlighter : public QSyntaxHighlighter {
 public:
  using QSyntaxHighlighter::QSyntaxHighlighter;

 protected:
  void highlightBlock(const QString &text) override {
    const int severity = RunTab::lineSeverity(text);
    if (severity == 0) return;
    const bool dark = document() && QPalette().color(QPalette::Base).lightness() < 128;
    QTextCharFormat f;
    f.setForeground(severity == 2 ? QColor(dark ? "#ff6b6b" : "#c62828") : QColor(dark ? "#ffb74d" : "#a35a00"));
    if (severity == 2) f.setFontWeight(QFont::Bold);
    setFormat(0, text.size(), f);
  }
};

}  // namespace

// ---------------------------------------------------------------------------------------------- ChiSquaredTrace

ChiSquaredTrace::ChiSquaredTrace(QWidget *parent) : QWidget(parent) {
  setObjectName("chiTrace");
  setToolTip(tr("The best chi-squared of the fit so far against the evaluation (log scale)."));
  setMinimumHeight(fontMetrics().height() * 3);
}

QSize ChiSquaredTrace::sizeHint() const { return QSize(fontMetrics().horizontalAdvance("M") * 24, fontMetrics().height() * 3); }

void ChiSquaredTrace::clear() {
  x_.clear();
  y_.clear();
  update();
}

void ChiSquaredTrace::add(double evaluation, double chiSquared) {
  if (!(chiSquared > 0) || !std::isfinite(chiSquared)) return;
  if (!y_.isEmpty() && chiSquared >= y_.last() && evaluation <= x_.last()) return;
  x_ << evaluation;
  y_ << chiSquared;
  update();
}

void ChiSquaredTrace::paintEvent(QPaintEvent *) {
  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing);
  const QRectF r = QRectF(rect()).adjusted(2, 2, -2, -2);
  p.setPen(QPen(palette().color(QPalette::Mid), 1));
  p.setBrush(palette().color(QPalette::Base));
  p.drawRoundedRect(r, 3, 3);
  if (x_.size() < 2) return;
  const double x0 = x_.first(), x1 = std::max(x_.last(), x0 + 1);
  double y0 = INFINITY, y1 = -INFINITY;
  for (double y : y_) {
    y0 = std::min(y0, std::log10(y));
    y1 = std::max(y1, std::log10(y));
  }
  if (y1 - y0 < 1e-9) y1 = y0 + 1;
  const QRectF in = r.adjusted(4, 4, -4, -4);
  QPainterPath path;
  for (int i = 0; i < x_.size(); i++) {
    const QPointF pt(in.left() + (x_[i] - x0) / (x1 - x0) * in.width(),
                     in.bottom() - (std::log10(y_[i]) - y0) / (y1 - y0) * in.height());
    if (i == 0) path.moveTo(pt);
    else path.lineTo(pt);
  }
  p.setPen(QPen(palette().color(QPalette::Highlight), 1.6));
  p.setBrush(Qt::NoBrush);
  p.drawPath(path);
}

// ------------------------------------------------------------------------------------------------------- RunTab

int RunTab::lineSeverity(const QString &line) {
  const QRegularExpression error("\\b(error|refused|aborted|failed|fatal)\\b",
                                        QRegularExpression::CaseInsensitiveOption);
  const QRegularExpression warning("\\bwarning\\b|^\\s*W_[A-Z_]+|\\bnot converged\\b",
                                          QRegularExpression::CaseInsensitiveOption);
  if (error.match(line).hasMatch()) return 2;
  if (warning.match(line).hasMatch()) return 1;
  return 0;
}

RunTab::RunTab(QWidget *parent) : QWidget(parent) {
  expert_ = expert::Enabled();
  command_ = new QComboBox;
  command_->setObjectName("calcType");
  for (const azure2::front::Command &c : azure2::front::Commands()) {
    if (!c.built) continue;
    command_->addItem(Title(QString::fromStdString(c.name)), QString::fromStdString(c.name));
    command_->setItemData(command_->count() - 1, QString::fromStdString(c.help), Qt::ToolTipRole);
  }
  connect(command_, &QComboBox::currentIndexChanged, this, &RunTab::commandChanged);
  commandHelp_ = new QLabel;
  commandHelp_->setObjectName("commandHelp");
  commandHelp_->setProperty("role", "hint");
  commandHelp_->setWordWrap(true);
  threads_ = new QSpinBox;
  threads_->setObjectName("threads");
  threads_->setRange(1, std::max(1, QThread::idealThreadCount()));
  threads_->setValue(std::min(4, threads_->maximum()));
  threads_->setToolTip(tr("The OpenMP threads of the run (--threads); results do not depend on it."));
  runButton_ = new QPushButton(tr("Save and &Run"));
  runButton_->setObjectName("runButton");
  runButton_->setProperty("primary", true);
  runButton_->setToolTip(tr("Save the project, then run the command (F5 runs calculate, F6 fit)."));
  stopButton_ = new QPushButton(QIcon(":/actions/stop.svg"), tr("Stop"));
  stopButton_->setObjectName("stopButton");
  stopButton_->setProperty("danger", true);
  stopButton_->setEnabled(false);
  stopButton_->setToolTip(tr("Stop the run: a fit ends at its best point and writes its outputs there (Shift+F5)."));
  connect(runButton_, &QPushButton::clicked, this, &RunTab::runRequested);
  connect(stopButton_, &QPushButton::clicked, this, &RunTab::stopRequested);

  QGridLayout *calcLayout = new QGridLayout;
  calcLayout->addWidget(new QLabel(tr("Calculation Type:")), 0, 0, Qt::AlignRight);
  calcLayout->addWidget(command_, 0, 1);
  calcLayout->addWidget(new QLabel(tr("Threads:")), 0, 2, Qt::AlignRight);
  calcLayout->addWidget(threads_, 0, 3);
  calcLayout->setColumnStretch(4, 1);
  calcLayout->addWidget(runButton_, 0, 5);
  calcLayout->addWidget(stopButton_, 0, 6);
  calcLayout->addWidget(commandHelp_, 1, 1, 1, 6);

  optionsBox_ = new QGroupBox(tr("Options of the calculation"));
  optionsGrid_ = new QGridLayout(optionsBox_);
  optionsGrid_->setColumnStretch(2, 1);
  optionsGrid_->setHorizontalSpacing(12);

  progress_ = new QLabel;
  progress_->setObjectName("runProgress");
  trace_ = new ChiSquaredTrace;
  trace_->hide();
  QHBoxLayout *progressRow = new QHBoxLayout;
  progressRow->addWidget(progress_, 1);
  progressRow->addWidget(trace_);

  log_ = new QPlainTextEdit;
  log_->setObjectName("runLog");
  log_->setReadOnly(true);
  log_->setMaximumBlockCount(20000);
  log_->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
  log_->setLineWrapMode(QPlainTextEdit::NoWrap);
  new LogHighlighter(log_->document());

  total_ = new QLabel(tr("No run yet."));
  total_->setObjectName("chiTotal");
  total_->setWordWrap(true);
  total_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  chi_ = new QTableWidget(0, 5);
  chi_->setObjectName("chiTable");
  chi_->setHorizontalHeaderLabels({tr("#"), tr("Data set"), tr("N"), tr("Chi-squared"), tr("Chi-squared/N")});
  chi_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  chi_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
  chi_->verticalHeader()->setVisible(false);
  chi_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  chi_->setSelectionBehavior(QAbstractItemView::SelectRows);
  chi_->setSortingEnabled(true);
  chi_->setToolTip(tr("The chi-squared of each data set of the last run (the program's \"Segment #\" lines); click a "
                      "column's title to sort by it."));
  plotButton_ = new QPushButton(QIcon(":/actions/plot.svg"), tr("Open &Plot"));
  plotButton_->setObjectName("openPlot");
  plotButton_->setEnabled(false);
  connect(plotButton_, &QPushButton::clicked, this, &RunTab::plotRequested);
  files_ = new QListWidget;
  files_->setObjectName("outputFiles");
  files_->setToolTip(tr("The files the run wrote in the output directory; double-click opens one."));
  connect(files_, &QListWidget::itemDoubleClicked, this, [this](QListWidgetItem *item) {
    QDesktopServices::openUrl(QUrl::fromLocalFile(QDir(outputDirectory_).filePath(item->text())));
  });
  QWidget *summary = new QWidget;
  QVBoxLayout *summaryLayout = new QVBoxLayout(summary);
  summaryLayout->setContentsMargins(0, 0, 0, 0);
  QHBoxLayout *totalRow = new QHBoxLayout;
  totalRow->addWidget(total_, 1);
  totalRow->addWidget(plotButton_);
  summaryLayout->addLayout(totalRow);
  summaryLayout->addWidget(chi_, 3);
  QLabel *filesTitle = new QLabel(tr("Files written:"));
  summaryLayout->addWidget(filesTitle);
  summaryLayout->addWidget(files_, 1);
  QSplitter *bottom = new QSplitter(Qt::Horizontal);
  bottom->addWidget(log_);
  bottom->addWidget(summary);
  bottom->setStretchFactor(0, 3);
  bottom->setStretchFactor(1, 2);

  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addLayout(calcLayout);
  layout->addWidget(optionsBox_);
  layout->addLayout(progressRow);
  layout->addWidget(bottom, 1);
  commandChanged();
}

QString RunTab::command() const { return command_->currentData().toString(); }

bool RunTab::setCommand(const QString &command) {
  const int i = command_->findData(command);
  if (i < 0) return false;
  command_->setCurrentIndex(i);
  return true;
}

int RunTab::threads() const { return threads_->value(); }

QString RunTab::valueOf(const QString &option) const {
  QWidget *w = editors_.value(option, nullptr);
  if (!w) return QString();
  if (QCheckBox *b = qobject_cast<QCheckBox *>(w)) return b->isChecked() ? "1" : "0";
  if (QComboBox *c = qobject_cast<QComboBox *>(w)) return c->currentData().toString();
  if (QSpinBox *s = qobject_cast<QSpinBox *>(w)) return s->value() == s->minimum() ? QString() : QString::number(s->value());
  if (QDoubleSpinBox *d = qobject_cast<QDoubleSpinBox *>(w)) {
    const double low = d->property("low").toDouble();
    const bool open = d->property("lowOpen").toBool();
    if (d->value() < low || (open && d->value() <= low)) return QString();
    return QString::number(d->value(), 'g', 10);
  }
  QLineEdit *e = qobject_cast<QLineEdit *>(w);
  if (!e) e = w->findChild<QLineEdit *>();
  return e ? e->text().trimmed() : QString();
}

void RunTab::commandChanged() {
  // keep what was typed, by option
  for (auto it = editors_.constBegin(); it != editors_.constEnd(); ++it) kept_[it.key()] = valueOf(it.key());
  for (QWidget *w : rowWidgets_) delete w;
  rowWidgets_.clear();
  rowOf_.clear();
  editors_.clear();
  needs_.clear();
  const QString command = this->command();
  commandHelp_->setText(command_->currentData(Qt::ToolTipRole).toString());
  int row = 0;
  for (const azure2::front::Option &o : azure2::front::Options()) {
    bool applies = o.commands.empty();
    for (const std::string &c : o.commands) applies = applies || QString::fromStdString(c) == command;
    if (!applies || o.option == "--max-calls=") continue;  // the call limit is the Fitting Settings tab's
    const QString option = QString::fromStdString(o.option);
    const QString unit = QString::fromStdString(o.unit);
    QString help = QString::fromStdString(o.help);
    QString name = option;
    name.remove('=');
    // the option's words ("--rate-entrance" -> "Rate entrance"); the option itself and its help in the tool tip
    QString label = name.mid(2).replace('-', ' ');
    label[0] = label[0].toUpper();
    const QString kept = kept_.value(option);
    QWidget *editor = nullptr;
    if (option == "--rate-entrance=" || option == "--rate-exit=") {
      // the pairs of the project in place of a typed key (the exit: -1, the total radiative capture, too)
      QComboBox *pairs = new QComboBox;
      pairs->addItem(tr("(choose a pair)"), QString());
      for (const QPair<int, QString> &p : pairs_) pairs->addItem(QString("%1: %2").arg(p.first).arg(p.second), QString::number(p.first));
      if (option == "--rate-exit=") pairs->addItem(tr("-1: total radiative capture of the entrance pair"), "-1");
      const int at = pairs->findData(kept);
      pairs->setCurrentIndex(at < 0 ? 0 : at);
      editor = pairs;
    } else if (o.kind == "flag") {
      QCheckBox *box = new QCheckBox;
      box->setChecked(kept == "1");
      connect(box, &QCheckBox::toggled, this, &RunTab::updateDependents);
      editor = box;
    } else if (o.kind == "word") {
      QComboBox *words = new QComboBox;
      if (!o.required) words->addItem(tr("(default)"), QString());
      for (const std::string &w : o.words) words->addItem(QString::fromStdString(w), QString::fromStdString(w));
      const int at = words->findData(kept);
      words->setCurrentIndex(at < 0 ? 0 : at);
      editor = words;
    } else if (o.kind == "whole") {
      QSpinBox *spin = new QSpinBox;
      const double low = std::max(Low(o), -1e9), high = std::min(High(o), 2e9);
      spin->setRange((int)low - 1, (int)high);
      spin->setSpecialValueText(o.required ? tr("required") : tr("default %1").arg(o.defaultValue));
      if (!unit.isEmpty()) spin->setSuffix(" " + unit);
      spin->setValue(kept.isEmpty() ? spin->minimum() : kept.toInt());
      editor = spin;
    } else if (o.kind == "number") {
      NumberBox *spin = new NumberBox;
      const double low = Low(o), high = High(o);
      spin->setDecimals(10);
      spin->setRange(low - 1, high);
      spin->setSingleStep(high - low > 100 ? 1 : 0.1);
      spin->setProperty("low", low);
      spin->setProperty("lowOpen", o.lowOpen);
      spin->setSpecialValueText(o.required ? tr("required") : tr("default %1").arg(o.defaultValue));
      if (!unit.isEmpty()) spin->setSuffix(" " + unit);
      spin->setValue(kept.isEmpty() ? spin->minimum() : kept.toDouble());
      spin->setMinimumWidth(fontMetrics().horizontalAdvance("M") * 9);
      editor = spin;
    } else {
      QLineEdit *edit = new QLineEdit(kept);
      edit->setPlaceholderText(o.required ? tr("required") : tr("none"));
      QWidget *host = new QWidget;
      QHBoxLayout *h = new QHBoxLayout(host);
      h->setContentsMargins(0, 0, 0, 0);
      h->addWidget(edit, 1);
      QToolButton *browse = new QToolButton;
      browse->setIcon(QIcon(":/actions/open.svg"));
      browse->setToolTip(tr("Choose the file"));
      connect(browse, &QToolButton::clicked, this, [this, edit]() {
        const QString f = QFileDialog::getOpenFileName(this, tr("Choose a file"));
        if (!f.isEmpty()) edit->setText(f);
      });
      h->addWidget(browse);
      host->setMinimumWidth(fontMetrics().horizontalAdvance("M") * 24);
      editor = host;
    }
    if (!o.when.empty() && !help.contains(QString::fromStdString(o.when))) help += " (" + QString::fromStdString(o.when) + ")";
    if (!o.instead.empty()) help += tr("; or %1").arg(QString::fromStdString(o.instead));
    editor->setObjectName("option_" + name.mid(2));
    editor->setToolTip(name + ": " + help);
    QLabel *labelWidget = new QLabel(label + ":");
    labelWidget->setToolTip(name + ": " + help);
    labelWidget->setBuddy(editor);
    QLabel *hint = new QLabel(help);
    hint->setProperty("role", "hint");
    hint->setWordWrap(true);
    hint->setObjectName("hint_" + name.mid(2));
    optionsGrid_->addWidget(labelWidget, row, 0, Qt::AlignRight | Qt::AlignVCenter);
    optionsGrid_->addWidget(editor, row, 1, Qt::AlignLeft | Qt::AlignVCenter);
    optionsGrid_->addWidget(hint, row, 2);
    rowWidgets_ << labelWidget << editor << hint;
    rowOf_.insert(option, {labelWidget, editor, hint});
    editors_.insert(option, editor);
    if (!o.needs.empty()) needs_.insert(option, QString::fromStdString(o.needs));
    row++;
  }
  optionsBox_->setVisible(!editors_.isEmpty());
  updateDependents();
  updateExpertRows();
}

namespace {

/// The options a casual user meets (qt_gui.md section 19.4): the start file, the band, the rate's pairs and
/// temperatures, MINOS's chi^2 rise and the sampler's walkers and steps.  The others are expert options.
bool StandardOption(const QString &option) {
  const char *const kStandard[] = {"--parameter-file=", "--covariance-band", "--rate-entrance=", "--rate-exit=",
                                          "--temperature-file=", "--temperature-min=", "--temperature-max=",
                                          "--temperature-step=", "--delta-chi2=", "--walkers=", "--steps="};
  for (const char *o : kStandard)
    if (option == o) return true;
  return false;
}

}  // namespace

void RunTab::setExpertMode(bool on) {
  expert_ = on;
  updateExpertRows();
}

bool RunTab::optionShown(const QString &option) const {
  const QList<QWidget *> row = rowOf_.value(option);
  return !row.isEmpty() && !row[1]->isHidden();
}

void RunTab::updateExpertRows() {
  // the commands of an expert (MINOS; the sampler, which its own tab runs) out of the drop-down in standard view,
  // unless chosen
  if (QListView *list = qobject_cast<QListView *>(command_->view()))
    for (int i = 0; i < command_->count(); i++) {
      const QString c = command_->itemData(i).toString();
      list->setRowHidden(i, !expert_ && i != command_->currentIndex() && (c == "minos" || c == "mcmc"));
    }
  for (auto it = rowOf_.constBegin(); it != rowOf_.constEnd(); ++it) {
    const QString value = valueOf(it.key());
    const bool set = !value.isEmpty() && !(qobject_cast<QCheckBox *>(it.value()[1]) && value == "0");
    const bool shown = expert_ || StandardOption(it.key()) || set;  // a set option is shown in either view
    it.value()[0]->setVisible(shown);
    it.value()[1]->setVisible(shown);
    it.value()[2]->setVisible(shown && expert_);  // the help line: the tool tip in standard view
    QWidget *editor = it.value()[1];
    if (!editor->property("baseTip").isValid()) editor->setProperty("baseTip", editor->toolTip());
    const bool marked = shown && !expert_ && !StandardOption(it.key());
    editor->setToolTip(editor->property("baseTip").toString() + (marked ? "\n" + expert::ShownBecauseSet() : QString()));
  }
}

void RunTab::updateDependents() {
  for (auto it = needs_.constBegin(); it != needs_.constEnd(); ++it) {
    QCheckBox *flag = qobject_cast<QCheckBox *>(editors_.value(it.value(), nullptr));
    if (QWidget *w = editors_.value(it.key(), nullptr)) w->setEnabled(!flag || flag->isChecked());
  }
}

QStringList RunTab::options() const {
  QStringList out;
  for (auto it = editors_.constBegin(); it != editors_.constEnd(); ++it) {
    const QString value = valueOf(it.key());
    if (qobject_cast<QCheckBox *>(it.value())) {
      if (value == "1") out << it.key();
      continue;
    }
    if (!value.isEmpty()) out << it.key() + value;
  }
  return out;
}

bool RunTab::setOption(const QString &option, const QString &value) {
  QWidget *w = editors_.value(option, nullptr);
  if (!w) return false;
  if (QCheckBox *b = qobject_cast<QCheckBox *>(w)) {
    b->setChecked(value == "1");
    return true;
  }
  if (QComboBox *c = qobject_cast<QComboBox *>(w)) {
    const int at = c->findData(value);
    if (at < 0) return false;
    c->setCurrentIndex(at);
    return true;
  }
  if (QSpinBox *s = qobject_cast<QSpinBox *>(w)) {
    bool ok = true;
    const int v = value.isEmpty() ? s->minimum() : value.toInt(&ok);
    if (!ok || v < s->minimum() || v > s->maximum()) return false;
    s->setValue(v);
    return true;
  }
  if (QDoubleSpinBox *d = qobject_cast<QDoubleSpinBox *>(w)) {
    bool ok = true;
    const double v = value.isEmpty() ? d->minimum() : value.toDouble(&ok);
    if (!ok || v < d->minimum() || v > d->maximum()) return false;
    d->setValue(v);
    return true;
  }
  QLineEdit *e = qobject_cast<QLineEdit *>(w);
  if (!e) e = w->findChild<QLineEdit *>();
  if (!e) return false;
  e->setText(value);
  return true;
}

void RunTab::setPairs(const QList<QPair<int, QString>> &pairs) {
  if (pairs == pairs_) return;
  pairs_ = pairs;
  commandChanged();  // the drop-downs again, the choices kept by key
}

void RunTab::setDataSetNames(const QStringList &names) { dataSetNames_ = names; }

void RunTab::presetOption(const QString &option, const QString &value) {
  kept_[option] = value;
  setOption(option, value);
}

void RunTab::setRunning(bool running) {
  runButton_->setEnabled(!running);
  stopButton_->setEnabled(running);
  command_->setEnabled(!running);
  optionsBox_->setEnabled(!running);
  threads_->setEnabled(!running);
  if (running) updateDependents();
}

QString RunTab::log() const { return log_->toPlainText(); }

void RunTab::clearLog() {
  log_->clear();
  progress_->clear();
  trace_->clear();
  trace_->hide();
  chi_->setRowCount(0);
  total_->setText(tr("No run yet."));
  files_->clear();
  plotButton_->setEnabled(false);
}

void RunTab::appendLine(const QString &line, bool transient) {
  if (transient) {
    progress_->setText(line.trimmed());
    return;
  }
  log_->appendPlainText(line);
}

void RunTab::showProgress(const QString &stage, qlonglong evaluation, qlonglong calls, double objective, double best) {
  if (stage == "mcmc") {
    progress_->setText(tr("mcmc: step %1 of %2, log posterior %3").arg(evaluation).arg(calls).arg(objective, 0, 'g', 10));
    return;
  }
  QString text = tr("%1: evaluation %2, chi-squared %3").arg(stage).arg(evaluation).arg(objective, 0, 'g', 10);
  if (!std::isnan(best)) text += tr(" (best %1)").arg(best, 0, 'g', 10);
  progress_->setText(text);
  if (stage == "fit" && std::isfinite(best)) {
    trace_->add((double)evaluation, best);
    trace_->setVisible(trace_->points() > 1);
  }
}

void RunTab::showOutcome(const RunOutcome &outcome) {
  outputDirectory_ = outcome.outputDirectory;
  chi_->setSortingEnabled(false);
  chi_->setRowCount(0);
  long totalPoints = 0;
  for (int i = 0; i < outcome.dataSetChiSquared.size(); i++) {
    const int row = chi_->rowCount();
    chi_->insertRow(row);
    const int segment = outcome.dataSetChiSquared[i].first;
    const double chi2 = outcome.dataSetChiSquared[i].second;
    const int n = outcome.dataSetPoints.value(i);
    totalPoints += std::max(0, n);
    const QString name = outcome.dataSetChiSquared.size() == dataSetNames_.size() ? dataSetNames_.value(i) : QString();
    chi_->setItem(row, 0, new NumberItem(QString::number(segment), segment));
    chi_->setItem(row, 1, new QTableWidgetItem(name));
    chi_->setItem(row, 2, new NumberItem(n > 0 ? QString::number(n) : QString(), n));
    chi_->setItem(row, 3, new NumberItem(QString::number(chi2, 'g', 10), chi2));
    chi_->setItem(row, 4, n > 0 ? new NumberItem(QString::number(chi2 / n, 'g', 6), chi2 / n) : new NumberItem(QString(), NAN));
  }
  chi_->setSortingEnabled(true);
  chi_->sortByColumn(0, Qt::AscendingOrder);
  if (outcome.haveChiSquared) {
    QString text = tr("Total chi-squared: <b>%1</b>").arg(outcome.totalChiSquared, 0, 'g', 12);
    if (totalPoints > 0)
      text += tr(" over %1 points (chi-squared/N = %2)").arg(totalPoints).arg(outcome.totalChiSquared / totalPoints, 0, 'g', 6);
    total_->setText(text);
  } else {
    total_->setText(outcome.exitCode == 0 ? tr("The run wrote no chi-squared.") : tr("No result: the run failed."));
  }
  files_->clear();
  files_->addItems(outcome.outputFiles);
  plotButton_->setEnabled(!outcome.outputFiles.isEmpty());
  QString end = outcome.exitCode == 0 ? tr("The run ended.") : tr("The run failed (exit %1).").arg(outcome.exitCode);
  if (outcome.stopped) end += tr("  It was stopped on request.");
  progress_->setText(end);
}
