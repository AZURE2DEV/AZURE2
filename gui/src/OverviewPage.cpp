// The Overview page (OverviewPage.h; docs/upgrade/qt_gui.md section 15).
#include "OverviewPage.h"

#include <QApplication>
#include <QDateTime>
#include <QDir>
#include <QEvent>
#include <QFileInfo>
#include <QFrame>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPainter>
#include <QPushButton>
#include <QScrollArea>
#include <QStackedWidget>
#include <QStyle>
#include <QStyledItemDelegate>
#include <QTableWidget>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <limits>

#include "AzureTheme.h"
#include "OverviewCharts.h"
#include "ProjectHolder.h"
#include "RunTab.h"

namespace front = azure2::front;

/// The parameters' worker's result, handed from its thread to the page's when the thread ends.
struct OverviewPage::Job {
  std::string projectFile, parameterFile;
  bool fit = false;
  QString key;
  bool ok = false;
  front::RunParameters out;
  std::string why;
};

namespace {

/// chi^2/N above it is named in "What to look at next" and drawn in the warning colour.
constexpr double kHighChiPerPoint = 3.0;
/// |pull| of a norm above it is named.
constexpr double kHighPull = 2.0;
constexpr int kColumnChiPerPoint = 4;

QFrame *Card(const QString &title, QVBoxLayout *&body, const QString &name) {
  QFrame *card = new QFrame;
  card->setObjectName(name);
  card->setProperty("role", "card");
  body = new QVBoxLayout(card);
  body->setSpacing(6);
  QLabel *heading = new QLabel(title);
  QFont f = heading->font();
  f.setPointSizeF(f.pointSizeF() * 1.12);
  f.setBold(true);
  heading->setFont(f);
  heading->setProperty("role", "title");
  body->addWidget(heading);
  return card;
}

QLabel *Text(const QString &name = QString()) {
  QLabel *l = new QLabel;
  l->setObjectName(name);
  l->setWordWrap(true);
  l->setTextFormat(Qt::RichText);
  l->setTextInteractionFlags(Qt::TextBrowserInteraction);
  return l;
}

QString Number(double x, int digits = 6) { return std::isfinite(x) ? QString::number(x, 'g', digits) : QString("--"); }

QString Muted(const QString &text) {
  return QString("<span style='color:%1'>%2</span>")
      .arg(QApplication::palette().color(QPalette::PlaceholderText).name(), text);
}

QString Warn(const QString &text) {
  return QString("<span style='color:%1; font-weight:bold'>%2</span>").arg(BarChart::warningColor().name(), text);
}

QString Link(const QString &href, const QString &text) {
  // the colour of the look in use: a rich text keeps the colour it was made with
  return QString("<a href=\"%1\" style='color:%2'>%3</a>")
      .arg(href, QApplication::palette().color(QPalette::Link).name(), text);
}

/// A table item that sorts by its number.
class NumberItem : public QTableWidgetItem {
 public:
  NumberItem(const QString &text, double value) : QTableWidgetItem(text) {
    setData(Qt::UserRole, value);
    setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
  }
  bool operator<(const QTableWidgetItem &other) const override {
    const double a = data(Qt::UserRole).toDouble(), b = other.data(Qt::UserRole).toDouble();
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && !std::isnan(b);
    return a < b;
  }
};

/// The chi^2/N column: a small bar left of the number, never under it, its length chi^2/N over the column's scale
/// (BarChart::ScaleTop: the largest, or ten times the median when one value is far above the others; a bar beyond
/// the scale is drawn full with a break).  The bar takes only the room the number leaves (qt_gui.md section 19.9).
class RatioBarDelegate : public QStyledItemDelegate {
 public:
  using QStyledItemDelegate::QStyledItemDelegate;
  double top = 1;
  void paint(QPainter *p, const QStyleOptionViewItem &option, const QModelIndex &index) const override {
    QStyledItemDelegate::paint(p, option, index);
    const double v = index.data(Qt::UserRole).toDouble();
    if (!std::isfinite(v) || top <= 0) return;
    const QRect cell = option.rect.adjusted(4, 0, -4, 0);
    const int text = option.fontMetrics.horizontalAdvance(index.data(Qt::DisplayRole).toString());
    const int room = cell.width() - text - 8;  // the bar's room: what the right-aligned number leaves
    if (room < 10) return;
    p->save();
    QColor c = v > kHighChiPerPoint ? BarChart::warningColor() : option.palette.color(QPalette::Highlight);
    c.setAlpha(170);
    const int h = std::max(4, cell.height() / 3);
    const QRect r(cell.left(), cell.center().y() - h / 2, std::max(2, (int)std::lround(room * std::min(1.0, v / top))), h);
    p->setPen(Qt::NoPen);
    p->setBrush(c);
    p->drawRect(r);
    if (v > top) {
      p->setPen(QPen(option.palette.color(QPalette::Base), 2));
      p->drawLine(QPoint(r.right() - 7, r.bottom() + 1), QPoint(r.right() - 3, r.top() - 1));
    }
    p->restore();
  }
};

QWidget *Tile(const QString &caption, QLabel *&value, const QString &name) {
  QFrame *box = new QFrame;
  box->setObjectName("tile_" + name);
  QVBoxLayout *l = new QVBoxLayout(box);
  l->setContentsMargins(6, 2, 6, 2);
  l->setSpacing(0);
  value = new QLabel("--");
  value->setObjectName("tileValue_" + name);
  QFont f = value->font();
  f.setPointSizeF(f.pointSizeF() * 1.35);
  f.setBold(true);
  value->setFont(f);
  value->setTextInteractionFlags(Qt::TextSelectableByMouse);
  QLabel *c = new QLabel(caption);
  c->setProperty("role", "hint");
  l->addWidget(value);
  l->addWidget(c);
  return box;
}

}  // namespace

OverviewPage::OverviewPage(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  setObjectName("overviewPage");
  refreshTimer_ = new QTimer(this);
  refreshTimer_->setSingleShot(true);
  refreshTimer_->setInterval(200);
  connect(refreshTimer_, &QTimer::timeout, this, &OverviewPage::refresh);
  build();
}

OverviewPage::~OverviewPage() {
  // a set-up in progress ends before the page goes (a running QThread may not be destroyed)
  if (worker_) worker_->wait();
}

void OverviewPage::build() {
  QVBoxLayout *outer = new QVBoxLayout(this);
  outer->setContentsMargins(0, 0, 0, 0);
  scroll_ = new QScrollArea;
  scroll_->setWidgetResizable(true);
  scroll_->setFrameShape(QFrame::NoFrame);
  QWidget *content = new QWidget;
  content->setObjectName("overviewContent");
  grid_ = new QGridLayout(content);
  grid_->setContentsMargins(10, 8, 10, 8);
  grid_->setSpacing(10);
  scroll_->setWidget(content);
  outer->addWidget(scroll_);
  auto linkTo = [this](QLabel *l) { connect(l, &QLabel::linkActivated, this, &OverviewPage::followLink); };

  // What to look at next, and the quick actions
  QVBoxLayout *b;
  stepsCard_ = Card(tr("What to look at next"), b, "overviewSteps");
  stepsText_ = Text("overviewStepsText");
  linkTo(stepsText_);
  b->addWidget(stepsText_);
  QHBoxLayout *actions = new QHBoxLayout;
  const QList<QPair<QString, QString>> buttons = {{"validate", tr("Validate")},
                                                  {"run:calculate", tr("Calculate")},
                                                  {"run:fit", tr("Fit")},
                                                  {"plot:", tr("Open Plot")}};
  for (const auto &button : buttons) {
    QPushButton *pb = new QPushButton(button.second);
    pb->setObjectName("overviewAction_" + QString(button.first).remove(':'));
    const QString link = button.first;
    connect(pb, &QPushButton::clicked, this, [this, link]() { followLink(link); });
    actions->addWidget(pb);
  }
  actions->addStretch(1);
  b->addLayout(actions);

  // The project
  projectCard_ = Card(tr("The project"), b, "overviewProject");
  projectHead_ = Text("overviewProjectHead");
  b->addWidget(projectHead_);
  levelsTitle_ = Text("overviewLevelsTitle");
  b->addWidget(levelsTitle_);
  levelBars_ = new BarChart(true);
  levelBars_->setObjectName("overviewLevelBars");
  levelBars_->setEmptyText(tr("No level yet"));
  b->addWidget(levelBars_);
  projectText_ = Text("overviewProjectText");
  linkTo(projectText_);
  b->addWidget(projectText_);
  validity_ = new QPushButton;
  validity_->setObjectName("overviewValidity");
  validity_->setFlat(true);
  validity_->setCursor(Qt::PointingHandCursor);
  validity_->setToolTip(tr("Show the Diagnostics"));
  connect(validity_, &QPushButton::clicked, this, [this]() { followLink("diagnostics"); });
  QHBoxLayout *validityRow = new QHBoxLayout;
  validityRow->addWidget(validity_);
  validityRow->addStretch(1);
  b->addLayout(validityRow);
  projectMore_ = Text("overviewProjectNotes");
  b->addWidget(projectMore_);
  b->addStretch(1);

  // The last run
  runCard_ = Card(tr("The last run"), b, "overviewRun");
  runStack_ = new QStackedWidget;
  QWidget *none = new QWidget;
  none->setObjectName("overviewNoRun");
  QVBoxLayout *nl = new QVBoxLayout(none);
  QLabel *noneText = Text();
  noneText->setText(tr("<b>No run yet.</b>"));
  noneText->setToolTip(tr("The output directory holds no results of this project: calculate the data sets from the "
                          "parameters as they stand, or fit them."));
  nl->addWidget(noneText);
  QHBoxLayout *nb = new QHBoxLayout;
  for (const auto &button : {qMakePair(QString("calculate"), tr("Calculate")), qMakePair(QString("fit"), tr("Fit"))}) {
    QPushButton *pb = new QPushButton(button.second);
    pb->setObjectName("overviewNoRun_" + button.first);
    pb->setProperty("primary", true);
    const QString command = button.first;
    connect(pb, &QPushButton::clicked, this, [this, command]() { emit commandRequested(command); });
    nb->addWidget(pb);
  }
  nb->addStretch(1);
  nl->addLayout(nb);
  nl->addStretch(1);
  QWidget *have = new QWidget;
  QVBoxLayout *hl = new QVBoxLayout(have);
  hl->setContentsMargins(0, 0, 0, 0);
  runHead_ = Text("overviewRunHead");
  hl->addWidget(runHead_);
  stale_ = new QLabel;
  stale_->setObjectName("overviewStale");
  stale_->setWordWrap(true);
  hl->addWidget(stale_);
  runStatus_ = Text("overviewRunStatus");
  hl->addWidget(runStatus_);
  QHBoxLayout *tiles = new QHBoxLayout;
  tiles->setSpacing(4);
  const QList<QPair<QString, QString>> tileNames = {{"chi2", tr("total chi² (data)")},
                                                    {"points", tr("points N")},
                                                    {"free", tr("free parameters p")},
                                                    {"dof", tr("chi² / (N - p)")}};
  for (const auto &t : tileNames) {
    QLabel *value = nullptr;
    tiles->addWidget(Tile(t.second, value, t.first));
    tiles_.insert(t.first, value);
  }
  tiles->addStretch(1);
  hl->addLayout(tiles);
  totalsText_ = Text("overviewTotals");
  hl->addWidget(totalsText_);
  table_ = new QTableWidget(0, 7);
  table_->setObjectName("overviewDataSets");
  table_->setHorizontalHeaderLabels(
      {tr("Data set"), tr("Reaction"), tr("N"), tr("chi²"), tr("chi²/N"), tr("Norm"), tr("Pull")});
  table_->horizontalHeaderItem(5)->setToolTip(tr("The norm the run used"));
  table_->horizontalHeaderItem(6)->setToolTip(
      tr("(norm - prior centre) / prior width, from the norm's prior chi-squared; empty without a prior"));
  table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table_->setSelectionBehavior(QAbstractItemView::SelectRows);
  table_->setSelectionMode(QAbstractItemView::SingleSelection);
  table_->setAlternatingRowColors(true);
  table_->verticalHeader()->hide();
  table_->verticalHeader()->setDefaultSectionSize(table_->fontMetrics().height() + 8);
  // the reaction takes what is left (elided when long); the bar column a fixed width
  table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  table_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
  table_->horizontalHeader()->setSectionResizeMode(kColumnChiPerPoint, QHeaderView::Fixed);
  table_->horizontalHeader()->resizeSection(kColumnChiPerPoint, table_->fontMetrics().horizontalAdvance("M") * 11);
  table_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Interactive);  // the names: at most 18 characters, elided
  table_->setTextElideMode(Qt::ElideRight);
  table_->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  table_->setItemDelegateForColumn(kColumnChiPerPoint, new RatioBarDelegate(table_));
  table_->setToolTip(tr("Click a data set to open it in the Plot tab; click a header to sort."));
  table_->setMinimumHeight(table_->fontMetrics().height() * 8);
  connect(table_, &QTableWidget::cellClicked, this, [this](int row, int) {
    if (QTableWidgetItem *item = table_->item(row, 0)) emit plotRequested(item->data(Qt::UserRole + 1).toString());
  });
  hl->addWidget(table_, 1);
  flagsText_ = Text("overviewFlags");
  linkTo(flagsText_);
  hl->addWidget(flagsText_);
  notesText_ = Text("overviewRunNotes");
  hl->addWidget(notesText_);
  runStack_->addWidget(none);
  runStack_->addWidget(have);
  b->addWidget(runStack_, 1);

  // The charts
  chartsCard_ = Card(tr("Charts"), b, "overviewCharts");
  QLabel *barsTitle = new QLabel(tr("chi²/N per data set"));
  barsTitle->setProperty("role", "hint");
  b->addWidget(barsTitle);
  chiBars_ = new BarChart(false);
  chiBars_->setObjectName("overviewChiBars");
  chiBars_->setEmptyText(tr("No run yet"));
  b->addWidget(chiBars_, 1);
  traceTitle_ = new QLabel(tr("Best chi² against the evaluation (fit)"));
  QLabel *traceTitle = traceTitle_;
  traceTitle->setProperty("role", "hint");
  b->addWidget(traceTitle);
  trace_ = new ChiSquaredTrace;
  trace_->setObjectName("overviewTrace");
  trace_->setMinimumHeight(fontMetrics().height() * 6);
  b->addWidget(trace_, 1);
  traceNote_ = Text("overviewTraceNote");
  b->addWidget(traceNote_);

  reflow();
}

void OverviewPage::reflow() {
  const int w = scroll_->viewport()->width();
  const int columns = w >= 1250 ? 3 : w >= 820 ? 2 : 1;
  if (columns == columns_) return;
  columns_ = columns;
  for (QWidget *card : {stepsCard_, projectCard_, runCard_, chartsCard_}) grid_->removeWidget(card);
  for (int c = 0; c < 3; c++) grid_->setColumnStretch(c, 0);
  for (int r = 0; r < 4; r++) grid_->setRowStretch(r, 0);
  grid_->addWidget(stepsCard_, 0, 0, 1, columns);
  if (columns == 3) {
    grid_->addWidget(projectCard_, 1, 0);
    grid_->addWidget(runCard_, 1, 1);
    grid_->addWidget(chartsCard_, 1, 2);
    grid_->setColumnStretch(0, 4);
    grid_->setColumnStretch(1, 7);
    grid_->setColumnStretch(2, 3);
    grid_->setRowStretch(1, 1);
  } else if (columns == 2) {
    grid_->addWidget(projectCard_, 1, 0);
    grid_->addWidget(runCard_, 1, 1);
    grid_->addWidget(chartsCard_, 2, 0, 1, 2);
    grid_->setColumnStretch(0, 2);
    grid_->setColumnStretch(1, 3);
    grid_->setRowStretch(1, 1);
  } else {
    grid_->addWidget(projectCard_, 1, 0);
    grid_->addWidget(runCard_, 2, 0);
    grid_->addWidget(chartsCard_, 3, 0);
    grid_->setColumnStretch(0, 1);
  }
}

void OverviewPage::resizeEvent(QResizeEvent *e) {
  QWidget::resizeEvent(e);
  reflow();
}

void OverviewPage::showEvent(QShowEvent *e) {
  QWidget::showEvent(e);
  reflow();
  startParameters();
}

void OverviewPage::changeEvent(QEvent *e) {
  QWidget::changeEvent(e);
  // the theme changed: the rich texts carry colours of the palette
  if ((e->type() == QEvent::PaletteChange || e->type() == QEvent::ApplicationPaletteChange || e->type() == QEvent::StyleChange) &&
      !project_.file.empty()) {
    fillProject();
    fillRun();
    fillSteps();
  }
}

void OverviewPage::refreshLater() { refreshTimer_->start(); }

void OverviewPage::setRunning(bool running) {
  running_ = running;
  if (running) clearProgress();
  refresh();
}

void OverviewPage::setValidity(bool ok, int errors, int warnings) {
  validated_ = true;
  valid_ = ok;
  errors_ = errors;
  warnings_ = warnings;
  fillProject();
  fillSteps();
}

void OverviewPage::setFitterChoice(const QString &fitter) {
  fitter_ = fitter;
  fillProject();
}

void OverviewPage::addProgress(const QString &stage, qlonglong evaluation, double best) {
  if (stage != "fit" || !std::isfinite(best)) return;
  trace_->add((double)evaluation, best);
  trace_->show();
  traceTitle_->show();
  traceNote_->setText(tr("This window's fit: %n evaluation(s) recorded.", nullptr, trace_->points()));
}

void OverviewPage::clearProgress() {
  trace_->clear();
  trace_->hide();
  traceTitle_->hide();
  traceNote_->setText(Muted(tr("Drawn while a fit runs here.")));
  traceNote_->setToolTip(tr("Recorded while a fit runs in this window: the engine writes no file of a fit's progress."));
}

void OverviewPage::refresh() {
  refreshTimer_->stop();
  project_ = front::SummarizeProject(holder_->project());
  // the window's own run writes the files as it goes: they are read when it has ended
  if (!running_) run_ = front::ReadLastRun(holder_->project());
  fillProject();
  fillRun();
  fillCharts();
  fillSteps();
  if (isVisible()) startParameters();
}

void OverviewPage::startParameters() {
  if (running_ || project_.file.empty()) return;
  const QFileInfo file(holder_->path());
  if (!file.exists()) return;
  const QString parameterFile = QString::fromStdString(run_.parameterFile);
  const QString key = QString("%1|%2|%3|%4|%5")
                          .arg(file.absoluteFilePath())
                          .arg(file.lastModified().toMSecsSinceEpoch())
                          .arg(parameterFile)
                          .arg(parameterFile.isEmpty() ? 0 : QFileInfo(parameterFile).lastModified().toMSecsSinceEpoch())
                          .arg(QString::fromStdString(run_.command));
  if (key == parametersKey_) return;
  if (worker_) {
    pending_ = true;
    return;
  }
  job_ = std::make_shared<Job>();
  job_->projectFile = file.absoluteFilePath().toStdString();
  job_->parameterFile = run_.parameterFile;
  job_->fit = run_.command == "fit";
  job_->key = key;
  std::shared_ptr<Job> job = job_;
  worker_ = QThread::create([job]() {
    job->ok = front::ReadRunParameters(job->projectFile, job->parameterFile, job->fit, job->out, job->why);
  });
  worker_->setObjectName("overviewParameters");
  connect(worker_, &QThread::finished, this, &OverviewPage::parametersDone);
  if (!haveParameters_) tiles_["free"]->setText("...");
  worker_->start();
}

void OverviewPage::parametersDone() {
  worker_->deleteLater();
  worker_ = nullptr;
  parametersKey_ = job_->key;
  haveParameters_ = job_->ok;
  parameters_ = job_->ok ? job_->out : front::RunParameters();
  parametersWhy_ = QString::fromStdString(job_->why);
  job_.reset();
  fillProject();
  fillRun();
  fillSteps();
  if (pending_) {
    pending_ = false;
    startParameters();
  }
}

void OverviewPage::fillProject() {
  const QFileInfo file(QString::fromStdString(project_.file));
  // the file's name; its folder in the tip (qt_gui.md section 19.6)
  QString head = QString("<b>%1</b>").arg(file.fileName().toHtmlEscaped());
  projectHead_->setToolTip(QDir::toNativeSeparators(file.absoluteFilePath()));
  if (holder_->project().Modified()) head += "<br>" + Warn(tr("modified, not saved"));
  QStringList nuclei;
  for (const std::string &n : project_.nuclei) nuclei << QString::fromStdString(n);
  head += "<br>" + tr("Compound nucleus: <b>%1</b>").arg(nuclei.isEmpty() ? tr("none yet") : nuclei.join(", "));
  // the pairs by name; each threshold (E_x of the compound nucleus) in the tip
  QStringList pairs, thresholds;
  for (const front::PairChoice &p : project_.pairs) {
    pairs << QString::fromStdString(p.label).toHtmlEscaped();
    thresholds << tr("%1: threshold at E_x = %2 MeV").arg(QString::fromStdString(p.label)).arg(p.threshold, 0, 'f', 3);
  }
  head += "<br>" + tr("Pairs: %1").arg(pairs.isEmpty() ? tr("none yet") : pairs.join(" &middot; "));
  if (!thresholds.isEmpty()) projectHead_->setToolTip(projectHead_->toolTip() + "\n" + thresholds.join("\n"));
  projectHead_->setText(head);

  levelsTitle_->setText(tr("<b>%n level(s)</b>", nullptr, project_.levels) +
                        tr(" in %n J<sup>&pi;</sup> group(s)", nullptr, (int)project_.levelGroups.size()));
  QStringList labels;
  QVector<double> counts;
  const bool several = project_.nuclei.size() > 1;
  for (const front::LevelGroupCount &g : project_.levelGroups) {
    labels << (several ? QString::fromStdString(g.nucleus) + " " : QString()) + QString::fromStdString(g.label);
    counts << g.levels;
  }
  levelBars_->setBars(labels, counts);

  QStringList lines;
  QString data = tr("Data sets: <b>%1</b> of %2 active").arg(project_.activeDataSets).arg(project_.dataSets);
  if (run_.haveChiSquared && !running_) data += tr(", %n point(s) in the last run", nullptr, run_.points);
  if (project_.extrapolations > 0)
    data += tr("; extrapolations: %1 of %2 active").arg(project_.activeExtrapolations).arg(project_.extrapolations);
  lines << data;
  QStringList effects;
  for (const front::RecordUse &u : project_.effects)
    effects << tr("%1 %2 (%3 used)").arg(QString::fromStdString(u.kind)).arg(u.records).arg(u.used);
  lines << tr("Experimental effects: %1").arg(effects.isEmpty() ? tr("none") : effects.join(", "));
  QString free;
  if (haveParameters_) {
    QStringList kinds;
    for (const auto &k : parameters_.freeByKind) kinds << QString("%1 %2").arg(QString::fromStdString(k.first)).arg(k.second);
    free = tr("<b>%1</b> of %2").arg(parameters_.free).arg(parameters_.parameters) +
           (kinds.isEmpty() ? QString() : " (" + kinds.join(", ") + ")");
  } else if (!parametersWhy_.isEmpty()) {
    {
      // the first line of the reason; the whole text is in the Diagnostics panel and in this line's tool tip
      const QString first = parametersWhy_.section('\n', 0, 0).left(160);
      free = Muted(tr("not known: %1%2").arg(first.toHtmlEscaped(), first.size() < parametersWhy_.size() ? tr(" ... (see Validate)") : QString()));
    }
  } else {
    free = Muted(QFileInfo::exists(holder_->path()) ? tr("being counted...") : tr("save the project to count them"));
  }
  lines << tr("Free parameters: %1").arg(free);
  QStringList ec;
  for (int key : project_.externalCapturePairs) ec << QString::number(key);
  QString basis = tr("Basis: <b>%1</b>, %2; external capture %3")
                      .arg(QString::fromStdString(project_.basis), QString::fromStdString(project_.matrix),
                           ec.isEmpty() ? tr("off") : tr("on (pair %1)").arg(ec.join(", ")));
  if (!fitter_.isEmpty()) basis += tr("; fitter %1").arg(fitter_.toHtmlEscaped());
  lines << basis;
  QStringList set;
  for (const front::OptionState &o : project_.options)
    if (o.on != o.defaultOn)
      set << QString("%1: %2").arg(QString::fromStdString(o.title).toHtmlEscaped(), o.on ? tr("on") : tr("off"));
  if (!set.isEmpty()) lines << tr("Options set away from their default: %1").arg(set.join("; "));
  projectText_->setText(lines.join("<br>"));

  QString validity;
  if (!validated_) validity = tr("Not validated");
  else if (!valid_ || errors_ > 0) validity = tr("%n error(s)", nullptr, errors_);
  else validity = tr("Valid");
  if (validated_ && warnings_ > 0) validity += tr(", %n warning(s)", nullptr, warnings_);
  validity_->setText(validity);
  validity_->setIcon(style()->standardIcon(!validated_                    ? QStyle::SP_MessageBoxQuestion
                                           : (!valid_ || errors_ > 0)     ? QStyle::SP_MessageBoxCritical
                                           : warnings_ > 0                ? QStyle::SP_MessageBoxWarning
                                                                          : QStyle::SP_DialogApplyButton));
  QStringList notes;
  for (const std::string &n : project_.notes) notes << QString::fromStdString(n).toHtmlEscaped();
  projectMore_->setText(notes.isEmpty() ? QString() : Muted(notes.join("<br>")));
  projectMore_->setVisible(!notes.isEmpty());
}

void OverviewPage::fillRun() {
  if (!run_.present) {
    runStack_->setCurrentIndex(0);
    return;
  }
  runStack_->setCurrentIndex(1);
  const QString command = QString::fromStdString(run_.command);
  const QString when = QDateTime::fromSecsSinceEpoch(run_.runTime).toString("yyyy-MM-dd hh:mm:ss");
  QString head = tr("<b>%1</b> &middot; %2 %3")
                     .arg(command.left(1).toUpper() + command.mid(1), when, Muted(tr("(inferred from the files)")));
  if (running_) head += "<br>" + Muted(tr("A run is in progress: its results are shown when it ends."));
  runHead_->setText(head);
  stale_->setVisible(run_.stale);
  QStringList staleWhy;
  if (run_.projectNewer)
    staleWhy << tr("The project changed since this run (the file was saved at %1).")
                    .arg(QDateTime::fromSecsSinceEpoch(run_.projectTime).toString("yyyy-MM-dd hh:mm:ss"));
  if (!run_.changedInputs.empty()) {
    QStringList names;
    for (const std::string &n : run_.changedInputs) names << QString::fromStdString(n);
    staleWhy << tr("The data of %1 changed since this run.").arg(names.join(", "));
  }
  stale_->setText(QString::fromUtf8("⚠ ") + staleWhy.join(" ") + tr(" Calculate again."));
  stale_->setStyleSheet(QString("color: %1; font-weight: bold;").arg(BarChart::warningColor().name()));
  if (!run_.status.empty()) {
    QString status = tr("Fit status: <b>%1</b>").arg(QString::fromStdString(run_.status).toHtmlEscaped());
    if (!run_.minimizer.empty()) status += tr(" (%1, %n call(s))", nullptr, (int)run_.calls).arg(QString::fromStdString(run_.minimizer));
    runStatus_->setText(status);
  } else {
    runStatus_->setText(command == "calculate" ? Muted(tr("A calculation (no fit)."))
                                               : QString());
  }
  runStatus_->setVisible(!runStatus_->text().isEmpty());

  tiles_["chi2"]->setText(run_.haveChiSquared ? Number(run_.totalChiSquared, 8) : QString("--"));
  tiles_["points"]->setText(run_.haveChiSquared ? QString::number(run_.points) : QString("--"));
  tiles_["free"]->setText(haveParameters_ ? QString::number(parameters_.free) : QString("..."));
  const double dof = haveParameters_ ? front::ChiSquaredPerDegreeOfFreedom(run_, parameters_)
                                     : std::numeric_limits<double>::quiet_NaN();
  tiles_["dof"]->setText(run_.haveChiSquared && haveParameters_ ? Number(dof, 4) : QString("--"));
  tiles_["dof"]->setStyleSheet(std::isfinite(dof) && dof > kHighChiPerPoint
                                   ? QString("color: %1;").arg(BarChart::warningColor().name())
                                   : QString());
  QStringList totals, labels;
  for (const auto &t : run_.totals)
    if (t.first != "Total-Chi-Squared" && t.first != "Total-N" && t.second != 0) {
      // the file's label in words ("Total-Norm-Chi-Squared" -> "norm chi²"), the label in the tip; a term that adds
      // nothing is not listed (qt_gui.md section 19.8)
      QString name = QString::fromStdString(t.first);
      const QString label = name;
      if (name.startsWith("Total-")) name.remove(0, 6);
      if (name.endsWith("-Chi-Squared")) name = name.chopped(12) + QString::fromUtf8(" chi²");
      name.replace('-', ' ');
      totals << QString("%1 %2").arg(name.toHtmlEscaped(), Number(t.second));
      labels << label;
    }
  totalsText_->setText(totals.isEmpty() ? QString() : Muted(totals.join(" &middot; ")));
  totalsText_->setVisible(!totals.isEmpty());
  totalsText_->setToolTip(labels.isEmpty() ? QString() : tr("chiSquared.out: %1").arg(labels.join(", ")));

  // the table: one row per data set of chiSquared.out
  table_->setSortingEnabled(false);
  table_->setRowCount((int)run_.dataSets.size());
  double top = 0;
  QVector<double> ratios;
  for (const front::DataSetResult &d : run_.dataSets) {
    ratios << d.chiSquaredPerPoint;
    if (std::isfinite(d.chiSquaredPerPoint)) top = std::max(top, d.chiSquaredPerPoint);
  }
  QVector<double> sorted = ratios;
  std::sort(sorted.begin(), sorted.end(), [](double a, double b) { return a > b; });
  const double third = sorted.size() >= 3 ? sorted[2] : sorted.isEmpty() ? 0 : sorted.last();
  bool capped = false;
  const double scale = BarChart::ScaleTop(ratios, &capped);
  static_cast<RatioBarDelegate *>(table_->itemDelegateForColumn(kColumnChiPerPoint))->top = scale;
  if (QTableWidgetItem *head = table_->horizontalHeaderItem(kColumnChiPerPoint))
    head->setToolTip(capped ? tr("chi²/N; the bars on a scale cut at ten times the median (%1): a bar with a break is "
                                 "longer").arg(Number(scale, 3))
                            : tr("chi²/N; the bars on a linear scale up to the largest"));
  for (int i = 0; i < (int)run_.dataSets.size(); i++) {
    const front::DataSetResult &d = run_.dataSets[i];
    const QString name = d.record.empty() ? tr("data set %1").arg(d.key) : QString::fromStdString(d.record);
    QTableWidgetItem *nameItem = new QTableWidgetItem(name);
    nameItem->setData(Qt::UserRole + 1, QString::fromStdString(d.record));
    nameItem->setToolTip(tr("%1 (%2), key %3: click to plot").arg(name, QString::fromStdString(d.quantity)).arg(d.key));
    table_->setItem(i, 0, nameItem);
    table_->setItem(i, 1, new QTableWidgetItem(QString::fromStdString(d.reaction)));
    table_->setItem(i, 2, new NumberItem(QString::number(d.points), d.points));
    table_->setItem(i, 3, new NumberItem(Number(d.chiSquared), d.chiSquared));
    table_->setItem(i, 4, new NumberItem(Number(d.chiSquaredPerPoint, 4), d.chiSquaredPerPoint));
    table_->setItem(i, 5, new NumberItem(Number(d.norm, 5), d.norm));
    table_->setItem(i, 6, new NumberItem(d.normPrior ? QString::number(d.normPull, 'f', 2) : QString(),
                                         d.normPrior ? d.normPull : std::numeric_limits<double>::quiet_NaN()));
    const bool high = std::isfinite(d.chiSquaredPerPoint) && d.chiSquaredPerPoint > kHighChiPerPoint;
    const bool worst = std::isfinite(d.chiSquaredPerPoint) && d.chiSquaredPerPoint >= third && top > 0;
    const bool pull = d.normPrior && std::fabs(d.normPull) > kHighPull;
    for (int c = 0; c < table_->columnCount(); c++) {
      QTableWidgetItem *item = table_->item(i, c);
      QFont f = item->font();
      f.setBold(worst || high);
      item->setFont(f);
      if (high && c == kColumnChiPerPoint) item->setForeground(BarChart::warningColor());
      if (pull && c == 6) item->setForeground(BarChart::warningColor());
    }
  }
  table_->setSortingEnabled(true);
  table_->resizeColumnToContents(0);
  table_->setColumnWidth(0, std::min(table_->columnWidth(0), table_->fontMetrics().horizontalAdvance(QString(18, 'n')) + 12));
  const int rows = std::min(16, std::max(3, table_->rowCount()));
  table_->setMinimumHeight(table_->horizontalHeader()->sizeHint().height() + rows * table_->verticalHeader()->defaultSectionSize() + 4);

  QStringList flags;
  if (haveParameters_) {
    QStringList atLimit, noError;
    for (const front::ParameterFlag &f : parameters_.flagged) {
      const QString n = QString::fromStdString(f.name).toHtmlEscaped();
      if (f.atLimit) atLimit << QString("%1 = %2 %3").arg(n, Number(f.value), QString::fromStdString(f.unit));
      if (f.noError) noError << n;
    }
    if (!atLimit.isEmpty()) flags << Warn(tr("At a limit:")) + " " + atLimit.join(", ");
    if (!noError.isEmpty()) flags << Warn(tr("No error:")) + " " + noError.join(", ");
    if (!flags.isEmpty()) flags << Link("tab:fitting", tr("Fitting Settings"));
  }
  flagsText_->setText(flags.join("<br>"));
  flagsText_->setVisible(!flags.isEmpty());
  QStringList notes;
  for (const std::string &n : run_.leftOut) notes << tr("%1: active, but no point in the run").arg(QString::fromStdString(n).toHtmlEscaped());
  for (const std::string &n : run_.notes) notes << QString::fromStdString(n).toHtmlEscaped();
  notesText_->setText(notes.isEmpty() ? QString() : Muted(notes.join("<br>")));
  notesText_->setVisible(!notes.isEmpty());
}

void OverviewPage::fillCharts() {
  QStringList labels;
  QVector<double> values;
  QVector<bool> high;
  for (const front::DataSetResult &d : run_.present ? run_.dataSets : std::vector<front::DataSetResult>()) {
    labels << QString::fromStdString(d.record);
    values << d.chiSquaredPerPoint;
    high << (std::isfinite(d.chiSquaredPerPoint) && d.chiSquaredPerPoint > kHighChiPerPoint);
  }
  chiBars_->setBars(labels, values, high);
  chiBars_->setReference(1.0, "1");
  if (trace_->points() == 0) clearProgress();
}

void OverviewPage::fillSteps() {
  steps_.clear();
  if (validated_ && (!valid_ || errors_ > 0))
    steps_ << Step{tr("The project does not validate: %n error(s).", nullptr, errors_), "diagnostics"};
  else if (validated_ && warnings_ > 0)
    steps_ << Step{tr("The validation gives %n warning(s).", nullptr, warnings_), "diagnostics"};
  // a project not finished yet: the next thing to add, before any run (qt_gui.md section 19.8)
  const bool noPairs = project_.pairs.empty(), noLevels = project_.levels == 0,
             noData = project_.dataSets == 0 && project_.extrapolations == 0;
  if (noPairs) steps_ << Step{tr("No particle pair yet: add the entrance and exit pairs."), "tab:pairs"};
  else if (noLevels) steps_ << Step{tr("No level yet: add the levels of the compound nucleus."), "tab:levels"};
  if (noData) steps_ << Step{tr("No data set or extrapolation yet: add one."), "tab:segments"};
  if (!running_ && !run_.present && !noPairs && !noLevels && !noData)
    steps_ << Step{tr("There is no run yet: calculate or fit."), "tab:calculate"};
  if (run_.present && run_.projectNewer)
    steps_ << Step{tr("The outputs are older than the project file: calculate again."), "tab:calculate"};
  if (run_.present && !run_.changedInputs.empty())
    steps_ << Step{tr("The data of %n data set(s) changed since the run: calculate again.", nullptr,
                      (int)run_.changedInputs.size()),
                   "record:data:" + QString::fromStdString(run_.changedInputs.front())};
  if (run_.present && !run_.status.empty() && run_.status != "converged")
    steps_ << Step{tr("The fit ended with status \"%1\".").arg(QString::fromStdString(run_.status)), "tab:calculate"};
  QStringList high, pulls;
  QString worst, worstPull;
  double worstValue = 0, worstPullValue = 0;
  for (const front::DataSetResult &d : run_.dataSets) {
    if (std::isfinite(d.chiSquaredPerPoint) && d.chiSquaredPerPoint > kHighChiPerPoint) {
      high << QString::fromStdString(d.record);
      if (d.chiSquaredPerPoint > worstValue) worstValue = d.chiSquaredPerPoint, worst = QString::fromStdString(d.record);
    }
    if (d.normPrior && std::fabs(d.normPull) > kHighPull) {
      pulls << QString::fromStdString(d.record);
      if (std::fabs(d.normPull) > worstPullValue)
        worstPullValue = std::fabs(d.normPull), worstPull = QString::fromStdString(d.record);
    }
  }
  if (!high.isEmpty())
    steps_ << Step{tr("%n data set(s) have chi²/N above 3 (the largest: %1, %2).", nullptr, high.size())
                       .arg(worst, Number(worstValue, 3)),
                   "plot:" + worst};
  if (!pulls.isEmpty())
    steps_ << Step{tr("%n normalization(s) lie more than 2 sigma from their prior (the farthest: %1).", nullptr, pulls.size())
                       .arg(worstPull),
                   "record:data:" + worstPull};
  if (haveParameters_ && run_.present) {
    int atLimit = 0, noError = 0;
    for (const front::ParameterFlag &f : parameters_.flagged) atLimit += f.atLimit, noError += f.noError;
    if (atLimit) steps_ << Step{tr("%n parameter(s) ended at a limit.", nullptr, atLimit), "tab:fitting"};
    if (noError) steps_ << Step{tr("%n free parameter(s) have no error after the fit.", nullptr, noError), "tab:fitting"};
  }
  if (!run_.leftOut.empty())
    steps_ << Step{tr("%n active data set(s) were left out by the run (no point in the window).", nullptr, (int)run_.leftOut.size()),
                   "record:data:" + QString::fromStdString(run_.leftOut.front())};
  // Calculate and Fit once on the page: in "The last run" while there is none, here once there is one
  for (const char *name : {"overviewAction_runcalculate", "overviewAction_runfit"})
    if (QPushButton *b = findChild<QPushButton *>(name)) b->setVisible(run_.present || running_);
  QStringList lines;
  for (const Step &s : steps_) lines << "&bull; " + Link(s.link, s.text.toHtmlEscaped());
  if (lines.isEmpty())
    lines << (run_.present ? tr("Nothing stands out: the project validates and every data set has chi²/N of 3 or less.")
                           : tr("Nothing to report."));
  stepsText_->setText(lines.join("<br>"));
}

void OverviewPage::followLink(const QString &link) {
  if (link == "diagnostics") emit diagnosticsRequested();
  else if (link == "validate") emit validateRequested();
  else if (link.startsWith("run:")) emit commandRequested(link.mid(4));
  else if (link.startsWith("plot:")) emit plotRequested(link.mid(5));
  else if (link.startsWith("tab:")) emit tabRequested(link.mid(4));
  else if (link.startsWith("record:")) {
    const int colon = link.indexOf(':', 7);
    emit recordRequested(link.mid(7, colon - 7), link.mid(colon + 1));
  }
}

QString OverviewPage::tileText(const QString &name) const { return tiles_.contains(name) ? tiles_[name]->text() : QString(); }

bool OverviewPage::staleShown() const { return runStack_->currentIndex() == 1 && !stale_->isHidden(); }

bool OverviewPage::noRunShown() const { return runStack_->currentIndex() == 0; }
