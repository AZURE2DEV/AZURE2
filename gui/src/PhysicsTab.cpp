// The Physics tab: the tab, the project's session as the pages share it, and the session-free pages (Functions,
// Stopping).  docs/upgrade/qt_gui.md section 15; the engine's side docs/upgrade/physics_tab.md.
#include "PhysicsTab.h"

#include <cmath>
#include <map>

#include <QCheckBox>
#include <QComboBox>
#include <QDateTime>
#include <QFileDialog>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QSplitter>
#include <QStackedWidget>
#include <QVBoxLayout>

#include "PhysicsWorker.h"
#include "ProjectHolder.h"

namespace front = azure2::front;

// ---------------------------------------------------------------- PhysicsContext

PhysicsContext::PhysicsContext(ProjectHolder *holder, PhysicsWorker *worker, QObject *parent)
    : QObject(parent), holder_(holder), worker_(worker) {}

QString PhysicsContext::projectFile() const {
  const QString path = holder_->path();
  return !path.isEmpty() && QFileInfo(path).isFile() ? path : QString();
}

QString PhysicsContext::needs() const {
  if (projectFile().isEmpty())
    return tr("This page reads the saved project through the engine's session: open a project (File > Open), or save "
              "the new one (File > Save).  The pages Functions and Stopping need no project.");
  return QString();
}

void PhysicsContext::setParameters(Parameters choice, const QString &file) {
  choice_ = choice;
  if (choice == ChosenFile) chosenFile_ = file;
}

QString PhysicsContext::lastFit() const {
  return projectFile().isEmpty() ? QString() : QString::fromStdString(front::LastParameterFile(holder_->project()));
}

QString PhysicsContext::parameterFile() const {
  if (choice_ == LastFit) return lastFit();
  if (choice_ == ChosenFile) return chosenFile_;
  return QString();
}

std::string PhysicsContext::key() const {
  const QString file = projectFile();
  if (file.isEmpty()) return std::string();
  const QString param = parameterFile();
  const QString stamp = QFileInfo(file).lastModified().toString(Qt::ISODateWithMs) + "|" +
                        (param.isEmpty() ? QString() : QFileInfo(param).lastModified().toString(Qt::ISODateWithMs));
  return (file + "|" + param + "|" + stamp).toStdString();
}

void PhysicsContext::run(const QString &lane, Job job, std::function<void(bool, const QString &)> done) {
  const std::string file = projectFile().toStdString(), param = parameterFile().toStdString(), k = key();
  if (file.empty()) {
    done(false, needs());
    return;
  }
  auto result = std::make_shared<std::pair<bool, std::string>>(false, std::string());
  worker_->submit(
      lane,
      [file, param, k, job, result](PhysicsSessionState &state, const std::atomic<bool> &cancel) {
        std::string why;
        result->first = state.Ensure(k, file, param, why) && job(state.session, cancel, why);
        result->second = why;
      },
      [result, done]() { done(result->first, QString::fromStdString(result->second)); });
}

void PhysicsContext::refresh() {
  auto info = std::make_shared<Info>();
  info->key = key();
  run(
      "info",
      [info](front::PhysicsSession &s, const std::atomic<bool> &, std::string &) {
        info->parametrization = s.Parametrization();
        info->pairs = s.Pairs();
        info->groups = s.Groups();
        info->poles = s.Poles();
        info->experiments = s.ThmExperiments();
        return true;
      },
      [this, info](bool ok, const QString &why) {
        info->ok = ok;
        info->why = why;
        info_ = *info;
        emit infoChanged();
      });
}

void PhysicsContext::ensure() {
  if (!current() && !projectFile().isEmpty()) refresh();
}

QString PhysicsContext::pairLabel(int pairKey) const {
  for (const front::PairChoice &p : front::Pairs(holder_->project()))
    if (p.key == pairKey) return QString::fromStdString(p.label);
  return tr("pair %1").arg(pairKey);
}

namespace physicsqt {

QString JPi(double J, int parity) {
  const long twice = std::lround(2 * J);
  return (twice % 2 ? QString("%1/2").arg(twice) : QString::number(twice / 2)) + (parity < 0 ? "-" : "+");
}

// ---------------------------------------------------------------- ProjectPage

ProjectPage::ProjectPage(PhysicsContext *context, QWidget *parent) : QWidget(parent), context_(context) {
  stack_ = new QStackedWidget;
  needs_ = new QLabel;
  needs_->setObjectName("physicsNeeds");
  needs_->setWordWrap(true);
  needs_->setAlignment(Qt::AlignCenter);
  needs_->setMargin(24);
  stack_->addWidget(needs_);
  QVBoxLayout *v = new QVBoxLayout(this);
  v->setContentsMargins(0, 0, 0, 0);
  v->addWidget(stack_);
  connect(context_, &PhysicsContext::infoChanged, this, [this]() {
    const PhysicsContext::Info &i = context_->info();
    setNeeds(i.ok ? QString() : tr("The engine could not set up the project's session:\n%1").arg(i.why));
    // a hidden page fills its drop-downs now and draws when it is shown (activate)
    if (isVisible()) shownKey_ = i.key;
    infoChanged();
  });
}

void ProjectPage::setContent(QWidget *content) {
  stack_->addWidget(content);
  stack_->setCurrentIndex(0);
}

void ProjectPage::setNeeds(const QString &text) {
  needs_->setText(text);
  stack_->setCurrentIndex(text.isEmpty() && stack_->count() > 1 ? 1 : 0);
}

QString ProjectPage::needsText() const { return stack_->currentIndex() == 0 ? needs_->text() : QString(); }

bool ProjectPage::contentShown() const { return stack_->currentIndex() == 1; }

void ProjectPage::activate() {
  const QString needs = context_->needs();
  if (!needs.isEmpty()) {
    setNeeds(needs);
    return;
  }
  if (context_->current()) {
    setNeeds(context_->info().ok ? QString() : tr("The engine could not set up the project's session:\n%1").arg(context_->info().why));
    if (context_->info().ok && shownKey_ != context_->info().key) {
      shownKey_ = context_->info().key;
      infoChanged();
    }
    return;
  }
  setNeeds(tr("Setting up the project's session (the levels and the data, as a calculation would)..."));
  context_->ensure();
}

// ---------------------------------------------------------------- FunctionsPage

FunctionsPage::FunctionsPage(bool stopping, PhysicsWorker *worker, QWidget *parent)
    : QWidget(parent), stopping_(stopping), worker_(worker) {
  const QString name = stopping ? "stopping" : "functions";
  setObjectName(name + "Page");
  rows_ = stopping ? front::StoppingFunctions() : front::PhysicsFunctions();
  function_ = new QComboBox;
  function_->setObjectName(name + "Function");
  for (const front::PhysicsFunction &f : rows_)
    function_->addItem(QString("%1  --  %2").arg(QString::fromStdString(f.name), QString::fromStdString(f.summary)));
  for (int i = 0; i < (int)rows_.size(); i++)
    function_->setItemData(i, QString::fromStdString(rows_[i].group + ": " + rows_[i].summary), Qt::ToolTipRole);
  return_ = new QComboBox;
  return_->setObjectName(name + "Return");
  for (QComboBox *box : {function_, return_}) {
    box->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    box->setMinimumContentsLength(24);
  }
  doc_ = new QLabel;
  doc_->setObjectName(name + "Doc");
  doc_->setWordWrap(true);
  doc_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  doc_->setAlignment(Qt::AlignTop | Qt::AlignLeft);
  QScrollArea *docScroll = new QScrollArea;
  docScroll->setWidget(doc_);
  docScroll->setWidgetResizable(true);
  docScroll->setFrameShape(QFrame::NoFrame);
  docScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  docScroll->setMaximumHeight(fontMetrics().height() * 8);
  form_ = new ArgumentForm;
  form_->setObjectName(name + "Form");
  QPushButton *drawButton = new QPushButton(tr("Draw"));
  drawButton->setObjectName(name + "Draw");
  drawButton->setProperty("primary", true);
  drawButton->setToolTip(tr("Evaluate and draw this curve alone"));
  QPushButton *add = new QPushButton(tr("Add curve"));
  add->setObjectName(name + "Add");
  add->setToolTip(tr("Evaluate and draw this curve over the others (another l, radius, material ...)"));
  QPushButton *clearButton = new QPushButton(tr("Clear"));
  connect(drawButton, &QPushButton::clicked, this, [this]() { draw(false); });
  connect(add, &QPushButton::clicked, this, [this]() { draw(true); });
  connect(clearButton, &QPushButton::clicked, this, &FunctionsPage::clear);
  connect(function_, &QComboBox::currentIndexChanged, this, &FunctionsPage::functionChanged);
  connect(form_->sweep()->log, &QCheckBox::toggled, this, [this](bool on) { panel_->logX()->setChecked(on); });
  connect(return_, &QComboBox::currentIndexChanged, this, &FunctionsPage::redraw);
  panel_ = new PlotPanel(name + "Panel");
  if (stopping) panel_->setCitation(tr("Stopping powers: %1").arg(QString::fromStdString(front::StoppingCitation())));

  QWidget *controls = new QWidget;
  QVBoxLayout *c = new QVBoxLayout(controls);
  c->setContentsMargins(0, 0, 0, 0);
  QHBoxLayout *top = new QHBoxLayout;
  top->addWidget(new QLabel(stopping ? tr("Stopping function") : tr("Function")));
  top->addWidget(function_, 1);
  c->addLayout(top);
  c->addWidget(docScroll);
  c->addWidget(form_);
  QHBoxLayout *ret = new QHBoxLayout;
  ret->addWidget(new QLabel(tr("Draw")));
  ret->addWidget(return_, 1);
  c->addLayout(ret);
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(drawButton);
  buttons->addWidget(add);
  buttons->addWidget(clearButton);
  buttons->addStretch(1);
  c->addLayout(buttons);
  c->addStretch(1);
  QScrollArea *controlsScroll = new QScrollArea;
  controlsScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  controlsScroll->setWidget(controls);
  controlsScroll->setWidgetResizable(true);
  controlsScroll->setFrameShape(QFrame::NoFrame);
  controlsScroll->setMinimumWidth(fontMetrics().horizontalAdvance(QString(46, 'x')));
  QSplitter *split = new QSplitter;
  split->addWidget(controlsScroll);
  split->addWidget(panel_);
  split->setStretchFactor(1, 3);
  (new QVBoxLayout(this))->addWidget(split);
  layout()->setContentsMargins(0, 0, 0, 0);
  functionChanged();
}

void FunctionsPage::functionChanged() {
  const int i = function_->currentIndex();
  if (i < 0 || i >= (int)rows_.size()) return;
  const front::PhysicsFunction &f = rows_[i];
  form_->setFunction(f);
  {
    const QSignalBlocker block(return_);
    return_->clear();
    for (const front::PhysicsReturn &r : f.returns)
      return_->addItem(QString::fromStdString(r.name) + (r.unit.empty() ? QString() : " (" + QString::fromStdString(r.unit) + ")"),
                       QString::fromStdString(r.doc));
    for (int k = 0; k < return_->count(); k++) return_->setItemData(k, QString::fromStdString(f.returns[k].doc), Qt::ToolTipRole);
  }
  doc_->setText(QString("<b>%1</b> (%2)<br>%3<br><i>%4</i>")
                    .arg(QString::fromStdString(f.name).toHtmlEscaped(), QString::fromStdString(f.group).toHtmlEscaped(),
                         QString::fromStdString(f.doc).toHtmlEscaped(), QString::fromStdString(f.library).toHtmlEscaped()));
  panel_->logX()->setChecked(form_->swept() >= 0 && form_->sweep()->log->isChecked());
  curves_.clear();
  draw(false);
}

void FunctionsPage::clear() {
  curves_.clear();
  redraw();
}

void FunctionsPage::draw(bool add) {
  const int row = function_->currentIndex();
  if (row < 0) return;
  std::vector<std::vector<double>> columns;
  QVector<double> grid;
  QString why;
  if (!form_->columns(columns, grid, why, stopping_)) {
    panel_->setStatus(tr("Not drawn: %1").arg(why));
    return;
  }
  // a swept whole number (l, an index) takes whole values only
  const int s = form_->swept();
  if (s >= 0 && rows_[row].arguments[s].kind != "real") {
    QVector<double> whole;
    for (double v : grid)
      if (whole.isEmpty() || std::lround(v) != whole.back()) whole << (double)std::lround(v);
    grid = whole;
    std::size_t column = 0;
    for (int k = 0; k < s; k++) column += !(stopping_ && rows_[row].arguments[k].kind == "text");
    columns[column].assign(grid.begin(), grid.end());
  }
  const std::string name = rows_[row].name, formula = stopping_ ? form_->text(0).toStdString() : std::string();
  const bool stopping = stopping_;
  auto curve = std::make_shared<Curve>();
  curve->grid = grid;
  curve->parts = form_->fixedParts();
  auto result = std::make_shared<std::pair<bool, std::string>>(false, std::string());
  worker_->submit(
      objectName(),
      [=](PhysicsSessionState &, const std::atomic<bool> &) {
        result->first = stopping ? front::EvaluateStopping(name, formula, columns, curve->values, result->second)
                                 : front::EvaluatePhysics(name, columns, curve->values, result->second);
      },
      [this, curve, result, add]() {
        if (!result->first) {
          panel_->setStatus(tr("Not drawn: %1").arg(QString::fromStdString(result->second)));
          return;
        }
        if (!add) curves_.clear();
        curves_ << *curve;
        redraw();
      });
}

void FunctionsPage::redraw() {
  const int row = function_->currentIndex(), ret = return_->currentIndex();
  if (row < 0 || ret < 0) return;
  const front::PhysicsFunction &f = rows_[row];
  QList<PlotSeries> series;
  QStringList notes;
  // a curve's label: the fixed arguments in which the curves differ (the function's name for one curve)
  QList<int> differ;
  for (int k = 0; !curves_.isEmpty() && k < curves_[0].parts.size(); k++)
    for (const Curve &c : curves_)
      if (c.parts.value(k) != curves_[0].parts[k]) {
        differ << k;
        break;
      }
  QStringList common;
  for (int k = 0; !curves_.isEmpty() && k < curves_[0].parts.size(); k++)
    if (!differ.contains(k)) common << curves_[0].parts[k];
  for (int i = 0; i < curves_.size(); i++) {
    const Curve &c = curves_[i];
    PlotSeries ps;
    QStringList own;
    for (int k : differ) own << c.parts.value(k);
    ps.label = own.isEmpty() ? QString::fromStdString(f.name) : own.join(", ");
    ps.color = CurveColour(i);
    ps.modelX = c.grid;
    if (ret < (int)c.values.values.size()) ps.modelY = QVector<double>(c.values.values[ret].begin(), c.values.values[ret].end());
    // the points the row did not compute are gaps (NaN), counted by their status
    std::map<int, int> counts;
    for (int st : c.values.status)
      if (st != 0) counts[st]++;
    for (const auto &[st, n] : counts) {
      const std::string words = stopping_ ? front::StoppingStatusText(f.name, st) : front::PhysicsStatusText(f.name, st);
      notes << tr("%1: %2 of %3 points not computed: %4").arg(ps.label).arg(n).arg(c.grid.size()).arg(QString::fromStdString(words));
    }
    series << ps;
  }
  const int s = form_->swept();
  QString xTitle;
  if (s >= 0) {
    const front::PhysicsArgument &a = f.arguments[s];
    xTitle = QString::fromStdString(a.name) + (a.unit.empty() ? QString() : " (" + QString::fromStdString(a.unit) + ")");
  }
  panel_->show(series, xTitle, return_->currentText());
  if (!common.isEmpty()) notes.prepend(tr("At %1").arg(common.join(", ")));
  panel_->setStatus(notes.join("\n"));
}

}  // namespace physicsqt

// ---------------------------------------------------------------- PhysicsTab

PhysicsTab::PhysicsTab(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  setObjectName("physicsTab");
  sessionWorker_ = new PhysicsWorker(this);
  freeWorker_ = new PhysicsWorker(this);
  context_ = new PhysicsContext(holder, sessionWorker_, this);
  pages_ = new QListWidget;
  pages_->setObjectName("physicsPages");
  const QStringList names = {tr("Functions"),   tr("Channels of the project"), tr("R-matrix of the project"),
                             tr("Pole map"),    tr("THM"),                     tr("Stopping")};
  const QStringList tips = {
      tr("Every function of the physics table (Coulomb functions, P, S, phases, widths, kinematics, THM vertex) "
         "over one of its arguments; no project needed"),
      tr("The project's channel functions of a pair: P, S, phases, F and G at the channel radius"),
      tr("One J^pi group: phase shifts, inelasticities, |U|, Argand plots, |T|, P and S per channel, |A|, with the poles"),
      tr("A real-energy map of log10 |det A^-1|, |U| or |A| with the poles on it"),
      tr("The THM experiments' vertex, Coulomb line shape and distortion"),
      tr("SRIM stopping powers, ranges, energy loss and straggling of a material; no project needed")};
  for (int i = 0; i < names.size(); i++) {
    QListWidgetItem *item = new QListWidgetItem(names[i]);
    item->setToolTip(tips[i]);
    pages_->addItem(item);
  }
  pages_->setMaximumWidth(fontMetrics().horizontalAdvance(tr("R-matrix of the project")) + 48);

  functions_ = new physicsqt::FunctionsPage(false, freeWorker_);
  stopping_ = new physicsqt::FunctionsPage(true, freeWorker_);
  channels_ = new physicsqt::ChannelsPage(context_);
  rmatrix_ = new physicsqt::RMatrixPage(context_);
  poleMap_ = new physicsqt::PoleMapPage(context_);
  thm_ = new physicsqt::ThmPage(context_);
  stack_ = new QStackedWidget;
  for (QWidget *w : std::initializer_list<QWidget *>{functions_, channels_, rmatrix_, poleMap_, thm_, stopping_})
    stack_->addWidget(w);

  // the session bar of the project pages: which parameters, the session's state
  sessionBar_ = new QWidget;
  parameters_ = new QComboBox;
  parameters_->setObjectName("physicsParameters");
  parameters_->addItem(tr("Start values of the project"), PhysicsContext::StartValues);
  parameters_->addItem(tr("The last fit (param.sav of the output directory)"), PhysicsContext::LastFit);
  parameters_->addItem(tr("A parameter file..."), PhysicsContext::ChosenFile);
  parameters_->setToolTip(tr("The parameters the session's levels take: the project's own values, the parameter file "
                             "the last fit wrote, or a file read as AZURE2 --parameter-file reads it"));
  connect(parameters_, &QComboBox::activated, this, &PhysicsTab::parametersChosen);
  QPushButton *reread = new QPushButton(tr("Re-read"));
  reread->setToolTip(tr("Set the session up again from the saved project (after an edit, save first)"));
  connect(reread, &QPushButton::clicked, this, [this]() {
    context_->refresh();
  });
  sessionStatus_ = new QLabel;
  sessionStatus_->setObjectName("physicsSessionStatus");
  sessionStatus_->setProperty("role", "hint");
  sessionStatus_->setWordWrap(true);
  QHBoxLayout *bar = new QHBoxLayout(sessionBar_);
  bar->setContentsMargins(0, 0, 0, 0);
  bar->addWidget(new QLabel(tr("Parameters")));
  bar->addWidget(parameters_);
  bar->addWidget(reread);
  bar->addWidget(sessionStatus_, 1);
  connect(context_, &PhysicsContext::infoChanged, this, [this]() {
    const PhysicsContext::Info &i = context_->info();
    if (!i.ok) {
      sessionStatus_->setText(tr("No session: %1").arg(i.why.section('\n', 0, 0)));
      return;
    }
    const QString param = context_->parameterFile();
    sessionStatus_->setText(tr("Session of %1 at %2; %3 parametrization; %4 J^pi groups, %5 levels%6")
                                .arg(QFileInfo(context_->projectFile()).fileName(),
                                     param.isEmpty() ? tr("the start values") : QFileInfo(param).fileName(),
                                     QString::fromStdString(i.parametrization))
                                .arg(i.groups.size())
                                .arg(i.poles.size())
                                .arg(i.extrapolations ? tr(" (set up with its extrapolations: no data set)") : QString()));
  });

  // the busy strip: either worker running
  busyStrip_ = new QWidget;
  busyStrip_->setObjectName("physicsBusy");
  busyBar_ = new QProgressBar;
  busyBar_->setRange(0, 0);
  busyBar_->setMaximumWidth(220);
  QPushButton *cancelButton = new QPushButton(tr("Cancel"));
  cancelButton->setObjectName("physicsCancel");
  cancelButton->setToolTip(tr("Stop the computation; a session's set-up in progress is finished in the background "
                              "and its result dropped"));
  connect(cancelButton, &QPushButton::clicked, this, &PhysicsTab::cancel);
  QHBoxLayout *busy = new QHBoxLayout(busyStrip_);
  busy->setContentsMargins(0, 0, 0, 0);
  busy->addWidget(new QLabel(tr("Computing...")));
  busy->addWidget(busyBar_);
  busy->addWidget(cancelButton);
  busy->addStretch(1);
  busyStrip_->hide();
  for (PhysicsWorker *w : {sessionWorker_, freeWorker_}) connect(w, &PhysicsWorker::busyChanged, this, &PhysicsTab::busyChanged);

  QWidget *right = new QWidget;
  QVBoxLayout *r = new QVBoxLayout(right);
  r->setContentsMargins(0, 0, 0, 0);
  r->addWidget(sessionBar_);
  r->addWidget(stack_, 1);
  r->addWidget(busyStrip_);
  QHBoxLayout *h = new QHBoxLayout(this);
  h->addWidget(pages_);
  h->addWidget(right, 1);
  connect(pages_, &QListWidget::currentRowChanged, this, &PhysicsTab::showPage);
  pages_->setCurrentRow(Functions);
  showPage(Functions);
  connect(holder_, &ProjectHolder::opened, this, &PhysicsTab::projectChanged);
}

PhysicsTab::~PhysicsTab() = default;

void PhysicsTab::showPage(int page) {
  if (page < 0 || page >= stack_->count()) return;
  if (pages_->currentRow() != page) pages_->setCurrentRow(page);
  stack_->setCurrentIndex(page);
  const bool project = page != Functions && page != Stopping;
  sessionBar_->setVisible(project);
  if (project) static_cast<physicsqt::ProjectPage *>(stack_->widget(page))->activate();
}

bool PhysicsTab::busy() const { return sessionWorker_->busy() || freeWorker_->busy(); }

void PhysicsTab::cancel() {
  sessionWorker_->cancel();
  freeWorker_->cancel();
}

void PhysicsTab::busyChanged() { busyStrip_->setVisible(busy()); }

void PhysicsTab::projectChanged() {
  const int page = stack_->currentIndex();
  if (page != Functions && page != Stopping) static_cast<physicsqt::ProjectPage *>(stack_->widget(page))->activate();
}

void PhysicsTab::parametersChosen() {
  const auto choice = (PhysicsContext::Parameters)parameters_->currentData().toInt();
  if (choice == PhysicsContext::LastFit && context_->lastFit().isEmpty()) {
    sessionStatus_->setText(tr("No param.sav in the project's output directory yet (a fit writes it): the start values "
                               "are kept."));
    parameters_->setCurrentIndex(0);
    return;
  }
  QString file;
  if (choice == PhysicsContext::ChosenFile) {
    file = QFileDialog::getOpenFileName(this, tr("Parameter file"), holder_->directory(),
                                        tr("Parameter files (*.sav *.par *.fit);;All files (*)"));
    if (file.isEmpty()) {
      parameters_->setCurrentIndex(context_->parameters());
      return;
    }
    parameters_->setItemText(2, tr("File: %1").arg(QFileInfo(file).fileName()));
  }
  context_->setParameters(choice, file);
  projectChanged();
}
