// The Physics tab's project pages: Channels, R-matrix, Pole map, THM (PhysicsTab.h; docs/upgrade/qt_gui.md section
// 15).  Each draws a per-session call of the facade; what can be drawn from each call is the facade's
// (PhysicsPlots), the values are the engine's.
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QVBoxLayout>

#include "HeatmapWidget.h"
#include "PhysicsTab.h"

namespace front = azure2::front;

namespace physicsqt {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

/// The plots of a call that a line plot draws (or, `maps`, the maps).
std::vector<front::PhysicsPlot> PlotsOf(const std::string &call, bool maps) {
  std::vector<front::PhysicsPlot> out;
  for (const front::PhysicsPlot &p : front::PhysicsPlots(call))
    if (p.mapRows.empty() != maps) out.push_back(p);
  return out;
}

void FillPlots(QComboBox *box, const std::vector<front::PhysicsPlot> &plots) {
  const QSignalBlocker block(box);
  const QString was = box->currentData().toString();
  box->clear();
  for (const front::PhysicsPlot &p : plots) {
    box->addItem(QString::fromStdString(p.title), QString::fromStdString(p.name));
    box->setItemData(box->count() - 1, QString::fromStdString(p.doc), Qt::ToolTipRole);
  }
  const int i = box->findData(was);
  box->setCurrentIndex(i >= 0 ? i : 0);
}

/// The plot of that name (a copy: the list is often a temporary).
std::optional<front::PhysicsPlot> PlotNamed(const std::vector<front::PhysicsPlot> &plots, const QString &name) {
  for (const front::PhysicsPlot &p : plots)
    if (QString::fromStdString(p.name) == name) return p;
  return std::nullopt;
}

/// The words of a per-session call (its title, doc, axis title).
front::PhysicsCall CallNamed(const std::string &name) {
  for (const front::PhysicsCall &c : front::PhysicsCalls())
    if (c.name == name) return c;
  return front::PhysicsCall();
}

QVector<double> Values(const std::vector<double> &v) { return QVector<double>(v.begin(), v.end()); }

/// The rows of a combo of J^pi groups.
void FillGroups(QComboBox *box, const std::vector<front::PhysicsGroup> &groups) {
  const QSignalBlocker block(box);
  const int was = box->currentData().toInt();
  box->clear();
  for (const front::PhysicsGroup &g : groups)
    box->addItem(QObject::tr("%1  (group %2: %3 levels, %4 channels)")
                     .arg(JPi(g.J, g.parity))
                     .arg(g.jgroup)
                     .arg(g.levels)
                     .arg(g.channels.size()),
                 g.jgroup);
  const int i = box->findData(was);
  box->setCurrentIndex(i >= 0 ? i : 0);
}

/// A sensible first E_x range of a group: from the lowest open particle threshold of its channels (or 1 MeV below
/// its lowest level) to 1 MeV above its highest level.
std::pair<double, double> GroupRange(const PhysicsContext::Info &info, int jgroup) {
  double lo = INFINITY, hi = -INFINITY, threshold = INFINITY;
  for (const front::PhysicsPole &p : info.poles)
    if (jgroup <= 0 || p.jgroup == jgroup) lo = std::min(lo, p.energy), hi = std::max(hi, p.energy);
  for (const front::PhysicsGroup &g : info.groups) {
    if (jgroup > 0 && g.jgroup != jgroup) continue;
    for (const front::PhysicsChannel &c : g.channels) {
      if (c.radiation != "P") continue;
      for (const front::PhysicsPair &p : info.pairs)
        if (p.key == c.pair) threshold = std::min(threshold, p.separation + p.excitation);
    }
  }
  if (!std::isfinite(lo)) lo = std::isfinite(threshold) ? threshold : 0, hi = lo + 3;
  // above the threshold (a bound level below it has no phase shift to show), at least 2 MeV wide
  const double a = std::max(0.0, std::isfinite(threshold) ? std::max(threshold + 0.01, lo - 1) : lo - 1);
  return {a, std::max(hi + 1, a + 2)};
}

QString FactsText(const front::PhysicsCurves &c) {
  QStringList parts;
  for (const auto &f : c.facts) parts << QString::fromStdString(f.first + " = " + f.second);
  QString out = parts.join("; ");
  for (const std::string &n : c.notes) out += "\n" + QString::fromStdString(n);
  return out;
}

/// The facts on the page without the list of level energies (the pole table below lists them; the tip has all).
QString ShortFactsText(const front::PhysicsCurves &c) {
  QStringList parts;
  int levels = 0;
  for (const auto &f : c.facts) {
    if (f.first.compare(0, 8, "E_level ") == 0) {
      levels++;
      continue;
    }
    parts << QString::fromStdString(f.first + " = " + f.second);
  }
  if (levels > 0) parts << QObject::tr("%n level(s)", nullptr, levels);
  QString out = parts.join("; ");
  for (const std::string &n : c.notes) out += "\n" + QString::fromStdString(n);
  return out;
}

QWidget *Scrolled(QWidget *w) {
  // the form at its natural height, the rest of the column empty (a form layout would stretch its rows)
  QWidget *box = new QWidget;
  QVBoxLayout *v = new QVBoxLayout(box);
  v->addWidget(w);
  v->addStretch(1);
  QScrollArea *s = new QScrollArea;
  s->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  s->setWidget(box);
  s->setWidgetResizable(true);
  s->setFrameShape(QFrame::NoFrame);
  return s;
}

}  // namespace

/// A phase shift in degrees, defined modulo 180 (the engine gives arg(U_cc)/2 in (-90, 90], in [0, 180) for identical
/// particles), drawn continuous: a step of more than 90 degrees between neighbouring points is taken as the wrap and
/// removed, and the curve is moved by a multiple of 180 so that it starts in (-90, 90].  Display only: a phase near 0
/// that wraps to 180 and back (c12c12_tumino2018, 12C+12C) no longer draws as a comb (qt_gui.md section 19.8).
QVector<double> UnwrappedPhaseDegrees(const QVector<double> &deg) {
  QVector<double> out = deg;
  double offset = 0;
  int first = -1;
  for (int i = 0; i < out.size(); i++) {
    if (!std::isfinite(deg[i])) continue;
    if (first >= 0) {
      int prev = i - 1;
      while (prev >= 0 && !std::isfinite(deg[prev])) prev--;
      const double step = deg[i] - deg[prev];
      if (step > 90) offset -= 180;
      else if (step < -90) offset += 180;
    } else {
      first = i;
    }
    out[i] = deg[i] + offset;
  }
  if (first >= 0) {
    const double shift = -180.0 * std::round(out[first] / 180.0);
    for (double &v : out) v += shift;
  }
  return out;
}


// ---------------------------------------------------------------- Channels

ChannelsPage::ChannelsPage(PhysicsContext *context, QWidget *parent) : ProjectPage(context, parent) {
  setObjectName("channelsPage");
  pair_ = new QComboBox;
  pair_->setObjectName("channelsPair");
  l_ = new QSpinBox;
  l_->setObjectName("channelsL");
  l_->setRange(0, 12);
  radius_ = new QLineEdit;
  radius_->setObjectName("channelsRadius");
  radius_->setToolTip(tr("The channel radius a_c (fm); empty: the pair's own radius"));
  twoJ_ = new QSpinBox;
  twoJ_->setRange(-1, 60);
  twoJ_->setValue(-1);
  twoJ_->setSpecialValueText(tr("none"));
  twoJ_->setToolTip(tr("2J of the channel: needed only for a pair with a spin-orbit term in its nuclear potential"));
  energies_ = new SweepFields("MeV (c.m.)");
  energies_->setObjectName("channelsEnergies");
  energies_->set(0.05, 5, 300, false);
  quantity_ = new QComboBox;
  quantity_->setObjectName("channelsQuantity");
  FillPlots(quantity_, PlotsOf("ChannelFunctions", false));
  panel_ = new PlotPanel("channelsPanel");
  connect(energies_->log, &QCheckBox::toggled, panel_->logX(), &QCheckBox::setChecked);
  QPushButton *drawButton = new QPushButton(tr("Draw"));
  drawButton->setObjectName("channelsDraw");
  drawButton->setProperty("primary", true);
  QPushButton *add = new QPushButton(tr("Add curve"));
  add->setToolTip(tr("Draw over the others (another pair, l or radius)"));
  QPushButton *clearButton = new QPushButton(tr("Clear"));
  connect(drawButton, &QPushButton::clicked, this, [this]() { draw(false); });
  connect(add, &QPushButton::clicked, this, [this]() { draw(true); });
  connect(clearButton, &QPushButton::clicked, this, [this]() {
    curves_.clear();
    redraw();
  });
  connect(quantity_, &QComboBox::currentIndexChanged, this, [this]() {
    const auto p = PlotNamed(front::PhysicsPlots("ChannelFunctions"), quantity_->currentData().toString());
    if (p) panel_->logY()->setChecked(p->logY);
    redraw();
  });
  if (const auto p = PlotNamed(PlotsOf("ChannelFunctions", false), quantity_->currentData().toString()))
    panel_->logY()->setChecked(p->logY);

  QWidget *controls = new QWidget;
  QFormLayout *f = new QFormLayout(controls);
  f->addRow(tr("Pair"), pair_);
  f->addRow(tr("l"), l_);
  f->addRow(tr("Radius (fm)"), radius_);
  f->addRow(tr("2J"), twoJ_);
  f->addRow(tr("E_c"), energies_);
  f->addRow(tr("Draw"), quantity_);
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(drawButton);
  buttons->addWidget(add);
  buttons->addWidget(clearButton);
  buttons->addStretch(1);
  f->addRow(buttons);
  // what the page computes: the tip of the pair's drop-down, not a paragraph (qt_gui.md section 19.6)
  pair_->setToolTip(QString::fromStdString(CallNamed("ChannelFunctions").doc));
  QSplitter *split = new QSplitter;
  split->addWidget(Scrolled(controls));
  split->addWidget(panel_);
  split->setStretchFactor(1, 3);
  setContent(split);
}

void ChannelsPage::infoChanged() {
  const PhysicsContext::Info &info = context_->info();
  const QSignalBlocker block(pair_);
  const int was = pair_->currentData().toInt();
  pair_->clear();
  for (const front::PhysicsPair &p : info.pairs)
    if (p.kind == "particle")
      pair_->addItem(tr("%1  (pair %2, a = %3 fm)").arg(context_->pairLabel(p.key)).arg(p.key).arg(p.radius), p.key);
  const int i = pair_->findData(was);
  pair_->setCurrentIndex(i >= 0 ? i : 0);
  if (curves_.isEmpty() && pair_->count() && isVisible()) draw(false);
}

void ChannelsPage::draw(bool add) {
  QVector<double> grid;
  QString why;
  if (pair_->count() == 0) return;
  if (!energies_->grid(grid, why)) {
    panel_->setStatus(tr("Not drawn: %1").arg(why));
    return;
  }
  const int pairKey = pair_->currentData().toInt(), l = l_->value(), twoJ = twoJ_->value();
  bool ok = true;
  const double radius = radius_->text().trimmed().isEmpty() ? 0.0 : radius_->text().toDouble(&ok);
  if (!ok || radius < 0) {
    panel_->setStatus(tr("Not drawn: the radius must be a number >= 0 (empty: the pair's)"));
    return;
  }
  auto curve = std::make_shared<Curve>();
  curve->label = tr("%1, l = %2%3").arg(context_->pairLabel(pairKey)).arg(l).arg(radius > 0 ? tr(", a = %1 fm").arg(radius) : QString());
  const std::vector<double> e(grid.begin(), grid.end());
  context_->run(
      "channels",
      [=](front::PhysicsSession &s, const std::atomic<bool> &, std::string &w) {
        return s.ChannelFunctions(pairKey, l, radius, twoJ, e, curve->answer, w);
      },
      [this, curve, add](bool done, const QString &w) {
        if (!done) {
          panel_->setStatus(tr("Not drawn: %1").arg(w));
          return;
        }
        if (!add) curves_.clear();
        curves_ << *curve;
        redraw();
      });
}

void ChannelsPage::redraw() {
  const std::vector<front::PhysicsPlot> plots = front::PhysicsPlots("ChannelFunctions");
  const auto p = PlotNamed(plots, quantity_->currentData().toString());
  if (!p) return;
  QList<PlotSeries> series;
  QStringList notes;
  for (int i = 0; i < curves_.size(); i++) {
    const front::PhysicsSeries *s = curves_[i].answer.Find(p->y);
    if (!s) continue;
    PlotSeries ps;
    ps.label = curves_[i].label;
    ps.color = CurveColour(i);
    ps.modelX = Values(curves_[i].answer.axisValues);
    ps.modelY = Values(s->values);
    if (const front::PhysicsSeries *open = curves_[i].answer.Find("open")) {
      const int closed = (int)std::count(open->values.begin(), open->values.end(), 0.0);
      if (closed) notes << tr("%1: %2 points below the threshold (E_c <= 0): closed channel").arg(ps.label).arg(closed);
    }
    series << ps;
  }
  panel_->show(series, QString::fromStdString(CallNamed("ChannelFunctions").axisTitle), QString::fromStdString(p->yTitle));
  panel_->setStatus(notes.join("\n"));
}

// ---------------------------------------------------------------- R-matrix

RMatrixPage::RMatrixPage(PhysicsContext *context, QWidget *parent) : ProjectPage(context, parent) {
  setObjectName("rmatrixPage");
  group_ = new QComboBox;
  group_->setObjectName("rmatrixGroup");
  energies_ = new SweepFields("MeV (E_x)");
  energies_->setObjectName("rmatrixEnergies");
  energies_->set(1, 5, 400, false);
  quantity_ = new QComboBox;
  quantity_->setObjectName("rmatrixQuantity");
  FillPlots(quantity_, PlotsOf("Block", false));
  series_ = new QListWidget;
  series_->setObjectName("rmatrixSeries");
  series_->setToolTip(tr("The lines of this quantity: channels, channel pairs or level pairs"));
  series_->setMaximumHeight(fontMetrics().height() * 7);
  poles_ = new QTableWidget;
  poles_->setObjectName("rmatrixPoles");
  poles_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  poles_->setSelectionBehavior(QAbstractItemView::SelectRows);
  poles_->setSelectionMode(QAbstractItemView::SingleSelection);
  poles_->verticalHeader()->hide();
  poles_->setToolTip(tr("The group's levels as the engine reports them (parameters.out): the physical E_x and the "
                        "observed widths.  Click a level to centre the plot on it."));
  connect(poles_, &QTableWidget::cellClicked, this, [this](int row, int) { centreOnPole(row); });
  facts_ = new QLabel;
  facts_->setObjectName("rmatrixFacts");
  facts_->setWordWrap(true);
  facts_->setProperty("role", "hint");
  panel_ = new PlotPanel("rmatrixPanel");
  QPushButton *drawButton = new QPushButton(tr("Calculate"));
  drawButton->setObjectName("rmatrixDraw");
  drawButton->setProperty("primary", true);
  connect(drawButton, &QPushButton::clicked, this, &RMatrixPage::compute);
  connect(group_, &QComboBox::currentIndexChanged, this, &RMatrixPage::groupChanged);
  connect(energies_, &SweepFields::changed, this, &RMatrixPage::compute);
  connect(quantity_, &QComboBox::currentIndexChanged, this, [this]() {
    const auto p = PlotNamed(PlotsOf("Block", false), quantity_->currentData().toString());
    if (p) panel_->logY()->setChecked(p->logY);
    fillSeries();
    redraw();
  });
  connect(series_, &QListWidget::itemChanged, this, &RMatrixPage::redraw);

  QWidget *controls = new QWidget;
  QFormLayout *f = new QFormLayout(controls);
  f->addRow(tr("J^pi group"), group_);
  f->addRow(tr("E_x"), energies_);
  f->addRow(tr("Draw"), quantity_);
  f->addRow(tr("Lines"), series_);
  f->addRow(drawButton);
  f->addRow(facts_);
  QWidget *poleBox = new QWidget;
  QVBoxLayout *pv = new QVBoxLayout(poleBox);
  pv->setContentsMargins(0, 0, 0, 0);
  QLabel *poleTitle = new QLabel(tr("Poles of the group (click: centre on it)"));
  poleTitle->setProperty("role", "hint");
  pv->addWidget(poleTitle);
  pv->addWidget(poles_);
  QSplitter *right = new QSplitter(Qt::Vertical);
  right->addWidget(panel_);
  right->addWidget(poleBox);
  right->setStretchFactor(0, 3);
  right->setStretchFactor(1, 1);
  QSplitter *split = new QSplitter;
  split->addWidget(Scrolled(controls));
  split->addWidget(right);
  split->setStretchFactor(1, 3);
  setContent(split);
}

void RMatrixPage::infoChanged() {
  FillGroups(group_, context_->info().groups);
  answer_ = front::PhysicsCurves();
  answerGroup_ = 0;
  if (group_->count()) groupChanged();
}

void RMatrixPage::groupChanged() {
  const auto range = GroupRange(context_->info(), group_->currentData().toInt());
  energies_->set(range.first, range.second, energies_->points->value(), false);
  fillPoles();
  compute();
}

void RMatrixPage::fillPoles() {
  const int jgroup = group_->currentData().toInt();
  const PhysicsContext::Info &info = context_->info();
  const front::PhysicsGroup *g = nullptr;
  for (const front::PhysicsGroup &x : info.groups)
    if (x.jgroup == jgroup) g = &x;
  poles_->clear();
  QStringList headers = {tr("J^pi"), tr("level"), tr("E_x (MeV)"), tr("Gamma (MeV)")};
  if (g)
    for (const front::PhysicsChannel &c : g->channels)
      headers << tr("c%1 %2 l=%3").arg(c.channel).arg(c.radiation == "P" ? context_->pairLabel(c.pair) : QString::fromStdString(c.radiation)).arg(c.l);
  poles_->setColumnCount(headers.size());
  poles_->setHorizontalHeaderLabels(headers);
  poles_->setRowCount(0);
  for (const front::PhysicsPole &p : info.poles) {
    if (p.jgroup != jgroup) continue;
    const int row = poles_->rowCount();
    poles_->insertRow(row);
    auto set = [&](int col, const QString &text, double value) {
      QTableWidgetItem *item = new QTableWidgetItem(text);
      item->setData(Qt::UserRole, value);
      poles_->setItem(row, col, item);
    };
    set(0, JPi(p.J, p.parity), p.J);
    set(1, QString::number(p.level), p.level);
    set(2, Number(p.energy), p.energy);
    set(3, Number(p.totalWidth), p.totalWidth);
    for (const front::PhysicsPole::Width &w : p.widths) {
      const int col = 3 + w.channel;
      if (col >= poles_->columnCount()) continue;
      const QString text = w.kind == "width" ? Number(w.width)
                           : w.kind == "anc" ? tr("ANC %1 fm^-1/2").arg(Number(w.value))
                                             : QString("%1 (%2)").arg(Number(w.value), QString::fromStdString(w.kind));
      set(col, text, w.kind == "width" ? w.width : w.value);
      poles_->item(row, col)->setToolTip(tr("%1; value as parameters.out writes it: %2; P = %3")
                                             .arg(QString::fromStdString(w.kind), Number(w.value), Number(w.penetrability)));
    }
  }
  poles_->resizeColumnsToContents();
}

void RMatrixPage::compute() {
  const int jgroup = group_->currentData().toInt();
  if (jgroup <= 0) return;
  QVector<double> grid;
  QString why;
  if (!energies_->grid(grid, why)) {
    panel_->setStatus(tr("Not drawn: %1").arg(why));
    return;
  }
  auto answer = std::make_shared<front::PhysicsCurves>();
  const std::vector<double> e(grid.begin(), grid.end());
  context_->run(
      "rmatrix",
      [=](front::PhysicsSession &s, const std::atomic<bool> &, std::string &w) { return s.Block(jgroup, e, *answer, w); },
      [this, answer, jgroup](bool done, const QString &w) {
        if (!done) {
          panel_->setStatus(tr("Not drawn: %1").arg(w));
          return;
        }
        answer_ = *answer;
        answerGroup_ = jgroup;
        facts_->setText(ShortFactsText(answer_));
        facts_->setToolTip(FactsText(answer_));
        fillSeries();
        redraw();
      });
}

void RMatrixPage::fillSeries() {
  const auto p = PlotNamed(PlotsOf("Block", false), quantity_->currentData().toString());
  QStringList unchecked;
  for (int i = 0; i < series_->count(); i++)
    if (series_->item(i)->checkState() != Qt::Checked) unchecked << series_->item(i)->text();
  const QSignalBlocker block(series_);
  series_->clear();
  if (!p) return;
  for (int k = 0; k < (int)answer_.series.size(); k++) {
    const front::PhysicsSeries &s = answer_.series[k];
    if (s.name != p->y) continue;
    if (p->diagonalOnly && s.Key("channel") != s.Key("channel_out")) continue;
    QListWidgetItem *item = new QListWidgetItem(QString::fromStdString(s.label));
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    item->setCheckState(unchecked.contains(item->text()) ? Qt::Unchecked : Qt::Checked);
    item->setData(Qt::UserRole, k);
    series_->addItem(item);
  }
}

void RMatrixPage::redraw() {
  const auto p = PlotNamed(PlotsOf("Block", false), quantity_->currentData().toString());
  if (!p) return;
  QList<PlotSeries> curves;
  // the group's poles: dashed lines on an energy axis, markers on the Argand curves
  QVector<double> poleE;
  QStringList poleLabels;
  for (const front::PhysicsPole &pole : context_->info().poles)
    if (pole.jgroup == answerGroup_) poleE << pole.energy, poleLabels << tr("%1 MeV").arg(Number(pole.energy));
  const QVector<double> x = Values(answer_.axisValues);
  for (int i = 0; i < series_->count(); i++) {
    if (series_->item(i)->checkState() != Qt::Checked) continue;
    const front::PhysicsSeries &s = answer_.series[series_->item(i)->data(Qt::UserRole).toInt()];
    PlotSeries ps;
    ps.label = QString::fromStdString(s.label);
    ps.color = CurveColour(i);
    ps.marker = i % 5;
    if (p->x.empty()) {
      ps.modelX = x;
      ps.modelY = Values(s.values);
      if (p->y == "delta") ps.modelY = UnwrappedPhaseDegrees(ps.modelY);
    } else {
      const front::PhysicsSeries *sx = answer_.Find(p->x, s.keys);
      if (!sx) continue;
      ps.modelX = Values(sx->values);
      ps.modelY = Values(s.values);
      // a marker at each pole: the grid point nearest its energy
      for (double e : poleE) {
        if (x.isEmpty()) break;
        int best = 0;
        for (int k = 1; k < x.size(); k++)
          if (std::fabs(x[k] - e) < std::fabs(x[best] - e)) best = k;
        if (e < x.front() || e > x.back()) continue;
        ps.dataX << ps.modelX[best];
        ps.dataY << ps.modelY[best];
      }
    }
    curves << ps;
  }
  const QString xTitle = p->x.empty() ? QString::fromStdString(CallNamed("Block").axisTitle) : QString::fromStdString(p->xTitle);
  panel_->show(curves, xTitle, QString::fromStdString(p->yTitle));
  if (p->x.empty()) panel_->plot()->setMarkers(poleE, poleLabels);
  else panel_->plot()->setMarkers({});
  panel_->setStatus(p->x.empty() ? QString() : tr("The markers sit at the grid energy nearest each level of the group."));
}

void RMatrixPage::centreOnPole(int row) {
  QTableWidgetItem *e = poles_->item(row, 2), *g = poles_->item(row, 3);
  if (!e || answer_.axisValues.empty()) return;
  const auto p = PlotNamed(PlotsOf("Block", false), quantity_->currentData().toString());
  if (!p || !p->x.empty()) return;  // the Argand plot has no energy axis
  const double energy = e->data(Qt::UserRole).toDouble(), width = g ? g->data(Qt::UserRole).toDouble() : 0;
  const double span = answer_.axisValues.back() - answer_.axisValues.front();
  const double half = std::max(3 * width, 0.02 * std::fabs(span));
  const QRectF v = panel_->plot()->view();
  panel_->plot()->zoomTo(QRectF(QPointF(energy - half, v.top()), QPointF(energy + half, v.bottom())));
}

// ---------------------------------------------------------------- Pole map

PoleMapPage::PoleMapPage(PhysicsContext *context, QWidget *parent) : ProjectPage(context, parent) {
  setObjectName("poleMapPage");
  // one line on the page, the reading of the map in its tip (qt_gui.md section 19.6)
  explanation_ = new QLabel(tr("A real-energy map along E_x (the complex plane is not available)."));
  explanation_->setToolTip(
      tr("A real-energy map: the excitation energy E_x runs along the real axis, each row is a J^pi group (log10 "
         "|det A^-1|), a pair of channels (|U_cc'|) or a pair of levels (|A|) of one group.  A minimum of |det A^-1| "
         "(a dark valley of log10 |det A^-1|) is where the inverse level matrix is nearly singular on the real axis: "
         "the energy of a narrow resonance, close to the level's physical energy (the white points).  |U| and |A| peak "
         "there.  The complex-energy plane, where the poles of U are the zeros of det A^-1 continued in E, is not "
         "available: the engine's Coulomb functions and level matrix take real energies only (docs/upgrade/"
         "physics_tab.md, section 5.3)."));
  explanation_->setObjectName("poleMapExplanation");
  explanation_->setWordWrap(true);
  explanation_->setProperty("role", "hint");
  kind_ = new QComboBox;
  kind_->setObjectName("poleMapKind");
  FillPlots(kind_, PlotsOf("Block", true));
  group_ = new QComboBox;
  group_->setObjectName("poleMapGroup");
  log_ = new QCheckBox(tr("log colour"));
  log_->setObjectName("poleMapLog");
  log_->setToolTip(tr("Colour by log10 of the value"));
  energies_ = new SweepFields("MeV (E_x)");
  energies_->setObjectName("poleMapEnergies");
  energies_->set(1, 5, 400, false);
  QPushButton *drawButton = new QPushButton(tr("Calculate"));
  drawButton->setObjectName("poleMapDraw");
  drawButton->setProperty("primary", true);
  map_ = new HeatmapWidget;
  map_->setObjectName("poleMap");
  cut_ = new PlotPanel("poleMapCut");
  readout_ = new QLabel;
  readout_->setObjectName("poleMapReadout");
  readout_->setProperty("role", "hint");
  connect(map_, &HeatmapWidget::cursorMoved, readout_, &QLabel::setText);
  connect(map_, &HeatmapWidget::rowClicked, this, &PoleMapPage::showCut);
  connect(kind_, &QComboBox::currentIndexChanged, this, &PoleMapPage::kindChanged);
  connect(group_, &QComboBox::currentIndexChanged, this, &PoleMapPage::compute);
  connect(log_, &QCheckBox::toggled, map_, &HeatmapWidget::setLogColour);
  connect(drawButton, &QPushButton::clicked, this, &PoleMapPage::compute);
  QPushButton *reset = new QPushButton(tr("Reset zoom"));
  connect(reset, &QPushButton::clicked, map_, &HeatmapWidget::resetZoom);
  QPushButton *exportButton = new QPushButton(tr("Export map..."));
  connect(exportButton, &QPushButton::clicked, this, [this]() {
    const QString file = QFileDialog::getSaveFileName(this, tr("Export the map"), QString(),
                                                      tr("PNG image (*.png);;SVG drawing (*.svg);;PDF (*.pdf);;Values as CSV (*.csv)"));
    QString why;
    if (!file.isEmpty() && !map_->exportTo(file, why)) QMessageBox::warning(this, tr("Export"), why);
  });

  QWidget *top = new QWidget;
  QVBoxLayout *t = new QVBoxLayout(top);
  t->setContentsMargins(0, 0, 0, 0);
  QHBoxLayout *row1 = new QHBoxLayout;
  row1->addWidget(new QLabel(tr("Map")));
  row1->addWidget(kind_, 1);
  row1->addWidget(new QLabel(tr("J^pi group")));
  row1->addWidget(group_, 1);
  row1->addWidget(log_);
  QHBoxLayout *row2 = new QHBoxLayout;
  row2->addWidget(new QLabel(tr("E_x")));
  row2->addWidget(energies_, 1);
  row2->addWidget(drawButton);
  row2->addWidget(reset);
  row2->addWidget(exportButton);
  t->addWidget(explanation_);
  t->addLayout(row1);
  t->addLayout(row2);
  QWidget *mapBox = new QWidget;
  QVBoxLayout *mv = new QVBoxLayout(mapBox);
  mv->setContentsMargins(0, 0, 0, 0);
  mv->addWidget(map_, 1);
  mv->addWidget(readout_);
  QSplitter *split = new QSplitter(Qt::Vertical);
  split->addWidget(mapBox);
  split->addWidget(cut_);
  split->setStretchFactor(0, 3);
  split->setStretchFactor(1, 2);
  QWidget *content = new QWidget;
  QVBoxLayout *c = new QVBoxLayout(content);
  c->setContentsMargins(0, 0, 0, 0);
  c->addWidget(top);
  c->addWidget(split, 1);
  setContent(content);
  kindChanged();
}

void PoleMapPage::infoChanged() {
  FillGroups(group_, context_->info().groups);
  const auto range = GroupRange(context_->info(), 0);
  energies_->set(range.first, range.second, energies_->points->value(), false);
  compute();  // a few Block calls: milliseconds; the map then matches the fields whether shown or not
}

void PoleMapPage::kindChanged() {
  const auto p = PlotNamed(PlotsOf("Block", true), kind_->currentData().toString());
  group_->setEnabled(p && p->mapRows == "series");
  {
    const QSignalBlocker block(log_);
    log_->setChecked(p && p->logY);
  }
  if (isVisible()) compute();
}

void PoleMapPage::compute() {
  const auto pp = PlotNamed(PlotsOf("Block", true), kind_->currentData().toString());
  if (!pp || context_->info().groups.empty()) return;
  const front::PhysicsPlot plot = *pp;
  QVector<double> grid;
  QString why;
  if (!energies_->grid(grid, why)) {
    readout_->setText(tr("Not drawn: %1").arg(why));
    return;
  }
  std::vector<int> groups;
  if (plot.mapRows == "jgroup")
    for (const front::PhysicsGroup &g : context_->info().groups) groups.push_back(g.jgroup);
  else
    groups.push_back(group_->currentData().toInt());
  auto answers = std::make_shared<std::vector<front::PhysicsCurves>>();
  const std::vector<double> e(grid.begin(), grid.end());
  context_->run(
      "polemap",
      [=](front::PhysicsSession &s, const std::atomic<bool> &cancel, std::string &w) {
        for (int g : groups) {
          if (cancel) return false;
          front::PhysicsCurves c;
          if (!s.Block(g, e, c, w)) return false;
          answers->push_back(std::move(c));
        }
        return true;
      },
      [this, answers, plot, groups, grid](bool done, const QString &w) {
        if (!done) {
          readout_->setText(tr("Not drawn: %1").arg(w));
          return;
        }
        const PhysicsContext::Info &info = context_->info();
        QStringList rows, tips;
        QVector<QVector<double>> values;
        QVector<QPointF> points;
        QVector<double> vertical;
        for (std::size_t a = 0; a < answers->size(); a++)
          for (const front::PhysicsSeries &s : (*answers)[a].series) {
            if (s.name != plot.y) continue;
            if (plot.mapRows == "jgroup") {
              QString jpi;
              for (const front::PhysicsGroup &g : info.groups)
                if (g.jgroup == groups[a]) jpi = JPi(g.J, g.parity);
              rows << jpi;
              tips << tr("%1 (group %2)").arg(jpi).arg(groups[a]);
              for (const front::PhysicsPole &p : info.poles)
                if (p.jgroup == groups[a]) points << QPointF(p.energy, rows.size() - 1);
            } else if (s.Key("channel") >= 0) {
              rows << QString("c%1-c%2").arg(s.Key("channel")).arg(s.Key("channel_out"));
              tips << QString::fromStdString(s.label);
            } else {
              rows << QString("L%1-L%2").arg(s.Key("level")).arg(s.Key("level_out"));
              tips << QString::fromStdString(s.label);
            }
            values << Values(s.values);
          }
        if (plot.mapRows != "jgroup")
          for (const front::PhysicsPole &p : info.poles)
            if (p.jgroup == groups[0]) vertical << p.energy;
        map_->setMap(grid, rows, values, QString::fromStdString(CallNamed("Block").axisTitle),
                     QString::fromStdString(plot.yTitle));
        map_->setRowTips(tips);
        map_->setLogColour(log_->isChecked());
        map_->setPoints(points);
        map_->setVerticalMarkers(vertical);
        cutMarkers_ = vertical;
        readout_->setText(rows.isEmpty() ? tr("No row: the group has no such series (no open particle channel?)")
                                         : tr("Click a row for its cut below."));
        showCut(0);
      });
}

void PoleMapPage::showCut(int row) {
  if (row < 0 || row >= map_->rows().size()) {
    cut_->show({}, QString(), QString());
    return;
  }
  PlotSeries s;
  s.label = map_->rows()[row];
  s.color = CurveColour(0);
  s.modelX = map_->x();
  s.modelY = map_->values()[row];
  cut_->logY()->setChecked(log_->isChecked());
  cut_->show({s}, QString::fromStdString(CallNamed("Block").axisTitle), map_->valueTitle());
  // the poles on the cut: the row's own points (the det map), or the group's lines (a map of one group)
  QVector<double> marks = cutMarkers_;
  for (const QPointF &pt : map_->points())
    if ((int)std::lround(pt.y()) == row) marks << pt.x();
  cut_->plot()->setMarkers(marks);
  cut_->setStatus(tr("Cut through row %1").arg(s.label));
}

// ---------------------------------------------------------------- THM

ThmPage::ThmPage(PhysicsContext *context, QWidget *parent) : ProjectPage(context, parent) {
  setObjectName("thmPage");
  experiment_ = new QComboBox;
  experiment_->setObjectName("thmExperiment");
  part_ = new QComboBox;
  part_->setObjectName("thmPart");
  for (const front::PhysicsCall &c : front::PhysicsCalls())
    if (!front::PhysicsPlots(c.name).empty() && c.name != "Block" && c.name != "ChannelFunctions") {
      part_->addItem(QString::fromStdString(c.title), QString::fromStdString(c.name));
      part_->setItemData(part_->count() - 1, QString::fromStdString(c.doc), Qt::ToolTipRole);
    }
  quantity_ = new QComboBox;
  quantity_->setObjectName("thmQuantity");
  energies_ = new SweepFields("MeV (c.m. x + A)");
  energies_->setObjectName("thmEnergies");
  energies_->set(0.05, 1.5, 100, false);
  facts_ = new QLabel;
  facts_->setObjectName("thmFacts");
  facts_->setWordWrap(true);
  facts_->setProperty("role", "hint");
  facts_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  panel_ = new PlotPanel("thmPanel");
  connect(energies_->log, &QCheckBox::toggled, panel_->logX(), &QCheckBox::setChecked);
  QPushButton *drawButton = new QPushButton(tr("Draw"));
  drawButton->setObjectName("thmDraw");
  drawButton->setProperty("primary", true);
  QPushButton *add = new QPushButton(tr("Add curve"));
  QPushButton *clearButton = new QPushButton(tr("Clear"));
  connect(drawButton, &QPushButton::clicked, this, [this]() { draw(false); });
  connect(add, &QPushButton::clicked, this, [this]() { draw(true); });
  connect(clearButton, &QPushButton::clicked, this, [this]() {
    curves_.clear();
    redraw();
  });
  connect(part_, &QComboBox::currentIndexChanged, this, &ThmPage::partChanged);
  connect(quantity_, &QComboBox::currentIndexChanged, this, [this]() {
    const auto p = PlotNamed(front::PhysicsPlots(part_->currentData().toString().toStdString()),
                                            quantity_->currentData().toString());
    if (p) panel_->logY()->setChecked(p->logY);
    redraw();
  });
  QWidget *controls = new QWidget;
  QFormLayout *f = new QFormLayout(controls);
  f->addRow(tr("Experiment"), experiment_);
  f->addRow(tr("Part"), part_);
  f->addRow(tr("E"), energies_);
  f->addRow(tr("Draw"), quantity_);
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(drawButton);
  buttons->addWidget(add);
  buttons->addWidget(clearButton);
  buttons->addStretch(1);
  f->addRow(buttons);
  f->addRow(facts_);
  QSplitter *split = new QSplitter;
  split->addWidget(Scrolled(controls));
  split->addWidget(panel_);
  split->setStretchFactor(1, 3);
  setContent(split);
  partChanged();
}

void ThmPage::infoChanged() {
  const PhysicsContext::Info &info = context_->info();
  if (info.ok && info.experiments.empty()) {
    setNeeds(tr("The project has no THM experiment (<thm> experiment[...]): the THM tab adds one."));
    return;
  }
  const QSignalBlocker block(experiment_);
  const QString was = experiment_->currentText();
  experiment_->clear();
  for (const std::string &e : info.experiments) experiment_->addItem(QString::fromStdString(e));
  const int i = experiment_->findText(was);
  experiment_->setCurrentIndex(i >= 0 ? i : 0);
}

void ThmPage::partChanged() {
  const std::string call = part_->currentData().toString().toStdString();
  FillPlots(quantity_, front::PhysicsPlots(call));
  if (const auto p = PlotNamed(front::PhysicsPlots(call), quantity_->currentData().toString()))
    panel_->logY()->setChecked(p->logY);
  for (const front::PhysicsCall &c : front::PhysicsCalls())
    if (c.name == call) energies_->points->setValue(c.slow && call == "ThmDistortion" ? 30 : 100);
  curves_.clear();
  redraw();
}

void ThmPage::draw(bool add) {
  QVector<double> grid;
  QString why;
  if (experiment_->count() == 0) return;
  if (!energies_->grid(grid, why)) {
    panel_->setStatus(tr("Not drawn: %1").arg(why));
    return;
  }
  const std::string experiment = experiment_->currentText().toStdString(), call = part_->currentData().toString().toStdString();
  auto curve = std::make_shared<Curve>();
  curve->label = experiment_->currentText();
  const std::vector<double> e(grid.begin(), grid.end());
  context_->run(
      "thm",
      [=](front::PhysicsSession &s, const std::atomic<bool> &cancel, std::string &w) {
        if (call == "ThmVertex") return s.ThmVertex(experiment, e, curve->answer, w);
        if (call == "ThmLineshape") return s.ThmLineshape(experiment, e, curve->answer, w);
        // the distortion: a DWBA per energy; in chunks, so that Cancel acts between them
        const std::size_t chunk = 10;
        for (std::size_t i = 0; i < e.size(); i += chunk) {
          if (cancel) return false;
          front::PhysicsCurves part;
          const std::vector<double> sub(e.begin() + i, e.begin() + std::min(e.size(), i + chunk));
          if (!s.ThmDistortion(experiment, sub, part, w)) return false;
          if (i == 0) {
            curve->answer = part;
            continue;
          }
          curve->answer.axisValues.insert(curve->answer.axisValues.end(), part.axisValues.begin(), part.axisValues.end());
          for (std::size_t k = 0; k < curve->answer.series.size() && k < part.series.size(); k++)
            curve->answer.series[k].values.insert(curve->answer.series[k].values.end(), part.series[k].values.begin(),
                                                  part.series[k].values.end());
        }
        return true;
      },
      [this, curve, add](bool done, const QString &w) {
        if (!done) {
          panel_->setStatus(tr("Not drawn: %1").arg(w));
          return;
        }
        if (!add) curves_.clear();
        curves_ << *curve;
        facts_->setText(FactsText(curve->answer));
        redraw();
      });
}

void ThmPage::redraw() {
  const std::string call = part_->currentData().toString().toStdString();
  const std::vector<front::PhysicsPlot> plots = front::PhysicsPlots(call);
  const auto p = PlotNamed(plots, quantity_->currentData().toString());
  QList<PlotSeries> series;
  int n = 0;
  if (p)
    for (const Curve &c : curves_)
      for (const front::PhysicsSeries &s : c.answer.series) {
        if (s.name != p->y) continue;
        PlotSeries ps;
        ps.label = curves_.size() > 1 ? c.label + ": " + QString::fromStdString(s.label) : QString::fromStdString(s.label);
        ps.color = CurveColour(n++);
        ps.modelX = Values(c.answer.axisValues);
        ps.modelY = Values(s.values);
        series << ps;
      }
  QString xTitle;
  for (const front::PhysicsCall &c : front::PhysicsCalls())
    if (c.name == call) xTitle = QString::fromStdString(c.axisTitle);
  panel_->show(series, xTitle, p ? QString::fromStdString(p->yTitle) : QString());
  panel_->setStatus(p && series.isEmpty() && !curves_.isEmpty() ? tr("This experiment's answer has no %1 (another vertex model?)").arg(QString::fromStdString(p->y)) : QString());
}

}  // namespace physicsqt
