#include "ThmDiagnosticsPage.h"

#include <QComboBox>
#include <QFileInfo>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QFrame>
#include <QHash>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QVBoxLayout>
#include <cmath>

#include "ThmPlotWidget.h"

namespace {

ThmPlotWidget *makePlot(const QString &x, const QString &y) {
  ThmPlotWidget *p = new ThmPlotWidget;
  p->setAxisLabels(x, y);
  return p;
}

const char *kEcm = "E_cm (MeV)";

}  // namespace

ThmDiagnosticsPage::ThmDiagnosticsPage(std::function<QString(int, ThmDiagnosticsRequest &)> prepare,
                                       std::function<QList<Target>()> targets, QWidget *parent) :
  QWidget(parent),
  prepare_(prepare),
  targets_(targets) {
  segmentCombo = new QComboBox;
  segmentCombo->setToolTip(tr("The THM data segment to look at; its experiment, if any, in front."));
  segmentCombo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
  segmentCombo->setMinimumContentsLength(24);
  computeButton = new QPushButton(tr("Compute"));
  computeButton->setToolTip(
      tr("Runs AZURE2 on the project as this workspace would leave it (a temporary copy; nothing is written) at the "
         "current parameters, on a grid over the segment's data."));
  connect(computeButton, &QPushButton::clicked, this, [this]() { compute(); });
  busyBar = new QProgressBar;
  busyBar->setRange(0, 0);  // busy indicator
  busyBar->setMaximumWidth(100);
  busyBar->setTextVisible(false);
  busyBar->hide();
  vertexGroupCombo = new QComboBox;
  vertexGroupCombo->setToolTip(tr("The J^pi group whose entrance channels the vertex panel shows."));
  vertexGroupCombo->setSizeAdjustPolicy(QComboBox::AdjustToContents);
  connect(vertexGroupCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(drawVertex()));
  statusLabel = new QLabel(tr("Press Compute to evaluate the selected segment."));
  statusLabel->setWordWrap(true);
  statusLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);

  // One toolbar row: what to compute, and the vertex group to show.
  QHBoxLayout *top = new QHBoxLayout;
  top->setSpacing(8);
  top->addWidget(new QLabel(tr("Segment:")));
  top->addWidget(segmentCombo, 1);
  top->addSpacing(4);
  top->addWidget(new QLabel(QString::fromUtf8("J<sup>π</sup>:")));
  top->addWidget(vertexGroupCombo);
  top->addSpacing(4);
  top->addWidget(computeButton);
  top->addWidget(busyBar);

  vertexPlot = makePlot(tr(kEcm), QString::fromUtf8("|M_l|²"));
  vertexPlot->setToolTip(
      tr("|M_l(E)|^2 = |(B_c - 1) j_l(pa) - pa j_l'(pa) [+ C_l]|^2 at the quasi-free point; dashed vertical lines: "
         "its nodes, where a resonance is suppressed. With a spectator window (ps=), dashed curves: <|M_l|^2> over "
         "the window, as AZURE2 uses it."));
  hoesPlot = makePlot(tr(kEcm), QString::fromUtf8("σ (b)"));
  hoesPlot->setLogY(true);
  hoesPlot->setToolTip(tr("The HOES cross section of the segment's channel (quasi-free point; no resolution, weight "
                          "or line shape; scaled to the other curve) and the on-shell angle-integrated one."));
  lineshapePlot = makePlot(tr(kEcm), QString::fromUtf8("|N_C|²"));
  lineshapePlot->setToolTip(tr("|N_C|^2 = exp[2 zeta arctan(2 (E_lambda - E)/Gamma_lambda)] of the levels with a "
                               "pole near the data (lineshape=on); 1 at the pole."));
  zetaPlot = makePlot(tr(kEcm), QString::fromUtf8("ζ"));
  zetaPlot->setToolTip(tr("zeta(E) = eta_sB - eta_0 of the segment's exit pair (lineshape=on)."));
  weightPlot = makePlot(tr(kEcm), tr("w"));
  weightPlot->setToolTip(tr("The weight table of the segment (weight[k]), as the engine interpolates it."));
  windowPlot = makePlot(tr("p_s (MeV/c)"), tr("w (per MeV/c)"));
  windowPlot->setToolTip(tr("w(p) = |phi(p)|^2 p^2 over [p_min, p_max] (a table: as given), unit area; dots: the "
                            "Gauss-Legendre nodes at which AZURE2 evaluates the vertex."));

  // Each plot in a card: a framed panel with a short bold title.
  auto card = [this](ThmPlotWidget *plot, const QString &title) {
    QFrame *f = new QFrame;
    f->setFrameShape(QFrame::StyledPanel);
    f->setFrameShadow(QFrame::Plain);
    f->setAutoFillBackground(true);
    f->setBackgroundRole(QPalette::Base);
    QLabel *t = new QLabel(title);
    QFont bold = t->font();
    bold.setBold(true);
    t->setFont(bold);
    t->setToolTip(plot->toolTip());
    QVBoxLayout *l = new QVBoxLayout;
    l->setContentsMargins(8, 6, 8, 6);
    l->setSpacing(2);
    l->addWidget(t);
    l->addWidget(plot, 1);
    f->setLayout(l);
    f->setMinimumHeight(220);
    cards_[plot] = f;
    return f;
  };
  vertexPanel = card(vertexPlot, QString::fromUtf8("Entrance vertex |M<sub>l</sub>|²"));
  card(hoesPlot, tr("HOES vs on-shell"));
  card(lineshapePlot, QString::fromUtf8("Line shape |N<sub>C</sub>|²"));
  card(zetaPlot, QString::fromUtf8("Line shape ζ(E)"));
  card(weightPlot, tr("Weight w(E)"));
  card(windowPlot, QString::fromUtf8("Spectator window w(p<sub>s</sub>)"));

  QWidget *panels = new QWidget;
  grid_ = new QGridLayout;
  grid_->setContentsMargins(0, 0, 0, 0);
  grid_->setSpacing(8);
  panels->setLayout(grid_);
  QScrollArea *scroll = new QScrollArea;
  scroll->setWidget(panels);
  scroll->setWidgetResizable(true);
  scroll->setFrameShape(QFrame::NoFrame);
  scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

  QVBoxLayout *main = new QVBoxLayout;
  main->setSpacing(8);
  main->addLayout(top);
  main->addWidget(statusLabel);
  main->addWidget(scroll, 1);
  setLayout(main);

  clearPlots(tr("not computed yet"));
  layoutPanels();
  refreshTargets();
}

ThmDiagnosticsPage::~ThmDiagnosticsPage() {
  if (thread_) {
    thread_->wait();  // the engine cannot be interrupted; it holds no GUI object
    delete thread_;
  }
}

void ThmDiagnosticsPage::refreshTargets() {
  const int current = segmentCombo->currentData().toInt();
  list_ = targets_();
  segmentCombo->blockSignals(true);
  segmentCombo->clear();
  for (const Target &t : list_) segmentCombo->addItem(t.text, t.segment);
  int at = segmentCombo->findData(current);
  segmentCombo->setCurrentIndex(at >= 0 ? at : (list_.isEmpty() ? -1 : 0));
  segmentCombo->blockSignals(false);
  computeButton->setEnabled(!list_.isEmpty() && !busy());
  if (list_.isEmpty()) statusLabel->setText(tr("No active THM data segment."));
  // Results of settings that have changed since are marked as such.
  if (!computedText_.isEmpty() && !busy()) {
    ThmDiagnosticsRequest now;
    if (prepare_(computedSegment_, now).isEmpty() && now.projectText != computedText_)
      statusLabel->setText(tr("Changed since the last Compute: press Compute again."));
  }
}

QString ThmDiagnosticsPage::reactionOf(int segment) const {
  for (const Target &t : list_)
    if (t.segment == segment) return t.reaction;
  return QString();
}

void ThmDiagnosticsPage::clearPlots(const QString &message) {
  for (ThmPlotWidget *p : {vertexPlot, hoesPlot, lineshapePlot, zetaPlot, weightPlot, windowPlot}) {
    p->clear();
    p->setMessage(message);
  }
  for (ThmPlotWidget *p : {lineshapePlot, zetaPlot, weightPlot, windowPlot}) cards_[p]->hide();
  vertexGroupCombo->clear();
  statusLabel->setToolTip(QString());
}

void ThmDiagnosticsPage::compute() {
  if (busy() || segmentCombo->currentIndex() < 0) return;
  ThmDiagnosticsRequest request;
  const QString why = prepare_(segmentCombo->currentData().toInt(), request);
  if (!why.isEmpty()) {
    statusLabel->setText(why);
    return;
  }
  computedText_ = request.projectText;
  computedSegment_ = request.segment;
  thread_ = new ThmDiagnosticsThread(request);
  connect(thread_, SIGNAL(finished()), this, SLOT(threadFinished()));
  computeButton->setEnabled(false);
  segmentCombo->setEnabled(false);
  busyBar->show();
  statusLabel->setText(tr("Computing segment %1 ...").arg(request.segment));
  thread_->start();
}

bool ThmDiagnosticsPage::computeNow() {
  if (busy() || segmentCombo->currentIndex() < 0) return false;
  ThmDiagnosticsRequest request;
  const QString why = prepare_(segmentCombo->currentData().toInt(), request);
  if (!why.isEmpty()) {
    statusLabel->setText(why);
    return false;
  }
  computedText_ = request.projectText;
  computedSegment_ = request.segment;
  showResult(ComputeThmDiagnostics(request));
  return result_.error.isEmpty();
}

void ThmDiagnosticsPage::threadFinished() {
  ThmDiagnosticsThread *t = thread_;
  thread_ = nullptr;
  busyBar->hide();
  segmentCombo->setEnabled(true);
  computeButton->setEnabled(!list_.isEmpty());
  showResult(t->result());
  t->deleteLater();
  emit computed();
}

void ThmDiagnosticsPage::showResult(const ThmDiagnosticsResult &result) {
  result_ = result;
  if (!result.error.isEmpty()) {
    clearPlots(tr("no result"));
    statusLabel->setText(result.error);
    return;
  }
  clearPlots(QString());
  const QString reaction = reactionOf(result.segment);
  // One short line; the details in its tooltip.
  statusLabel->setText(tr("Segment %1%2 \u00b7 %3 points \u00b7 %4 \u2026 %5 MeV")
                           .arg(result.segment)
                           .arg(result.experiment.isEmpty() ? QString() : tr(" (%1)").arg(result.experiment))
                           .arg(result.energy.size())
                           .arg(result.eLo, 0, 'g', 4)
                           .arg(result.eHi, 0, 'g', 4));
  statusLabel->setToolTip(
      tr("Segment %1%2: %3 points from %4 to %5 MeV; vertex %6, B + T_s = %7 MeV.")
          .arg(result.segment)
          .arg(result.experiment.isEmpty() ? QString() : tr(" (experiment %1)").arg(result.experiment))
          .arg(result.energy.size())
          .arg(result.eLo, 0, 'g', 4)
          .arg(result.eHi, 0, 'g', 4)
          .arg(result.vertexMode)
          .arg(result.binding, 0, 'g', 6) +
      (result.window ? tr(" Spectator window: %1 nodes, <T_s> = %2 MeV.")
                           .arg(result.nodeP.size())
                           .arg(result.meanTs, 0, 'g', 6)
                     : QString()) +
      (result.lineshape && result.nc2Hidden
           ? tr(" Line shape: %1 level(s) with a pole outside the data range, or beyond the first four, not drawn.")
                 .arg(result.nc2Hidden)
           : QString()));

  // Entrance vertex: the J^pi groups in the combo; drawVertex draws one.
  vertexGroupCombo->blockSignals(true);
  for (const ThmDiagnosticsResult::VertexGroup &g : result.vertex) vertexGroupCombo->addItem(g.jpi);
  vertexGroupCombo->blockSignals(false);
  vertexPlot->setTitle(reaction);
  drawVertex();

  // HOES and on-shell.
  hoesPlot->setTitle(reaction);
  ThmPlotWidget::Series hoes;
  hoes.x = result.hoesEnergy;
  hoes.y = result.hoes;
  hoes.color = thmPlotColor(0);
  hoes.label = result.onShell.isEmpty() ? tr("HOES") : tr("HOES (scaled)");
  hoesPlot->addSeries(hoes);
  if (!result.onShell.isEmpty()) {
    ThmPlotWidget::Series on;
    on.x = result.onShellEnergy;
    on.y = result.onShell;
    on.color = thmPlotColor(1);
    on.style = Qt::DashLine;
    on.label = tr("on-shell");
    hoesPlot->addSeries(on);
  }

  // Line shape.
  if (result.lineshape) {
    cards_[lineshapePlot]->show();
    cards_[zetaPlot]->show();
    lineshapePlot->setTitle(reaction);
    zetaPlot->setTitle(reaction);
    for (int i = 0; i < result.nc2.size(); i++) {
      ThmPlotWidget::Series s;
      s.x = result.energy;
      s.y = result.nc2[i].y;
      s.color = thmPlotColor(i);
      s.label = result.nc2[i].label;
      lineshapePlot->addSeries(s);
    }
    if (result.nc2.isEmpty()) lineshapePlot->setMessage(tr("no level with a pole inside the data range"));
    ThmPlotWidget::Series z;
    z.x = result.energy;
    z.y = result.zeta;
    z.color = thmPlotColor(0);
    zetaPlot->addSeries(z);
  }

  // Weight.
  if (!result.weight.isEmpty()) {
    cards_[weightPlot]->show();
    weightPlot->setTitle(reaction);
    ThmPlotWidget::Series w;
    w.x = result.energy;
    w.y = result.weight;
    w.color = thmPlotColor(2);
    w.label = QFileInfo(result.weightFile).fileName();
    weightPlot->addSeries(w);
    bool positive = true;
    double lo = 1e300, hi = 0.0;
    for (double v : result.weight) {
      positive = positive && v > 0.0;
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
    weightPlot->setLogY(positive && hi > 20.0 * lo);  // a factor that varies by decades
  }

  // Spectator-momentum window.
  if (result.window) {
    cards_[windowPlot]->show();
    windowPlot->setTitle(reaction);
    if (!result.windowP.isEmpty()) {
      ThmPlotWidget::Series w;
      w.x = result.windowP;
      w.y = result.windowW;
      w.color = thmPlotColor(0);
      w.label = QString::fromUtf8("|φ|² p²");
      windowPlot->addSeries(w);
      ThmPlotWidget::Series n;
      n.x = result.nodeP;
      n.y = result.nodeW;
      n.color = thmPlotColor(1);
      n.symbols = true;
      n.label = tr("%1 nodes").arg(result.nodeP.size());
      windowPlot->addSeries(n);
    } else {
      windowPlot->setMessage(tr("one node, p_s = %1 MeV/c").arg(result.nodeP.value(0), 0, 'g', 6));
    }
  }
  layoutPanels();
}

bool ThmDiagnosticsPage::panelShown(ThmPlotWidget *plot) const {
  // The first two always; the optional ones are hidden and shown explicitly.
  return plot == vertexPlot || plot == hoesPlot || !cards_.value(plot)->isHidden();
}

int ThmDiagnosticsPage::columnsFor(int shown) const {
  // Three columns when there are more panels than two rows of two hold (or
  // three), and the page is wide enough for panels of ~280 px; else two.
  const bool wide = width() >= 3 * 280;
  return wide && (shown == 3 || shown > 4) ? 3 : 2;
}

void ThmDiagnosticsPage::layoutPanels() {
  // The panels that are shown fill a grid of equal cells in reading order,
  // so none leaves a hole; the plots of a row share their frame top and those
  // of a column their left edge (ThmPlotWidget::setAlignedWith).
  QList<ThmPlotWidget *> shown;
  for (ThmPlotWidget *p : {vertexPlot, hoesPlot, lineshapePlot, zetaPlot, weightPlot, windowPlot}) {
    grid_->removeWidget(cards_[p]);
    if (panelShown(p)) shown << p;
  }
  columns_ = columnsFor(shown.size());
  const int rows = (shown.size() + columns_ - 1) / columns_;
  for (int c = 0; c < 3; c++) grid_->setColumnStretch(c, c < columns_ ? 1 : 0);
  for (int r = 0; r < 3; r++) grid_->setRowStretch(r, r < rows ? 1 : 0);
  for (int i = 0; i < shown.size(); i++) grid_->addWidget(cards_[shown[i]], i / columns_, i % columns_);
  for (ThmPlotWidget *p : {lineshapePlot, zetaPlot, weightPlot, windowPlot})
    if (!panelShown(p)) grid_->addWidget(cards_[p], 3, 0);  // parked
  for (int i = 0; i < shown.size(); i++) {
    QList<ThmPlotWidget *> row, column;
    for (int j = 0; j < shown.size(); j++) {
      if (j / columns_ == i / columns_) row << shown[j];
      if (j % columns_ == i % columns_) column << shown[j];
    }
    shown[i]->setAlignedWith(row, column);
  }
}

void ThmDiagnosticsPage::resizeEvent(QResizeEvent *event) {
  QWidget::resizeEvent(event);
  int shown = 0;
  for (ThmPlotWidget *p : {vertexPlot, hoesPlot, lineshapePlot, zetaPlot, weightPlot, windowPlot})
    if (panelShown(p)) shown++;
  if (columnsFor(shown) != columns_) layoutPanels();
}

void ThmDiagnosticsPage::drawVertex() {
  vertexPlot->clear();
  const int g = vertexGroupCombo->currentIndex();
  if (g < 0 || g >= result_.vertex.size()) {
    vertexPlot->setMessage(result_.error.isEmpty() ? tr("no entrance channel") : tr("no result"));
    return;
  }
  const ThmDiagnosticsResult::VertexGroup &group = result_.vertex[g];
  bool window = false;
  for (int i = 0; i < group.curves.size(); i++) {
    ThmPlotWidget::Series s;
    s.x = result_.energy;
    s.y = group.curves[i].y;
    s.color = thmPlotColor(i);
    s.label = group.curves[i].label;
    vertexPlot->addSeries(s);
    if (!group.curves[i].yWindow.isEmpty()) {
      ThmPlotWidget::Series a;
      a.x = result_.energy;
      a.y = group.curves[i].yWindow;
      a.color = thmPlotColor(i);
      a.style = Qt::DashLine;  // one legend entry for all of them, below
      vertexPlot->addSeries(a);
      window = true;
    }
    for (double e : group.curves[i].nodes) {
      ThmPlotWidget::Marker m;
      m.x = e;
      m.color = thmPlotColor(i);
      m.label = QString::number(e, 'f', 3);
      vertexPlot->addMarker(m);
    }
  }
  if (window) vertexPlot->addLegendEntry(tr("p_s window"), Qt::DashLine);
}
