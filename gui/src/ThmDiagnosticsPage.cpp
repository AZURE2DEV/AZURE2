#include "ThmDiagnosticsPage.h"

#include <QComboBox>
#include <QFileInfo>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
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
  p->setMinimumHeight(250);
  return p;
}

const char *kEcm = "E (c.m. of the THM entrance pair, MeV)";

}  // namespace

ThmDiagnosticsPage::ThmDiagnosticsPage(std::function<QString(int, ThmDiagnosticsRequest &)> prepare,
                                       std::function<QList<Target>()> targets, QWidget *parent) :
  QWidget(parent),
  prepare_(prepare),
  targets_(targets) {
  segmentCombo = new QComboBox;
  segmentCombo->setToolTip(tr("The THM data segment to look at; its experiment, if any, in front."));
  computeButton = new QPushButton(tr("Compute"));
  computeButton->setToolTip(
      tr("Runs AZURE2 on the project as this workspace would leave it (a temporary copy; nothing is written into "
         "the project) at the current parameters, on a grid over the segment's data."));
  connect(computeButton, &QPushButton::clicked, this, [this]() { compute(); });
  busyBar = new QProgressBar;
  busyBar->setRange(0, 0);  // busy indicator
  busyBar->setMaximumWidth(120);
  busyBar->setTextVisible(false);
  busyBar->hide();
  statusLabel = new QLabel(tr("Read-only. Press Compute to evaluate the selected segment."));
  statusLabel->setWordWrap(true);

  QHBoxLayout *top = new QHBoxLayout;
  top->addWidget(new QLabel(tr("Segment:")));
  top->addWidget(segmentCombo, 1);
  top->addWidget(computeButton);
  top->addWidget(busyBar);

  vertexGroupCombo = new QComboBox;
  vertexGroupCombo->setToolTip(tr("The J^pi group whose entrance channels are shown."));
  connect(vertexGroupCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(drawVertex()));
  vertexPlot = makePlot(tr(kEcm), QString::fromUtf8("|M_l|\u00b2"));
  vertexPlot->setToolTip(
      tr("Entrance vertex |M_l(E)|^2 = |(B_c - 1) j_l(pa) - pa j_l'(pa) [+ C_l]|^2 at the quasi-free point "
         "(spectator at rest), with the boundary B_c of the <thm> vertex option; dashed lines: its nodes. A "
         "resonance near a node is suppressed in the HOES cross section."));
  vertexPanel = new QWidget;
  QVBoxLayout *vl = new QVBoxLayout;
  vl->setContentsMargins(0, 0, 0, 0);
  QHBoxLayout *vh = new QHBoxLayout;
  vh->addWidget(new QLabel(tr("Entrance J^pi:")));
  vh->addWidget(vertexGroupCombo);
  vh->addStretch(1);
  vl->addLayout(vh);
  vl->addWidget(vertexPlot, 1);
  vertexPanel->setLayout(vl);

  hoesPlot = makePlot(tr(kEcm), tr("cross section (b)"));
  hoesPlot->setLogY(true);
  hoesPlot->setToolTip(tr("The HOES cross section of the segment's channel (the options of the <thm> block, at the "
                          "quasi-free point, without resolution, weight or line shape; its arbitrary scale matched "
                          "to the other curve) and the on-shell angle-integrated cross section of the same channel."));
  lineshapePlot = makePlot(tr(kEcm), QString::fromUtf8("|N_C|\u00b2"));
  lineshapePlot->setToolTip(tr("|N_C|^2 = exp[2 zeta arctan(2 (E_lambda - E)/Gamma_lambda)] of the levels with a pole "
                               "near the data (lineshape=on); 1 at the pole."));
  zetaPlot = makePlot(tr(kEcm), QString::fromUtf8("ζ"));
  zetaPlot->setToolTip(tr("zeta(E) = eta_sB - eta_0 of the segment's exit pair (lineshape=on)."));
  weightPlot = makePlot(tr(kEcm), tr("w(E)"));
  weightPlot->setToolTip(tr("The weight table of the segment (weight[k]), as the engine interpolates it."));

  QWidget *panels = new QWidget;
  QGridLayout *grid = new QGridLayout;
  grid->addWidget(vertexPanel, 0, 0);
  grid->addWidget(hoesPlot, 0, 1);
  grid->addWidget(lineshapePlot, 1, 0);
  grid->addWidget(zetaPlot, 1, 1);
  grid->addWidget(weightPlot, 2, 0);
  grid->setColumnStretch(0, 1);
  grid->setColumnStretch(1, 1);
  panels->setLayout(grid);
  QScrollArea *scroll = new QScrollArea;
  scroll->setWidget(panels);
  scroll->setWidgetResizable(true);
  scroll->setFrameShape(QFrame::NoFrame);

  QVBoxLayout *main = new QVBoxLayout;
  main->addLayout(top);
  main->addWidget(statusLabel);
  main->addWidget(scroll, 1);
  setLayout(main);

  clearPlots(tr("not computed yet"));
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
  if (list_.isEmpty()) statusLabel->setText(tr("No active THM data segment to look at."));
  // Results of settings that have changed since are marked as such.
  if (!computedText_.isEmpty() && !busy()) {
    ThmDiagnosticsRequest now;
    if (prepare_(computedSegment_, now).isEmpty() && now.projectText != computedText_)
      statusLabel->setText(tr("The project or this workspace changed since the last Compute; press Compute again."));
  }
}

QString ThmDiagnosticsPage::reactionOf(int segment) const {
  for (const Target &t : list_)
    if (t.segment == segment) return t.reaction;
  return QString();
}

void ThmDiagnosticsPage::clearPlots(const QString &message) {
  for (ThmPlotWidget *p : {vertexPlot, hoesPlot, lineshapePlot, zetaPlot, weightPlot}) {
    p->clear();
    p->setMessage(message);
  }
  lineshapePlot->hide();
  zetaPlot->hide();
  weightPlot->hide();
  vertexGroupCombo->clear();
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
  statusLabel->setText(
      tr("Segment %1%2: %3 points from %4 to %5 MeV; vertex %6, B + T_s = %7 MeV.")
          .arg(result.segment)
          .arg(result.experiment.isEmpty() ? QString() : tr(" (experiment %1)").arg(result.experiment))
          .arg(result.energy.size())
          .arg(result.eLo, 0, 'g', 4)
          .arg(result.eHi, 0, 'g', 4)
          .arg(result.vertexMode)
          .arg(result.binding, 0, 'g', 6) +
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
    lineshapePlot->show();
    zetaPlot->show();
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
    weightPlot->show();
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
}

void ThmDiagnosticsPage::drawVertex() {
  vertexPlot->clear();
  const int g = vertexGroupCombo->currentIndex();
  if (g < 0 || g >= result_.vertex.size()) {
    vertexPlot->setMessage(result_.error.isEmpty() ? tr("no entrance channel") : tr("no result"));
    return;
  }
  const ThmDiagnosticsResult::VertexGroup &group = result_.vertex[g];
  for (int i = 0; i < group.curves.size(); i++) {
    ThmPlotWidget::Series s;
    s.x = result_.energy;
    s.y = group.curves[i].y;
    s.color = thmPlotColor(i);
    s.label = group.curves[i].label;
    vertexPlot->addSeries(s);
    for (double e : group.curves[i].nodes) {
      ThmPlotWidget::Marker m;
      m.x = e;
      m.color = thmPlotColor(i);
      m.label = QString::number(e, 'f', 3);
      vertexPlot->addMarker(m);
    }
  }
}
