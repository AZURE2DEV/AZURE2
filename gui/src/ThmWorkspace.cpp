#include "ThmWorkspace.h"

#include <QHBoxLayout>
#include <QTextStream>
#include <QRegularExpression>
#include <QTextDocument>
#include <QLabel>
#include <QMessageBox>
#include <QPointer>
#include <QPushButton>
#include <QTabWidget>
#include <QVBoxLayout>

#include "AZURESetup.h"
#include "ChannelsModel.h"
#include "Config.h"
#include "ThmChannelsPage.h"
#ifdef AZURE2_THM_DIAGNOSTICS
#include "ThmDiagnosticsPage.h"
#endif
#include "ThmExperimentsPage.h"
#include "ThmModelPage.h"

#ifdef AZURE2_THM_DIAGNOSTICS
namespace {
// Active THM data segments, with the experiment they are in on the Experiments page.
QList<ThmDiagnosticsPage::Target> diagnosticsTargets(AZURESetup *setup, ThmExperimentsPage *experimentsPage) {
  QList<ThmDiagnosticsPage::Target> out;
  PairsModel *pairs = setup->getPairsTab()->getPairsModel();
  const QList<SegmentsDataData> lines = setup->getSegmentsTab()->getSegmentsDataModel()->getLines();
  const QList<PairsData> pairList = pairs->getPairs();
  for (int k = 1; k <= lines.size(); k++) {
    const SegmentsDataData &s = lines.at(k - 1);
    if (!s.isActive || !s.isTHM) continue;
    ThmDiagnosticsPage::Target t;
    t.segment = k;
    if (s.entrancePairIndex >= 1 && s.entrancePairIndex <= pairList.size() && s.exitPairIndex >= 1 &&
        s.exitPairIndex <= pairList.size())
      t.reaction = pairs->getReactionLabel(pairList.at(s.entrancePairIndex - 1), pairList.at(s.exitPairIndex - 1))
                       .remove(QRegularExpression("\\s*\\[0(\\.0*)? MeV\\]"));  // ground states: no "[0.000 MeV]"
    QTextDocument d;
    d.setHtml(t.reaction);
    QString experiment;
    for (const ThmExperimentRecord &r : experimentsPage->records())
      if (r.segments.contains(k)) experiment = r.name;
    t.text = (experiment.isEmpty() ? QString() : ThmWorkspace::tr("experiment %1: ").arg(experiment)) +
             ThmWorkspace::tr("segment %1, %2 (%3)").arg(k).arg(d.toPlainText().simplified(), s.dataFile);
    out << t;
  }
  return out;
}
}  // namespace
#endif

ThmWorkspace::ThmWorkspace(AZURESetup *setup, const ThmSettings &settings, QWidget *parent) :
  QDialog(parent),
  setup_(setup) {
  PairsModel *pairs = setup->getPairsTab()->getPairsModel();
  SegmentsDataModel *data = setup->getSegmentsTab()->getSegmentsDataModel();
  SegmentsTestModel *test = setup->getSegmentsTab()->getSegmentsTestModel();
  hasThm_ = !ThmChannelsPage::thmEntrancePairs(data, test).isEmpty();

  noThmLabel = new QLabel(
      tr("This project has no THM segment. The THM workspace edits the Trojan Horse (HOES) observable of THM "
         "segments: tick \"THM\" on a segment in the Segments tab (isDiff >= 10), then open it again."));
  noThmLabel->setWordWrap(true);
  noThmLabel->setVisible(!hasThm_);

  modelPage = new ThmModelPage(settings, setup->projectDirectory());
  experimentsPage = new ThmExperimentsPage(settings.experimentLines, data, pairs, setup->projectDirectory(),
                                           !!(setup->GetConfig().paramMask & Config::USE_BRUNE_FORMALISM));
  // A ps window is refused with a spectator energy for its pair: the Model page's values.
  // (Guarded: the pages go one by one when the dialog is destroyed.)
  QPointer<ThmModelPage> model(modelPage);
  experimentsPage->setSpectatorEnergy(
      [model](int pairKey) { return model ? model->settings().spectatorEnergyOf(pairKey) : 0.0; });
  experimentsPage->refreshDerived();
  channelsPage = new ThmChannelsPage(pairs, setup->getLevelsTab()->getLevelsModel(),
                                     setup->getLevelsTab()->getChannelsModel(), data, test);
  pages = new QTabWidget;
  pages->addTab(modelPage, tr("Model"));
  pages->addTab(experimentsPage, tr("Experiments"));
  pages->addTab(channelsPage, tr("Channels"));
#ifdef AZURE2_THM_DIAGNOSTICS
  diagnosticsPage = new ThmDiagnosticsPage(
      [this](int segment, ThmDiagnosticsRequest &request) {
        QString why;
        if (!projectSnapshot(request.projectText, &why)) return why;
        request.projectDir = setup_->projectDirectory();
        request.paramMask = setup_->GetConfig().paramMask;
        request.segment = segment;
        return QString();
      },
      [this]() { return diagnosticsTargets(setup_, experimentsPage); });
  pages->addTab(diagnosticsPage, tr("Diagnostics"));
  pages->setTabToolTip(3, tr("Read-only plots computed by AZURE2 on request: entrance vertex, HOES and on-shell "
                             "cross sections, line shape, weight, spectator-momentum window"));
  connect(pages, &QTabWidget::currentChanged, this, [this](int) {
    if (pages->currentWidget() == diagnosticsPage) diagnosticsPage->refreshTargets();
  });
#endif
  connect(pages, &QTabWidget::currentChanged, experimentsPage, [this](int) {
    if (pages->currentWidget() == experimentsPage) experimentsPage->refreshDerived();
  });
  pages->setTabToolTip(0, tr("Options of the THM observable: the <thm> block"));
  pages->setTabToolTip(1, tr("Segments sharing one profiled norm and a background: experiment[...] lines"));
  pages->setTabToolTip(2, tr("Binding energy and width input flag of the THM entrance channels"));
  pages->setEnabled(hasThm_);

  QPushButton *cancelButton = new QPushButton(tr("Cancel"));
  acceptButton = new QPushButton(tr("Accept"));
  acceptButton->setDefault(true);
  acceptButton->setEnabled(hasThm_);
  connect(acceptButton, SIGNAL(clicked()), this, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(acceptButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(noThmLabel);
  mainLayout->addWidget(pages, 1);
  mainLayout->addLayout(buttonBox);
  setLayout(mainLayout);
  setWindowTitle(tr("THM Workspace"));
  resize(760, 720);
}


bool ThmWorkspace::projectSnapshot(QString &text, QString *error) {
  if (setup_->GetConfig().configfile.empty()) {
    if (error) *error = tr("Save the project first: the engine runs on its file.");
    return false;
  }
  // The <thm> block as Accept would install it (AZURESetup::setThmSettings).
  ThmSettings wanted = modelPage->settings();
  wanted.experimentLines = experimentsPage->experimentLines();
  ThmSettings current;
  const bool readable = setup_->thmSettings(current);
  bool present = setup_->hasThmOptionsBlock();
  QStringList block = setup_->thmOptionsLines();
  if (!(readable && current == wanted)) {
    present = !wanted.isDefault();
    block = present ? wanted.compose(setup_->hasThmOptionsBlock() ? setup_->thmOptionsLines() : QStringList())
                    : QStringList();
  }
  // The Channels page's values, installed for the time of writing.
  PairsModel *pairs = setup_->getPairsTab()->getPairsModel();
  ChannelsModel *channels = setup_->getLevelsTab()->getChannelsModel();
  const QList<PairsData> pairsBefore = pairs->getPairs();
  const QList<ChannelsData> channelsBefore = channels->getChannels();
  channelsPage->apply();
  text.clear();
  QTextStream out(&text);
  const bool written = setup_->writeProject(out, setup_->projectDirectory());
  out.flush();
  const QList<PairsData> pairsAfter = pairs->getPairs();
  for (int i = 0; i < pairsBefore.size() && i < pairsAfter.size(); i++)
    if (pairsAfter[i].bindingEnergy != pairsBefore[i].bindingEnergy)
      pairs->setData(pairs->index(i, 15), pairsBefore[i].bindingEnergy, Qt::EditRole);
  const QList<ChannelsData> channelsAfter = channels->getChannels();
  for (int i = 0; i < channelsBefore.size() && i < channelsAfter.size(); i++)
    if (channelsAfter[i].gammaIsRWA != channelsBefore[i].gammaIsRWA)
      channels->setData(channels->index(i, 7), channelsBefore[i].gammaIsRWA, Qt::EditRole);
  if (!written) {
    if (error) *error = tr("The project could not be written.");
    return false;
  }
  // Replace the block that was written (right after </targetInt>) by the wanted one.
  const QString anchor = "</targetInt>\n";
  int at = text.indexOf(anchor);
  if (at < 0) {
    if (error) *error = tr("The project has no <targetInt> block.");
    return false;
  }
  at += anchor.size();
  if (setup_->hasThmOptionsBlock()) {
    const int end = text.indexOf("</thm>\n", at);
    if (end >= 0) text.remove(at, end + 7 - at);
  }
  if (present) {
    QString b;
    QTextStream bs(&b);
    AZURESetup::writeThmBlock(bs, block);
    bs.flush();
    text.insert(at, b);
  }
  return true;
}

QString ThmWorkspace::validate() {
  if (!hasThm_) return QString();
  QString why = modelPage->check();
  if (!why.isEmpty()) {
    pages->setCurrentWidget(modelPage);
    return why;
  }
  why = experimentsPage->check();
  if (!why.isEmpty()) {
    pages->setCurrentWidget(experimentsPage);
    return why;
  }
  why = channelsPage->check();
  if (!why.isEmpty()) pages->setCurrentWidget(channelsPage);
  return why;
}

void ThmWorkspace::apply() {
  if (!hasThm_) return;
  ThmSettings s = modelPage->settings();
  s.experimentLines = experimentsPage->experimentLines();
  setup_->setThmSettings(s);
  channelsPage->apply();
}

void ThmWorkspace::accept() {
  const QString why = validate();
  if (!why.isEmpty()) {
    QMessageBox::warning(this, tr("THM Workspace"), why);
    return;
  }
  apply();
  QDialog::accept();
}
