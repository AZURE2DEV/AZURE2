#include "ThmWorkspace.h"

#include <QHBoxLayout>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QTabWidget>
#include <QVBoxLayout>

#include "AZURESetup.h"
#include "ThmChannelsPage.h"
#include "ThmExperimentsPage.h"
#include "ThmModelPage.h"

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
  experimentsPage = new ThmExperimentsPage(settings.experimentLines, data, pairs);
  channelsPage = new ThmChannelsPage(pairs, setup->getLevelsTab()->getLevelsModel(),
                                     setup->getLevelsTab()->getChannelsModel(), data, test);
  pages = new QTabWidget;
  pages->addTab(modelPage, tr("Model"));
  pages->addTab(experimentsPage, tr("Experiments"));
  pages->addTab(channelsPage, tr("Channels"));
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
  resize(720, 680);
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
