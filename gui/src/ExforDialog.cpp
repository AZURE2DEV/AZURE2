#include "ExforDialog.h"
#include "PairsModel.h"
#include "ElementMap.h"

#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QHeaderView>
#include <QPlainTextEdit>
#include <QLabel>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QSplitter>
#include <QFormLayout>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QFile>
#include <QDir>
#include <QRegExp>
#include <QTextDocument>
#include <cmath>

ExforDialog::ExforDialog(PairsModel *pairsModel, QWidget *parent) :
  QDialog(parent),
  pairsModel_(pairsModel),
  exfor_(new ExforData(this)),
  currentDatasetId_(),
  haveData_(false) {
  setWindowTitle(tr("Retrieve Experimental Data from EXFOR"));
  resize(820, 640);

  // ---- Query group ---------------------------------------------------------
  QGroupBox *queryBox = new QGroupBox(tr("Reaction"));
  QGridLayout *q = new QGridLayout(queryBox);

  entranceCombo_ = new QComboBox;
  exitCombo_ = new QComboBox;
  targetEdit_ = new QLineEdit;
  reactionEdit_ = new QLineEdit;
  quantityCombo_ = new QComboBox;
  quantityCombo_->setEditable(true);
  quantityCombo_->addItem("SIG");  // cross section
  quantityCombo_->addItem("DA");   // differential (angular distribution)
  searchButton_ = new QPushButton(tr("Search EXFOR"));
  statusLabel_ = new QLabel;

  q->addWidget(new QLabel(tr("Entrance pair:")), 0, 0);
  q->addWidget(entranceCombo_, 0, 1);
  q->addWidget(new QLabel(tr("Exit pair:")), 0, 2);
  q->addWidget(exitCombo_, 0, 3);
  q->addWidget(new QLabel(tr("Target:")), 1, 0);
  q->addWidget(targetEdit_, 1, 1);
  q->addWidget(new QLabel(tr("Reaction:")), 1, 2);
  q->addWidget(reactionEdit_, 1, 3);
  q->addWidget(new QLabel(tr("Quantity:")), 2, 0);
  q->addWidget(quantityCombo_, 2, 1);
  q->addWidget(searchButton_, 2, 3);
  q->addWidget(statusLabel_, 3, 0, 1, 4);

  // ---- Results table -------------------------------------------------------
  resultsTable_ = new QTableWidget;
  resultsTable_->setColumnCount(6);
  resultsTable_->setHorizontalHeaderLabels(
      QStringList() << tr("Author / Year") << tr("Reaction") << tr("Points")
                    << tr("E range (MeV, lab)") << tr("DatasetID") << tr("Reference"));
  resultsTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
  resultsTable_->setSelectionMode(QAbstractItemView::SingleSelection);
  resultsTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  resultsTable_->horizontalHeader()->setStretchLastSection(true);
  resultsTable_->verticalHeader()->setVisible(false);

  datasetInfoLabel_ = new QLabel(tr("<i>Select a dataset to view details...</i>"));
  datasetInfoLabel_->setWordWrap(true);
  datasetInfoLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::LinksAccessibleByMouse);
  datasetInfoLabel_->setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
  datasetInfoLabel_->setContentsMargins(5, 5, 5, 5);

  downloadButton_ = new QPushButton(tr("Show / Preview Dataset"));
  downloadButton_->setEnabled(false);

  // ---- Preview -------------------------------------------------------------
  rawView_ = new QPlainTextEdit;
  rawView_->setReadOnly(true);
  rawView_->setLineWrapMode(QPlainTextEdit::NoWrap);
  azureView_ = new QPlainTextEdit;
  azureView_->setReadOnly(true);
  azureView_->setLineWrapMode(QPlainTextEdit::NoWrap);

  QGroupBox *rawBox = new QGroupBox(tr("Full dataset (EXFOR computational CSV)"));
  QVBoxLayout *rawLay = new QVBoxLayout(rawBox);
  rawLay->addWidget(rawView_);
  QGroupBox *azBox = new QGroupBox(
      tr("AZURE2 format  (E[MeV,lab]  angle[deg]  cross[b or b/sr]  error)"));
  QVBoxLayout *azLay = new QVBoxLayout(azBox);
  azLay->addWidget(azureView_);

  QSplitter *previewSplit = new QSplitter(Qt::Horizontal);
  previewSplit->addWidget(rawBox);
  previewSplit->addWidget(azBox);

  saveButton_ = new QPushButton(tr("Save to data folder..."));
  saveButton_->setEnabled(false);
  QPushButton *closeButton = new QPushButton(tr("Close"));

  QHBoxLayout *bottom = new QHBoxLayout;
  bottom->addWidget(saveButton_);
  bottom->addStretch();
  bottom->addWidget(closeButton);

  // ---- Assemble ------------------------------------------------------------
  QVBoxLayout *main = new QVBoxLayout(this);
  main->addWidget(queryBox);
  main->addWidget(resultsTable_, 1);
  main->addWidget(datasetInfoLabel_);
  main->addWidget(downloadButton_);
  main->addWidget(previewSplit, 2);
  main->addLayout(bottom);

  populatePairCombos();
  updateQueryFromPairs();

  connect(entranceCombo_, SIGNAL(currentIndexChanged(int)),
          this, SLOT(updateQueryFromPairs()));
  connect(exitCombo_, SIGNAL(currentIndexChanged(int)),
          this, SLOT(updateQueryFromPairs()));
  connect(searchButton_, SIGNAL(clicked()), this, SLOT(doSearch()));
  connect(downloadButton_, SIGNAL(clicked()), this, SLOT(downloadSelected()));
  connect(resultsTable_, SIGNAL(itemSelectionChanged()),
          this, SLOT(onDatasetSelectionChanged()));
  connect(saveButton_, SIGNAL(clicked()), this, SLOT(saveToDataFolder()));
  connect(closeButton, SIGNAL(clicked()), this, SLOT(reject()));

  connect(exfor_, SIGNAL(searchFinished(const QList<ExforDataset> &)),
          this, SLOT(onSearchFinished(const QList<ExforDataset> &)));
  connect(exfor_,
          SIGNAL(downloadFinished(const QString &, const QList<ExforPoint> &, bool)),
          this,
          SLOT(onDownloadFinished(const QString &, const QList<ExforPoint> &, bool)));
  connect(exfor_, SIGNAL(errorOccurred(const QString &)),
          this, SLOT(onError(const QString &)));
}

void ExforDialog::populatePairCombos() {
  entranceCombo_->clear();
  exitCombo_->clear();
  if (!pairsModel_) return;
  QList<PairsData> pairs = pairsModel_->getPairs();
  for (int i = 0; i < pairs.size(); ++i) {
    QString label = pairsModel_->getParticleLabel(pairs[i]);
    // Strip the HTML markup used in the table for a plain combo entry.
    QTextDocument doc;
    doc.setHtml(label);
    QString plain = doc.toPlainText().trimmed();
    entranceCombo_->addItem(plain, i);
    exitCombo_->addItem(plain, i);
  }
  if (pairs.size() > 1) exitCombo_->setCurrentIndex(1);
}

QString ExforDialog::nucleusCode(int z, double mass) {
  int a = qRound(mass);
  std::map<int, QString>::const_iterator it = elementMap.find(z);
  QString sym = (it != elementMap.end()) ? it->second : QString::number(z);
  return QString("%1-%2").arg(sym).arg(a);
}

QString ExforDialog::lightParticleCode(int z, double mass, int pairType) {
  if (pairType == 10) return "G";  // photon
  int a = qRound(mass);
  if (z == 0) return "N";  // neutron
  if (z == 1 && a == 1) return "P";
  if (z == 1 && a == 2) return "D";
  if (z == 1 && a == 3) return "T";
  if (z == 2 && a == 3) return "HE3";
  if (z == 2 && a == 4) return "A";
  return nucleusCode(z, mass);  // heavy ion
}

void ExforDialog::updateQueryFromPairs() {
  if (!pairsModel_) return;
  QList<PairsData> pairs = pairsModel_->getPairs();
  int ei = entranceCombo_->currentIndex();
  int xi = exitCombo_->currentIndex();
  if (ei < 0 || ei >= pairs.size() || xi < 0 || xi >= pairs.size()) return;

  const PairsData &en = pairs[ei];
  const PairsData &ex = pairs[xi];

  targetEdit_->setText(nucleusCode(en.heavyZ, en.heavyM));

  QString proj = lightParticleCode(en.lightZ, en.lightM, en.pairType);
  QString emission;
  bool elastic = (en.lightZ == ex.lightZ && qRound(en.lightM) == qRound(ex.lightM) &&
                  en.heavyZ == ex.heavyZ && qRound(en.heavyM) == qRound(ex.heavyM) &&
                  en.pairType == ex.pairType);
  if (elastic) {
    emission = "EL";
    quantityCombo_->setCurrentText("DA");
  } else {
    emission = lightParticleCode(ex.lightZ, ex.lightM, ex.pairType);
    quantityCombo_->setCurrentText(emission == "G" ? "SIG" : "SIG");
  }
  reactionEdit_->setText(QString("%1,%2").arg(proj).arg(emission));
}

void ExforDialog::doSearch() {
  resultsTable_->setRowCount(0);
  rawView_->clear();
  azureView_->clear();
  haveData_ = false;
  saveButton_->setEnabled(false);
  downloadButton_->setEnabled(false);

  statusLabel_->setText(tr("Searching EXFOR..."));
  searchButton_->setEnabled(false);
  exfor_->searchDatasets(targetEdit_->text().trimmed(),
                         reactionEdit_->text().trimmed(),
                         quantityCombo_->currentText().trimmed());
}

void ExforDialog::onSearchFinished(const QList<ExforDataset> &datasets) {
  searchButton_->setEnabled(true);
  resultsTable_->setRowCount(datasets.size());
  for (int i = 0; i < datasets.size(); ++i) {
    const ExforDataset &d = datasets[i];
    QString erange = QString("%1 - %2")
                         .arg(d.enMin / 1.0e6, 0, 'g', 4)
                         .arg(d.enMax / 1.0e6, 0, 'g', 4);
    resultsTable_->setItem(i, 0, new QTableWidgetItem(d.author));
    resultsTable_->setItem(i, 1, new QTableWidgetItem(d.reactionCode));
    resultsTable_->setItem(i, 2, new QTableWidgetItem(QString::number(d.npts)));
    resultsTable_->setItem(i, 3, new QTableWidgetItem(erange));
    resultsTable_->setItem(i, 4, new QTableWidgetItem(d.id));
    resultsTable_->setItem(i, 5, new QTableWidgetItem(d.reference));
  }
  resultsTable_->resizeColumnsToContents();
  statusLabel_->setText(tr("Found %1 dataset(s).").arg(datasets.size()));
}

void ExforDialog::onDatasetSelectionChanged() {
  QList<QTableWidgetItem *> selected = resultsTable_->selectedItems();
  downloadButton_->setEnabled(!selected.isEmpty());
  if (!selected.isEmpty()) {
    int row = selected.first()->row();
    if (row >= 0 && row < exfor_->datasets().size()) {
      const ExforDataset &d = exfor_->datasets().at(row);
      QString info = QString("<b>Dataset ID:</b> %1 &nbsp;&nbsp;|&nbsp;&nbsp; "
                             "<b>Reaction:</b> %2 &nbsp;&nbsp;|&nbsp;&nbsp; "
                             "<b>Points:</b> %3<br>"
                             "<b>Author / Year:</b> %4<br>"
                             "<b>Reference:</b> %5")
                         .arg(d.id)
                         .arg(d.reactionCode)
                         .arg(d.npts)
                         .arg(d.author)
                         .arg(d.reference);
      datasetInfoLabel_->setText(info);
    }
  } else {
    datasetInfoLabel_->setText(tr("<i>Select a dataset to view details...</i>"));
  }
}

void ExforDialog::downloadSelected() {
  int row = resultsTable_->currentRow();
  if (row < 0) return;
  QTableWidgetItem *idItem = resultsTable_->item(row, 4);
  if (!idItem) return;
  currentDatasetId_ = idItem->text();
  haveData_ = false;
  saveButton_->setEnabled(false);
  statusLabel_->setText(tr("Downloading dataset %1...").arg(currentDatasetId_));
  downloadButton_->setEnabled(false);

  // Pass the entrance-channel charges/masses so S-factor data (b*MeV, b*keV,
  // ...) can be converted back to cross section.
  double z1 = 0.0, z2 = 0.0, m1 = 0.0, m2 = 0.0;
  if (pairsModel_) {
    QList<PairsData> pairs = pairsModel_->getPairs();
    int ei = entranceCombo_->currentIndex();
    if (ei >= 0 && ei < pairs.size()) {
      z1 = pairs[ei].lightZ;
      m1 = pairs[ei].lightM;
      z2 = pairs[ei].heavyZ;
      m2 = pairs[ei].heavyM;
    }
  }
  exfor_->downloadDataset(currentDatasetId_, z1, z2, m1, m2);
}

void ExforDialog::onDownloadFinished(const QString &rawCsv,
                                     const QList<ExforPoint> &points,
                                     bool differential) {
  downloadButton_->setEnabled(true);
  rawView_->setPlainText(rawCsv);
  azureView_->setPlainText(ExforData::toAzureText(points));
  haveData_ = true;
  saveButton_->setEnabled(true);
  statusLabel_->setText(
      tr("Dataset %1: %2 point(s), %3.")
          .arg(currentDatasetId_)
          .arg(points.size())
          .arg(differential ? tr("differential (b/sr)")
                            : tr("angle-integrated (barn)")));
}

void ExforDialog::onError(const QString &message) {
  searchButton_->setEnabled(true);
  downloadButton_->setEnabled(!resultsTable_->selectedItems().isEmpty());
  statusLabel_->setText(tr("Error: %1").arg(message));
}

QString ExforDialog::suggestedFileName() const {
  QString rc = reactionEdit_->text();
  rc.replace(QRegExp("[^A-Za-z0-9]"), "_");
  return QString("%1_%2.dat").arg(currentDatasetId_).arg(rc);
}

void ExforDialog::saveToDataFolder() {
  if (!haveData_) return;

  QString startDir = QDir::current().filePath("data");
  if (!QDir(startDir).exists()) startDir = QDir::currentPath();
  QString defaultPath = QDir(startDir).filePath(suggestedFileName());

  QString fileName = QFileDialog::getSaveFileName(
      this, tr("Save EXFOR data as AZURE2 data file"), defaultPath,
      tr("Data files (*.dat *.txt);;All files (*)"));
  if (fileName.isEmpty()) return;

  QFile file(fileName);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QMessageBox::warning(this, tr("Save Failed"),
                         tr("Could not open file for writing:\n%1").arg(fileName));
    return;
  }
  QTextStream out(&file);
  out << ExforData::toAzureText(exfor_->lastPoints());
  file.close();

  statusLabel_->setText(tr("Saved %1 point(s) to %2")
                            .arg(exfor_->lastPoints().size())
                            .arg(QDir::toNativeSeparators(fileName)));
}
