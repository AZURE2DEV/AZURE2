#include "ThmModelPage.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QPushButton>
#include <QSet>
#include <QSpinBox>
#include <QTableWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <functional>
#include <sstream>

#include "ChooseFileButton.h"
#include "Config.h"

namespace {

// Shortest text that reads back as the same double.
QString numberText(double x) { return QString::number(x, 'g', QLocale::FloatingPointShortest); }

// The engine reads numbers with operator>> on an istringstream; so does this,
// so that "0.5 " or "5" are taken exactly as AZURE2 takes them.
bool readDouble(const QString &text, double &x) {
  std::istringstream s(text.toStdString());
  return !!(s >> x);
}
const char *kDocs =
    "docs/source/theory/thm_implementation.rst, section \"Options (&lt;thm&gt; block)\"; "
    "online: <a href=\"https://rdeboer1.github.io/AZURE2/\">rdeboer1.github.io/AZURE2</a> "
    "(Theory &gt; Trojan Horse (HOES) Observable).";

}  // namespace

// ---------------------------------------------------------------------------
// ThmModelPage

ThmModelPage::ThmModelPage(const ThmSettings &settings, const QString &projectDir, QWidget *parent) :
  QWidget(parent),
  projectDir_(projectDir) {
  entranceLCombo = new QComboBox;
  entranceLCombo->addItems(QStringList() << "incoherent" << "coherent");
  entranceLCombo->setToolTip(
      tr("entranceL: how the entrance partial waves l of one channel spin add. incoherent "
         "(default): the l cross terms vanish once the exit direction is integrated; coherent: "
         "one amplitude, as mrmpy and AZURE2 before September 2026."));
  vertexCombo = new QComboBox;
  vertexCombo->addItems(QStringList() << "constant" << "perlevel" << "onshell");
  vertexCombo->setToolTip(
      tr("vertex: boundary constant B in the transfer vertex (B-1) j_l(pa) - pa j_l'(pa). "
         "constant (default): shift function at the lowest level of the J group (La Cognata "
         "2010); perlevel: S_c(E_lambda) of each level under Brune, as mrmpy (alias real); "
         "onshell: L_c(E) = S_c + iP_c, the outgoing-wave log-derivative (Tribble 2014 eq. 2.76)."));
  kinematicsCombo = new QComboBox;
  kinematicsCombo->addItems(QStringList() << "lacognata" << "triple" << "kf3body" << "lambda32");
  kinematicsCombo->setToolTip(
      tr("kinematics: the kinematic factors the model carries, matching how the HOES data were "
         "divided out of the triple cross section. lacognata (default): exit k_f/mu_f (Tumino "
         "2021 eq. 50); triple: raw d3sigma/|phi|^2; kf3body: divided by the full three-body KF "
         "(Typel & Baur 2003); lambda32: divided by lambda3/lambda2 (Pizzone 2011)."));
  coulombIntegralCheck = new QCheckBox(tr("External Coulomb term in the vertex"));
  coulombIntegralCheck->setToolTip(
      tr("coulombIntegral: add 2 eta k Int_a^inf O_l(kr)/O_l(ka) j_l(pr) dr to the vertex "
         "(Tribble 2014 eq. 2.79; Typel & Baur eq. A.4). Mostly a normalization for small eta, "
         "a shape change for large eta."));
  spectatorEnergySpin = new QDoubleSpinBox;
  spectatorEnergySpin->setRange(0.0, 1000.0);
  spectatorEnergySpin->setDecimals(6);
  spectatorEnergySpin->setSuffix(" MeV");
  spectatorEnergySpin->setToolTip(
      tr("spectatorEnergy: mean spectator kinetic energy <p_sx^2>/2mu_sx added to E + B in the "
         "half-off-shell momentum (Typel & Baur 2003 eq. 11), for every THM entrance pair. "
         "0 = the spectator at rest (default)."));

  QGroupBox *globalBox = new QGroupBox(tr("THM observable"));
  QGridLayout *g = new QGridLayout;
  g->addWidget(new QLabel(tr("Entrance partial waves:")), 0, 0, Qt::AlignRight);
  g->addWidget(entranceLCombo, 0, 1);
  g->addWidget(new QLabel(tr("Vertex boundary:")), 1, 0, Qt::AlignRight);
  g->addWidget(vertexCombo, 1, 1);
  g->addWidget(new QLabel(tr("Kinematic factors:")), 2, 0, Qt::AlignRight);
  g->addWidget(kinematicsCombo, 2, 1);
  g->addWidget(coulombIntegralCheck, 3, 1);
  g->addWidget(new QLabel(tr("Spectator energy:")), 4, 0, Qt::AlignRight);
  g->addWidget(spectatorEnergySpin, 4, 1);
  g->setColumnStretch(1, 1);
  globalBox->setLayout(g);

  spectatorTable = new QTableWidget(0, 2);
  spectatorTable->setHorizontalHeaderLabels(QStringList() << tr("Pair") << tr("Energy (MeV)"));
  spectatorTable->horizontalHeader()->setStretchLastSection(true);
  spectatorTable->verticalHeader()->hide();
  spectatorTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  spectatorTable->setToolTip(
      tr("spectatorEnergy[<pair>]: the spectator energy for one THM entrance pair (its key in "
         "the Particle Pairs tab), overriding the global value."));
  QPushButton *addSpectator = new QPushButton(tr("Add"));
  QPushButton *removeSpectator = new QPushButton(tr("Remove"));
  connect(addSpectator, &QPushButton::clicked, this, [this]() { addSpectatorRow(1, 0.0); });
  connect(removeSpectator, SIGNAL(clicked()), this, SLOT(removeSelectedSpectatorRows()));
  QGroupBox *spectatorBox = new QGroupBox(tr("Spectator energy per pair"));
  QGridLayout *sl = new QGridLayout;
  sl->addWidget(spectatorTable, 0, 0, 1, 3);
  sl->addWidget(addSpectator, 1, 1);
  sl->addWidget(removeSpectator, 1, 2);
  sl->setColumnStretch(0, 1);
  spectatorBox->setLayout(sl);

  weightTable = new QTableWidget(0, 4);
  weightTable->setHorizontalHeaderLabels(QStringList() << tr("Segments") << tr("Segment")
                                                       << tr("Weight file (E_cm MeV, w)") << QString());
  weightTable->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
  weightTable->verticalHeader()->hide();
  weightTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  weightTable->setToolTip(
      tr("weight[<segment>]: w(E) multiplying the THM model of one segment before the resolution "
         "folding, e.g. Mukhamedzhanov's Coulomb-distortion factor R(E) (PRC 99 (2019) 064618). "
         "Two columns, E_cm of the THM entrance pair (MeV) and w > 0, strictly increasing E, "
         "'#' comments; a relative path is taken from the directory of the .azr. The segment "
         "number counts every line of the Data (or Test) segments table, active or not."));
  QPushButton *addWeight = new QPushButton(tr("Add"));
  QPushButton *removeWeight = new QPushButton(tr("Remove"));
  connect(addWeight, &QPushButton::clicked, this, [this]() { addWeightRow(false, 1, QString()); });
  connect(removeWeight, SIGNAL(clicked()), this, SLOT(removeSelectedWeightRows()));
  QGroupBox *weightBox = new QGroupBox(tr("Energy-dependent weight per segment"));
  QGridLayout *wl = new QGridLayout;
  wl->addWidget(weightTable, 0, 0, 1, 3);
  wl->addWidget(addWeight, 1, 1);
  wl->addWidget(removeWeight, 1, 2);
  wl->setColumnStretch(0, 1);
  weightBox->setLayout(wl);

  QLabel *docs = new QLabel(tr("These options apply to THM segments (isDiff >= 10) only; "
                               "defaults write no &lt;thm&gt; block. Documentation: ") +
                            QString(kDocs));
  docs->setWordWrap(true);
  docs->setOpenExternalLinks(true);
  docs->setTextFormat(Qt::RichText);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(globalBox);
  mainLayout->addWidget(spectatorBox);
  mainLayout->addWidget(weightBox);
  mainLayout->addWidget(docs);
  setLayout(mainLayout);

  setSettings(settings);
}

void ThmModelPage::clearRows() {
  spectatorTable->setRowCount(0);
  weightTable->setRowCount(0);
}

void ThmModelPage::addSpectatorRow(int pair, double energy) {
  int row = spectatorTable->rowCount();
  spectatorTable->insertRow(row);
  QSpinBox *pairSpin = new QSpinBox;
  pairSpin->setRange(1, 9999);
  pairSpin->setValue(pair);
  QLineEdit *energyEdit = new QLineEdit(numberText(energy));
  QDoubleValidator *v = new QDoubleValidator(0.0, 1.0e6, 12, energyEdit);
  v->setLocale(QLocale::c());
  energyEdit->setValidator(v);
  spectatorTable->setCellWidget(row, 0, pairSpin);
  spectatorTable->setCellWidget(row, 1, energyEdit);
}

void ThmModelPage::addWeightRow(bool test, int segment, const QString &file) {
  int row = weightTable->rowCount();
  weightTable->insertRow(row);
  QComboBox *block = new QComboBox;
  block->addItems(QStringList() << tr("Data") << tr("Test"));
  block->setCurrentIndex(test ? 1 : 0);
  QSpinBox *segSpin = new QSpinBox;
  segSpin->setRange(1, 9999);
  segSpin->setValue(segment);
  QLineEdit *fileEdit = new QLineEdit(file);
  ChooseFileButton *choose = new ChooseFileButton(tr("..."));
  choose->setLineEdit(fileEdit);
  connect(choose, SIGNAL(clicked(QLineEdit *)), this, SLOT(chooseWeightFile(QLineEdit *)));
  weightTable->setCellWidget(row, 0, block);
  weightTable->setCellWidget(row, 1, segSpin);
  weightTable->setCellWidget(row, 2, fileEdit);
  weightTable->setCellWidget(row, 3, choose);
}

void ThmModelPage::setSettings(const ThmSettings &s) {
  entranceLCombo->setCurrentText(s.entranceL);
  vertexCombo->setCurrentText(s.vertex);
  kinematicsCombo->setCurrentText(s.kinematics);
  coulombIntegralCheck->setChecked(s.coulombIntegral);
  spectatorEnergySpin->setValue(s.spectatorEnergy);
  clearRows();
  for (auto it = s.spectatorByPair.begin(); it != s.spectatorByPair.end(); ++it)
    addSpectatorRow(it.key(), it.value());
  for (auto it = s.weight.begin(); it != s.weight.end(); ++it) addWeightRow(false, it.key(), it.value());
  for (auto it = s.weightTest.begin(); it != s.weightTest.end(); ++it)
    addWeightRow(true, it.key(), it.value());
  // The spin box rounds to its decimals; keep the exact value that was read so
  // that an untouched dialog gives back exactly the settings it was given.
  spectatorEnergySpin->setProperty("exactValue", s.spectatorEnergy);
  experimentLines_ = s.experimentLines;
}

ThmSettings ThmModelPage::settings() const {
  ThmSettings s;
  s.entranceL = entranceLCombo->currentText();
  s.vertex = vertexCombo->currentText();
  s.kinematics = kinematicsCombo->currentText();
  s.coulombIntegral = coulombIntegralCheck->isChecked();
  double exact = spectatorEnergySpin->property("exactValue").toDouble();
  s.spectatorEnergy = spectatorEnergySpin->value() == spectatorEnergySpin->valueFromText(
                                                          spectatorEnergySpin->textFromValue(exact))
                          ? exact
                          : spectatorEnergySpin->value();
  for (int r = 0; r < spectatorTable->rowCount(); r++) {
    QSpinBox *pair = qobject_cast<QSpinBox *>(spectatorTable->cellWidget(r, 0));
    QLineEdit *energy = qobject_cast<QLineEdit *>(spectatorTable->cellWidget(r, 1));
    double x = -1.0;
    if (!readDouble(energy->text().trimmed(), x)) x = -1.0;  // refused by validate
    s.spectatorByPair[pair->value()] = x;
  }
  for (int r = 0; r < weightTable->rowCount(); r++) {
    QComboBox *block = qobject_cast<QComboBox *>(weightTable->cellWidget(r, 0));
    QSpinBox *segment = qobject_cast<QSpinBox *>(weightTable->cellWidget(r, 1));
    QLineEdit *file = qobject_cast<QLineEdit *>(weightTable->cellWidget(r, 2));
    (block->currentIndex() == 1 ? s.weightTest : s.weight)[segment->value()] = file->text().trimmed();
  }
  s.experimentLines = experimentLines_;
  return s;
}

QString ThmModelPage::check() const {
  // Two rows for the same key would silently keep only the last one.
  QSet<int> pairs, data, test;
  for (int r = 0; r < spectatorTable->rowCount(); r++) {
    int key = qobject_cast<QSpinBox *>(spectatorTable->cellWidget(r, 0))->value();
    if (pairs.contains(key)) return tr("Pair %1 is listed twice.").arg(key);
    pairs.insert(key);
  }
  for (int r = 0; r < weightTable->rowCount(); r++) {
    bool isTest = qobject_cast<QComboBox *>(weightTable->cellWidget(r, 0))->currentIndex() == 1;
    int key = qobject_cast<QSpinBox *>(weightTable->cellWidget(r, 1))->value();
    QSet<int> &seen = isTest ? test : data;
    if (seen.contains(key)) return tr("%1 segment %2 has two weight files.").arg(isTest ? "Test" : "Data").arg(key);
    seen.insert(key);
  }
  return settings().validate(projectDir_);
}

void ThmModelPage::chooseWeightFile(QLineEdit *lineEdit) {
  QString start = projectDir_.isEmpty() ? QDir::currentPath() : projectDir_;
  QString file = QFileDialog::getOpenFileName(this, tr("THM weight table"), start,
                                              tr("Weight tables (*.dat *.txt *.w);;All files (*)"));
  if (file.isEmpty()) return;
  // Inside the project directory: store it relative, so the project can move.
  if (!projectDir_.isEmpty()) {
    QString rel = QDir(projectDir_).relativeFilePath(file);
    if (!rel.startsWith("..")) file = rel;
  }
  lineEdit->setText(file);
}

void ThmModelPage::removeSelectedSpectatorRows() {
  QList<int> rows;
  for (const QModelIndex &i : spectatorTable->selectionModel()->selectedRows()) rows << i.row();
  if (rows.isEmpty() && spectatorTable->rowCount() > 0)
    rows << (spectatorTable->currentRow() >= 0 ? spectatorTable->currentRow() : spectatorTable->rowCount() - 1);
  std::sort(rows.begin(), rows.end(), std::greater<int>());
  for (int r : rows) spectatorTable->removeRow(r);
}

void ThmModelPage::removeSelectedWeightRows() {
  QList<int> rows;
  for (const QModelIndex &i : weightTable->selectionModel()->selectedRows()) rows << i.row();
  if (rows.isEmpty() && weightTable->rowCount() > 0)
    rows << (weightTable->currentRow() >= 0 ? weightTable->currentRow() : weightTable->rowCount() - 1);
  std::sort(rows.begin(), rows.end(), std::greater<int>());
  for (int r : rows) weightTable->removeRow(r);
}
