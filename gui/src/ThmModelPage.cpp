#include "ThmModelPage.h"
#include "ThmNumberText.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
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

// ---------------------------------------------------------------------------
// ThmModelPage

ThmModelPage::ThmModelPage(const ThmSettings &settings, const QString &projectDir, QWidget *parent) :
  QWidget(parent),
  projectDir_(projectDir) {
  entranceLCombo = new QComboBox;
  entranceLCombo->addItems(QStringList() << "incoherent" << "coherent");
  entranceLCombo->setToolTip(
      tr("entranceL: incoherent (default): the l cross terms vanish once the exit direction is integrated; "
         "coherent: one amplitude, as mrmpy and AZURE2 before September 2026; not with an exit-angle window "
         "(Experiments page)."));
  vertexCombo = new QComboBox;
  vertexCombo->addItems(QStringList() << "constant" << "perlevel" << "onshell");
  vertexCombo->setToolTip(
      tr("vertex: boundary B in (B-1) j_l(pa) - pa j_l'(pa). constant (default): shift function at the lowest "
         "level of the J group; perlevel: S_c(E_lambda) of each level (Brune); onshell: L_c(E) = S_c + iP_c."));
  kinematicsCombo = new QComboBox;
  kinematicsCombo->addItems(QStringList() << "lacognata" << "triple" << "kf3body" << "lambda32");
  kinematicsCombo->setToolTip(
      tr("kinematics: as the HOES data were extracted. lacognata (default): exit k_f/mu_f; triple: raw "
         "d3sigma/|phi|^2; kf3body: divided by the three-body KF; lambda32: divided by lambda3/lambda2."));
  coulombIntegralCheck = new QCheckBox(tr("External Coulomb term in the vertex"));
  coulombIntegralCheck->setToolTip(
      tr("coulombIntegral: add 2 eta k Int_a^inf O_l(kr)/O_l(ka) j_l(pr) dr to the vertex (Tribble 2014 eq. "
         "2.79)."));
  spectatorEnergySpin = new QDoubleSpinBox;
  spectatorEnergySpin->setRange(0.0, 1000.0);
  spectatorEnergySpin->setDecimals(6);
  spectatorEnergySpin->setSuffix(" MeV");
  spectatorEnergySpin->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  spectatorEnergySpin->setToolTip(
      tr("spectatorEnergy: mean spectator kinetic energy added to E + B in the half-off-shell momentum, for "
         "every THM entrance pair; 0 = at rest (default)."));

  QGroupBox *globalBox = new QGroupBox(tr("THM observable"));
  QFormLayout *g = new QFormLayout;
  g->setLabelAlignment(Qt::AlignRight | Qt::AlignVCenter);
  g->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  g->setHorizontalSpacing(8);
  g->setVerticalSpacing(6);
  g->addRow(tr("Entrance partial waves:"), entranceLCombo);
  g->addRow(tr("Vertex boundary:"), vertexCombo);
  g->addRow(tr("Kinematic factors:"), kinematicsCombo);
  g->addRow(tr("Spectator energy:"), spectatorEnergySpin);
  g->addRow(QString(), coulombIntegralCheck);
  globalBox->setLayout(g);

  auto tableBox = [](const QString &title, QTableWidget *table, QPushButton *add, QPushButton *remove) {
    QGroupBox *box = new QGroupBox(title);
    QGridLayout *l = new QGridLayout;
    l->setVerticalSpacing(6);
    l->addWidget(table, 0, 0, 1, 3);
    l->addWidget(add, 1, 1);
    l->addWidget(remove, 1, 2);
    l->setColumnStretch(0, 1);
    box->setLayout(l);
    return box;
  };

  spectatorTable = new QTableWidget(0, 2);
  spectatorTable->setHorizontalHeaderLabels(QStringList() << tr("Pair") << tr("Energy (MeV)"));
  spectatorTable->horizontalHeader()->setStretchLastSection(true);
  spectatorTable->verticalHeader()->hide();
  spectatorTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  spectatorTable->setToolTip(tr("spectatorEnergy[<pair>]: overrides the global value for one THM entrance pair "
                                "(its key in the Particle Pairs tab)."));
  QPushButton *addSpectator = new QPushButton(tr("Add"));
  QPushButton *removeSpectator = new QPushButton(tr("Remove"));
  connect(addSpectator, &QPushButton::clicked, this, [this]() { addSpectatorRow(1, 0.0); });
  connect(removeSpectator, SIGNAL(clicked()), this, SLOT(removeSelectedSpectatorRows()));
  QGroupBox *spectatorBox = tableBox(tr("Spectator energy per pair"), spectatorTable, addSpectator, removeSpectator);

  weightTable = new QTableWidget(0, 4);
  weightTable->setHorizontalHeaderLabels(QStringList() << tr("Segments") << tr("Segment")
                                                       << tr("Weight file (E_cm MeV, w)") << QString());
  weightTable->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
  weightTable->horizontalHeader()->setSectionResizeMode(3, QHeaderView::ResizeToContents);
  weightTable->verticalHeader()->hide();
  weightTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  weightTable->setToolTip(
      tr("weight[<segment>]: w(E) multiplying the THM model of one segment before the resolution folding, e.g. a "
         "Coulomb-distortion factor. Two columns, E_cm (MeV) and w > 0; relative to the project directory. The "
         "segment number counts every line of the Data (or Test) segments table."));
  QPushButton *addWeight = new QPushButton(tr("Add"));
  QPushButton *removeWeight = new QPushButton(tr("Remove"));
  connect(addWeight, &QPushButton::clicked, this, [this]() { addWeightRow(false, 1, QString()); });
  connect(removeWeight, SIGNAL(clicked()), this, SLOT(removeSelectedWeightRows()));
  QGroupBox *weightBox = tableBox(tr("Energy-dependent weight per segment"), weightTable, addWeight, removeWeight);

  QHBoxLayout *top = new QHBoxLayout;
  top->setSpacing(12);
  top->addWidget(globalBox, 3);
  top->addWidget(spectatorBox, 2);
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->setSpacing(12);
  mainLayout->addLayout(top);
  mainLayout->addWidget(weightBox, 1);
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
  pairSpin->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  pairSpin->setFrame(false);
  pairSpin->setRange(1, 9999);
  pairSpin->setValue(pair);
  QLineEdit *energyEdit = new QLineEdit(ThmText::number(energy));
  energyEdit->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  energyEdit->setFrame(false);
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
  segSpin->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  segSpin->setFrame(false);
  segSpin->setRange(1, 9999);
  segSpin->setValue(segment);
  QLineEdit *fileEdit = new QLineEdit(file);
  fileEdit->setFrame(false);
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
    if (!ThmText::readDouble(energy->text().trimmed(), x)) x = -1.0;  // refused by validate
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
