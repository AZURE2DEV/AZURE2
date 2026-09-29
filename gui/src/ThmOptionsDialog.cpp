#include "ThmOptionsDialog.h"

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
#include <QMessageBox>
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
bool readInt(const QString &text, int &x, bool strict) {
  std::istringstream s(text.toStdString());
  if (!(s >> x)) return false;
  std::string rest;
  return !strict || !(s >> rest);
}

const char *kDocs =
    "docs/source/theory/thm_implementation.rst, section \"Options (&lt;thm&gt; block)\"; "
    "online: <a href=\"https://rdeboer1.github.io/AZURE2/\">rdeboer1.github.io/AZURE2</a> "
    "(Theory &gt; Trojan Horse (HOES) Observable).";

}  // namespace

// ---------------------------------------------------------------------------
// ThmSettings

bool ThmSettings::operator==(const ThmSettings &o) const {
  return entranceL == o.entranceL && vertex == o.vertex && kinematics == o.kinematics &&
         coulombIntegral == o.coulombIntegral && spectatorEnergy == o.spectatorEnergy &&
         spectatorByPair == o.spectatorByPair && weight == o.weight && weightTest == o.weightTest &&
         experimentLines == o.experimentLines;
}

QList<QPair<QString, QString>> ThmSettings::keyValues() const {
  const ThmSettings d;
  QList<QPair<QString, QString>> kv;
  if (entranceL != d.entranceL) kv << qMakePair(QString("entranceL"), entranceL);
  if (vertex != d.vertex) kv << qMakePair(QString("vertex"), vertex);
  if (kinematics != d.kinematics) kv << qMakePair(QString("kinematics"), kinematics);
  if (coulombIntegral != d.coulombIntegral) kv << qMakePair(QString("coulombIntegral"), QString("1"));
  if (spectatorEnergy != d.spectatorEnergy)
    kv << qMakePair(QString("spectatorEnergy"), numberText(spectatorEnergy));
  for (auto it = spectatorByPair.begin(); it != spectatorByPair.end(); ++it)
    kv << qMakePair(QString("spectatorEnergy[%1]").arg(it.key()), numberText(it.value()));
  for (auto it = weight.begin(); it != weight.end(); ++it)
    kv << qMakePair(QString("weight[%1]").arg(it.key()), it.value());
  for (auto it = weightTest.begin(); it != weightTest.end(); ++it)
    kv << qMakePair(QString("weightTest[%1]").arg(it.key()), it.value());
  return kv;
}

bool ThmSettings::parseLine(const QString &rawLine, QString &key, QString &value, ThmSettings &s) {
  key.clear();
  value.clear();
  QString line = rawLine;
  int hash = line.indexOf('#');
  if (hash >= 0) line.truncate(hash);
  line = line.trimmed();
  if (line.isEmpty()) return true;
  if (line.startsWith("experiment[")) {
    // A THM experiment record: not an option of this dialog.  Kept verbatim
    // (an empty key, as a comment); Config::ReadThmBlock validates it.
    s.experimentLines << rawLine;
    return true;
  }
  int eq = line.indexOf('=');
  if (eq < 0) return false;
  QString k = line.left(eq);
  while (!k.isEmpty() && (k.endsWith(' ') || k.endsWith('\t'))) k.chop(1);
  QString v = line.mid(eq + 1);
  while (!v.isEmpty() && (v.startsWith(' ') || v.startsWith('\t'))) v.remove(0, 1);

  if (k == "vertex") {
    if (v == "real") v = "perlevel";
    if (v != "onshell" && v != "constant" && v != "perlevel") return false;
    s.vertex = v;
  } else if (k == "kinematics") {
    if (v != "lacognata" && v != "triple" && v != "kf3body" && v != "lambda32") return false;
    s.kinematics = v;
  } else if (k == "entranceL") {
    if (v != "coherent" && v != "incoherent") return false;
    s.entranceL = v;
  } else if (k == "coulombIntegral") {
    if (v == "1" || v == "true" || v == "on")
      s.coulombIntegral = true;
    else if (v == "0" || v == "false" || v == "off")
      s.coulombIntegral = false;
    else
      return false;
    v = s.coulombIntegral ? "1" : "0";
  } else if (k.startsWith("spectatorEnergy")) {
    double x;
    if (!readDouble(v, x) || !(x >= 0.0)) return false;
    v = numberText(x);
    if (k == "spectatorEnergy") {
      s.spectatorEnergy = x;
    } else if (k.size() > 17 && k[15] == '[' && k.endsWith(']')) {
      int pair;
      if (!readInt(k.mid(16, k.size() - 17), pair, false)) return false;
      s.spectatorByPair[pair] = x;
      k = QString("spectatorEnergy[%1]").arg(pair);
    } else {
      return false;
    }
  } else if (k.startsWith("weight")) {
    bool test = k.startsWith("weightTest[");
    int open = test ? 10 : 6;
    int segment;
    if (!(k.size() > open + 2 && k[open] == '[' && k.endsWith(']') && !v.isEmpty())) return false;
    if (!readInt(k.mid(open + 1, k.size() - open - 2), segment, true) || segment < 1) return false;
    (test ? s.weightTest : s.weight)[segment] = v;
    k = QString(test ? "weightTest[%1]" : "weight[%1]").arg(segment);
  } else {
    return false;
  }
  key = k;
  value = v;
  return true;
}

bool ThmSettings::parse(const QStringList &lines, ThmSettings &out, QString *error) {
  ThmSettings s;
  for (const QString &line : lines) {
    QString key, value;
    if (!parseLine(line, key, value, s)) {
      if (error) *error = QString("<thm> line not understood: '%1'").arg(line.trimmed());
      return false;
    }
  }
  out = s;
  return true;
}

QString ThmSettings::validate(const QString &projectDir) const {
  const QMap<int, QString> *maps[2] = {&weight, &weightTest};
  const char *names[2] = {"weight", "weightTest"};
  for (int m = 0; m < 2; m++) {
    for (auto it = maps[m]->begin(); it != maps[m]->end(); ++it) {
      QString name = it.value().trimmed();
      if (name.isEmpty()) return QString("%1[%2]: no file given.").arg(names[m]).arg(it.key());
      if (name.contains('#'))
        return QString("%1[%2]: the path contains '#', which starts a comment in the .azr.")
            .arg(names[m])
            .arg(it.key());
      QString path = QDir::isAbsolutePath(name) || projectDir.isEmpty()
                         ? name
                         : QDir(projectDir).filePath(name);
      ThmWeightTable table;
      std::string why = table.Read(path.toStdString());
      if (!why.empty())
        return QString("%1[%2]: %3").arg(names[m]).arg(it.key()).arg(QString::fromStdString(why));
    }
  }
  for (auto it = spectatorByPair.begin(); it != spectatorByPair.end(); ++it)
    if (!(it.value() >= 0.0)) return QString("spectatorEnergy[%1] must be >= 0.").arg(it.key());
  if (!(spectatorEnergy >= 0.0)) return QString("spectatorEnergy must be >= 0.");
  return QString();
}

QStringList ThmSettings::compose(const QStringList &oldLines) const {
  const QList<QPair<QString, QString>> kv = keyValues();
  QMap<QString, QString> wanted;
  for (const auto &p : kv) wanted[p.first] = p.second;
  QSet<QString> written;
  QStringList out;
  for (const QString &line : oldLines) {
    ThmSettings scratch;
    QString key, value;
    if (!parseLine(line, key, value, scratch) || key.isEmpty()) {
      out << line;  // comment, blank (or a line the engine would refuse: left alone)
      continue;
    }
    if (!wanted.contains(key) || written.contains(key)) continue;  // back to default, or repeated
    written.insert(key);
    if (value == wanted[key]) {
      out << line;
      continue;
    }
    // Rewrite the value, keeping the indentation and any inline comment.
    int hash = line.indexOf('#');
    QString code = hash >= 0 ? line.left(hash) : line;
    QString comment = hash >= 0 ? line.mid(hash) : QString();
    int start = 0;
    while (start < code.size() && code[start].isSpace()) start++;
    int end = code.size();
    while (end > start && code[end - 1].isSpace()) end--;
    out << code.left(start) + key + "=" + wanted[key] + code.mid(end) + comment;
  }
  for (const auto &p : kv)
    if (!written.contains(p.first)) out << p.first + "=" + p.second;
  return out;
}

// ---------------------------------------------------------------------------
// ThmOptionsDialog

ThmOptionsDialog::ThmOptionsDialog(const ThmSettings &settings, const QString &projectDir,
                                   QWidget *parent) :
  QDialog(parent),
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

  QPushButton *cancelButton = new QPushButton(tr("Cancel"));
  QPushButton *okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);
  connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(globalBox);
  mainLayout->addWidget(spectatorBox);
  mainLayout->addWidget(weightBox);
  mainLayout->addWidget(docs);
  mainLayout->addLayout(buttonBox);
  setLayout(mainLayout);
  setWindowTitle(tr("THM Options"));
  resize(640, 560);

  setSettings(settings);
}

void ThmOptionsDialog::clearRows() {
  spectatorTable->setRowCount(0);
  weightTable->setRowCount(0);
}

void ThmOptionsDialog::addSpectatorRow(int pair, double energy) {
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

void ThmOptionsDialog::addWeightRow(bool test, int segment, const QString &file) {
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

void ThmOptionsDialog::setSettings(const ThmSettings &s) {
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

ThmSettings ThmOptionsDialog::settings() const {
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

void ThmOptionsDialog::accept() {
  // Two rows for the same key would silently keep only the last one.
  QSet<int> pairs, data, test;
  for (int r = 0; r < spectatorTable->rowCount(); r++) {
    int key = qobject_cast<QSpinBox *>(spectatorTable->cellWidget(r, 0))->value();
    if (pairs.contains(key)) {
      QMessageBox::warning(this, tr("THM Options"), tr("Pair %1 is listed twice.").arg(key));
      return;
    }
    pairs.insert(key);
  }
  for (int r = 0; r < weightTable->rowCount(); r++) {
    bool isTest = qobject_cast<QComboBox *>(weightTable->cellWidget(r, 0))->currentIndex() == 1;
    int key = qobject_cast<QSpinBox *>(weightTable->cellWidget(r, 1))->value();
    QSet<int> &seen = isTest ? test : data;
    if (seen.contains(key)) {
      QMessageBox::warning(this, tr("THM Options"),
                           tr("%1 segment %2 has two weight files.").arg(isTest ? "Test" : "Data").arg(key));
      return;
    }
    seen.insert(key);
  }
  QString error = settings().validate(projectDir_);
  if (!error.isEmpty()) {
    QMessageBox::warning(this, tr("THM Options"), error);
    return;
  }
  QDialog::accept();
}

void ThmOptionsDialog::chooseWeightFile(QLineEdit *lineEdit) {
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

void ThmOptionsDialog::removeSelectedSpectatorRows() {
  QList<int> rows;
  for (const QModelIndex &i : spectatorTable->selectionModel()->selectedRows()) rows << i.row();
  if (rows.isEmpty() && spectatorTable->rowCount() > 0)
    rows << (spectatorTable->currentRow() >= 0 ? spectatorTable->currentRow() : spectatorTable->rowCount() - 1);
  std::sort(rows.begin(), rows.end(), std::greater<int>());
  for (int r : rows) spectatorTable->removeRow(r);
}

void ThmOptionsDialog::removeSelectedWeightRows() {
  QList<int> rows;
  for (const QModelIndex &i : weightTable->selectionModel()->selectedRows()) rows << i.row();
  if (rows.isEmpty() && weightTable->rowCount() > 0)
    rows << (weightTable->currentRow() >= 0 ? weightTable->currentRow() : weightTable->rowCount() - 1);
  std::sort(rows.begin(), rows.end(), std::greater<int>());
  for (int r : rows) weightTable->removeRow(r);
}
