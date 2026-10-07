#include "ThmExperimentsPage.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QScrollArea>
#include <QSet>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QStandardItemModel>
#include <QTimer>
#include <QStyle>
#include <QTableWidget>
#include <QTextDocument>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "PairsModel.h"
#include "SegmentsDataModel.h"
#include "DataLine.h"
#include "ThmExperiment.h"
#include "ThmLineshape.h"
#include "ThmNumberSpin.h"
#include "Config.h"
#include "Constants.h"
#include "ThmDistortion.h"
#include "ThmOptical.h"

namespace {

// MeV/u, the value EData::SetupThmExperiments converts the binding with.
const double kAmu = 931.49410242;

// A whole token as a number, as the engine reads Ebeam.
bool readWholeDouble(const QString &text, double &x) {
  std::istringstream s(text.toStdString());
  std::string rest;
  return !!(s >> x) && !(s >> rest) && std::isfinite(x);
}

void splitWindow(const QString &text, QString &lo, QString &hi);

QString plain(const QString &html) {
  QTextDocument d;
  d.setHtml(html);
  return d.toPlainText().simplified();
}

}  // namespace

QStringList ThmExperimentsPage::nuclideNames() {
  // The engine's table, looked up rather than copied (ThmNuclide::Find).
  QStringList names;
  for (int Z = 0; Z <= 118; Z++)
    for (int A = std::max(Z, 1); A <= 300; A++)
      if (const ThmNuclide *n = ThmNuclide::Find(Z, A)) names << QString::fromStdString(n->name);
  return names;
}

ThmExperimentsPage::ThmExperimentsPage(const QStringList &experimentLines, SegmentsDataModel *segments,
                                       PairsModel *pairs, const QString &projectDir, bool brune, QWidget *parent) :
  QWidget(parent),
  segments_(segments),
  pairs_(pairs),
  projectDir_(projectDir),
  brune_(brune),
  oldLines_(experimentLines) {
  oldRecords_ = ThmExperimentRecord::read(experimentLines);
  records_ = oldRecords_;
  pairBinding_ = [this](int pairKey) {
    const QList<PairsData> pairs = pairs_->getPairs();
    return pairKey >= 1 && pairKey <= pairs.size() ? pairs.at(pairKey - 1).bindingEnergy : 0.0;
  };

  experimentTable = new QTableWidget(0, 4);
  experimentTable->setHorizontalHeaderLabels(QStringList() << tr("Name") << tr("Segments") << tr("Background")
                                                           << tr("Three-body reaction"));
  // A compact list: name and segments; the rest is in the row's tooltip.
  experimentTable->setColumnHidden(2, true);
  experimentTable->setColumnHidden(3, true);
  experimentTable->horizontalHeader()->setStretchLastSection(false);
  experimentTable->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
  experimentTable->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
  experimentTable->verticalHeader()->hide();
  experimentTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  experimentTable->setSelectionMode(QAbstractItemView::SingleSelection);
  experimentTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
  experimentTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  experimentTable->setToolTip(
      tr("experiment[<name>]: THM data segments measured together share one profiled norm and, optionally, a "
         "background. Segments in no experiment keep their own norm."));
  addButton = new QPushButton(tr("Add"));
  removeButton = new QPushButton(tr("Remove"));
  connect(addButton, &QPushButton::clicked, this, [this]() { addExperiment(); });
  connect(removeButton, &QPushButton::clicked, this, [this]() { removeCurrent(); });
  connect(experimentTable->selectionModel(), SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)),
          this, SLOT(tableSelectionChanged()));

  QGroupBox *listBox = new QGroupBox(tr("Experiments"));
  listBox->setMaximumWidth(fontMetrics().horizontalAdvance("M") * 12);
  QGridLayout *ll = new QGridLayout;
  ll->addWidget(experimentTable, 0, 0, 1, 2);
  ll->addWidget(addButton, 1, 0);
  ll->addWidget(removeButton, 1, 1);
  listBox->setLayout(ll);

  nameEdit = new QLineEdit;
  nameEdit->setToolTip(tr("Letters, digits and _ - . + only."));
  connect(nameEdit, SIGNAL(textEdited(const QString &)), this, SLOT(nameEdited(const QString &)));
  segmentList = new QListWidget;
  segmentList->setToolTip(
      tr("segments=: THM data segments (isDiff >= 10) with a free norm that are in no other experiment. Numbers "
         "count every line of the Data segments table, active or not."));
  // About three rows: the sections below keep their room.
  segmentList->setMinimumHeight(2 * segmentList->fontMetrics().lineSpacing() + 10);
  segmentList->setMaximumHeight(4 * segmentList->fontMetrics().lineSpacing() + 10);
  segmentList->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
  connect(segmentList, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(segmentItemChanged(QListWidgetItem *)));
  backgroundCombo = new QComboBox;
  backgroundCombo->addItems(QStringList() << "none" << "const" << "linear" << "quadratic");
  backgroundCombo->setToolTip(
      tr("background=: b0 + b1 E + b2 E^2 (E: c.m. energy of the THM entrance pair) added to the folded model, "
         "profiled with the norm."));
  connect(backgroundCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(backgroundChanged(int)));
  // The exit angle (theta=): angle-integrated, or dsigma/dOmega averaged over a window.
  const QString thetaTip =
      tr("theta=: with a window the observable is dsigma/dOmega averaged over the window of the c.m. angle of exit "
         "particle 1 relative to the x-A direction (entrance particle 1 relative to 2); min = max: one angle. If "
         "the paper quotes the other exit particle's angle, give the supplementary window, 180 - max to 180 - min. "
         "Not with entranceL=coherent (Model page).");
  thetaCombo = new QComboBox;
  thetaCombo->addItem(tr("all (angle-integrated)"), "all");
  thetaCombo->addItem(tr("window"), "window");
  thetaCombo->setToolTip(thetaTip);
  connect(thetaCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(thetaEdited()));
  thetaMinEdit = new ThmNumberSpin(QString::fromUtf8("°"), 0.0, 180.0, 5.0);
  thetaMaxEdit = new ThmNumberSpin(QString::fromUtf8("°"), 0.0, 180.0, 5.0);
  thetaMinEdit->setToolTip(tr("theta_min, c.m. (0-180)."));
  thetaMaxEdit->setToolTip(tr("theta_max, c.m. (theta_min-180)."));
  thetaMinEdit->setWrittenText("0");
  thetaMaxEdit->setWrittenText("180");
  connect(thetaMinEdit, SIGNAL(valueChanged(double)), this, SLOT(thetaEdited()));
  connect(thetaMaxEdit, SIGNAL(valueChanged(double)), this, SLOT(thetaEdited()));

  const QStringList nuclides = nuclideNames();
  auto nuclideCombo = [&](const QString &tip) {
    QComboBox *c = new QComboBox;
    c->setEditable(true);
    c->addItems(nuclides);
    c->setCurrentIndex(-1);
    c->setInsertPolicy(QComboBox::NoInsert);
    c->setToolTip(tip + tr(" A nuclide of the built-in table (AME2020) or Z,A,mass (nuclear mass in u)."));
    connect(c, SIGNAL(editTextChanged(const QString &)), this, SLOT(kinematicsEdited()));
    return c;
  };
  beamCombo = nuclideCombo(tr("beam=: the projectile."));
  targetCombo = nuclideCombo(tr("target=: the target."));
  spectatorCombo = nuclideCombo(tr("spectator=: the spectator s of the Trojan horse a = x + s."));
  beamEnergyEdit = new ThmNumberSpin(" MeV", 0.0, 1.0e5, 1.0);
  beamEnergyEdit->setSpecialValueText(QString::fromUtf8("—"));  // 0: not given
  beamEnergyEdit->setToolTip(tr("Ebeam=: the lab beam energy (> 0)."));
  connect(beamEnergyEdit, SIGNAL(valueChanged(double)), this, SLOT(kinematicsEdited()));
  lineshapeCheck = new QCheckBox(tr("Coulomb (lineshape=on)"));
  lineshapeCheck->setToolTip(
      tr("lineshape=on: the Coulomb field of the charged spectator skews and shifts each resonance (factor N_C per "
         "level; Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233). Needs the three-body reaction and the Brune "
         "parameterization."));
  connect(lineshapeCheck, SIGNAL(toggled(bool)), this, SLOT(lineshapeToggled(bool)));

  // Spectator-momentum window (ps=, psNodes=).
  const QString psPhysics =
      tr("At fixed E the spectator direction fixes |p_s|, and the three-body phase space is d cos(theta_cm) = "
         "|p_s| d|p_s|/(beta k_sF k_aA): the vertex is averaged over the accepted directions -- |p_s| inside "
         "[p_min, p_max] and reachable at E, and the Directions below -- with the weight |phi(p_s)|^2 "
         "d cos(theta_cm), p_xA^2/2mu_xA = E + B + p_s^2/2mu_sx. Not together with a non-zero spectator energy "
         "(Model page) for the same entrance pair.");
  psKindCombo = new QComboBox;
  psKindCombo->addItem(tr("point (quasi-free)"), "delta");
  psKindCombo->addItem(QString::fromUtf8("Hulthén"), "hulthen");
  psKindCombo->addItem(tr("Gaussian"), "gauss");
  psKindCombo->addItem(tr("table"), "table");
  psKindCombo->setToolTip(tr("ps=: the momentum distribution |phi(p_s)|^2 of the spectator. Point: p_s = 0, "
                             "nothing written. ") +
                          psPhysics);
  connect(psKindCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(psEdited()));
  auto psSpin = [&](const QString &suffix, double hi, double step, const QString &tip) {
    ThmNumberSpin *e = new ThmNumberSpin(suffix, 0.0, hi, step);
    e->setToolTip(tip);
    connect(e, SIGNAL(valueChanged(double)), this, SLOT(psEdited()));
    return e;
  };
  const QString perFm = QString::fromUtf8(" fm⁻¹");
  psMinEdit = psSpin(" MeV/c", 1.0e4, 1.0, tr("p_min of the accepted window (>= 0)."));
  psMaxEdit = psSpin(" MeV/c", 1.0e4, 1.0, tr("p_max of the accepted window (>= p_min; = p_min: one point)."));
  psCustomCheck = new QCheckBox(tr("custom a, b"));
  psCustomCheck->setToolTip(QString::fromUtf8("Hulthén phi(p) ~ 1/(a^2+q^2) - 1/(b^2+q^2). Unticked: the "
                                              "deuteron, a = 0.2317, b = 1.202 fm^-1 (Tribble 2014 eq. 4.4)."));
  connect(psCustomCheck, SIGNAL(toggled(bool)), this, SLOT(psEdited()));
  psAEdit = psSpin(perFm, 100.0, 0.01, tr("Hulthen a (> 0)."));
  psBEdit = psSpin(perFm, 100.0, 0.01, tr("Hulthen b (> a)."));
  psAEdit->setWrittenText("0.2317");
  psBEdit->setWrittenText("1.202");
  psFwhmEdit = psSpin(" MeV/c", 1.0e4, 1.0, tr("FWHM of |phi(p_s)|^2 = exp(-4 ln2 p^2/FWHM^2) (> 0)."));
  psTableEdit = new QLineEdit;
  psTableEdit->setPlaceholderText(tr("file"));
  psTableEdit->setToolTip(tr("ps=table:<file>: two columns, p_s (MeV/c, increasing) and the momentum distribution "
                             "|phi(p_s)|^2 (any scale; e.g. the measured yield over the kinematic factor); the "
                             "window is its range. Relative to the project directory."));
  connect(psTableEdit, SIGNAL(textEdited(const QString &)), this, SLOT(psEdited()));
  psTableButton = new QPushButton("...");
  psTableButton->setToolTip(tr("Choose the table; a file inside the project directory is stored relative to it."));
  connect(psTableButton, SIGNAL(clicked()), this, SLOT(chooseTable()));
  psNodesSpin = new QSpinBox;
  psNodesSpin->setRange(1, 64);
  psNodesSpin->setValue(16);
  psNodesSpin->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  psNodesSpin->setToolTip(tr("psNodes=: Gauss-Legendre nodes in cos(theta_cm) on the directions inside the "
                             "window (1-64, default 16); with a direction window its own nodes are used."));
  connect(psNodesSpin, SIGNAL(valueChanged(int)), this, SLOT(psNodesChanged(int)));

  // Spectator directions (spectatorAngles=, spectatorAngleNodes=): with a
  // computed distortion, R(E) or the DW vertex averaged over the accepted
  // directions instead of the one spectatorAngle direction.
  directionCombo = new QComboBox;
  directionCombo->addItem(tr("one (Distortion: angle)"), "one");
  directionCombo->addItem(tr("lab window"), "lab");
  directionCombo->addItem(tr("c.m. window"), "cm");
  directionCombo->addItem(tr("lab table"), "labtable");
  directionCombo->addItem(tr("c.m. table"), "cmtable");
  directionCombo->setToolTip(
      tr("spectatorAngles=: the accepted spectator directions (polar angle to the beam, lab or c.m.) over which "
         "the vertex (with a distribution above) and R(E) or the DW vertex are averaged, weight d cos(theta_cm) x "
         "acceptance x |phi|^2. At fixed E the direction fixes |p_s|, so the momentum window above cuts the same "
         "directions. A table gives the acceptance per angle. All: every direction inside the momentum window; "
         "one: the single direction of the Distortion section's angle."));
  connect(directionCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(directionEdited()));
  auto dirSpin = [&](const QString &tip) {
    ThmNumberSpin *e = new ThmNumberSpin(QString::fromUtf8("°"), 0.0, 180.0, 1.0);
    e->setToolTip(tip);
    connect(e, SIGNAL(valueChanged(double)), this, SLOT(directionEdited()));
    return e;
  };
  directionMinEdit = dirSpin(tr("The smallest accepted polar angle of the spectator to the beam (0-180)."));
  directionMaxEdit = dirSpin(tr("The largest accepted polar angle (>= the smallest; equal: one direction)."));
  directionTableEdit = new QLineEdit;
  directionTableEdit->setPlaceholderText(tr("file"));
  directionTableEdit->setToolTip(tr("spectatorAngles=[cm:]table:<file>: two columns, the angle (deg, increasing) "
                                    "and the acceptance (>= 0); the window is its range. Relative to the project "
                                    "directory."));
  connect(directionTableEdit, SIGNAL(textEdited(const QString &)), this, SLOT(directionEdited()));
  directionTableButton = new QPushButton("...");
  directionTableButton->setToolTip(psTableButton->toolTip());
  connect(directionTableButton, SIGNAL(clicked()), this, SLOT(chooseDirectionTable()));
  directionNodesSpin = new QSpinBox;
  directionNodesSpin->setRange(1, 64);
  directionNodesSpin->setValue(8);
  directionNodesSpin->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  directionNodesSpin->setToolTip(tr("spectatorAngleNodes=: Gauss-Legendre nodes in cos(theta_cm) per interval "
                                    "(1-64, default 8; a lab window can map to two c.m. intervals)."));
  connect(directionNodesSpin, SIGNAL(valueChanged(int)), this, SLOT(directionNodesChanged(int)));

  // Distortion factor (distortion= and its keys).
  distortionCombo = new QComboBox;
  distortionCombo->addItem(tr("None"), "none");
  distortionCombo->addItem(tr("Coulomb"), "coulomb");
  distortionCombo->addItem(tr("Optical"), "optical");
  distortionCombo->addItem(tr("Table"), "table");
  distortionCombo->setToolTip(
      tr("distortion=: R(E) multiplies the model of every segment before the folding (zero-range DWBA; "
         "Mukhamedzhanov & Pang PRC 99 (2019) 064618). Coulomb: point-Coulomb waves in a + A and s + F; optical: "
         "per channel; table: w(E) from a file. Coulomb and optical need the three-body reaction."));
  connect(distortionCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(distortionKindChanged()));
  angleKindCombo = new QComboBox;
  angleKindCombo->addItem(tr("quasi-free"), "qf");
  angleKindCombo->addItem(tr("lab"), "lab");
  angleKindCombo->addItem(tr("c.m."), "cm");
  angleKindCombo->setToolTip(tr("spectatorAngle=: the direction of the spectator. Quasi-free: k_sF along k_aA "
                                "(default); lab: an angle to the beam converted at every E; c.m.: fixed."));
  connect(angleKindCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(distortionEdited()));
  angleEdit = new ThmNumberSpin(QString::fromUtf8("°"), 0.0, 180.0, 1.0);
  angleEdit->setToolTip(tr("The spectator angle to the beam (0-180)."));
  connect(angleEdit, SIGNAL(valueChanged(double)), this, SLOT(distortionEdited()));
  distortionRefEdit = new ThmNumberSpin(" MeV", 0.0, 1.0e3, 0.1);
  distortionRefEdit->setSpecialValueText(tr("auto"));
  distortionRefEdit->setToolTip(tr("distortionRef=: E_ref, where R = 1 (c.m. of x + A). Auto: the middle of the "
                                   "data. Only the scale, which the profiled norm absorbs."));
  connect(distortionRefEdit, SIGNAL(valueChanged(double)), this, SLOT(distortionEdited()));
  ratioCombo = new QComboBox;
  ratioCombo->addItem(tr("DWBA/PWBA"), "dwpw");
  ratioCombo->addItem(tr("DWBA"), "dw");
  ratioCombo->setToolTip(tr("distortionRatio=: rho = |M/M_PW|^2 (dwpw, default: the correction to data divided by "
                            "the momentum distribution) or |M|^2 (dw, the papers' ratio); R = rho(E)/rho(E_ref)."));
  connect(ratioCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(distortionEdited()));
  boundCombo = new QComboBox;
  boundCombo->addItem(tr("Whittaker"), "whittaker");
  boundCombo->addItem(tr("Yukawa"), "yukawa");
  boundCombo->setToolTip(tr("boundState=: the s-x bound state tail, W(2 kappa r)/r (default) or exp(-kappa r)/r."));
  connect(boundCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(distortionEdited()));
  rminEdit = new ThmNumberSpin(" fm", 0.0, 50.0, 0.5);
  rminEdit->setSpecialValueText(QString::fromUtf8("—"));
  rminEdit->setToolTip(tr("boundState=...:rmin: the bound state is zero below r_min (0-50 fm)."));
  connect(rminEdit, SIGNAL(valueChanged(double)), this, SLOT(distortionEdited()));
  const char *channelKeys[2] = {"opticalAA", "opticalSF"};
  for (int c = 0; c < 2; c++) {
    opticalCombo[c] = new QComboBox;
    opticalCombo[c]->addItem(tr("plane"), "plane");
    opticalCombo[c]->addItem(tr("Coulomb"), "coulomb");
    opticalCombo[c]->addItem(tr("global"), "global");
    opticalCombo[c]->addItem(QString::fromUtf8("Woods–Saxon"), "ws");
    opticalCombo[c]->setToolTip(tr("%1=: plane (no distortion), point Coulomb (default), a global optical "
                                   "potential evaluated at the channel's energy, or a Woods-Saxon optical potential "
                                   "plus Coulomb.")
                                    .arg(channelKeys[c]));
    connect(opticalCombo[c], SIGNAL(currentIndexChanged(int)), this, SLOT(opticalKindChanged()));
    globalCombo[c] = new QComboBox;
    for (const ThmGlobalOptical &g : ThmGlobalOpticals()) {
      globalCombo[c]->addItem(g.name, g.name);
      globalCombo[c]->setItemData(globalCombo[c]->count() - 1,
                                  QString::fromUtf8("%1 — %2; %3 on A = %4–%5, E_lab = %6–%7 MeV")
                                      .arg(g.name, g.reference, g.projectiles)
                                      .arg(g.aMin)
                                      .arg(g.aMax)
                                      .arg(g.eMin)
                                      .arg(g.eMax),
                                  Qt::ToolTipRole);
    }
    globalCombo[c]->setToolTip(tr("The global optical potential (spin-orbit dropped), evaluated at the projectile's "
                                  "lab energy: fixed for a + A, following E_sF for s + F."));
    connect(globalCombo[c], SIGNAL(currentIndexChanged(int)), this, SLOT(globalNameChanged()));
    opticalButton[c] = new QPushButton(QString::fromUtf8("Edit…"));
    opticalButton[c]->setAutoDefault(false);
    connect(opticalButton[c], &QPushButton::clicked, this, [this, c]() { editOptical(c); });
  }
  distortionTableEdit = new QLineEdit;
  distortionTableEdit->setPlaceholderText(tr("file"));
  distortionTableEdit->setToolTip(tr("distortion=table:<file>: two columns, E (c.m. of x + A, MeV, increasing) and "
                                     "w > 0, the weight[k] format; every point inside it. Relative to the project "
                                     "directory."));
  connect(distortionTableEdit, SIGNAL(textEdited(const QString &)), this, SLOT(distortionEdited()));
  distortionTableButton = new QPushButton("...");
  distortionTableButton->setToolTip(psTableButton->toolTip());
  connect(distortionTableButton, SIGNAL(clicked()), this, SLOT(chooseDistortionTable()));

  // Derived values, as AZURE2 prints them (the whole text: derivedText, in the tooltips).
  auto value = [this]() {
    QLabel *l = new QLabel(QString::fromUtf8("—"));
    l->setTextInteractionFlags(Qt::TextSelectableByMouse);
    return l;
  };
  bindingValue = value();
  qfValue = value();
  zetaValue = value();
  meanTsValue = value();
  distortionValue = value();
  bindingWarningIcon = new QLabel;
  bindingWarningIcon->setPixmap(
      style()->standardIcon(QStyle::SP_MessageBoxWarning)
          .pixmap(style()->pixelMetric(QStyle::PM_SmallIconSize), style()->pixelMetric(QStyle::PM_SmallIconSize)));
  bindingWarningLabel = new QLabel;
  bindingWarningLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
  messageIcon = new QLabel;
  const int icon = style()->pixelMetric(QStyle::PM_SmallIconSize);
  messageIcon->setPixmap(style()->standardIcon(QStyle::SP_MessageBoxWarning).pixmap(icon, icon));
  messageIcon->setAlignment(Qt::AlignTop);
  messageLabel = new QLabel;
  messageLabel->setWordWrap(true);
  messageLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);

  // Sections: four-column forms (label, field, label, field) whose label
  // columns have one width in every section, so that the fields line up.
  QList<QGridLayout *> forms;
  QList<QLabel *> leftLabels, rightLabels;
  auto label = [&](const QString &text, bool left, const QString &tip = QString()) {
    QLabel *l = new QLabel(text);
    l->setToolTip(tip);
    (left ? leftLabels : rightLabels) << l;
    return l;
  };
  auto form = [&]() {
    QGridLayout *g = new QGridLayout;
    g->setHorizontalSpacing(8);
    g->setVerticalSpacing(6);
    g->setColumnStretch(1, 1);
    g->setColumnStretch(3, 1);
    forms << g;
    return g;
  };
  const Qt::Alignment right = Qt::AlignRight | Qt::AlignVCenter;

  QGroupBox *experimentBox = new QGroupBox(tr("Experiment"));
  QGridLayout *xl = form();
  xl->addWidget(label(tr("Name:"), true), 0, 0, right);
  xl->addWidget(nameEdit, 0, 1);
  xl->addWidget(label(tr("Background:"), false), 0, 2, right);
  xl->addWidget(backgroundCombo, 0, 3);
  xl->addWidget(label(tr("Segments:"), true), 1, 0, Qt::AlignRight | Qt::AlignTop);
  xl->addWidget(segmentList, 1, 1, 1, 3);
  xl->addWidget(label(tr("Exit angle:"), true, thetaTip), 2, 0, right);
  xl->addWidget(thetaCombo, 2, 1);
  QLabel *thetaLabel = label(QString::fromUtf8("θ<sub>cm</sub>:"), false, thetaTip);
  QHBoxLayout *thl = new QHBoxLayout;
  thl->setContentsMargins(0, 0, 0, 0);
  thl->addWidget(thetaMinEdit, 1);
  thl->addWidget(new QLabel(QString::fromUtf8("–")));
  thl->addWidget(thetaMaxEdit, 1);
  QWidget *thetaBox = new QWidget;
  thetaBox->setLayout(thl);
  xl->addWidget(thetaLabel, 2, 2, right);
  xl->addWidget(thetaBox, 2, 3);
  thetaWindowRow_ = {thetaLabel, thetaBox};
  experimentBox->setLayout(xl);

  kinematicsBox = new QGroupBox(tr("Three-body reaction"));
  kinematicsBox->setCheckable(true);
  kinematicsBox->setChecked(false);
  kinematicsBox->setToolTip(
      tr("beam, target, spectator and Ebeam: all four or none. One of beam and target is a nucleus of the "
         "segments' entrance pair, the other the second nucleus plus the spectator."));
  connect(kinematicsBox, SIGNAL(toggled(bool)), this, SLOT(kinematicsEdited()));
  QGridLayout *kl = form();
  kl->addWidget(label(tr("Beam:"), true), 0, 0, right);
  kl->addWidget(beamCombo, 0, 1);
  kl->addWidget(label(tr("Target:"), false), 0, 2, right);
  kl->addWidget(targetCombo, 0, 3);
  kl->addWidget(label(tr("Spectator:"), true), 1, 0, right);
  kl->addWidget(spectatorCombo, 1, 1);
  kl->addWidget(label(tr("Beam energy:"), false, tr("Lab frame")), 1, 2, right);
  kl->addWidget(beamEnergyEdit, 1, 3);
  kl->addWidget(label("B(x+s):", true, tr("Binding energy of the Trojan horse a = x + s")), 2, 0, right);
  kl->addWidget(bindingValue, 2, 1);
  kl->addWidget(label(QString::fromUtf8("E<sub>qf</sub>:"), false, tr("Quasi-free energy E(x+A) - B")), 2, 2, right);
  kl->addWidget(qfValue, 2, 3);
  // The line shape needs the reaction: in its section.
  kl->addWidget(label(tr("Line shape:"), true), 3, 0, right);
  kl->addWidget(lineshapeCheck, 3, 1);
  kl->addWidget(label(QString::fromUtf8("\u03b6:"), false, tr("zeta at the lowest and highest data point")), 3, 2,
                right);
  kl->addWidget(zetaValue, 3, 3);
  // B from the masses against the pair's B (field 32): one line when they differ.
  QHBoxLayout *bw = new QHBoxLayout;
  bw->setContentsMargins(0, 0, 0, 0);
  bw->addWidget(bindingWarningIcon);
  bw->addWidget(bindingWarningLabel, 1);
  kl->addLayout(bw, 4, 1, 1, 3);
  kinematicsBox->setLayout(kl);

  psBox = new QGroupBox(tr("Spectator acceptance"));
  psBox->setToolTip(psPhysics);
  // Only the fields of the chosen distribution are shown (showPsRows).
  QGridLayout *pl = form();
  pl->addWidget(label(tr("Distribution:"), true), 0, 0, right);
  pl->addWidget(psKindCombo, 0, 1);
  pl->addWidget(psCustomCheck, 0, 3);
  QLabel *minLabel = label(QString::fromUtf8("p<sub>min</sub>:"), true);
  QLabel *maxLabel = label(QString::fromUtf8("p<sub>max</sub>:"), false);
  pl->addWidget(minLabel, 1, 0, right);
  pl->addWidget(psMinEdit, 1, 1);
  pl->addWidget(maxLabel, 1, 2, right);
  pl->addWidget(psMaxEdit, 1, 3);
  QLabel *tableLabel = label(tr("Table:"), true);
  QHBoxLayout *tl = new QHBoxLayout;
  tl->setContentsMargins(0, 0, 0, 0);
  tl->addWidget(psTableEdit, 1);
  tl->addWidget(psTableButton);
  pl->addWidget(tableLabel, 1, 0, right);
  pl->addLayout(tl, 1, 1, 1, 3);
  QLabel *aLabel = label("a:", true), *bLabel = label("b:", false);
  pl->addWidget(aLabel, 2, 0, right);
  pl->addWidget(psAEdit, 2, 1);
  pl->addWidget(bLabel, 2, 2, right);
  pl->addWidget(psBEdit, 2, 3);
  QLabel *fwhmLabel = label(tr("FWHM:"), true);
  pl->addWidget(fwhmLabel, 2, 0, right);
  pl->addWidget(psFwhmEdit, 2, 1);
  QLabel *meanLabel = label(QString::fromUtf8("\u27e8T<sub>s</sub>\u27e9:"), true, tr("Mean spectator energy"));
  QLabel *nodesLabel = label(tr("Nodes:"), false);
  pl->addWidget(meanLabel, 3, 0, right);
  pl->addWidget(meanTsValue, 3, 1);
  pl->addWidget(nodesLabel, 3, 2, right);
  pl->addWidget(psNodesSpin, 3, 3);
  psWindowRow_ = {minLabel, psMinEdit, maxLabel, psMaxEdit};
  psHulthenRow_ = {psCustomCheck, aLabel, psAEdit, bLabel, psBEdit};
  psGaussRow_ = {fwhmLabel, psFwhmEdit};
  psTableRow_ = {tableLabel, psTableEdit, psTableButton};
  psNodesRow_ = {nodesLabel, psNodesSpin, meanLabel, meanTsValue};
  QLabel *directionLabel = label(tr("Directions:"), true, tr("Spectator directions (spectatorAngles=)"));
  QLabel *directionNodesLabel = label(tr("Nodes:"), false);
  pl->addWidget(directionLabel, 4, 0, right);
  pl->addWidget(directionCombo, 4, 1);
  pl->addWidget(directionNodesLabel, 4, 2, right);
  pl->addWidget(directionNodesSpin, 4, 3);
  QLabel *thMinLabel = label(QString::fromUtf8("θ<sub>min</sub>:"), true);
  QLabel *thMaxLabel = label(QString::fromUtf8("θ<sub>max</sub>:"), false);
  pl->addWidget(thMinLabel, 5, 0, right);
  pl->addWidget(directionMinEdit, 5, 1);
  pl->addWidget(thMaxLabel, 5, 2, right);
  pl->addWidget(directionMaxEdit, 5, 3);
  QLabel *directionTableLabel = label(tr("Table:"), true);
  QHBoxLayout *dirl = new QHBoxLayout;
  dirl->setContentsMargins(0, 0, 0, 0);
  dirl->addWidget(directionTableEdit, 1);
  dirl->addWidget(directionTableButton);
  pl->addWidget(directionTableLabel, 5, 0, right);
  pl->addLayout(dirl, 5, 1, 1, 3);
  directionRow_ = {directionLabel, directionCombo};
  directionWindowRow_ = {thMinLabel, directionMinEdit, thMaxLabel, directionMaxEdit, directionNodesLabel,
                         directionNodesSpin};
  directionTableRow_ = {directionTableLabel, directionTableEdit, directionTableButton, directionNodesLabel,
                        directionNodesSpin};
  psBox->setLayout(pl);

  distortionBox = new QGroupBox(tr("Distortion"));
  distortionBox->setToolTip(distortionCombo->toolTip());
  // Only the fields of the chosen kind are shown (showDistortionRows).
  QGridLayout *dl = form();
  dl->addWidget(label(tr("Distortion:"), true), 0, 0, right);
  dl->addWidget(distortionCombo, 0, 1);
  QLabel *ratioLabel = label(tr("Ratio:"), false);
  dl->addWidget(ratioLabel, 0, 2, right);
  dl->addWidget(ratioCombo, 0, 3);
  QLabel *angleLabel = label(tr("Angle:"), true, tr("Spectator direction"));
  QHBoxLayout *al = new QHBoxLayout;
  al->setContentsMargins(0, 0, 0, 0);
  al->addWidget(angleKindCombo, 1);
  al->addWidget(angleEdit, 1);
  QWidget *angleBox = new QWidget;
  angleBox->setLayout(al);
  dl->addWidget(angleLabel, 1, 0, right);
  dl->addWidget(angleBox, 1, 1);
  QLabel *refLabel = label(QString::fromUtf8("E<sub>ref</sub>:"), false, tr("Reference energy, R(E_ref) = 1"));
  dl->addWidget(refLabel, 1, 2, right);
  dl->addWidget(distortionRefEdit, 1, 3);
  QLabel *boundLabel = label(tr("Bound state:"), true);
  dl->addWidget(boundLabel, 2, 0, right);
  dl->addWidget(boundCombo, 2, 1);
  QLabel *rminLabel = label(QString::fromUtf8("r<sub>min</sub>:"), false);
  dl->addWidget(rminLabel, 2, 2, right);
  dl->addWidget(rminEdit, 2, 3);
  QWidget *opticalBox[2];
  QLabel *opticalLabel[2] = {label("a + A:", true, tr("opticalAA=: the entrance channel")),
                             label("s + F:", false, tr("opticalSF=: the spectator's exit channel"))};
  for (int c = 0; c < 2; c++) {
    QHBoxLayout *ol = new QHBoxLayout;
    ol->setContentsMargins(0, 0, 0, 0);
    ol->addWidget(opticalCombo[c], 1);
    ol->addWidget(globalCombo[c], 1);
    ol->addWidget(opticalButton[c]);
    opticalBox[c] = new QWidget;
    opticalBox[c]->setLayout(ol);
    dl->addWidget(opticalLabel[c], 3, 2 * c, right);
    dl->addWidget(opticalBox[c], 3, 2 * c + 1);
  }
  QLabel *distortionTableLabel = label(tr("Table:"), true);
  QHBoxLayout *dtl = new QHBoxLayout;
  dtl->setContentsMargins(0, 0, 0, 0);
  dtl->addWidget(distortionTableEdit, 1);
  dtl->addWidget(distortionTableButton);
  QWidget *distortionTableBox = new QWidget;
  distortionTableBox->setLayout(dtl);
  dl->addWidget(distortionTableLabel, 4, 0, right);
  dl->addWidget(distortionTableBox, 4, 1, 1, 3);
  QLabel *rLabel = label("R(E):", true, tr("R at the lowest and highest data point"));
  dl->addWidget(rLabel, 5, 0, right);
  dl->addWidget(distortionValue, 5, 1, 1, 3);
  distortionComputedRows_ = {ratioLabel, ratioCombo, angleLabel, angleBox, refLabel, distortionRefEdit,
                             boundLabel, boundCombo, rminLabel, rminEdit};
  distortionOpticalRow_ = {opticalLabel[0], opticalBox[0], opticalLabel[1], opticalBox[1]};
  distortionTableRow_ = {distortionTableLabel, distortionTableBox};
  distortionValueRow_ = {rLabel, distortionValue};
  distortionBox->setLayout(dl);

  int leftWidth = 0, rightWidth = 0;
  for (QLabel *l : leftLabels) leftWidth = std::max(leftWidth, l->sizeHint().width());
  for (QLabel *l : rightLabels) rightWidth = std::max(rightWidth, l->sizeHint().width());
  for (QGridLayout *g : forms) {
    g->setColumnMinimumWidth(0, leftWidth);
    g->setColumnMinimumWidth(2, rightWidth);
  }

  QHBoxLayout *message = new QHBoxLayout;
  message->setContentsMargins(0, 0, 0, 0);
  message->addWidget(messageIcon);
  message->addWidget(messageLabel, 1);

  // Coherent background (cbackground=): a THM-only amplitude that interferes
  // with the resonances; its Re and Im are fit parameters.  Unchecked: no key.
  coherentBox = new QGroupBox(tr("Coherent background"));
  coherentBox->setCheckable(true);
  coherentBox->setChecked(false);
  coherentBox->setToolTip(
      tr("cbackground=: c(E) = c0 + c1 E times the entrance vertex, added to the resonant HOES amplitude of one "
         "J^pi, entrance (s,l) and exit (s',l') before squaring; THM only (not in direct data). Re and Im of c "
         "are fit parameters; a ticked value is fixed."));
  connect(coherentBox, SIGNAL(toggled(bool)), this, SLOT(coherentToggled(bool)));
  coherentTable = new QTableWidget(0, 8);
  coherentTable->setHorizontalHeaderLabels(QStringList() << QString::fromUtf8("Jπ") << tr("Exit") << tr("Channels")
                                                         << tr("Form") << "Re c0" << "Im c0" << "Re c1" << "Im c1");
  coherentTable->horizontalHeaderItem(1)->setToolTip(tr("Exit pair key (a segment of the experiment)"));
  coherentTable->horizontalHeaderItem(2)->setToolTip(
      tr("all: every entrance (s,l) and exit (s',l') of the J^pi group, each its own amplitude; or one: "
         "s,l,s',l' (e.g. 1/2,0,1/2,1)"));
  coherentTable->horizontalHeaderItem(3)->setToolTip(tr("const: c0; linear: c0 + c1 E (E: c.m. energy, MeV)"));
  for (int k = 4; k < 8; k++) coherentTable->horizontalHeaderItem(k)->setToolTip(tr("Start value; ticked: fixed"));
  coherentTable->verticalHeader()->hide();
  coherentTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  coherentTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  coherentTable->setSelectionMode(QAbstractItemView::SingleSelection);
  coherentTable->setMaximumHeight(coherentTable->horizontalHeader()->sizeHint().height() + 4 * 30);
  connect(coherentTable, SIGNAL(itemChanged(QTableWidgetItem *)), this, SLOT(coherentEdited()));
  coherentAddButton = new QPushButton("+");
  coherentRemoveButton = new QPushButton(QString::fromUtf8("−"));
  coherentAddButton->setToolTip(tr("Add a term"));
  coherentRemoveButton->setToolTip(tr("Remove the selected term"));
  for (QPushButton *b : {coherentAddButton, coherentRemoveButton}) b->setFixedWidth(32);
  connect(coherentAddButton, &QPushButton::clicked, this, [this]() { addCoherentTerm(); });
  connect(coherentRemoveButton, SIGNAL(clicked()), this, SLOT(removeCoherentTerm()));
  QHBoxLayout *cbButtons = new QHBoxLayout;
  cbButtons->setContentsMargins(0, 0, 0, 0);
  cbButtons->addWidget(coherentAddButton);
  cbButtons->addWidget(coherentRemoveButton);
  cbButtons->addStretch(1);
  QVBoxLayout *cl = new QVBoxLayout;
  cl->addWidget(coherentTable);
  cl->addLayout(cbButtons);
  coherentBox->setLayout(cl);

  editorBox = new QWidget;
  QVBoxLayout *el = new QVBoxLayout;
  el->setContentsMargins(0, 0, 0, 0);
  el->setSpacing(8);
  el->addWidget(experimentBox);
  el->addWidget(coherentBox);
  el->addWidget(kinematicsBox);
  el->addWidget(psBox);
  el->addWidget(distortionBox);
  el->addLayout(message);
  el->addStretch(1);
  editorBox->setLayout(el);

  // The editor scrolls rather than squeezing its sections when space is short.
  QScrollArea *editorScroll = new QScrollArea;
  editorScroll->setWidget(editorBox);
  editorScroll->setWidgetResizable(true);
  editorScroll->setFrameShape(QFrame::NoFrame);
  QHBoxLayout *mainLayout = new QHBoxLayout;
  mainLayout->setSpacing(12);
  mainLayout->addWidget(listBox, 0);
  mainLayout->addWidget(editorScroll, 1);
  setLayout(mainLayout);

  refreshTable();
  selectExperiment(records_.isEmpty() ? -1 : 0);
}

QStringList ThmExperimentsPage::experimentLines() const {
  return ThmExperimentRecord::compose(oldLines_, oldRecords_, records_);
}

// ---------------------------------------------------------------------------
// Table and selection

void ThmExperimentsPage::refreshRow(int row) {
  const ThmExperimentRecord &r = records_.at(row);
  QString reaction = QString::fromUtf8("—");
  if (r.hasKinematics())
    reaction = tr("%1 + %2 at %3 MeV, spectator %4")
                   .arg(r.beam.isEmpty() ? "?" : r.beam, r.target.isEmpty() ? "?" : r.target,
                        r.beamEnergy.isEmpty() ? "?" : r.beamEnergy, r.spectator.isEmpty() ? "?" : r.spectator) +
               (r.lineshape ? tr(", line shape") : QString()) + (r.hasWindow() ? tr(", p_s window") : QString());
  if (r.hasDistortion()) reaction += tr(", distortion %1").arg(r.distortion);
  if (r.hasTheta()) reaction += QString::fromUtf8(", θ %1°").arg(r.theta);
  const QString cells[4] = {r.name, ThmExperimentRecord::segmentsListText(r.segments), r.background, reaction};
  for (int c = 0; c < 4; c++) {
    QTableWidgetItem *item = experimentTable->item(row, c);
    if (!item) {
      item = new QTableWidgetItem;
      experimentTable->setItem(row, c, item);
    }
    item->setText(cells[c]);
    QString tip = tr("Background: %1\nReaction: %2").arg(r.background, reaction);
    if (!r.extraTokens.isEmpty()) tip += "\n" + tr("Also kept as written: %1").arg(r.extraTokens.join(' '));
    item->setToolTip(tip);
  }
}

void ThmExperimentsPage::refreshTable() {
  experimentTable->blockSignals(true);
  experimentTable->setRowCount(records_.size());
  for (int row = 0; row < records_.size(); row++) refreshRow(row);
  experimentTable->blockSignals(false);
}

void ThmExperimentsPage::selectExperiment(int row) {
  current_ = (row >= 0 && row < records_.size()) ? row : -1;
  loading_ = true;
  if (current_ >= 0)
    experimentTable->selectRow(current_);
  else
    experimentTable->clearSelection();
  loading_ = false;
  loadEditor();
}

void ThmExperimentsPage::tableSelectionChanged() {
  if (loading_) return;
  QList<QModelIndex> rows = experimentTable->selectionModel()->selectedRows();
  current_ = rows.isEmpty() ? -1 : rows.first().row();
  loadEditor();
}

void ThmExperimentsPage::addExperiment() {
  QSet<QString> names;
  for (const ThmExperimentRecord &r : records_) names.insert(r.name);
  int n = 1;
  while (names.contains(QString("E%1").arg(n))) n++;
  ThmExperimentRecord r;
  r.name = QString("E%1").arg(n);
  records_ << r;
  refreshTable();
  selectExperiment(records_.size() - 1);
}

void ThmExperimentsPage::removeCurrent() {
  if (current_ < 0) return;
  int row = current_;
  records_.removeAt(row);
  refreshTable();
  selectExperiment(std::min(row, (int)records_.size() - 1));
}

// ---------------------------------------------------------------------------
// Editor

QString ThmExperimentsPage::segmentLabel(int key) const {
  const QList<SegmentsDataData> lines = segments_->getLines();
  if (key < 1 || key > lines.size()) return tr("%1: no such line in the Data segments table").arg(key);
  const SegmentsDataData &s = lines.at(key - 1);
  const QList<PairsData> pairs = pairs_->getPairs();
  auto pairText = [&](int pairKey) {
    if (pairKey < 1 || pairKey > pairs.size()) return QString("?");
    const PairsData &p = pairs.at(pairKey - 1);
    return plain(pairs_->getParticleLabel(p, 0)) + "+" + plain(pairs_->getParticleLabel(p, 1));
  };
  QString text = QString("%1: %2 -> %3, %4-%5 MeV, %6")
                     .arg(key)
                     .arg(pairText(s.entrancePairIndex), pairText(s.exitPairIndex))
                     .arg(s.lowEnergy)
                     .arg(s.highEnergy)
                     .arg(s.dataFile);
  QStringList notes;
  if (!s.isActive) notes << tr("inactive");
  if (!s.isTHM) notes << tr("not THM");
  if (!s.varyNorm) notes << tr("fixed norm");
  if (!notes.isEmpty()) text += " (" + notes.join(", ") + ")";
  return text;
}

void ThmExperimentsPage::fillSegmentList() {
  segmentList->blockSignals(true);
  segmentList->clear();
  if (current_ >= 0) {
    const ThmExperimentRecord &r = records_.at(current_);
    QSet<int> others;
    for (int i = 0; i < records_.size(); i++)
      if (i != current_)
        for (int k : records_.at(i).segments) others.insert(k);
    const QList<SegmentsDataData> lines = segments_->getLines();
    QList<int> keys;
    for (int k = 1; k <= lines.size(); k++)
      if ((lines.at(k - 1).isTHM && lines.at(k - 1).varyNorm && !others.contains(k)) || r.segments.contains(k))
        keys << k;
    for (int k : r.segments)
      if (!keys.contains(k)) keys << k;  // beyond the table: shown so that it can be unticked
    for (int k : keys) {
      QListWidgetItem *item = new QListWidgetItem(segmentLabel(k), segmentList);
      item->setData(Qt::UserRole, k);
      item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
      item->setCheckState(r.segments.contains(k) ? Qt::Checked : Qt::Unchecked);
    }
  }
  segmentList->blockSignals(false);
}

void ThmExperimentsPage::loadEditor() {
  loading_ = true;
  editorBox->setEnabled(current_ >= 0);
  removeButton->setEnabled(current_ >= 0);
  ThmExperimentRecord r = current_ >= 0 ? records_.at(current_) : ThmExperimentRecord();
  nameEdit->setText(r.name);
  fillSegmentList();
  int bg = backgroundCombo->findText(r.background);
  backgroundCombo->setCurrentIndex(bg >= 0 ? bg : 0);
  loadTheta(r);
  kinematicsBox->setChecked(r.hasKinematics());
  beamCombo->setEditText(r.beam);
  targetCombo->setEditText(r.target);
  spectatorCombo->setEditText(r.spectator);
  beamEnergyEdit->setWrittenText(r.beamEnergy);
  lineshapeCheck->setChecked(r.lineshape);
  const bool complete = !r.beam.isEmpty() && !r.target.isEmpty() && !r.spectator.isEmpty() && !r.beamEnergy.isEmpty();
  lineshapeCheck->setEnabled(complete);
  loadPs(r);
  psBox->setEnabled(complete);
  loadDistortion(r);
  loadDirections(r);
  updateDistortionItems(complete);
  loadCoherent(r);
  if (derivedTimer_) derivedTimer_->stop();
  showDerivedNow(r);
  loading_ = false;
}

int ThmExperimentsPage::derivedDelayMs_ = 200;

void ThmExperimentsPage::showDerived(const ThmExperimentRecord &x) {
  if (derivedDelayMs_ <= 0) {
    showDerivedNow(x);
    return;
  }
  if (!derivedTimer_) {
    derivedTimer_ = new QTimer(this);
    derivedTimer_->setSingleShot(true);
    connect(derivedTimer_, &QTimer::timeout, this, [this]() {
      if (current_ >= 0 && current_ < records_.size()) showDerivedNow(records_.at(current_));
    });
  }
  derivedTimer_->start(derivedDelayMs_);
}

void ThmExperimentsPage::settleDerived() {
  if (!derivedTimer_ || !derivedTimer_->isActive()) return;
  derivedTimer_->stop();
  if (current_ >= 0 && current_ < records_.size()) showDerivedNow(records_.at(current_));
}

void ThmExperimentsPage::showDerivedNow(const ThmExperimentRecord &x) {
  derivedCount_++;
  QString why, info = derivedInfo(x, &why);
  // R (or the table's w) at the ends of the data.
  double rLo = 0.0, rHi = 0.0;
  QString distortionWhy;
  const QString distortion = x.hasDistortion() ? distortionInfo(x, &distortionWhy, nullptr, nullptr, &rLo, &rHi)
                                               : QString();
  if (!distortion.isEmpty()) info += (info.isEmpty() ? "" : "\n") + distortion;
  if (why.isEmpty()) why = distortionWhy;
  if (why.isEmpty() && !x.segments.isEmpty()) {
    // The engine's refusal of a distortion key (a malformed value, a key
    // without its kind), at once rather than on Accept.
    const QString parse = ThmSettings::checkExperimentLines(QStringList() << x.line());
    for (const char *key : {"distortion", "optical", "spectatorAngle", "boundState", "theta", "cbackground"})
      if (parse.contains(key)) why = parse.mid(parse.indexOf("]: ") + 3);
  }
  if (why.isEmpty() && !x.cbackground.isEmpty()) why = coherentCheck(x);
  // A theta window with entranceL=coherent (the Model page), as the engine refuses it.
  if (why.isEmpty() && x.hasTheta() && entranceL_() == "coherent") why = coherentRefusal();
  updateThetaItems();
  derivedText_ = info.isEmpty() ? why : why.isEmpty() ? info : info + "\n" + why;
  // The values in the sections, compact; the whole text in their tooltips.
  const QString none = QString::fromUtf8("\u2014");
  auto mev = [](double v) { return QString::number(v, 'g', 6) + " MeV"; };
  Reaction r;
  QString ignored;
  const bool haveReaction = reaction(x, r, &ignored);
  bindingValue->setText(haveReaction ? mev(r.bind) : none);
  qfValue->setText(haveReaction ? mev(r.exa - r.bind) : none);
  qfValue->setToolTip(haveReaction ? tr("E(x+A) = %1, E_qf = E(x+A) - B").arg(mev(r.exa)) : QString());
  QString zeta;
  if (haveReaction && x.lineshape) lineshapeInfo(x, &ignored, &zeta);
  zetaValue->setText(zeta.isEmpty() ? none : zeta);
  ThmSpectatorWindow window;
  const bool haveWindow = haveReaction && x.hasWindow() && !windowInfo(x, &ignored, &window).isEmpty();
  if (haveWindow) {
    const double a = window.MeanEs(window.dataE.front()), b = window.MeanEs(window.dataE.back());
    meanTsValue->setText(QString::fromUtf8("%1 … %2 MeV").arg(QString::number(a, 'g', 3), QString::number(b, 'g', 3)));
  } else {
    meanTsValue->setText(none);
  }
  for (QLabel *l : {bindingValue, zetaValue, meanTsValue}) l->setToolTip(info);
  distortionValue->setText(distortion.isEmpty() ? none
                                                : QString::fromUtf8("%1 … %2")
                                                      .arg(QString::number(rLo, 'g', 3), QString::number(rHi, 'g', 3)));
  distortionValue->setToolTip(distortion);
  // B(x+s) from the masses against the pair's B.
  QString detail;
  const QString mismatch = haveReaction ? bindingMismatch(x, &detail) : QString();
  bindingWarningLabel->setText(mismatch);
  bindingWarningLabel->setToolTip(detail);
  bindingWarningIcon->setToolTip(detail);
  bindingWarningLabel->setVisible(!mismatch.isEmpty());
  bindingWarningIcon->setVisible(!mismatch.isEmpty());
  messageLabel->setText(why);
  messageLabel->setVisible(!why.isEmpty());
  messageIcon->setVisible(!why.isEmpty());
}

void ThmExperimentsPage::nameEdited(const QString &text) {
  if (loading_ || current_ < 0) return;
  records_[current_].name = text.trimmed();
  refreshRow(current_);
}

void ThmExperimentsPage::setSegmentsOfCurrent(const QList<int> &segments) {
  if (current_ < 0) return;
  storeSegments(segments);
  loadEditor();
}

void ThmExperimentsPage::storeSegments(const QList<int> &segments) {
  ThmExperimentRecord &r = records_[current_];
  QList<int> sorted = segments;
  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  QList<int> written;
  if (!ThmExperimentRecord::expandSegments(r.segmentsText, written) || written != sorted) r.segmentsText.clear();
  r.segments = sorted;
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::segmentItemChanged(QListWidgetItem *) {
  if (loading_ || current_ < 0) return;
  QList<int> chosen;
  for (int i = 0; i < segmentList->count(); i++)
    if (segmentList->item(i)->checkState() == Qt::Checked) chosen << segmentList->item(i)->data(Qt::UserRole).toInt();
  storeSegments(chosen);  // the list stays as it is while one of its items signals
}

void ThmExperimentsPage::backgroundChanged(int) {
  if (loading_ || current_ < 0) return;
  records_[current_].background = backgroundCombo->currentText();
  refreshRow(current_);
}

// ---------------------------------------------------------------------------
// Exit angle (theta=)

QString ThmExperimentsPage::coherentRefusal() {
  // EData::BuildThmGroups' words.
  return tr("theta= computes the interference of the entrance partial waves exactly (at fixed angle they "
            "interfere); entranceL=coherent is an approximation of the angle-integrated observable and cannot be "
            "combined with it.");
}

QString ThmExperimentsPage::thetaText() const {
  if (thetaCombo->currentData().toString() != "window") return "all";
  return thetaMinEdit->writtenText() + "-" + thetaMaxEdit->writtenText();
}

void ThmExperimentsPage::loadTheta(const ThmExperimentRecord &r) {
  const bool was = loading_;
  loading_ = true;
  thetaCombo->setCurrentIndex(r.hasTheta() ? 1 : 0);
  QString lo = "0", hi = "180";
  if (r.hasTheta()) splitWindow(r.theta, lo, hi);
  thetaMinEdit->setWrittenText(lo);
  thetaMaxEdit->setWrittenText(hi);
  loading_ = was;
  showThetaRows();
  updateThetaItems();
}

void ThmExperimentsPage::showThetaRows() {
  for (QWidget *w : thetaWindowRow_) w->setVisible(thetaCombo->currentData().toString() == "window");
}

void ThmExperimentsPage::updateThetaItems() {
  QStandardItemModel *m = qobject_cast<QStandardItemModel *>(thetaCombo->model());
  if (!m) return;
  const bool coherent = entranceL_() == "coherent";
  QStandardItem *item = m->item(1);
  const bool enabled = !coherent || thetaCombo->currentIndex() == 1;
  item->setFlags(enabled ? item->flags() | Qt::ItemIsEnabled : item->flags() & ~Qt::ItemIsEnabled);
  item->setToolTip(coherent ? coherentRefusal() : QString());
}

void ThmExperimentsPage::thetaEdited() {
  showThetaRows();
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  const QString text = thetaText();
  // "all" is the default: written only if the file wrote it.
  r.theta = text == "all" ? (r.theta == "all" ? r.theta : QString()) : text;
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::kinematicsEdited() {
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  const bool on = kinematicsBox->isChecked();
  r.beam = on ? beamCombo->currentText().trimmed() : QString();
  r.target = on ? targetCombo->currentText().trimmed() : QString();
  r.spectator = on ? spectatorCombo->currentText().trimmed() : QString();
  r.beamEnergy = on ? beamEnergyEdit->writtenText() : QString();
  // The line shape needs the reaction: offered once all four keys are given,
  // and dropped with them.
  lineshapeCheck->setEnabled(on && !r.beam.isEmpty() && !r.target.isEmpty() && !r.spectator.isEmpty() &&
                             !r.beamEnergy.isEmpty());
  if (!on && lineshapeCheck->isChecked()) {
    lineshapeCheck->blockSignals(true);
    lineshapeCheck->setChecked(false);
    lineshapeCheck->blockSignals(false);
  }
  r.lineshape = on && lineshapeCheck->isChecked();
  // So does the spectator-momentum window.
  psBox->setEnabled(lineshapeCheck->isEnabled());
  if (!on && (!r.ps.isEmpty() || !r.psNodes.isEmpty())) {
    r.ps.clear();
    r.psNodes.clear();
    loading_ = true;
    loadPs(r);
    loading_ = false;
  }
  // Coulomb and optical distortion need it too (a table does not).
  if (!on && r.hasComputedDistortion()) {
    r.distortion.clear();
    r.opticalAA.clear();
    r.opticalSF.clear();
    r.spectatorAngle.clear();
    r.distortionRef.clear();
    r.distortionRatio.clear();
    r.boundState.clear();
    r.spectatorAngles.clear();
    r.spectatorAngleNodes.clear();
    loading_ = true;
    loadDistortion(r);
    loadDirections(r);
    loading_ = false;
  }
  updateDistortionItems(lineshapeCheck->isEnabled());
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::lineshapeToggled(bool on) {
  if (loading_ || current_ < 0) return;
  records_[current_].lineshape = on && kinematicsBox->isChecked();
  refreshRow(current_);
  showDerived(records_.at(current_));
}

// ---------------------------------------------------------------------------
// The engine's rules

bool ThmExperimentsPage::reaction(const ThmExperimentRecord &x, Reaction &out, QString *error) const {
  if (error) error->clear();
  if (x.beam.isEmpty() || x.target.isEmpty() || x.spectator.isEmpty() || x.beamEnergy.isEmpty()) return false;
  ThmNuclide b, t, sp;
  std::string why;
  if (!(why = ThmNuclide::Parse(x.beam.toStdString(), b)).empty() ||
      !(why = ThmNuclide::Parse(x.target.toStdString(), t)).empty() ||
      !(why = ThmNuclide::Parse(x.spectator.toStdString(), sp)).empty()) {
    if (error) *error = QString::fromStdString(why);
    return false;
  }
  double ebeam;
  if (!readWholeDouble(x.beamEnergy, ebeam) || !(ebeam > 0.0)) {
    if (error) *error = tr("Ebeam='%1': expected the lab beam energy in MeV, > 0").arg(x.beamEnergy);
    return false;
  }
  // The entrance pair x + A of the segments in use (EData::SetupThmExperiments).
  const QList<SegmentsDataData> lines = segments_->getLines();
  const QList<PairsData> pairs = pairs_->getPairs();
  int pairKey = -1;
  for (int k : x.segments) {
    if (k < 1 || k > lines.size() || !lines.at(k - 1).isActive) continue;
    int e = lines.at(k - 1).entrancePairIndex;
    if (pairKey < 0)
      pairKey = e;
    else if (e != pairKey) {
      if (error)
        *error = tr("beam/target/spectator describe one reaction, but its segments have different entrance pairs.");
      return false;
    }
  }
  if (pairKey < 1 || pairKey > pairs.size()) return false;  // nothing to compare with yet
  const PairsData &pair = pairs.at(pairKey - 1);
  const int Z[2] = {pair.lightZ, pair.heavyZ};
  const int A[2] = {(int)std::lround(pair.lightM), (int)std::lround(pair.heavyM)};
  const double M[2] = {pair.lightM, pair.heavyM};
  int horse = -1, other = -1;
  for (int h = 0; h < 2 && horse < 0; h++) {
    const ThmNuclide &th = h == 0 ? b : t, &tg = h == 0 ? t : b;
    for (int k = 0; k < 2; k++)
      if (tg.Z == Z[k] && tg.A == A[k] && th.Z - sp.Z == Z[1 - k] && th.A - sp.A == A[1 - k]) {
        horse = h;
        other = k;
        break;
      }
  }
  if (horse < 0) {
    if (error)
      *error = tr("beam %1 + target %2 with spectator %3 does not give the entrance pair of its segments (Z,A) = "
                  "(%4,%5) + (%6,%7): one of beam/target must be a nucleus of the pair and the other the second "
                  "nucleus plus the spectator.")
                   .arg(QString::fromStdString(b.name), QString::fromStdString(t.name),
                        QString::fromStdString(sp.name))
                   .arg(Z[0])
                   .arg(A[0])
                   .arg(Z[1])
                   .arg(A[1]);
    return false;
  }
  const ThmNuclide &th = horse == 0 ? b : t, &nA = horse == 0 ? t : b;
  const ThmNuclide *tabX = ThmNuclide::Find(th.Z - sp.Z, th.A - sp.A);
  const double mX = tabX ? tabX->mass : M[1 - other];
  const double mA = nA.mass;
  out.beam = b;
  out.target = t;
  out.spectator = sp;
  out.horse = th;
  out.pairKey = pairKey;
  out.beamEnergy = ebeam;
  out.mX = mX;
  out.horseIsBeam = horse == 0;
  out.other = nA;
  out.bind = (mX + sp.mass - th.mass) * kAmu;
  out.exa = horse == 0 ? ebeam * mX / th.mass * mA / (mX + mA) : ebeam * mX / (mA + mX);
  return true;
}

QString ThmExperimentsPage::derivedInfo(const ThmExperimentRecord &x, QString *error) const {
  Reaction r;
  if (!reaction(x, r, error)) return QString();
  QString text =
      tr("B(x+s) = %3 MeV (Trojan horse %1 = x + %2)\nquasi-free E(x+A) = %4 MeV, E_qf = E(x+A) - B = %5 MeV")
          .arg(QString::fromStdString(r.horse.name), QString::fromStdString(r.spectator.name))
          .arg(QString::number(r.bind, 'g', 6))
          .arg(QString::number(r.exa, 'g', 6))
          .arg(QString::number(r.exa - r.bind, 'g', 6));
  QString shapeError, windowError;
  if (x.lineshape) {
    const QString shape = lineshapeInfo(x, &shapeError);
    if (!shape.isEmpty()) text += "\n" + shape;
  }
  if (x.hasWindow()) {
    const QString window = windowInfo(x, &windowError);
    if (!window.isEmpty()) text += "\n" + window;
  }
  if (error) *error = !shapeError.isEmpty() ? shapeError : windowError;
  return text;
}

bool ThmExperimentsPage::pointRange(const QList<int> &segments, double &lo, double &hi) const {
  QVector<double> energies;
  lo = 1.0e300;
  hi = -1.0e300;
  if (!pointEnergies(segments, energies)) return false;
  for (double e : energies) {
    lo = std::min(lo, e);
    hi = std::max(hi, e);
  }
  return lo <= hi;
}

bool ThmExperimentsPage::pointEnergies(const QList<int> &segments, QVector<double> &energies) const {
  // The c.m. energies of the points AZURE2 reads (ESegment::FillData): rows
  // of the data file inside the segment's lab energy range (and angle range
  // if differential), converted with the entrance pair's masses.
  const QList<SegmentsDataData> lines = segments_->getLines();
  const QList<PairsData> pairs = pairs_->getPairs();
  energies.clear();
  for (int k : segments) {
    if (k < 1 || k > lines.size() || !lines.at(k - 1).isActive) continue;
    const SegmentsDataData &s = lines.at(k - 1);
    if (s.entrancePairIndex < 1 || s.entrancePairIndex > pairs.size()) return false;
    const PairsData &pair = pairs.at(s.entrancePairIndex - 1);
    const double factor = pair.heavyM / (pair.lightM + pair.heavyM);
    const bool differential = s.dataType == 1 || s.dataType == 4 || s.dataType == 7 || s.dataType == 8;
    QString path = s.dataFile;
    if (QFileInfo(path).isRelative()) path = QDir(projectDir_).filePath(path);
    std::ifstream in(QFile::encodeName(path).constData());
    if (!in) return false;
    while (true) {
      DataLine line(in);
      if (line.atEnd()) break;
      if (!line.valid()) return false;
      if (line.energy() < s.lowEnergy || line.energy() > s.highEnergy) continue;
      if (differential && (line.angle() < s.lowAngle || line.angle() > s.highAngle)) continue;
      energies << line.energy() * factor;
    }
  }
  return !energies.isEmpty();
}

QString ThmExperimentsPage::lineshapeInfo(const ThmExperimentRecord &x, QString *error, QString *compact) const {
  Reaction r;
  if (!reaction(x, r, error)) return QString();
  if (!brune_) {
    if (error)
      *error = tr("lineshape=on uses the observed level energies and widths as the resonance poles; it needs the "
                  "Brune parameterization.");
    return QString();
  }
  const QList<PairsData> pairs = pairs_->getPairs();
  const PairsData &entrance = pairs.at(r.pairKey - 1);
  ThmLineshape shape;  // as EData::BuildThmGroups sets it up
  shape.Zs = r.spectator.Z;
  shape.ms = r.spectator.mass;
  shape.ZF = entrance.lightZ + entrance.heavyZ;
  shape.mF = entrance.lightM + entrance.heavyM;
  shape.eAA = r.beamEnergy * r.target.mass / (r.beam.mass + r.target.mass);
  shape.bind = r.bind;
  double lo, hi;
  if (!pointRange(x.segments, lo, hi)) return QString();  // the engine reports an unreadable file itself
  if (!(shape.EsF(hi) > 0.0)) {
    if (error)
      *error = tr("lineshape=on: at E = %1 MeV the spectator has no energy left (E_sF = E_aA - B - E = %2 - %3 - "
                  "%1 MeV <= 0); check Ebeam.")
                   .arg(QString::number(hi, 'g', 6), QString::number(shape.eAA, 'g', 6),
                        QString::number(shape.bind, 'g', 6));
    return QString();
  }
  // zeta per exit pair of the segments, at the lowest and highest point energy.
  const QList<SegmentsDataData> lines = segments_->getLines();
  QList<int> exits;
  QStringList parts;
  for (int k : x.segments) {
    if (k < 1 || k > lines.size() || !lines.at(k - 1).isActive) continue;
    const int key = lines.at(k - 1).exitPairIndex;
    if (exits.contains(key) || key < 1 || key > pairs.size()) continue;
    exits << key;
    const PairsData &p = pairs.at(key - 1);
    const bool lightFirst = p.lightM <= p.heavyM;
    const int ZB = lightFirst ? p.heavyZ : p.lightZ;
    const double mB = lightFirst ? p.heavyM : p.lightM;
    parts << tr("%1 into %2")
                 .arg(QString::fromUtf8("%1 \u2026 %2")
                          .arg(QString::number(shape.Zeta(lo, ZB, mB), 'g', 3),
                               QString::number(shape.Zeta(hi, ZB, mB), 'g', 3)),
                      plain(pairs_->getParticleLabel(p, 0)) + "+" + plain(pairs_->getParticleLabel(p, 1)));
  }
  if (compact) *compact = parts.join("; ").replace(" into ", " (") .replace("; ", "); ") + ")";
  return QString::fromUtf8("\u03b6 = ") + parts.join("; ") +
         tr(" at E = %1 \u2026 %2 MeV (the ends of the data)")
             .arg(QString::number(lo, 'g', 4), QString::number(hi, 'g', 4));
}

QString ThmExperimentsPage::check() const {
  QString why = ThmSettings::checkExperimentLines(experimentLines(), coulombIntegral_());
  if (!why.isEmpty()) return why;
  // The startup checks of EData::SetupThmExperiments that need the project.
  const QList<SegmentsDataData> lines = segments_->getLines();
  for (const ThmExperimentRecord &x : records_) {
    const QString where = QString("<thm> experiment[%1]: ").arg(x.name);
    for (int k : x.segments) {
      if (k > lines.size())
        return where + tr("segment %1: <segmentsData> has only %2 line(s).").arg(k).arg(lines.size());
      const SegmentsDataData &s = lines.at(k - 1);
      if (!s.isActive) continue;  // left out with a warning
      if (!s.isTHM) return where + tr("segment %1 is not a THM segment (isDiff < 10).").arg(k);
      if (!s.varyNorm)
        return where + tr("segment %1 has a fixed norm; the segments of an experiment share one free "
                          "(profiled) norm, so free it.")
                           .arg(k);
    }
    if (x.hasTheta() && entranceL_() == "coherent") return where + coherentRefusal();
    if (!x.cbackground.isEmpty()) {
      const QString c = coherentCheck(x);
      if (!c.isEmpty()) return where + c;
    }
    QString error;
    derivedInfo(x, &error);
    if (!error.isEmpty()) return where + error;
    if (x.hasDistortion()) {
      distortionInfo(x, &error);
      if (!error.isEmpty()) return where + error;
    }
  }
  return QString();
}

// ---------------------------------------------------------------------------
// Spectator-momentum window

namespace {

// "pmin-pmax" as written: split at the first '-' with a number on both sides
// (the engine's ReadWindow); otherwise at the first '-' after the start.
void splitWindow(const QString &text, QString &lo, QString &hi) {
  for (int k = 1; k + 1 < text.size(); k++) {
    if (text[k] != '-') continue;
    double a, b;
    if (readWholeDouble(text.left(k), a) && readWholeDouble(text.mid(k + 1), b)) {
      lo = text.left(k);
      hi = text.mid(k + 1);
      return;
    }
  }
  const int dash = text.indexOf('-', 1);
  lo = dash < 0 ? text : text.left(dash);
  hi = dash < 0 ? QString() : text.mid(dash + 1);
}

}  // namespace

QString ThmExperimentsPage::projectRelative(const QString &file, const QString &projectDir) {
  if (projectDir.isEmpty() || file.isEmpty()) return file;
  const QString rel = QDir(projectDir).relativeFilePath(file);
  return rel.startsWith("..") || QDir::isAbsolutePath(rel) ? file : rel;
}

void ThmExperimentsPage::loadPs(const ThmExperimentRecord &r) {
  // The controls show the value as written, so that composing them again
  // (psText) gives it back unchanged.
  QString kind = "delta", lo, hi, a = "0.2317", b = "1.202", fwhm, table;
  bool custom = false;
  if (r.ps.startsWith("table:")) {
    kind = "table";
    table = r.ps.mid(6);
  } else if (!r.ps.isEmpty() && r.ps != "delta") {
    const QStringList f = r.ps.split(':');
    kind = f[0];
    if (kind == "hulthen" && f.size() == 3) {
      custom = true;
      const QStringList ab = f[1].split(',');
      a = ab.value(0);
      b = ab.value(1);
    } else if (kind == "gauss" && f.size() == 3) {
      fwhm = f[1];
    }
    splitWindow(f.last(), lo, hi);
    if (f.size() < 2) lo = hi = QString();
  }
  const bool was = loading_;
  loading_ = true;
  int at = psKindCombo->findData(kind);
  psKindCombo->setCurrentIndex(at >= 0 ? at : 0);
  psMinEdit->setWrittenText(lo);
  psMaxEdit->setWrittenText(hi);
  psCustomCheck->setChecked(custom);
  psAEdit->setWrittenText(a);
  psBEdit->setWrittenText(b);
  psFwhmEdit->setWrittenText(fwhm);
  psTableEdit->setText(table);
  int nodes = 16;
  if (!r.psNodes.isEmpty()) {
    bool ok = false;
    const int n = r.psNodes.trimmed().toInt(&ok);
    if (ok) nodes = n;
  }
  psNodesSpin->setValue(std::min(64, std::max(1, nodes)));
  loading_ = was;
  showPsRows();
}

void ThmExperimentsPage::showPsRows() {
  // Only the fields of the chosen distribution are shown.
  const QString kind = psKindCombo->currentData().toString();
  for (QWidget *w : psWindowRow_) w->setVisible(kind == "hulthen" || kind == "gauss");
  for (QWidget *w : psHulthenRow_) w->setVisible(kind == "hulthen");
  for (QWidget *w : psGaussRow_) w->setVisible(kind == "gauss");
  for (QWidget *w : psTableRow_) w->setVisible(kind == "table");
  for (QWidget *w : psNodesRow_) w->setVisible(kind != "delta");
  psAEdit->setEnabled(psCustomCheck->isChecked());
  psBEdit->setEnabled(psCustomCheck->isChecked());
}

QString ThmExperimentsPage::psText() const {
  const QString kind = psKindCombo->currentData().toString();
  const QString window = psMinEdit->writtenText() + "-" + psMaxEdit->writtenText();
  if (kind == "hulthen")
    return psCustomCheck->isChecked()
               ? "hulthen:" + psAEdit->writtenText() + "," + psBEdit->writtenText() + ":" + window
               : "hulthen:" + window;
  if (kind == "gauss") return "gauss:" + psFwhmEdit->writtenText() + ":" + window;
  if (kind == "table") return "table:" + psTableEdit->text().trimmed();
  return QString();
}

void ThmExperimentsPage::psEdited() {
  showPsRows();
  // The Directions rows depend on the distribution (not while the page is
  // still being built: the spin boxes of the window signal as they are set up).
  if (directionCombo && distortionCombo && angleKindCombo && angleEdit) showDirectionRows();
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  if (!psCustomCheck->isChecked()) {
    // Back to the deuteron's values, shown as the engine's defaults.
    const QSignalBlocker blockA(psAEdit), blockB(psBEdit);
    psAEdit->setWrittenText("0.2317");
    psBEdit->setWrittenText("1.202");
  }
  r.ps = psText();
  if (!r.ps.isEmpty() && !r.spectatorAngle.isEmpty()) {
    // A momentum window is a set of directions: the one direction is dropped
    // (the engine refuses both).
    r.spectatorAngle.clear();
    loading_ = true;
    angleKindCombo->setCurrentIndex(std::max(0, angleKindCombo->findData("qf")));
    angleEdit->setWrittenText(QString());
    loading_ = false;
  }
  if (r.ps.isEmpty() && !r.psNodes.isEmpty()) {
    // psNodes needs a window: dropped with it.
    r.psNodes.clear();
    loading_ = true;
    psNodesSpin->setValue(16);
    loading_ = false;
  }
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::psNodesChanged(int n) {
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  r.psNodes = n == 16 ? QString() : QString::number(n);  // 16 is the engine's default
  showDerived(r);
}

void ThmExperimentsPage::chooseTable() {
  const QString start = projectDir_.isEmpty() ? QDir::currentPath() : projectDir_;
  const QString file = QFileDialog::getOpenFileName(this, tr("Spectator-momentum table"), start,
                                                    tr("Tables (*.dat *.txt);;All files (*)"));
  if (file.isEmpty()) return;
  psTableEdit->setText(projectRelative(file, projectDir_));
  psEdited();
}

void ThmExperimentsPage::refreshDerived() {
  if (current_ >= 0) showDerived(records_.at(current_));
}

// ---------------------------------------------------------------------------
// Spectator directions

void ThmExperimentsPage::loadDirections(const ThmExperimentRecord &r) {
  // As written, so that directionText gives it back unchanged.
  QString kind = "one", lo, hi, table, v = r.spectatorAngles;
  const bool cm = v.startsWith("cm:");
  if (cm) v = v.mid(3);
  if (v.startsWith("table:")) {
    kind = cm ? "cmtable" : "labtable";
    table = v.mid(6);
  } else if (!v.isEmpty()) {
    kind = cm ? "cm" : "lab";
    splitWindow(v, lo, hi);
  }
  const bool was = loading_;
  loading_ = true;
  directionCombo->setCurrentIndex(std::max(0, directionCombo->findData(kind)));
  directionMinEdit->setWrittenText(lo);
  directionMaxEdit->setWrittenText(hi);
  directionTableEdit->setText(table);
  int nodes = 8;
  if (!r.spectatorAngleNodes.isEmpty()) {
    bool ok = false;
    const int n = r.spectatorAngleNodes.trimmed().toInt(&ok);
    if (ok) nodes = n;
  }
  directionNodesSpin->setValue(std::min(64, std::max(1, nodes)));
  loading_ = was;
  showDirectionRows();
}

void ThmExperimentsPage::showDirectionRows() {
  // With a computed distortion or a momentum distribution (the engine
  // refuses the window with neither); the window's or the table's fields by
  // the kind.
  const QString distortion = distortionCombo->currentData().toString();
  const bool psWindow = psKindCombo->currentData().toString() != "delta";
  const bool computed = distortion == "coulomb" || distortion == "optical" || psWindow;
  const QString kind = directionCombo->currentData().toString();
  // Without a direction window, a momentum window accepts every direction
  // inside it; otherwise the one direction of the Distortion section.
  directionCombo->setItemText(0, psWindow ? tr("all inside the p_s window") : tr("one (Distortion: angle)"));
  // psNodes= counts the nodes of the momentum window alone.
  psNodesSpin->setEnabled(kind == "one");
  for (QWidget *w : directionWindowRow_) w->setVisible(false);
  for (QWidget *w : directionTableRow_) w->setVisible(false);
  for (QWidget *w : directionRow_) w->setVisible(computed);
  if (computed && (kind == "lab" || kind == "cm"))
    for (QWidget *w : directionWindowRow_) w->setVisible(true);
  if (computed && kind.endsWith("table"))
    for (QWidget *w : directionTableRow_) w->setVisible(true);
  // A window sets the directions: the Distortion section's one angle is off.
  const bool window = (computed && kind != "one") || psWindow;
  angleKindCombo->setEnabled(!window);
  angleEdit->setEnabled(!window && angleKindCombo->currentData().toString() != "qf");
  angleKindCombo->setToolTip(window ? tr("The spectator directions are set by the window in Spectator acceptance.")
                                    : tr("spectatorAngle=: the direction of the spectator. Quasi-free: k_sF along "
                                         "k_aA (default); lab: an angle to the beam converted at every E; c.m.: "
                                         "fixed."));
}

QString ThmExperimentsPage::directionText() const {
  const QString kind = directionCombo->currentData().toString();
  const QString window = directionMinEdit->writtenText() + "-" + directionMaxEdit->writtenText();
  if (kind == "lab") return window;
  if (kind == "cm") return "cm:" + window;
  if (kind == "labtable") return "table:" + directionTableEdit->text().trimmed();
  if (kind == "cmtable") return "cm:table:" + directionTableEdit->text().trimmed();
  return QString();
}

void ThmExperimentsPage::directionEdited() {
  showDirectionRows();
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  r.spectatorAngles = directionText();
  if (!r.spectatorAngles.isEmpty() && !r.psNodes.isEmpty()) {
    // The direction window has its own nodes (the engine refuses psNodes= with it).
    r.psNodes.clear();
    loading_ = true;
    psNodesSpin->setValue(16);
    loading_ = false;
  }
  if (r.spectatorAngles.isEmpty() && !r.spectatorAngleNodes.isEmpty()) {
    // The nodes need a window: dropped with it.
    r.spectatorAngleNodes.clear();
    loading_ = true;
    directionNodesSpin->setValue(8);
    loading_ = false;
  }
  if (!r.spectatorAngles.isEmpty() && !r.spectatorAngle.isEmpty()) {
    // One direction and a window exclude each other (the engine refuses both).
    r.spectatorAngle.clear();
    loading_ = true;
    angleKindCombo->setCurrentIndex(std::max(0, angleKindCombo->findData("qf")));
    angleEdit->setWrittenText(QString());
    loading_ = false;
    showDirectionRows();
  }
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::directionNodesChanged(int n) {
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  r.spectatorAngleNodes = n == 8 ? QString() : QString::number(n);  // 8 is the engine's default
  showDerived(r);
}

void ThmExperimentsPage::chooseDirectionTable() {
  const QString start = projectDir_.isEmpty() ? QDir::currentPath() : projectDir_;
  const QString file = QFileDialog::getOpenFileName(this, tr("Spectator acceptance table"), start,
                                                    tr("Tables (*.dat *.txt);;All files (*)"));
  if (file.isEmpty()) return;
  directionTableEdit->setText(projectRelative(file, projectDir_));
  directionEdited();
}

QString ThmExperimentsPage::windowInfo(const ThmExperimentRecord &x, QString *error, ThmSpectatorWindow *out) const {
  Reaction r;
  if (!reaction(x, r, error)) return QString();
  // The engine's parse of the line, and its table reader (Config::ReadThmBlock).
  std::vector<ThmExperiment> parsed;
  std::string why = ParseThmExperimentLine(x.line().toStdString(), parsed);
  if (!why.empty() || parsed.empty()) {
    if (error) *error = "<thm> " + QString::fromStdString(why);
    return QString();
  }
  ThmExperiment &e = parsed.front();
  if (e.psKind == ThmExperiment::PS_DELTA) return QString();
  if (spectatorEnergy_(r.pairKey) != 0.0) {
    if (error)
      *error = tr("a ps window and spectatorEnergy both set the spectator motion of entrance pair %1; use one "
                  "(ps=delta keeps spectatorEnergy).")
                   .arg(r.pairKey);
    return QString();
  }
  if (e.psKind == ThmExperiment::PS_TABLE) {
    QString path = QString::fromStdString(e.psTable);
    if (QFileInfo(path).isRelative() && !projectDir_.isEmpty()) path = QDir(projectDir_).filePath(path);
    why = ReadThmPsTable(QFile::encodeName(path).toStdString(), e.psTableP, e.psTableW);
    if (!why.empty()) {
      if (error) *error = "ps: " + QString::fromStdString(why);
      return QString();
    }
    e.psMin = e.psTableP.front();
    e.psMax = e.psTableP.back();
  }
  if (e.angleWindow == 2) {
    // spectatorAngles=[cm:]table:<file>, relative to the .azr (Config::ReadThmBlock).
    QString path = QString::fromStdString(e.angleTable);
    if (QFileInfo(path).isRelative() && !projectDir_.isEmpty()) path = QDir(projectDir_).filePath(path);
    why = ReadThmAngleTable(QFile::encodeName(path).toStdString(), e.angleTableT, e.angleTableW);
    if (!why.empty()) {
      if (error) *error = "spectatorAngles: " + QString::fromStdString(why);
      return QString();
    }
    e.angleMin = e.angleTableT.front();
    e.angleMax = e.angleTableT.back();
  }
  // mu_sx, the kinematics and the data energies as EData::BuildThmGroups takes them.
  const double muSx = r.mX * r.spectator.mass / (r.mX + r.spectator.mass) * kAmu;
  ThmDistortion d;
  QVector<double> energies;
  if (!fillKinematics(x, d) || !pointEnergies(x.segments, energies) || energies.isEmpty()) return QString();
  ThmSpectatorWindow window;
  why = BuildThmSpectatorWindow(e, muSx, d.kin, std::vector<double>(energies.begin(), energies.end()), window);
  if (!why.empty()) {
    if (error) *error = "ps: " + QString::fromStdString(why);
    return QString();
  }
  if (out) *out = window;
  const double lo = window.dataE.front(), hi = window.dataE.back();
  auto number = [](double v) { return QString::number(v, 'g', 6); };
  QString text = tr("p_s window: %1; mu_sx = %2 MeV, <T_s> = %3 MeV at E = %4 MeV, %5 MeV at E = %6 MeV")
                     .arg(QString::fromStdString(window.description), number(muSx), number(window.MeanEs(lo)),
                          number(lo), number(window.MeanEs(hi)), number(hi));
  if (e.vertexDW)
    text += tr(" (the DW vertex weights these directions with the |phi~(q)|^2 of its own bound state; the "
               "distribution above sets only the cut)");
  return text;
}

// ---------------------------------------------------------------------------
// Distortion factor

namespace {

// A key the page sets to `value`: "" for the engine's default, unless the
// file wrote the default out ("qf", "dwpw", ...), which is then kept.
QString keyValue(const QString &old, const QString &value, const QString &byDefault) {
  if (value != byDefault) return value;
  return old == byDefault ? old : QString();
}

const char *kTenZeros = "0,0,0,0,0,0,0,0,0,0";

// opticalAA/SF as written -> combo data.
QString opticalKind(const QString &value) {
  if (value == "plane") return "plane";
  if (value.isEmpty() || value == "coulomb") return "coulomb";
  if (!value.isEmpty() && value.at(0).isLetter()) return "global";
  return "ws";
}

// The first global potential that describes a partner (Z, A) of a channel, or "".
QString defaultGlobal(int z1, int a1, int z2, int a2) {
  for (int m = 0; m < (int)ThmGlobalOpticals().size(); m++)
    if (ThmGlobalOpticalFor(m, z1, a1) || ThmGlobalOpticalFor(m, z2, a2)) return ThmGlobalOpticals()[m].name;
  return QString();
}

}  // namespace

void ThmExperimentsPage::loadDistortion(const ThmExperimentRecord &r) {
  const bool was = loading_;
  loading_ = true;
  QString kind = "none", table;
  if (r.distortion == "coulomb" || r.distortion == "optical") kind = r.distortion;
  if (r.distortion.startsWith("table:")) {
    kind = "table";
    table = r.distortion.mid(6);
  }
  int at = distortionCombo->findData(kind);
  distortionCombo->setCurrentIndex(at >= 0 ? at : 0);
  distortionTableEdit->setText(table);
  QString angleKind = "qf", angle;
  if (r.spectatorAngle.startsWith("cm:")) {
    angleKind = "cm";
    angle = r.spectatorAngle.mid(3);
  } else if (!r.spectatorAngle.isEmpty() && r.spectatorAngle != "qf") {
    angleKind = "lab";
    angle = r.spectatorAngle;
  }
  angleKindCombo->setCurrentIndex(std::max(0, angleKindCombo->findData(angleKind)));
  angleEdit->setWrittenText(angle);
  distortionRefEdit->setWrittenText(r.distortionRef);
  ratioCombo->setCurrentIndex(r.distortionRatio == "dw" ? 1 : 0);
  const QStringList bound = r.boundState.split(':');
  boundCombo->setCurrentIndex(bound.value(0) == "yukawa" ? 1 : 0);
  rminEdit->setWrittenText(bound.size() > 1 ? bound.value(1) : QString());
  const QString *optical[2] = {&r.opticalAA, &r.opticalSF};
  for (int c = 0; c < 2; c++) {
    const QString k = opticalKind(*optical[c]);
    opticalCombo[c]->setCurrentIndex(std::max(0, opticalCombo[c]->findData(k)));
    lastOptical_[c] = k == "ws" ? *optical[c] : QString(kTenZeros);
    if (k == "global") {
      lastGlobal_[c] = *optical[c];
      const int at = globalCombo[c]->findData(optical[c]->section(':', 0, 0));
      if (at >= 0) globalCombo[c]->setCurrentIndex(at);
    }
  }
  loading_ = was;
  showDistortionRows();
}

void ThmExperimentsPage::showDistortionRows() {
  // Only the fields of the chosen kind; none: the combo alone.
  const QString kind = distortionCombo->currentData().toString();
  const bool computed = kind == "coulomb" || kind == "optical";
  for (QWidget *w : distortionComputedRows_) w->setVisible(computed);
  for (QWidget *w : distortionOpticalRow_) w->setVisible(kind == "optical");
  for (QWidget *w : distortionTableRow_) w->setVisible(kind == "table");
  for (QWidget *w : distortionValueRow_) w->setVisible(kind != "none");
  angleEdit->setEnabled(angleKindCombo->currentData().toString() != "qf");
  for (int c = 0; c < 2; c++) {
    const QString mode = opticalCombo[c]->currentData().toString();
    const bool ws = mode == "ws", global = mode == "global";
    globalCombo[c]->setVisible(global && kind == "optical");
    opticalButton[c]->setEnabled(ws || global);
    QString tip = tr("The ten Woods-Saxon numbers");
    if (ws) tip = tr("V,R,a,W,RW,aW,WD,RD,aD,RC = %1").arg(lastOptical_[c]);
    if (global && current_ >= 0 && current_ < records_.size()) tip = globalSummary(records_.at(current_), c);
    opticalButton[c]->setToolTip(tip);
  }
  showDirectionRows();
}

void ThmExperimentsPage::updateDistortionItems(bool complete) {
  // Coulomb and optical need the three-body reaction (all four keys), as the
  // engine's CheckThmExperiments; the current kind stays selectable, so that
  // its refusal can be seen.
  QStandardItemModel *m = qobject_cast<QStandardItemModel *>(distortionCombo->model());
  if (!m) return;
  for (int i = 0; i < distortionCombo->count(); i++) {
    const QString kind = distortionCombo->itemData(i).toString();
    const bool needs = kind == "coulomb" || kind == "optical";
    QStandardItem *item = m->item(i);
    const bool enabled = !needs || complete || distortionCombo->currentIndex() == i;
    item->setFlags(enabled ? item->flags() | Qt::ItemIsEnabled : item->flags() & ~Qt::ItemIsEnabled);
    item->setToolTip(needs && !complete ? tr("Needs the three-body reaction (beam, target, spectator, Ebeam).")
                                        : QString());
  }
}

void ThmExperimentsPage::distortionKindChanged() {
  showDistortionRows();
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  const QString kind = distortionCombo->currentData().toString();
  if (kind == "table")
    r.distortion = "table:" + distortionTableEdit->text().trimmed();
  else
    r.distortion = keyValue(r.distortion, kind, "none");
  // Keys the new kind does not take (the engine refuses them) go.
  if (kind != "optical") {
    r.opticalAA.clear();
    r.opticalSF.clear();
  }
  if (kind != "coulomb" && kind != "optical") {
    r.spectatorAngle.clear();
    r.distortionRef.clear();
    r.distortionRatio.clear();
    r.boundState.clear();
    r.spectatorAngles.clear();
    r.spectatorAngleNodes.clear();
    loadDirections(r);
  }
  loadDistortion(r);
  updateDistortionItems(lineshapeCheck->isEnabled());
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::distortionEdited() {
  showDistortionRows();
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  // Only the key of the control that changed is rewritten; the others stay as written.
  QObject *from = sender();
  if (from == distortionTableEdit) {
    r.distortion = "table:" + distortionTableEdit->text().trimmed();
  } else if (from == angleKindCombo || from == angleEdit) {
    const QString kind = angleKindCombo->currentData().toString();
    r.spectatorAngle = kind == "qf"    ? keyValue(r.spectatorAngle, "qf", "qf")
                       : kind == "cm" ? "cm:" + angleEdit->writtenText()
                                      : angleEdit->writtenText();
  } else if (from == distortionRefEdit) {
    r.distortionRef = distortionRefEdit->writtenText();
  } else if (from == ratioCombo) {
    r.distortionRatio = keyValue(r.distortionRatio, ratioCombo->currentData().toString(), "dwpw");
  } else if (from == boundCombo || from == rminEdit) {
    const QString rmin = rminEdit->writtenText();
    r.boundState = keyValue(r.boundState, boundCombo->currentData().toString() + (rmin.isEmpty() ? "" : ":" + rmin),
                            "whittaker");
  }
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::opticalKindChanged() {
  showDistortionRows();
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  QString *optical[2] = {&r.opticalAA, &r.opticalSF};
  for (int c = 0; c < 2; c++) {
    if (sender() != opticalCombo[c]) continue;
    const QString kind = opticalCombo[c]->currentData().toString();
    if (kind == "global") {
      if (lastGlobal_[c].isEmpty()) {
        // The first model for the channel's light partner.
        Reaction re;
        QString ignored;
        if (reaction(r, re, &ignored)) {
          const ThmNuclide &s = re.spectator, &h = re.horse, &o = re.other;
          lastGlobal_[c] = c == 0 ? defaultGlobal(h.Z, h.A, o.Z, o.A)
                                  : defaultGlobal(s.Z, s.A, h.Z - s.Z + o.Z, h.A - s.A + o.A);
        }
        if (lastGlobal_[c].isEmpty()) lastGlobal_[c] = globalCombo[c]->currentData().toString();
      }
      *optical[c] = lastGlobal_[c];
      loading_ = true;
      const int at = globalCombo[c]->findData(lastGlobal_[c].section(':', 0, 0));
      if (at >= 0) globalCombo[c]->setCurrentIndex(at);
      loading_ = false;
    } else {
      *optical[c] = kind == "ws" ? lastOptical_[c] : keyValue(*optical[c], kind, "coulomb");
    }
  }
  showDistortionRows();
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::globalNameChanged() {
  if (loading_ || current_ < 0) return;
  for (int c = 0; c < 2; c++) {
    if (sender() != globalCombo[c] || opticalCombo[c]->currentData().toString() != "global") continue;
    const bool extrapolate = lastGlobal_[c].endsWith(":extrapolate");
    setGlobalText(c, globalCombo[c]->currentData().toString() + (extrapolate ? ":extrapolate" : ""));
  }
}

void ThmExperimentsPage::setGlobalText(int channel, const QString &value) {
  if (current_ < 0 || channel < 0 || channel > 1) return;
  ThmExperimentRecord &r = records_[current_];
  lastGlobal_[channel] = value;
  (channel == 0 ? r.opticalAA : r.opticalSF) = value;
  loadDistortion(r);
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::setOpticalText(int channel, const QString &tenNumbers) {
  if (current_ < 0 || channel < 0 || channel > 1) return;
  ThmExperimentRecord &r = records_[current_];
  lastOptical_[channel] = tenNumbers;
  (channel == 0 ? r.opticalAA : r.opticalSF) = tenNumbers;
  loadDistortion(r);
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::editOptical(int channel) {
  if (opticalCombo[channel]->currentData().toString() == "global") {
    if (current_ < 0) return;
    double elab[2] = {0.0, 0.0}, p[2][10], lo = 0.0, hi = 0.0;
    const bool have = globalEnds(records_.at(current_), channel, elab, p, &lo, &hi);
    ThmGlobalOpticalDialog dialog(channel == 0 ? tr("a + A optical potential (opticalAA)")
                                               : tr("s + F optical potential (opticalSF)"),
                                  lastGlobal_[channel], have, elab, p, lo, hi, this);
    if (dialog.exec() == QDialog::Accepted && dialog.text() != lastGlobal_[channel])
      setGlobalText(channel, dialog.text());
    return;
  }
  ThmOpticalDialog dialog(channel == 0 ? tr("a + A optical potential (opticalAA)")
                                       : tr("s + F optical potential (opticalSF)"),
                          lastOptical_[channel], this);
  if (dialog.exec() == QDialog::Accepted && dialog.text() != lastOptical_[channel])
    setOpticalText(channel, dialog.text());
}

void ThmExperimentsPage::chooseDistortionTable() {
  const QString start = projectDir_.isEmpty() ? QDir::currentPath() : projectDir_;
  const QString file = QFileDialog::getOpenFileName(this, tr("Distortion table"), start,
                                                    tr("Tables (*.dat *.txt);;All files (*)"));
  if (file.isEmpty() || current_ < 0) return;
  distortionTableEdit->setText(projectRelative(file, projectDir_));
  records_[current_].distortion = "table:" + distortionTableEdit->text();
  refreshRow(current_);
  showDerived(records_.at(current_));
}

QString ThmExperimentsPage::distortionInfo(const ThmExperimentRecord &x, QString *error, double *loOut,
                                           double *hiOut, double *rLo, double *rHi) const {
  if (error) error->clear();
  if (!x.hasDistortion()) return QString();
  // The engine's parse of the line (the keys' own checks are
  // ThmSettings::checkExperimentLines').
  std::vector<ThmExperiment> parsed;
  if (!ParseThmExperimentLine(x.line().toStdString(), parsed).empty() || parsed.empty()) return QString();
  ThmExperiment &e = parsed.front();
  QVector<double> energies;
  if (!pointEnergies(x.segments, energies)) return QString();  // the engine reports an unreadable file itself
  double lo = 1.0e300, hi = -1.0e300;
  for (double v : energies) {
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  if (loOut) *loOut = lo;
  if (hiOut) *hiOut = hi;
  auto number = [](double v) { return QString::number(v, 'g', 6); };
  const QString key = x.line() + "|" + projectDir_ + "|" + number(lo) + "|" + number(hi);
  if (key == cacheKey_) {
    if (error) *error = cacheError_;
    if (rLo) *rLo = cache_[0];
    if (rHi) *rHi = cache_[1];
    return cacheText_;
  }
  QString text, why;
  double r0 = 0.0, r1 = 0.0;
  if (e.distortion == ThmExperiment::DIST_TABLE) {
    // Config::ReadThmBlock reads the table relative to the .azr;
    // EData::BuildThmGroups wants every point inside it.
    QString path = QString::fromStdString(e.distortionTable);
    if (QFileInfo(path).isRelative() && !projectDir_.isEmpty()) path = QDir(projectDir_).filePath(path);
    ThmWeightTable table;
    table.name = e.distortionTable;
    const std::string bad = table.Read(QFile::encodeName(path).toStdString());
    if (!bad.empty()) {
      why = "distortion: " + QString::fromStdString(bad);
    } else if (!table.Covers(lo) || !table.Covers(hi)) {
      why = tr("distortion: the points span E_cm = %1 to %2 MeV, beyond the table '%3' [%4, %5] MeV.")
                .arg(number(lo), number(hi), QString::fromStdString(table.name), number(table.e.front()),
                     number(table.e.back()));
    } else {
      r0 = table(lo);
      r1 = table(hi);
      text = tr("w(E) from the table '%1': %2 at E = %3 MeV, %4 at E = %5 MeV")
                 .arg(QString::fromStdString(table.name), number(r0), number(lo), number(r1), number(hi));
    }
  } else {
    // As EData::BuildThmGroups sets it up.
    ThmDistortion d;
    if (!fillKinematics(x, d)) return QString();
    const ThmDistortion::Kinematics dk = d.kin;
    d.experiment = e.name;
    d.angleKind = e.angleKind == 1 ? ThmDistortion::LAB : e.angleKind == 2 ? ThmDistortion::CM : ThmDistortion::QF;
    d.angle = e.angle;
    d.sf.kind = e.distortion == ThmExperiment::DIST_OPTICAL && e.opticalSF.kind == 0 ? ThmDistortion::Channel::PLANE
                                                                                      : ThmDistortion::Channel::POINT_COULOMB;
    d.sf.mu = dk.ms * (dk.mx + dk.mA) / (dk.ms + dk.mx + dk.mA) * uconv;
    d.vcm = std::sqrt(2.0 * dk.mBeam * uconv * dk.beamEnergy) / ((dk.mBeam + dk.mTarget) * uconv);
    for (double v : energies) {
      const std::string bad = d.CheckEnergy(v);
      if (!bad.empty()) {
        why = "distortion: " + QString::fromStdString(bad) + ".";
        break;
      }
    }
    if (why.isEmpty() && e.angleWindow == 2) {
      // spectatorAngles=[cm:]table:<file>, relative to the .azr (Config::ReadThmBlock).
      QString path = QString::fromStdString(e.angleTable);
      if (QFileInfo(path).isRelative() && !projectDir_.isEmpty()) path = QDir(projectDir_).filePath(path);
      const std::string bad = ReadThmAngleTable(QFile::encodeName(path).toStdString(), e.angleTableT, e.angleTableW);
      if (!bad.empty())
        why = "spectatorAngles: " + QString::fromStdString(bad);
      else
        e.angleMin = e.angleTableT.front(), e.angleMax = e.angleTableT.back();
    }
    if (why.isEmpty()) {
      // The same radial grid and waves as the engine's (they depend on the
      // lower end of its ln R grid), without its grid of R: R is evaluated
      // directly at the ends, as AZURE2 prints it.
      d.dataLo = lo;
      d.dataHi = hi;
      const std::string bad = d.Build(e, dk, lo - 0.5, lo - 0.5, 0.5 * (lo + hi));
      if (!bad.empty()) why = "distortion: " + QString::fromStdString(bad) + ".";
    }
    if (why.isEmpty() && d.angWindow) {
      // The branches up to the highest point, and every point with an accepted
      // direction (EData::BuildThmGroups).
      d.SetAngleSlots(hi + 0.5 < d.eAA - d.kin.bind ? hi + 0.5 : hi);
      for (double v : energies) {
        const std::string bad = d.CheckWindow(v);
        if (!bad.empty()) {
          why = "distortion: " + QString::fromStdString(bad) + ".";
          break;
        }
      }
    }
    if (why.isEmpty()) {
      const ThmDistortion::Point a = d.Evaluate(lo), b = d.Evaluate(hi);
      r0 = a.ok ? d.R(a) : 0.0;
      r1 = b.ok ? d.R(b) : 0.0;
      auto point = [&](const ThmDistortion::Point &p) {
        return tr("(E_sF = %1, eta_sF = %2, theta_cm = %3 deg)")
            .arg(number(p.esf), number(p.etasf), number(p.thetaCm));
      };
      text = tr("Distortion factor R(E), zero-range DWBA: %1.\n").arg(QString::fromStdString(d.description)) +
             tr("k_aA = %1 fm^-1, eta_aA = %2, kappa = %3 fm^-1, eta_b = %4, beta = m_s/m_a = %5; E_ref = %6 MeV\n")
                 .arg(number(d.aa.k), number(d.aa.eta), number(d.kappa), number(d.etaB), number(d.beta),
                      number(d.eRef)) +
             tr("R = %1 at E = %2 MeV %3, %4 at E = %5 MeV %6.")
                 .arg(number(r0), number(lo), point(a), number(r1), number(hi), point(b));
    }
  }
  cacheKey_ = key;
  cacheText_ = text;
  cacheError_ = why;
  cache_[0] = r0;
  cache_[1] = r1;
  if (error) *error = why;
  if (rLo) *rLo = r0;
  if (rHi) *rHi = r1;
  return text;
}

bool ThmExperimentsPage::fillKinematics(const ThmExperimentRecord &x, ThmDistortion &d) const {
  Reaction r;
  QString ignored;
  if (!reaction(x, r, &ignored)) return false;
  ThmDistortion::Kinematics &dk = d.kin;
  dk.Za = r.horse.Z;
  dk.ZA = r.other.Z;
  dk.Zs = r.spectator.Z;
  dk.Zx = r.horse.Z - r.spectator.Z;
  dk.Aa = r.horse.A;
  dk.AA = r.other.A;
  dk.As = r.spectator.A;
  dk.Ax = r.horse.A - r.spectator.A;
  dk.ma = r.horse.mass;
  dk.mA = r.other.mass;
  dk.ms = r.spectator.mass;
  dk.mx = r.mX;
  dk.horseIsBeam = r.horseIsBeam;
  dk.mBeam = r.beam.mass;
  dk.mTarget = r.target.mass;
  dk.beamEnergy = r.beamEnergy;
  dk.bind = r.bind;
  d.eAA = dk.beamEnergy * dk.mTarget / (dk.mBeam + dk.mTarget);
  return true;
}

bool ThmExperimentsPage::globalEnds(const ThmExperimentRecord &x, int channel, double elab[2], double p[2][10],
                                    double *loOut, double *hiOut) const {
  std::vector<ThmExperiment> parsed;
  if (!ParseThmExperimentLine(x.line().toStdString(), parsed).empty() || parsed.empty()) return false;
  const ThmExperiment &e = parsed.front();
  if (e.distortion != ThmExperiment::DIST_OPTICAL) return false;
  QVector<double> energies;
  if (!pointEnergies(x.segments, energies) || energies.isEmpty()) return false;
  const double lo = *std::min_element(energies.begin(), energies.end());
  const double hi = *std::max_element(energies.begin(), energies.end());
  ThmDistortion d;
  if (!fillKinematics(x, d)) return false;
  d.dataLo = lo;
  d.dataHi = hi;
  d.Setup(e, d.kin, lo);  // a range refusal still leaves the channels set up
  if (!d.GlobalEnds(channel, lo, hi, elab, p)) return false;
  if (loOut) *loOut = lo;
  if (hiOut) *hiOut = hi;
  return true;
}

QString ThmExperimentsPage::globalSummary(const ThmExperimentRecord &x, int channel) const {
  double elab[2], p[2][10], lo = 0.0, hi = 0.0;
  const QString name = (channel == 0 ? x.opticalAA : x.opticalSF).section(':', 0, 0);
  if (!globalEnds(x, channel, elab, p, &lo, &hi)) return name;
  auto pair = [](double a, double b, const QString &unit) {
    return std::fabs(a - b) <= 5.0e-4 * std::max(std::fabs(a), 1.0)
               ? QString::number(a, 'g', 4) + unit
               : QString::fromUtf8("%1–%2%3").arg(QString::number(a, 'g', 4), QString::number(b, 'g', 4), unit);
  };
  return tr("%1 at E_lab = %2: V = %3, W = %4, WD = %5 (E = %6 to %7 MeV)")
      .arg(name, pair(elab[0], elab[1], " MeV"), pair(p[0][0], p[1][0], " MeV"), pair(p[0][3], p[1][3], " MeV"),
           pair(p[0][6], p[1][6], " MeV"), QString::number(lo, 'g', 4), QString::number(hi, 'g', 4));
}

QString ThmExperimentsPage::bindingMismatch(const ThmExperimentRecord &x, QString *detail) const {
  if (detail) detail->clear();
  Reaction r;
  QString ignored;
  if (!reaction(x, r, &ignored)) return QString();
  const double pairB = pairBinding_(r.pairKey);
  if (!(std::fabs(pairB - r.bind) > 1.0e-3)) return QString();
  if (detail)
    *detail = tr("B(x+s) from the masses of %1 = x + %2 is %3 MeV, but the entrance pair %4 carries B = %5 MeV "
                 "(field 32 of its channel lines); the THM vertex uses field 32, the kinematics of this experiment "
                 "(the quasi-free energy above, E_sF of the line shape and of the distortion factor) use the masses.")
                  .arg(QString::fromStdString(r.horse.name), QString::fromStdString(r.spectator.name))
                  .arg(QString::number(r.bind, 'g', 6))
                  .arg(r.pairKey)
                  .arg(QString::number(pairB, 'g', 6));
  return QString::fromUtf8("B from masses %1 MeV ≠ pair B %2 MeV (vertex uses the pair value)")
      .arg(QString::number(r.bind, 'g', 4), QString::number(pairB, 'g', 4));
}

QString ThmExperimentsPage::bindingMismatchOfPair(int pairKey, double pairB) const {
  QStringList out;
  for (const ThmExperimentRecord &x : records_) {
    Reaction r;
    QString ignored;
    if (!reaction(x, r, &ignored) || r.pairKey != pairKey) continue;
    if (std::fabs(pairB - r.bind) > 1.0e-3)
      out << QString::fromUtf8("experiment[%1]: B from masses %2 MeV ≠ pair B %3 MeV (vertex uses the pair value)")
                 .arg(x.name, QString::number(r.bind, 'g', 4), QString::number(pairB, 'g', 4));
  }
  return out.join('\n');
}

// ---------------------------------------------------------------------------
// Woods-Saxon dialog

ThmOpticalDialog::ThmOpticalDialog(const QString &title, const QString &tenNumbers, QWidget *parent) :
  QDialog(parent) {
  setWindowTitle(title);
  const QStringList values = tenNumbers.split(',');
  const char *names[10] = {"V", "R", "a", "W", "R<sub>W</sub>", "a<sub>W</sub>", "W<sub>D</sub>", "R<sub>D</sub>",
                           "a<sub>D</sub>", "R<sub>C</sub>"};
  const QString rows[4] = {tr("Real volume"), tr("Imaginary volume"), tr("Imaginary surface"), tr("Coulomb")};
  QGridLayout *g = new QGridLayout;
  g->setHorizontalSpacing(8);
  g->setVerticalSpacing(6);
  for (int k = 0; k < 10; k++) {
    const bool depth = k < 9 && k % 3 == 0;
    fields[k] = depth ? new ThmNumberSpin(" MeV", -1.0e4, 1.0e4, 1.0) : new ThmNumberSpin(" fm", 0.0, 100.0, 0.1);
    fields[k]->setWrittenText(values.value(k, "0"));
    const int row = k / 3, col = k % 3;
    if (col == 0) {
      QLabel *r = new QLabel(rows[row]);
      g->addWidget(r, row, 0);
    }
    QLabel *l = new QLabel(QString(names[k]) + ":");
    g->addWidget(l, row, 1 + 2 * col, Qt::AlignRight | Qt::AlignVCenter);
    g->addWidget(fields[k], row, 2 + 2 * col);
  }
  fields[0]->setToolTip(tr("Real depth (> 0 attractive); 0 switches the term off."));
  fields[3]->setToolTip(tr("Imaginary volume depth (> 0 absorptive); 0 switches the term off."));
  fields[6]->setToolTip(tr("Imaginary surface depth, 4 W_D e^x/(1+e^x)^2; 0 switches the term off."));
  fields[9]->setToolTip(tr("Uniform-sphere Coulomb radius; 0: a point charge."));
  QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  QVBoxLayout *l = new QVBoxLayout;
  l->addLayout(g);
  l->addWidget(buttons);
  setLayout(l);
}

ThmGlobalOpticalDialog::ThmGlobalOpticalDialog(const QString &title, const QString &value, bool have,
                                               const double elab[2], const double p[2][10], double lo, double hi,
                                               QWidget *parent) :
  QDialog(parent), name_(value.section(':', 0, 0)) {
  setWindowTitle(title);
  const int m = ThmGlobalOpticalIndex(name_.toStdString());
  QGridLayout *g = new QGridLayout;
  g->setHorizontalSpacing(12);
  g->setVerticalSpacing(4);
  int row = 0;
  if (m >= 0) {
    const ThmGlobalOptical &model = ThmGlobalOpticals()[m];
    QLabel *ref = new QLabel(QString::fromUtf8("<b>%1</b> — %2").arg(name_, model.reference));
    g->addWidget(ref, row++, 0, 1, 3);
    QLabel *range = new QLabel(QString::fromUtf8("Valid for %1 on A = %2–%3, E<sub>lab</sub> = %4–%5 MeV")
                                   .arg(model.projectiles)
                                   .arg(model.aMin)
                                   .arg(model.aMax)
                                   .arg(model.eMin)
                                   .arg(model.eMax));
    g->addWidget(range, row++, 0, 1, 3);
  }
  const char *names[10] = {"V", "R", "a", "W", "R<sub>W</sub>", "a<sub>W</sub>", "W<sub>D</sub>", "R<sub>D</sub>",
                           "a<sub>D</sub>", "R<sub>C</sub>"};
  for (int end = 0; end < 2; end++) {
    QLabel *h = new QLabel(have ? QString::fromUtf8("E = %1 MeV<br>E<sub>lab</sub> = %2 MeV")
                                      .arg(QString::number(end == 0 ? lo : hi, 'g', 4),
                                           QString::number(elab[end], 'g', 4))
                                : QString::fromUtf8("—"));
    h->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    g->addWidget(h, row, 1 + end);
  }
  row++;
  for (int k = 0; k < 10; k++) {
    const bool depth = k < 9 && k % 3 == 0;
    g->addWidget(new QLabel(QString(names[k]) + ":"), row + k, 0, Qt::AlignRight | Qt::AlignVCenter);
    for (int end = 0; end < 2; end++) {
      valueLabels[end][k] = new QLabel(have ? QString::number(p[end][k], 'f', depth ? 3 : 4) + (depth ? " MeV" : " fm")
                                            : QString::fromUtf8("—"));
      valueLabels[end][k]->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
      valueLabels[end][k]->setTextInteractionFlags(Qt::TextSelectableByMouse);
      g->addWidget(valueLabels[end][k], row + k, 1 + end);
    }
  }
  row += 10;
  extrapolateCheck = new QCheckBox(tr("Use outside the validity range (:extrapolate)"));
  extrapolateCheck->setChecked(value.endsWith(":extrapolate"));
  extrapolateCheck->setToolTip(tr("Without it AZURE2 refuses the potential where the target mass or the "
                                  "projectile's lab energy over the data lies outside the range above; with it, a "
                                  "WARNING."));
  g->addWidget(extrapolateCheck, row++, 0, 1, 3);
  QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  QVBoxLayout *l = new QVBoxLayout;
  l->addLayout(g);
  l->addWidget(buttons);
  setLayout(l);
}

QString ThmGlobalOpticalDialog::text() const {
  return name_ + (extrapolateCheck->isChecked() ? ":extrapolate" : "");
}

QString ThmOpticalDialog::text() const {
  QStringList out;
  for (ThmNumberSpin *f : fields) out << f->writtenText();
  return out.join(',');
}

// ---------------------------------------------------------------------------
// Coherent background (cbackground=)

namespace {
const char *const kCoherentForms[2] = {"const", "linear"};
}

QString ThmExperimentsPage::coherentCheck(const ThmExperimentRecord &x) const {
  std::vector<ThmExperiment::CoherentTerm> terms;
  const std::string why = ParseThmCoherentBackground(x.cbackground.toStdString(), terms);
  if (!why.empty()) return QString::fromStdString(why);
  // What EData::BuildThmGroups checks without the compound nucleus.
  if (entranceL_() == "coherent")
    return tr("cbackground= adds an amplitude per entrance bucket (s, l); entranceL=coherent merges the l of a "
              "channel spin, so it cannot be combined with it.");
  const QList<SegmentsDataData> lines = segments_->getLines();
  QSet<int> exits;
  for (int k : x.segments)
    if (k >= 1 && k <= lines.size()) exits << lines.at(k - 1).exitPairIndex;
  QStringList jpis = jpiChoices_();
  for (const ThmExperiment::CoherentTerm &t : terms) {
    const QString jpi = QString::fromStdString(ThmSpinText(t.J)) + (t.parity > 0 ? "+" : "-");
    if (!exits.contains(t.exitKey))
      return tr("cbackground %1:%2: no segment of the experiment has exit pair %2.").arg(jpi).arg(t.exitKey);
    if (!jpis.isEmpty() && !jpis.contains(jpi))
      return tr("cbackground %1:%2: the model has no J^pi = %1 group.").arg(jpi).arg(t.exitKey);
  }
  return QString();
}

void ThmExperimentsPage::addCoherentRow(const ThmExperiment::CoherentTerm &t, const QStringList &valueText) {
  const QSignalBlocker blocker(coherentTable);
  const int row = coherentTable->rowCount();
  coherentTable->insertRow(row);
  QComboBox *jpi = new QComboBox;
  jpi->setEditable(true);
  jpi->addItems(jpiChoices_());
  jpi->setEditText(QString::fromStdString(ThmSpinText(t.J)) + (t.parity > 0 ? "+" : "-"));
  QComboBox *exit = new QComboBox;
  exit->setEditable(true);
  QStringList keys;
  if (current_ >= 0) {
    const QList<SegmentsDataData> lines = segments_->getLines();
    for (int k : records_.at(current_).segments)
      if (k >= 1 && k <= lines.size() && !keys.contains(QString::number(lines.at(k - 1).exitPairIndex)))
        keys << QString::number(lines.at(k - 1).exitPairIndex);
  }
  exit->addItems(keys);
  exit->setEditText(QString::number(t.exitKey));
  QComboBox *form = new QComboBox;
  form->addItems(QStringList() << kCoherentForms[0] << kCoherentForms[1]);
  form->setCurrentIndex(t.form == 2 ? 1 : 0);
  coherentTable->setCellWidget(row, 0, jpi);
  coherentTable->setCellWidget(row, 1, exit);
  coherentTable->setCellWidget(row, 3, form);
  coherentTable->setItem(row, 2,
                         new QTableWidgetItem(t.hasChannels ? QString::fromStdString(ThmSpinText(t.s)) + "," +
                                                                  QString::number(t.l) + "," +
                                                                  QString::fromStdString(ThmSpinText(t.sp)) + "," +
                                                                  QString::number(t.lp)
                                                            : QString("all")));
  for (int k = 0; k < 4; k++) {
    QTableWidgetItem *v = new QTableWidgetItem(k < valueText.size() ? valueText.at(k) : QString("0"));
    v->setCheckState(t.fixed[k] ? Qt::Checked : Qt::Unchecked);
    coherentTable->setItem(row, 4 + k, v);
  }
  connect(jpi, &QComboBox::editTextChanged, this, [this]() { coherentEdited(); });
  connect(exit, &QComboBox::editTextChanged, this, [this]() { coherentEdited(); });
  connect(form, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this]() {
    showCoherentRows();
    coherentEdited();
  });
}

// c1 is shown for a linear term only.
void ThmExperimentsPage::showCoherentRows() {
  const QSignalBlocker blocker(coherentTable);
  for (int row = 0; row < coherentTable->rowCount(); row++) {
    QComboBox *form = qobject_cast<QComboBox *>(coherentTable->cellWidget(row, 3));
    const bool linear = form && form->currentIndex() == 1;
    for (int k = 6; k < 8; k++) {
      QTableWidgetItem *v = coherentTable->item(row, k);
      if (!v) continue;
      v->setFlags(linear ? (Qt::ItemIsEnabled | Qt::ItemIsEditable | Qt::ItemIsUserCheckable | Qt::ItemIsSelectable)
                         : Qt::ItemFlags());
    }
  }
  coherentTable->setVisible(coherentBox->isChecked());
  coherentAddButton->setVisible(coherentBox->isChecked());
  coherentRemoveButton->setVisible(coherentBox->isChecked());
  coherentRemoveButton->setEnabled(coherentTable->rowCount() > 0);
}

void ThmExperimentsPage::loadCoherent(const ThmExperimentRecord &r) {
  const QSignalBlocker blocker(coherentBox);
  coherentTable->setRowCount(0);
  coherentBox->setChecked(!r.cbackground.isEmpty());
  std::vector<ThmExperiment::CoherentTerm> terms;
  if (!r.cbackground.isEmpty() && ParseThmCoherentBackground(r.cbackground.toStdString(), terms).empty()) {
    // The values as written (the engine's parse gives the numbers).
    const QStringList items = r.cbackground.split(';');
    for (size_t i = 0; i < terms.size(); i++) {
      QStringList values;
      const QString item = i < (size_t)items.size() ? items.at(i) : QString();
      const int eq = item.indexOf('=');
      if (eq >= 0)
        for (QString v : item.mid(eq + 1).split(',')) values << (v.endsWith('f') ? v.left(v.size() - 1) : v);
      addCoherentRow(terms[i], values);
    }
  }
  showCoherentRows();
}

QString ThmExperimentsPage::coherentText() const {
  QStringList terms;
  for (int row = 0; row < coherentTable->rowCount(); row++) {
    QComboBox *jpi = qobject_cast<QComboBox *>(coherentTable->cellWidget(row, 0));
    QComboBox *exit = qobject_cast<QComboBox *>(coherentTable->cellWidget(row, 1));
    QComboBox *form = qobject_cast<QComboBox *>(coherentTable->cellWidget(row, 3));
    const QTableWidgetItem *channels = coherentTable->item(row, 2);
    if (!jpi || !exit || !form || !channels) continue;
    const bool linear = form->currentIndex() == 1;
    QString term = jpi->currentText().trimmed() + ":" + exit->currentText().trimmed();
    const QString ch = channels->text().trimmed();
    if (!ch.isEmpty() && ch != "all") term += ":" + ch;
    if (linear) term += ":linear";
    QStringList values;
    bool needed = false;
    for (int k = 0; k < (linear ? 4 : 2); k++) {
      const QTableWidgetItem *v = coherentTable->item(row, 4 + k);
      QString text = v ? v->text().trimmed() : QString("0");
      if (text.isEmpty()) text = "0";
      const bool fixed = v && v->checkState() == Qt::Checked;
      double x = 0.0;
      needed = needed || fixed || !readWholeDouble(text, x) || x != 0.0;
      values << text + (fixed ? "f" : "");
    }
    if (needed) term += "=" + values.join(',');
    terms << term;
  }
  return terms.join(';');
}

void ThmExperimentsPage::addCoherentTerm() {
  if (current_ < 0) return;
  ThmExperiment::CoherentTerm t;
  const QStringList jpis = jpiChoices_();
  if (!jpis.isEmpty()) {
    std::vector<ThmExperiment::CoherentTerm> parsed;
    if (ParseThmCoherentBackground((jpis.first() + ":1").toStdString(), parsed).empty()) t = parsed.front();
  }
  const QList<SegmentsDataData> lines = segments_->getLines();
  for (int k : records_.at(current_).segments)
    if (k >= 1 && k <= lines.size()) {
      t.exitKey = lines.at(k - 1).exitPairIndex;
      break;
    }
  {
    const QSignalBlocker blocker(coherentBox);
    coherentBox->setChecked(true);
  }
  addCoherentRow(t, QStringList());
  showCoherentRows();
  coherentEdited();
}

void ThmExperimentsPage::removeCoherentTerm() {
  int row = coherentTable->currentRow();
  if (row < 0) row = coherentTable->rowCount() - 1;
  if (row < 0) return;
  coherentTable->removeRow(row);
  if (coherentTable->rowCount() == 0) {
    const QSignalBlocker blocker(coherentBox);
    coherentBox->setChecked(false);
  }
  showCoherentRows();
  coherentEdited();
}

void ThmExperimentsPage::coherentToggled(bool on) {
  if (on && coherentTable->rowCount() == 0) {
    addCoherentTerm();  // a checked box always has a term to edit
    return;
  }
  if (!on) coherentTable->setRowCount(0);
  showCoherentRows();
  coherentEdited();
}

void ThmExperimentsPage::coherentEdited() {
  if (loading_ || current_ < 0) return;
  ThmExperimentRecord &r = records_[current_];
  r.cbackground = coherentBox->isChecked() ? coherentText() : QString();
  refreshRow(current_);
  showDerived(r);
}

void ThmExperimentsPage::setJpiChoices(std::function<QStringList()> choices) {
  jpiChoices_ = choices;
  // Rows already shown (the page loads its first experiment on construction).
  for (int row = 0; row < coherentTable->rowCount(); row++)
    if (QComboBox *jpi = qobject_cast<QComboBox *>(coherentTable->cellWidget(row, 0))) {
      const QSignalBlocker blocker(jpi);
      const QString text = jpi->currentText();
      jpi->clear();
      jpi->addItems(jpiChoices_());
      jpi->setEditText(text);
    }
}
