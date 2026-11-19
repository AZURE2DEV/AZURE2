#include "PlotTab.h"

#include <algorithm>
#include <cmath>

#include <QCheckBox>
#include <QColorDialog>
#include <QComboBox>
#include <QMenu>
#include <QFileDialog>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QRadioButton>
#include <QSet>
#include <QSplitter>
#include <QVBoxLayout>

#include "PlotWidget.h"

namespace {

const QColor kPalette[] = {QColor(0x1f, 0x77, 0xb4), QColor(0xd6, 0x27, 0x28), QColor(0x2c, 0xa0, 0x2c),
                           QColor(0xff, 0x7f, 0x0e), QColor(0x94, 0x67, 0xbd), QColor(0x8c, 0x56, 0x4b),
                           QColor(0xe3, 0x77, 0xc2), QColor(0x7f, 0x7f, 0x7f), QColor(0xbc, 0xbd, 0x22),
                           QColor(0x17, 0xbe, 0xcf)};

QListWidget *CurveList() {
  QListWidget *l = new QListWidget;
  l->setSelectionMode(QAbstractItemView::SingleSelection);
  return l;
}

/// The list's curves to draw: the checked ones.
QList<QListWidgetItem *> Checked(QListWidget *l) {
  QList<QListWidgetItem *> out;
  for (int i = 0; i < l->count(); i++)
    if (l->item(i)->checkState() == Qt::Checked) out << l->item(i);
  return out;
}

}  // namespace

PlotTab::PlotTab(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  dataList_ = CurveList();
  dataList_->setObjectName("plotDataList");
  extrapList_ = CurveList();
  extrapList_->setObjectName("plotExtrapList");
  connect(dataList_, &QListWidget::itemChanged, this, &PlotTab::draw);
  connect(extrapList_, &QListWidget::itemChanged, this, &PlotTab::draw);
  // a curve's colour and marker from the list's context menu
  for (QListWidget *l : {dataList_, extrapList_}) {
    l->setContextMenuPolicy(Qt::CustomContextMenu);
    l->setToolTip(tr("Check the curves to draw; right-click one for its colour and marker."));
    connect(l, &QListWidget::customContextMenuRequested, this, [this, l](const QPoint &at) {
      QListWidgetItem *item = l->itemAt(at);
      if (!item) return;
      const QString record = item->data(Qt::UserRole).toString().section('/', 1);
      QMenu menu;
      menu.addAction(tr("Colour..."), this, [this, record]() {
        const QColor c = QColorDialog::getColor(styles_.value(record, qMakePair(QColor(), 0)).first, this,
                                                tr("The colour of %1").arg(record));
        if (c.isValid()) setCurveStyle(record, c, styles_.value(record, qMakePair(QColor(), -1)).second);
      });
      QMenu *markers = menu.addMenu(tr("Marker"));
      const QStringList names = {tr("Circle"), tr("Square"), tr("Diamond"), tr("Triangle"), tr("Inverted triangle")};
      for (int m = 0; m < names.size(); m++)
        markers->addAction(names[m], this, [this, record, m]() {
          setCurveStyle(record, styles_.value(record, qMakePair(QColor(), 0)).first, m);
        });
      menu.addAction(tr("Default style"), this, [this, record]() { setCurveStyle(record, QColor(), -1); });
      menu.exec(l->viewport()->mapToGlobal(at));
    });
  }

  xAxis_ = new QComboBox;
  xAxis_->setObjectName("plotXAxis");
  xAxis_->addItem(tr("CoM Energy"), 0);
  xAxis_->addItem(tr("Excitation Energy"), 1);
  xAxis_->addItem(tr("CoM Angle"), 2);
  ySigma_ = new QRadioButton(tr("Cross Section"));
  ySigma_->setToolTip(tr("The calculated column of the output file: the cross section, or the data set's observable "
                         "(A_y, a phase shift) as written."));
  ySigma_->setChecked(true);
  yS_ = new QRadioButton(tr("S-Factor"));
  yS_->setToolTip(tr("The S-factor column, for a cross section; another observable stays as it is."));
  logX_ = new QCheckBox(tr("Log x"));
  logY_ = new QCheckBox(tr("Log y"));
  logY_->setToolTip(tr("Points that are not positive are left out on a log axis."));
  band_ = new QCheckBox(tr("Uncertainty band"));
  band_->setChecked(true);
  band_->setToolTip(tr("The 1-sigma band of the calculation where the run wrote one (the Calculate tab's option "
                       "--covariance-band writes the .band files)."));
  residuals_ = new QCheckBox(tr("Residuals"));
  residuals_->setChecked(true);
  legend_ = new QCheckBox(tr("Legend"));
  legend_->setChecked(true);
  for (QCheckBox *b : {logX_, logY_, band_, residuals_, legend_}) connect(b, &QCheckBox::toggled, this, &PlotTab::draw);
  connect(xAxis_, &QComboBox::currentIndexChanged, this, &PlotTab::draw);
  connect(ySigma_, &QRadioButton::toggled, this, &PlotTab::draw);

  plot_ = new PlotWidget;
  cursor_ = new QLabel;
  cursor_->setObjectName("plotCursor");
  cursor_->setProperty("role", "hint");
  cursor_->setMinimumWidth(fontMetrics().horizontalAdvance("x = -0.000000e-00, y = -0.000000e-00"));
  connect(plot_, &PlotWidget::cursorMoved, this, [this](double x, double y, bool inside) {
    cursor_->setText(inside ? QString("x = %1, y = %2").arg(x, 0, 'g', 6).arg(y, 0, 'g', 6) : QString());
  });
  notes_ = new QLabel;
  notes_->setObjectName("plotNotes");
  notes_->setWordWrap(true);

  QPushButton *refresh = new QPushButton(tr("&Refresh"));
  refresh->setToolTip(tr("Read the output files of the last run again."));
  connect(refresh, &QPushButton::clicked, this, &PlotTab::refresh);
  QPushButton *reset = new QPushButton(tr("Reset &Zoom"));
  connect(reset, &QPushButton::clicked, plot_, &PlotWidget::resetZoom);
  QPushButton *exportButton = new QPushButton(tr("&Export..."));
  connect(exportButton, &QPushButton::clicked, this, &PlotTab::exportDialog);

  auto listBox = [this](const QString &title, QListWidget *list, bool extrapolations) {
    QGroupBox *box = new QGroupBox(title);
    QVBoxLayout *v = new QVBoxLayout(box);
    v->addWidget(list);
    QPushButton *all = new QPushButton(tr("All"));
    QPushButton *none = new QPushButton(tr("None"));
    all->setObjectName(extrapolations ? "plotExtrapAll" : "plotDataAll");
    none->setObjectName(extrapolations ? "plotExtrapNone" : "plotDataNone");
    connect(all, &QPushButton::clicked, this, [this, extrapolations]() { checkAll(extrapolations, true); });
    connect(none, &QPushButton::clicked, this, [this, extrapolations]() { checkAll(extrapolations, false); });
    QHBoxLayout *h = new QHBoxLayout;
    h->addWidget(all);
    h->addWidget(none);
    h->addStretch(1);
    v->addLayout(h);
    return box;
  };
  QGroupBox *dataBox = listBox(tr("Segments From Data"), dataList_, false);
  QGroupBox *extrapBox = listBox(tr("Segments Without Data"), extrapList_, true);
  QWidget *left = new QWidget;
  QVBoxLayout *leftLayout = new QVBoxLayout(left);
  leftLayout->setContentsMargins(0, 0, 0, 0);
  leftLayout->addWidget(dataBox, 2);
  leftLayout->addWidget(extrapBox, 1);
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(refresh);
  buttons->addWidget(reset);
  buttons->addWidget(exportButton);
  leftLayout->addLayout(buttons);

  QHBoxLayout *axes = new QHBoxLayout;
  axes->addWidget(new QLabel(tr("X:")));
  axes->addWidget(xAxis_);
  axes->addSpacing(12);
  axes->addWidget(new QLabel(tr("Y:")));
  axes->addWidget(ySigma_);
  axes->addWidget(yS_);
  axes->addSpacing(12);
  for (QCheckBox *b : {logX_, logY_, band_, residuals_, legend_}) axes->addWidget(b);
  axes->addStretch(1);
  axes->addWidget(cursor_);
  QWidget *right = new QWidget;
  QVBoxLayout *rightLayout = new QVBoxLayout(right);
  rightLayout->setContentsMargins(0, 0, 0, 0);
  rightLayout->addLayout(axes);
  rightLayout->addWidget(plot_, 1);
  rightLayout->addWidget(notes_);

  QSplitter *split = new QSplitter;
  split->addWidget(left);
  split->addWidget(right);
  split->setStretchFactor(1, 3);
  (new QVBoxLayout(this))->addWidget(split);
  connect(holder_, &ProjectHolder::opened, this, &PlotTab::refresh);
}

void PlotTab::refresh() {
  QStringList selected;
  for (QListWidget *l : {dataList_, extrapList_})
    for (QListWidgetItem *i : Checked(l)) selected << i->data(Qt::UserRole).toString();
  curves_ = azure2::front::ReadCurves(holder_->project());
  for (QListWidget *l : {dataList_, extrapList_}) {
    l->blockSignals(true);
    l->clear();
  }
  for (int i = 0; i < (int)curves_.curves.size(); i++) {
    const azure2::front::Curve &c = curves_.curves[i];
    const QString id = QString::fromStdString(c.section + "/" + c.record);
    QListWidgetItem *item = new QListWidgetItem(QString("#%1: %2  %3 (%4)")
                                                    .arg(c.key)
                                                    .arg(QString::fromStdString(c.record),
                                                         QString::fromStdString(c.reaction),
                                                         QString::fromStdString(c.quantity)));
    item->setData(Qt::UserRole, id);
    item->setData(Qt::UserRole + 1, i);
    item->setToolTip(tr("%1, block %2; %3 points%4")
                         .arg(QString::fromStdString(c.file))
                         .arg(c.block + 1)
                         .arg(c.points.empty() ? c.coefficients.size() : c.points.size())
                         .arg(c.haveBand ? tr(", with its uncertainty band") : QString()));
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    item->setCheckState(selected.contains(id) ? Qt::Checked : Qt::Unchecked);
    (c.section == "data" ? dataList_ : extrapList_)->addItem(item);
  }
  for (QListWidget *l : {dataList_, extrapList_}) l->blockSignals(false);
  QStringList notes;
  for (const std::string &n : curves_.notes) notes << QString::fromStdString(n);
  if (curves_.labFrame) notes << tr("kinflag.dat: the files hold lab energies and angles.");
  notes_->setText(notes.join("\n"));
  notes_->setVisible(!notes.isEmpty());
  draw();
}

void PlotTab::selectRecords(const QStringList &names) {
  for (QListWidget *l : {dataList_, extrapList_}) {
    l->blockSignals(true);
    for (int i = 0; i < l->count(); i++) {
      const QString id = l->item(i)->data(Qt::UserRole).toString();
      l->item(i)->setCheckState(names.contains(id.section('/', 1)) ? Qt::Checked : Qt::Unchecked);
    }
    l->blockSignals(false);
  }
  draw();
}

QStringList PlotTab::availableRecords() const {
  QStringList out;
  for (QListWidget *l : {dataList_, extrapList_})
    for (int i = 0; i < l->count(); i++) out << l->item(i)->data(Qt::UserRole).toString().section('/', 1);
  return out;
}

void PlotTab::setCurveStyle(const QString &record, const QColor &color, int marker) {
  if (!color.isValid() && marker < 0) styles_.remove(record);
  else styles_[record] = qMakePair(color, marker);
  draw();
}

void PlotTab::checkAll(bool extrapolations, bool on) {
  QListWidget *l = extrapolations ? extrapList_ : dataList_;
  l->blockSignals(true);
  for (int i = 0; i < l->count(); i++) l->item(i)->setCheckState(on ? Qt::Checked : Qt::Unchecked);
  l->blockSignals(false);
  draw();
}

QStringList PlotTab::checkedRecords() const {
  QStringList out;
  for (QListWidget *l : {dataList_, extrapList_})
    for (QListWidgetItem *i : Checked(l)) out << i->data(Qt::UserRole).toString().section('/', 1);
  return out;
}

QStringList PlotTab::dataRecords() const {
  QStringList out;
  for (const azure2::front::Curve &c : curves_.curves)
    if (c.section == "data") out << QString::fromStdString(c.record);
  return out;
}

QString PlotTab::notes() const { return notes_->text(); }

void PlotTab::draw() {
  const int x = xAxis_->currentData().toInt();
  const bool s = yS_->isChecked();
  QList<PlotSeries> series;
  bool anyS = false, anyOther = false;
  QSet<QString> observables;
  int n = 0;
  for (QListWidget *l : {dataList_, extrapList_})
    for (QListWidgetItem *item : Checked(l)) {
      const int index = item->data(Qt::UserRole + 1).toInt();
      if (index < 0 || index >= (int)curves_.curves.size()) continue;
      const azure2::front::Curve &c = curves_.curves[index];
      PlotSeries ps;
      ps.label = QString::fromStdString(c.record);
      ps.color = kPalette[n % 10];
      ps.marker = n % 5;
      n++;
      if (styles_.contains(ps.label)) {  // the user's own choice
        if (styles_[ps.label].first.isValid()) ps.color = styles_[ps.label].first;
        if (styles_[ps.label].second >= 0) ps.marker = styles_[ps.label].second;
      }
      const bool useS = s && c.sFactor;
      if (!useS) observables.insert(QString::fromStdString(c.observable));
      anyS = anyS || useS;
      anyOther = anyOther || !useS;
      if (s && !c.sFactor) ps.label += tr(" (no S factor)");
      if (c.observable == "legendre") {
        // the first coefficient against energy; the others in the file
        for (const std::vector<double> &row : c.coefficients)
          if (row.size() > 1) {
            ps.modelX << row[0];
            ps.modelY << row[1];
          }
        ps.label += tr(" (a_0)");
        series << ps;
        continue;
      }
      for (const azure2::front::CurvePoint &p : c.points) {
        const double xv = x == 0 ? p.energy : x == 1 ? p.excitation : p.angle;
        const double model = useS ? p.modelS : p.model;
        ps.modelX << xv;
        ps.modelY << model;
        if (c.haveBand) {
          const double b = useS ? p.bandS : p.band;
          ps.bandLow << model - b;
          ps.bandHigh << model + b;
        }
        if (c.section == "data") {
          const double d = useS ? p.dataS : p.data, e = useS ? p.dataSError : p.dataError;
          ps.dataX << xv;
          ps.dataY << d;
          ps.dataError << e;
          ps.residualX << xv;
          ps.residual << (e > 0 ? (d - model) / e : NAN);
        }
      }
      // the curve is drawn in the order of x (a data set's rows follow its file)
      QVector<int> order(ps.modelX.size());
      for (int i = 0; i < order.size(); i++) order[i] = i;
      std::sort(order.begin(), order.end(), [&ps](int a, int b) { return ps.modelX[a] < ps.modelX[b]; });
      PlotSeries sorted = ps;
      for (int i = 0; i < order.size(); i++) {
        sorted.modelX[i] = ps.modelX[order[i]];
        sorted.modelY[i] = ps.modelY[order[i]];
        if (!ps.bandLow.isEmpty()) {
          sorted.bandLow[i] = ps.bandLow[order[i]];
          sorted.bandHigh[i] = ps.bandHigh[order[i]];
        }
      }
      series << sorted;
    }
  const QString xTitle = x == 0 ? tr("E (c.m.) [MeV]") : x == 1 ? tr("E_x [MeV]") : tr("theta (c.m.) [deg]");
  // the output files' columns are in the engine's units: sigma (b), dsigma/dOmega (b/sr), S (MeV b)
  QString yTitle = tr("calculated column (engine units)");
  if (anyS && !anyOther) yTitle = tr("S factor [MeV b]");
  else if (observables.size() == 1) {
    const QString o = *observables.begin();
    yTitle = o == "cross_section" || o == "total_capture" ? tr("cross section [b]")
             : o == "differential"                        ? tr("dsigma/dOmega (c.m.) [b/sr]")
             : o == "analyzing_power"                     ? tr("analyzing power")
             : o == "phase_shift"                         ? tr("phase shift")
                                                          : tr("Legendre coefficient a_0");
  }
  plot_->setAxisTitles(xTitle, yTitle);
  plot_->setShowBand(band_->isChecked());
  plot_->setShowResiduals(residuals_->isChecked());
  plot_->setShowLegend(legend_->isChecked());
  plot_->blockSignals(true);
  plot_->setLogX(logX_->isChecked());
  plot_->setLogY(logY_->isChecked());
  plot_->blockSignals(false);
  plot_->setSeries(series);
}

bool PlotTab::exportPlot(const QString &file, QString &why) const { return plot_->exportTo(file, why); }

void PlotTab::exportDialog() {
  const QString file = QFileDialog::getSaveFileName(this, tr("Export the plot"), holder_->directory(),
                                                    tr("PNG image (*.png);;SVG drawing (*.svg);;PDF (*.pdf);;Points as CSV (*.csv)"));
  if (file.isEmpty()) return;
  QString why;
  if (!exportPlot(file, why)) QMessageBox::warning(this, tr("Export"), why);
}
