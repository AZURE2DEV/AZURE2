// The Physics tab's shared pieces (PhysicsForms.h, docs/upgrade/qt_gui.md section 15).
#include "PhysicsForms.h"

#include <cmath>
#include <limits>

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>

namespace physicsqt {

namespace front = azure2::front;

QColor CurveColour(int n) {
  static const QColor kColours[] = {QColor(0x1f, 0x77, 0xb4), QColor(0xd6, 0x27, 0x28), QColor(0x2c, 0xa0, 0x2c),
                                    QColor(0xff, 0x7f, 0x0e), QColor(0x94, 0x67, 0xbd), QColor(0x8c, 0x56, 0x4b),
                                    QColor(0xe3, 0x77, 0xc2), QColor(0x7f, 0x7f, 0x7f), QColor(0xbc, 0xbd, 0x22),
                                    QColor(0x17, 0xbe, 0xcf)};
  return kColours[((n % 10) + 10) % 10];
}

QVector<double> Grid(double a, double b, int n, bool log) {
  QVector<double> out;
  if (n < 1) return out;
  if (n == 1) return {a};
  for (int i = 0; i < n; i++) {
    const double t = double(i) / (n - 1);
    out << (log ? std::pow(10.0, std::log10(a) + t * (std::log10(b) - std::log10(a))) : a + t * (b - a));
  }
  out.back() = b;
  return out;
}

QString Number(double v) { return QString::number(v, 'g', 8); }

// ---------------------------------------------------------------- SweepFields

SweepFields::SweepFields(const QString &unit, QWidget *parent) : QWidget(parent) {
  from = new QLineEdit;
  to = new QLineEdit;
  from->setObjectName("sweepFrom");
  to->setObjectName("sweepTo");
  for (QLineEdit *e : {from, to}) {
    // wide enough for a value as the pages write it ("14.769448", "-0.0012345"): it was cut to "769448" (qt_gui.md
    // section 19.8)
    e->setMinimumWidth(fontMetrics().horizontalAdvance("-00.0000000") + 16);
    e->setMaximumWidth(fontMetrics().horizontalAdvance("-0000.0000000") + 16);
    connect(e, &QLineEdit::editingFinished, this, &SweepFields::changed);
  }
  points = new QSpinBox;
  points->setObjectName("sweepPoints");
  points->setRange(2, 5000);
  points->setValue(300);
  points->setToolTip(tr("The number of points of the grid"));
  log = new QCheckBox(tr("log"));
  log->setObjectName("sweepLog");
  log->setToolTip(tr("Points evenly spaced in log10 (both ends > 0)"));
  unit_ = new QLabel(unit);
  connect(points, &QSpinBox::valueChanged, this, &SweepFields::changed);
  connect(log, &QCheckBox::toggled, this, &SweepFields::changed);
  // two rows: the range with its unit, then the points and the log switch (one row did not fit the pages' form
  // column once the fields show their whole values; qt_gui.md section 19.8)
  QVBoxLayout *v = new QVBoxLayout(this);
  v->setContentsMargins(0, 0, 0, 0);
  v->setSpacing(2);
  QHBoxLayout *h = new QHBoxLayout;
  h->addWidget(from);
  h->addWidget(new QLabel(tr("to")));
  h->addWidget(to);
  h->addWidget(unit_);
  h->addStretch(1);
  QHBoxLayout *h2 = new QHBoxLayout;
  h2->addWidget(points);
  h2->addWidget(new QLabel(tr("points")));
  h2->addWidget(log);
  h2->addStretch(1);
  v->addLayout(h);
  v->addLayout(h2);
}

void SweepFields::set(double a, double b, int n, bool lg) {
  const QSignalBlocker b1(from), b2(to), b3(points), b4(log);
  from->setText(Number(a));
  to->setText(Number(b));
  points->setValue(n);
  log->setChecked(lg);
}

void SweepFields::setUnit(const QString &unit) { unit_->setText(unit); }

bool SweepFields::grid(QVector<double> &out, QString &why) const {
  bool okA = false, okB = false;
  const double a = from->text().toDouble(&okA), b = to->text().toDouble(&okB);
  if (!okA || !okB || !std::isfinite(a) || !std::isfinite(b)) {
    why = tr("the range's ends must be numbers");
    return false;
  }
  if (log->isChecked() && (a <= 0 || b <= 0)) {
    why = tr("log spacing needs both ends > 0");
    return false;
  }
  out = Grid(a, b, points->value(), log->isChecked());
  return true;
}

// ---------------------------------------------------------------- ArgumentForm

ArgumentForm::ArgumentForm(QWidget *parent) : QWidget(parent) {
  QVBoxLayout *v = new QVBoxLayout(this);
  v->setContentsMargins(0, 0, 0, 0);
  QWidget *fieldsBox = new QWidget;
  fields_ = new QFormLayout(fieldsBox);
  fields_->setContentsMargins(0, 0, 0, 0);
  fields_->setFieldGrowthPolicy(QFormLayout::FieldsStayAtSizeHint);
  v->addWidget(fieldsBox);
  sweepChoice_ = new QComboBox;
  sweepChoice_->setObjectName("sweepArgument");
  sweepChoice_->setToolTip(tr("The argument on the horizontal axis"));
  sweep_ = new SweepFields(QString());
  QHBoxLayout *h = new QHBoxLayout;
  h->addWidget(new QLabel(tr("Sweep")));
  h->addWidget(sweepChoice_);
  h->addStretch(1);
  v->addLayout(h);
  v->addWidget(sweep_);
  connect(sweepChoice_, &QComboBox::currentIndexChanged, this, [this]() {
    sweptChanged();
    emit changed();
  });
  connect(sweep_, &SweepFields::changed, this, &ArgumentForm::changed);
}

void ArgumentForm::setFunction(const front::PhysicsFunction &function) {
  function_ = function;
  while (fields_->rowCount() > 0) fields_->removeRow(0);
  editors_.clear();
  const QSignalBlocker block(sweepChoice_);
  sweepChoice_->clear();
  for (int k = 0; k < (int)function.arguments.size(); k++) {
    const front::PhysicsArgument &a = function.arguments[k];
    const double start = a.hasDefault ? a.defaultValue : a.hasExample ? a.example : 0.0;
    QWidget *editor = nullptr;
    if (a.kind == "flag") {
      QCheckBox *c = new QCheckBox;
      c->setChecked(start != 0);
      connect(c, &QCheckBox::toggled, this, &ArgumentForm::changed);
      editor = c;
    } else if (a.kind == "choice") {
      QComboBox *c = new QComboBox;
      for (const std::string &choice : a.choices) c->addItem(QString::fromStdString(choice));
      c->setCurrentIndex((int)start);
      connect(c, &QComboBox::currentIndexChanged, this, &ArgumentForm::changed);
      editor = c;
    } else if (a.kind == "integer") {
      QSpinBox *s = new QSpinBox;
      const double lo = std::isfinite(a.lower) ? a.lower + (a.lowerOpen ? 1 : 0) : -1000000;
      const double hi = std::isfinite(a.upper) ? a.upper - (a.upperOpen ? 1 : 0) : 1000000;
      s->setRange((int)std::max(-1e9, lo), (int)std::min(1e9, hi));
      s->setValue((int)start);
      connect(s, &QSpinBox::valueChanged, this, &ArgumentForm::changed);
      editor = s;
    } else {
      QLineEdit *e = new QLineEdit(a.kind == "text" ? QString::fromStdString(a.exampleText) : Number(start));
      e->setMinimumWidth(fontMetrics().horizontalAdvance("0.000000000e-00"));
      connect(e, &QLineEdit::editingFinished, this, &ArgumentForm::changed);
      editor = e;
    }
    editor->setObjectName("arg_" + QString::fromStdString(a.name));
    editor->setToolTip(QString::fromStdString(a.doc + " -- domain: " + a.domain));
    QWidget *row = new QWidget;
    QHBoxLayout *h = new QHBoxLayout(row);
    h->setContentsMargins(0, 0, 0, 0);
    h->addWidget(editor);
    if (!a.unit.empty()) h->addWidget(new QLabel(QString::fromStdString(a.unit)));
    h->addStretch(1);
    QLabel *label = new QLabel(QString::fromStdString(a.name));
    label->setToolTip(editor->toolTip());
    fields_->addRow(label, row);
    editors_ << editor;
    if (a.hasSweep) sweepChoice_->addItem(QString::fromStdString(a.name), k);
  }
  const int axis = sweepChoice_->findText(QString::fromStdString(function.axis));
  sweepChoice_->setCurrentIndex(axis >= 0 ? axis : 0);
  sweepChoice_->setEnabled(sweepChoice_->count() > 1);
  sweptChanged();
}

void ArgumentForm::sweptChanged() {
  const int s = swept();
  for (int k = 0; k < editors_.size(); k++) editors_[k]->setEnabled(k != s);
  sweep_->setVisible(s >= 0);
  if (s < 0) return;
  const front::PhysicsArgument &a = function_.arguments[s];
  sweep_->setUnit(QString::fromStdString(a.unit));
  sweep_->set(a.sweepLower, a.sweepUpper, a.kind == "integer" ? std::max(2, (int)std::min(a.sweepUpper - a.sweepLower + 1, 300.0)) : 300,
              a.sweepLog);
}

int ArgumentForm::swept() const { return sweepChoice_->count() ? sweepChoice_->currentData().toInt() : -1; }

void ArgumentForm::setSwept(int k) {
  const int i = sweepChoice_->findData(k);
  if (i >= 0) sweepChoice_->setCurrentIndex(i);
}

double ArgumentForm::value(int k) const {
  QWidget *e = editors_.value(k);
  if (QCheckBox *c = qobject_cast<QCheckBox *>(e)) return c->isChecked() ? 1 : 0;
  if (QComboBox *c = qobject_cast<QComboBox *>(e)) return c->currentIndex();
  if (QSpinBox *s = qobject_cast<QSpinBox *>(e)) return s->value();
  if (QLineEdit *l = qobject_cast<QLineEdit *>(e)) {
    bool ok = false;
    const double v = l->text().toDouble(&ok);
    return ok ? v : std::numeric_limits<double>::quiet_NaN();
  }
  return std::numeric_limits<double>::quiet_NaN();
}

QString ArgumentForm::text(int k) const {
  QLineEdit *l = qobject_cast<QLineEdit *>(editors_.value(k));
  return l ? l->text().trimmed() : QString();
}

void ArgumentForm::setValue(int k, double v) {
  QWidget *e = editors_.value(k);
  if (QCheckBox *c = qobject_cast<QCheckBox *>(e)) c->setChecked(v != 0);
  else if (QComboBox *c = qobject_cast<QComboBox *>(e)) c->setCurrentIndex((int)v);
  else if (QSpinBox *s = qobject_cast<QSpinBox *>(e)) s->setValue((int)v);
  else if (QLineEdit *l = qobject_cast<QLineEdit *>(e)) l->setText(Number(v));
}

void ArgumentForm::setText(int k, const QString &text) {
  if (QLineEdit *l = qobject_cast<QLineEdit *>(editors_.value(k))) l->setText(text);
}

bool ArgumentForm::columns(std::vector<std::vector<double>> &out, QVector<double> &grid, QString &why,
                           bool skipText) const {
  out.clear();
  grid.clear();
  const int s = swept();
  if (s >= 0 && !sweep_->grid(grid, why)) return false;
  for (int k = 0; k < count(); k++) {
    const front::PhysicsArgument &a = function_.arguments[k];
    if (a.kind == "text") {
      if (skipText) continue;
      why = tr("%1 is a text").arg(QString::fromStdString(a.name));
      return false;
    }
    if (k == s) {
      out.emplace_back(grid.begin(), grid.end());
      continue;
    }
    const double v = value(k);
    if (!std::isfinite(v)) {
      why = tr("%1: not a number").arg(QString::fromStdString(a.name));
      return false;
    }
    out.push_back({v});
  }
  if (s < 0) grid = {0.0};
  return true;
}

QString ArgumentForm::fixedText() const { return fixedParts().join(", "); }

QStringList ArgumentForm::fixedParts() const {
  QStringList parts;
  for (int k = 0; k < count(); k++) {
    const front::PhysicsArgument &a = function_.arguments[k];
    if (k == swept()) continue;
    QString v;
    if (a.kind == "text") v = text(k);
    else if (a.kind == "choice") v = QString::fromStdString(a.choices.at(std::max(0, (int)value(k))));
    else if (a.kind == "flag") v = value(k) != 0 ? tr("on") : tr("off");
    else v = Number(value(k)) + (a.unit.empty() ? QString() : " " + QString::fromStdString(a.unit));
    parts << QString("%1 = %2").arg(QString::fromStdString(a.name), v);
  }
  return parts;
}

// ---------------------------------------------------------------- PlotPanel

PlotPanel::PlotPanel(const QString &name, QWidget *parent) : QWidget(parent) {
  setObjectName(name);
  plot_ = new PlotWidget;
  plot_->setObjectName(name + "Plot");
  plot_->setShowResiduals(false);
  logX_ = new QCheckBox(tr("Log x"));
  logY_ = new QCheckBox(tr("Log y"));
  legend_ = new QCheckBox(tr("Legend"));
  legend_->setChecked(true);
  logX_->setObjectName(name + "LogX");
  logY_->setObjectName(name + "LogY");
  logY_->setToolTip(tr("Values that are not positive are left out on a log axis."));
  connect(logX_, &QCheckBox::toggled, plot_, &PlotWidget::setLogX);
  connect(logY_, &QCheckBox::toggled, plot_, &PlotWidget::setLogY);
  connect(legend_, &QCheckBox::toggled, plot_, &PlotWidget::setShowLegend);
  QPushButton *reset = new QPushButton(tr("Reset zoom"));
  connect(reset, &QPushButton::clicked, plot_, &PlotWidget::resetZoom);
  QPushButton *exportButton = new QPushButton(tr("Export..."));
  exportButton->setToolTip(tr("The plot as PNG, SVG or PDF, or its points as CSV"));
  connect(exportButton, &QPushButton::clicked, this, [this]() {
    const QString file = QFileDialog::getSaveFileName(
        this, tr("Export the plot"), QString(), tr("PNG image (*.png);;SVG drawing (*.svg);;PDF (*.pdf);;Points as CSV (*.csv)"));
    if (file.isEmpty()) return;
    QString why;
    if (!plot_->exportTo(file, why)) QMessageBox::warning(this, tr("Export"), why);
  });
  cursor_ = new QLabel;
  cursor_->setObjectName(name + "Cursor");
  cursor_->setProperty("role", "hint");
  cursor_->setMinimumWidth(fontMetrics().horizontalAdvance("x = -0.000000e-00, y = -0.000000e-00"));
  connect(plot_, &PlotWidget::cursorMoved, this, [this](double x, double y, bool inside) {
    cursor_->setText(inside ? QString("x = %1, y = %2").arg(x, 0, 'g', 6).arg(y, 0, 'g', 6) : QString());
  });
  status_ = new QLabel;
  status_->setObjectName(name + "Status");
  status_->setWordWrap(true);
  status_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  citation_ = new QLabel;
  citation_->setObjectName(name + "Citation");
  citation_->setWordWrap(true);
  citation_->setProperty("role", "hint");
  citation_->hide();
  QHBoxLayout *strip = new QHBoxLayout;
  for (QWidget *w : std::initializer_list<QWidget *>{logX_, logY_, legend_, reset, exportButton}) strip->addWidget(w);
  strip->addStretch(1);
  strip->addWidget(cursor_);
  QVBoxLayout *v = new QVBoxLayout(this);
  v->setContentsMargins(0, 0, 0, 0);
  v->addLayout(strip);
  v->addWidget(plot_, 1);
  v->addWidget(status_);
  v->addWidget(citation_);
}

void PlotPanel::show(const QList<PlotSeries> &curves, const QString &xTitle, const QString &yTitle) {
  plot_->setAxisTitles(xTitle, yTitle);
  plot_->setSeries(curves);
}

void PlotPanel::setStatus(const QString &text) {
  status_->setText(text);
  status_->setVisible(!text.isEmpty());
}

QString PlotPanel::status() const { return status_->text(); }

void PlotPanel::setCitation(const QString &text) {
  citation_->setText(text);
  citation_->setVisible(!text.isEmpty());
}

QString PlotPanel::cursor() const { return cursor_->text(); }

}  // namespace physicsqt
