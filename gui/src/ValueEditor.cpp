#include "ValueEditor.h"
#include "ExpertMode.h"

#include <QCheckBox>
#include <QColor>
#include <QComboBox>
#include <QIcon>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRegularExpression>
#include <QRegularExpressionValidator>
#include <QSpinBox>
#include <QStackedWidget>
#include <QToolButton>
#include <QVBoxLayout>
#include <cmath>
#include <functional>

#include "ProjectHolder.h"

namespace front = azure2::front;

namespace {

QString Q(const std::string &s) { return QString::fromStdString(s); }
std::string S(const QString &s) { return s.trimmed().toStdString(); }

// The value types the editors tell apart (front::ValueInfo::type, the schema's ValueTypeName).
bool Is(const front::ValueInfo &info, const char *type) { return info.type == type; }

/// A line for a number as written: its text is kept (no reformatting); the validator lets through what strtod reads.
QLineEdit *NumberLine(const QString &placeholder, bool integer = false) {
  QLineEdit *e = new QLineEdit;
  e->setPlaceholderText(placeholder);
  const QRegularExpression rx(integer ? "[-+]?\\d*" : "[-+]?(\\d+\\.?\\d*|\\.\\d*)?([eE][-+]?\\d*)?");
  e->setValidator(new QRegularExpressionValidator(rx, e));
  e->setMinimumWidth(e->fontMetrics().horizontalAdvance("0.000000e-00"));
  e->setMaximumWidth(e->fontMetrics().horizontalAdvance("0.0000000000e-000"));
  return e;
}

/// Put `word` in a drop-down as its current item (added, marked, when the list has it not).
void SelectWord(QComboBox *combo, const QString &word) {
  int at = combo->findData(word);
  if (at < 0) {
    combo->addItem(word.isEmpty() ? QObject::tr("(none)") : word + QObject::tr("  (not in the list)"), word);
    at = combo->count() - 1;
  }
  combo->setCurrentIndex(at);
}

QString Word(const QComboBox *combo) { return combo->currentIndex() >= 0 ? combo->currentData().toString() : QString(); }

/// The empty choice of an optional value: its default greyed ("yes (default)"), or "(none)" without one.
void AddEmptyChoice(QComboBox *combo, const QString &defaultText) {
  if (defaultText.isEmpty()) {
    combo->addItem(QObject::tr("(none)"), QString());
    return;
  }
  combo->addItem(QObject::tr("%1 (default)").arg(defaultText), QString());
  combo->setItemData(combo->count() - 1, QColor(Qt::gray), Qt::ForegroundRole);
  combo->setItemData(combo->count() - 1, QObject::tr("Left out: the default, %1, applies").arg(defaultText), Qt::ToolTipRole);
}

/// A drop-down of words: `labels` shown, `words` written; the empty choice first when `optional` ("(none)", or the
/// default greyed).
void FillCombo(QComboBox *combo, const QStringList &words, const QStringList &labels, bool optional,
               const QString &defaultText = QString()) {
  const QString was = Word(combo);
  const bool hadChoice = combo->count() > 0;
  combo->blockSignals(true);
  combo->clear();
  if (optional) AddEmptyChoice(combo, defaultText);
  for (int i = 0; i < words.size(); i++) combo->addItem(i < labels.size() ? labels[i] : words[i], words[i]);
  if (hadChoice) SelectWord(combo, was);
  combo->blockSignals(false);
}

QStringList Words(const std::vector<std::string> &v) {
  QStringList out;
  for (const std::string &s : v) out << Q(s);
  return out;
}

QHBoxLayout *Row(QWidget *host) {
  QHBoxLayout *h = new QHBoxLayout(host);
  h->setContentsMargins(0, 0, 0, 0);
  return h;
}

/// The units a value of `info` takes (its dimension's, or those of the record's quantity).
QStringList UnitsOf(const front::ValueInfo &info, const EditorContext &c) {
  if (info.unitOfQuantity) return Words(front::QuantityUnits(c.quantity ? S(c.quantity()) : std::string()));
  return Words(info.units);
}


}  // namespace

QString RoundedNumbers(const QString &text) {
  const QStringList words = text.split(' ');
  QStringList out;
  out.reserve(words.size());
  for (const QString &w : words) {
    int digits = 0;
    bool started = false, plain = !w.isEmpty();
    for (const QChar ch : w) {
      if (ch == 'e' || ch == 'E') break;
      if (ch.isDigit()) {
        started = started || ch != '0';
        if (started) digits++;
      } else if (ch != '.' && ch != '-' && ch != '+') {
        plain = false;
        break;
      }
    }
    bool ok = false;
    const double v = plain && digits > 12 ? w.toDouble(&ok) : 0.0;
    out << (ok ? QString::number(v, 'g', 12) : w);
  }
  return out.join(' ');
}

namespace {

/// A line showing numbers: the text rounded for display (RoundedNumbers), the exact text kept and given back while
/// the line shows what it was given, so that a field the user did not touch writes the file's text (19.9).
void ShowNumbers(QLineEdit *e, const QString &exact) {
  const QString shown = RoundedNumbers(exact);
  e->setText(shown);
  e->setProperty("exactText", exact);
  e->setProperty("shownText", shown);
}

QString NumbersOf(const QLineEdit *e) {
  const QVariant shown = e->property("shownText");
  return shown.isValid() && e->text() == shown.toString() ? e->property("exactText").toString() : e->text();
}

}  // namespace

QString PlainHint(const std::string &grammar) {
  const QString g = Q(grammar).trimmed();
  if (g.isEmpty()) return g;
  if (g.contains("<A><Symbol>")) return QObject::tr("e.g. 20Ne");
  if (g.contains('|') || g.contains('[') || g.count('<') > 1) return QString();
  QString plain = g;
  plain.remove('<').remove('>');
  return plain;
}

namespace {

// ---------------------------------------------------------------- one line
class WordEditor : public PieceEditor {
 public:
  WordEditor(const front::ValueInfo &info, const std::string &defaultText) : info_(info) {
    const bool number = Is(info, "number") || Is(info, "integer");
    edit_ = number ? NumberLine(PlainHint(info.grammar), Is(info, "integer")) : new QLineEdit;
    if (!number) edit_->setPlaceholderText(PlainHint(info.grammar));
    if (!defaultText.empty()) edit_->setPlaceholderText(tr("default: %1").arg(Q(defaultText)));
    Row(this)->addWidget(edit_, number ? 0 : 1);
    if (number) static_cast<QHBoxLayout *>(layout())->addStretch(1);
    connect(edit_, &QLineEdit::textEdited, this, &PieceEditor::changed);
  }
  bool set(const front::ValueParts &p) override {
    QString text = Q(p.text);
    if (Is(info_, "number_list")) text = Words(p.numbers).join(' ');
    if (Is(info_, "composition")) {
      QStringList parts;
      for (const auto &c : p.composition) parts << Q(c.first) + " " + Q(c.second);
      text = parts.join(' ');
    }
    ShowNumbers(edit_, text);
    return true;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    const QString whole = NumbersOf(edit_);
    const QStringList tokens = whole.split(' ', Qt::SkipEmptyParts);
    if (Is(info_, "number_list")) {
      for (const QString &t : tokens) p.numbers.push_back(S(t));
    } else if (Is(info_, "composition")) {
      for (int i = 0; i + 1 < tokens.size(); i += 2) p.composition.push_back({S(tokens[i]), S(tokens[i + 1])});
      if (tokens.size() % 2) p.composition.push_back({S(tokens.last()), ""});
    } else {
      p.text = S(whole);
    }
    return p;
  }

 private:
  front::ValueInfo info_;
  QLineEdit *edit_;
};

// ---------------------------------------------------------------- a drop-down: enum, unit, reference, quantity
class ChoiceEditor : public PieceEditor {
 public:
  ChoiceEditor(const front::ValueInfo &info, const EditorContext &c, bool optional, bool quantity)
      : info_(info), context_(c), optional_(optional), quantity_(quantity) {
    combo_ = new QComboBox;
    combo_->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    QHBoxLayout *h = Row(this);
    h->addWidget(combo_);
    if (info.epoch) {
      epoch_ = new QComboBox;
      epoch_->setToolTip(tr("The epoch of the target (an areal density of its own)"));
      h->addWidget(new QLabel(tr("epoch")));
      h->addWidget(epoch_);
      connect(epoch_, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
      connect(combo_, QOverload<int>::of(&QComboBox::activated), this, [this]() { fillEpochs(QString()); });
    }
    h->addStretch(1);
    refresh();
    connect(combo_, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
  }
  void refresh() override {
    QStringList words, labels;
    if (quantity_) {
      for (const front::QuantityInfo &q : front::Quantities(context_.kind)) {
        words << Q(q.name);
        labels << Q(q.title) + "  (" + Q(q.name) + ")";
      }
    } else if (Is(info_, "unit")) {
      words = UnitsOf(info_, context_);
    } else if (Is(info_, "reference")) {
      words = RecordNamesOfKinds(context_.holder, info_.values);
    } else {
      words = Words(info_.values);
    }
    FillCombo(combo_, words, labels, optional_ || combo_->count() == 0, optional_ ? Q(context_.defaultText) : QString());
    if (epoch_) fillEpochs(Word(epoch_));
  }
  bool set(const front::ValueParts &p) override {
    refresh();
    combo_->blockSignals(true);
    SelectWord(combo_, Q(p.text));
    combo_->blockSignals(false);
    if (epoch_) fillEpochs(Q(p.epoch));
    return true;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    p.text = S(Word(combo_));
    if (epoch_) p.epoch = S(Word(epoch_));
    return p;
  }
  void prefill() override {
    refresh();
    for (int i = 0; i < combo_->count(); i++)
      if (!combo_->itemData(i).toString().isEmpty()) {
        combo_->setCurrentIndex(i);
        break;
      }
  }

 private:
  /// The epochs of the chosen target: the names of its repeatable structured fields whose first part is a name.
  void fillEpochs(const QString &keep) {
    QStringList names;
    const QString chosen = Word(combo_);
    for (const std::string &kind : info_.values) {
      const front::KindInfo *k = context_.holder->kind(kind);
      if (!k) continue;
      for (const std::string &section : k->sections)
        for (const front::Record &r : context_.holder->records(Q(section))) {
          if (r.kind != kind || Q(r.name) != chosen) continue;
          for (const front::FieldInfo &f : k->fields) {
            if (!f.repeatable || f.value.alternatives.empty() || f.value.alternatives[0].parts.empty() ||
                f.value.alternatives[0].parts[0].value.type != "name")
              continue;
            for (const front::Field &x : r.fields) {
              front::ValueParts parts;
              std::string why;
              if (x.key == f.key && front::ParseFieldValue(r.kind, x.key, x.value, "", parts, why) && !parts.parts.empty())
                names << Q(parts.parts[0].value.text);
            }
          }
        }
    }
    epoch_->blockSignals(true);
    epoch_->clear();
    epoch_->addItem(tr("(none)"), QString());
    for (const QString &n : names) epoch_->addItem(n, n);
    SelectWord(epoch_, keep);
    epoch_->setEnabled(names.size() > 0 || !keep.isEmpty());
    epoch_->blockSignals(false);
  }

  front::ValueInfo info_;
  EditorContext context_;
  bool optional_, quantity_;
  QComboBox *combo_;
  QComboBox *epoch_ = nullptr;
};

/// A unit drop-down of a value's units (kept as written when the list has it not).
QComboBox *UnitCombo(const QStringList &units) {
  QComboBox *c = new QComboBox;
  c->setSizeAdjustPolicy(QComboBox::AdjustToContents);
  FillCombo(c, units, {}, false);
  return c;
}

// ---------------------------------------------------------------- number(s) and a unit
class NumberUnitEditor : public PieceEditor {
 public:
  NumberUnitEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info), context_(c) {
    QHBoxLayout *h = Row(this);
    const bool range = Is(info, "range");
    lo_ = NumberLine(range ? tr("from") : tr("value"));
    h->addWidget(lo_);
    if (range) {
      h->addWidget(new QLabel(tr("to")));
      hi_ = NumberLine(tr("to"));
      h->addWidget(hi_);
      connect(hi_, &QLineEdit::textEdited, this, &PieceEditor::changed);
    }
    unit_ = UnitCombo(UnitsOf(info, c));
    h->addWidget(unit_);
    h->addStretch(1);
    connect(lo_, &QLineEdit::textEdited, this, &PieceEditor::changed);
    connect(unit_, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
  }
  bool set(const front::ValueParts &p) override {
    ShowNumbers(lo_, Q(hi_ ? p.lo : p.number));
    if (hi_) ShowNumbers(hi_, Q(p.hi));
    SelectWord(unit_, Q(p.unit));
    return true;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    (hi_ ? p.lo : p.number) = S(NumbersOf(lo_));
    if (hi_) p.hi = S(NumbersOf(hi_));
    p.unit = S(Word(unit_));
    return p;
  }
  void refresh() override { FillCombo(unit_, UnitsOf(info_, context_), {}, false); }
  void prefill() override {
    refresh();
    if (unit_->currentIndex() < 0 && unit_->count()) unit_->setCurrentIndex(0);
  }

 private:
  front::ValueInfo info_;
  EditorContext context_;
  QLineEdit *lo_, *hi_ = nullptr;
  QComboBox *unit_;
};

// ---------------------------------------------------------------- a nuisance: value, unit, Fit, prior, grid_max
class NuisanceEditor : public PieceEditor {
 public:
  NuisanceEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info), context_(c) {
    QVBoxLayout *v = new QVBoxLayout(this);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(2);
    QWidget *first = new QWidget;
    QHBoxLayout *h = Row(first);
    if (!info.refKinds.empty()) {
      source_ = new QComboBox;
      source_->setToolTip(tr("A value of this record, or a shared record of %1 that stands for it")
                              .arg(Words(info.refKinds).join(", ")));
      h->addWidget(source_);
      connect(source_, QOverload<int>::of(&QComboBox::activated), this, [this]() {
        showParts();
        emit changed();
      });
    }
    value_ = NumberLine(tr("value"));
    h->addWidget(value_);
    if (!info.units.empty()) {
      unit_ = UnitCombo(Words(info.units));
      h->addWidget(unit_);
      connect(unit_, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
    }
    fit_ = new QCheckBox(tr("Fit"));
    fit_->setToolTip(tr("Vary this value in a fit (else it is fixed)"));
    h->addWidget(fit_);
    if (info.profiled) {
      profiled_ = new QCheckBox(tr("Profiled"));
      profiled_->setToolTip(tr("Profiled out at every step of the fit (THM data sets of an experiment); no value"));
      h->addWidget(profiled_);
      connect(profiled_, &QCheckBox::toggled, this, [this]() {
        showParts();
        emit changed();
      });
    }
    priorButton_ = new QToolButton;
    priorButton_->setText(tr("Prior..."));
    priorButton_->setCheckable(true);
    priorButton_->setToolTip(tr("A Gaussian prior on the value, and the largest value a fold's grid must span"));
    h->addWidget(priorButton_);
    h->addStretch(1);
    v->addWidget(first);

    priorBox_ = new QWidget;
    QHBoxLayout *ph = Row(priorBox_);
    ph->setContentsMargins(16, 0, 0, 0);
    usePrior_ = new QCheckBox(tr("Gaussian prior: centre"));
    centre_ = NumberLine(tr("centre"));
    width_ = NumberLine(tr("width"));
    percent_ = new QCheckBox(tr("% of the centre"));
    gridMax_ = NumberLine(tr("none"));
    gridMax_->setToolTip(tr("grid_max: the largest value the sub-point grids of a fold must span (in the value's unit)"));
    ph->addWidget(usePrior_);
    ph->addWidget(centre_);
    ph->addWidget(new QLabel(tr("width")));
    ph->addWidget(width_);
    ph->addWidget(percent_);
    ph->addSpacing(12);
    ph->addWidget(new QLabel(tr("grid max")));
    ph->addWidget(gridMax_);
    ph->addStretch(1);
    v->addWidget(priorBox_);
    priorBox_->setVisible(false);
    connect(priorButton_, &QToolButton::toggled, priorBox_, &QWidget::setVisible);
    for (QLineEdit *e : {value_, centre_, width_, gridMax_}) connect(e, &QLineEdit::textEdited, this, &PieceEditor::changed);
    for (QCheckBox *b : {fit_, usePrior_, percent_}) connect(b, &QCheckBox::toggled, this, &PieceEditor::changed);
    connect(usePrior_, &QCheckBox::toggled, this, [this](bool on) {
      centre_->setEnabled(on);
      width_->setEnabled(on);
      percent_->setEnabled(on);
    });
    refresh();
  }
  void refresh() override {
    if (!source_) return;
    QStringList words = {QString()}, labels = {tr("value:")};
    for (const QString &n : RecordNamesOfKinds(context_.holder, info_.refKinds)) {
      words << n;
      labels << tr("shared: %1").arg(n);
    }
    FillCombo(source_, words, labels, false);
  }
  bool set(const front::ValueParts &p) override {
    if (source_) SelectWord(source_, p.reference ? Q(p.text) : QString());
    ShowNumbers(value_, Q(p.number));
    if (unit_) SelectWord(unit_, Q(p.unit));
    fit_->setChecked(p.fit);
    if (profiled_) profiled_->setChecked(p.profiled);
    usePrior_->setChecked(p.hasPrior);
    ShowNumbers(centre_, Q(p.priorCentre));
    ShowNumbers(width_, Q(p.priorWidth));
    percent_->setChecked(p.priorPercent);
    ShowNumbers(gridMax_, Q(p.gridMax));
    centre_->setEnabled(p.hasPrior);
    width_->setEnabled(p.hasPrior);
    percent_->setEnabled(p.hasPrior);
    priorButton_->setChecked(p.hasPrior || !p.gridMax.empty());
    showParts();
    return true;
  }
  void prefill() override {
    if (info_.units.empty() && value_->text().isEmpty()) value_->setText("1");
    if (unit_ && unit_->currentIndex() < 0 && unit_->count()) unit_->setCurrentIndex(0);
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    if (source_ && !Word(source_).isEmpty()) {
      p.reference = true;
      p.text = S(Word(source_));
      return p;
    }
    p.profiled = profiled_ && profiled_->isChecked();
    p.number = S(NumbersOf(value_));
    if (unit_) p.unit = S(Word(unit_));
    p.fit = fit_->isChecked();
    p.hasPrior = usePrior_->isChecked();
    p.priorCentre = S(NumbersOf(centre_));
    p.priorWidth = S(NumbersOf(width_));
    p.priorPercent = percent_->isChecked();
    p.gridMax = S(NumbersOf(gridMax_));
    return p;
  }

 private:
  void showParts() {
    const bool own = !source_ || Word(source_).isEmpty();
    const bool profiled = profiled_ && profiled_->isChecked();
    value_->setVisible(own && !profiled);
    if (unit_) unit_->setVisible(own && !profiled);
    fit_->setVisible(own && !profiled);
    if (profiled_) profiled_->setVisible(own);
    priorButton_->setVisible(own && !profiled);
    priorBox_->setVisible(own && !profiled && priorButton_->isChecked());
  }

  front::ValueInfo info_;
  EditorContext context_;
  QComboBox *source_ = nullptr;
  QLineEdit *value_, *centre_, *width_, *gridMax_;
  QComboBox *unit_ = nullptr;
  QCheckBox *fit_, *profiled_ = nullptr, *usePrior_, *percent_;
  QToolButton *priorButton_;
  QWidget *priorBox_;
};

// ---------------------------------------------------------------- a path with Browse
class PathEditor : public PieceEditor {
 public:
  PathEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info), context_(c) {
    QHBoxLayout *h = Row(this);
    edit_ = new QLineEdit;
    edit_->setPlaceholderText(tr("a file, relative to the project"));
    QPushButton *browse = new QPushButton(tr("Browse..."));
    h->addWidget(edit_, 1);
    h->addWidget(browse);
    connect(edit_, &QLineEdit::textEdited, this, &PieceEditor::changed);
    connect(browse, &QPushButton::clicked, this, [this]() {
      const QString dir = context_.holder->directory();
      const QString file = QFileDialog::getOpenFileName(this, tr("Choose a file"), dir);
      if (file.isEmpty()) return;
      edit_->setText(QDir(dir).relativeFilePath(file));
      emit changed();
    });
  }
  bool set(const front::ValueParts &p) override {
    edit_->setText(Q(p.text));
    quoted_ = p.quoted;
    return true;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    p.text = S(edit_->text());
    p.quoted = quoted_;
    return p;
  }

 private:
  front::ValueInfo info_;
  EditorContext context_;
  QLineEdit *edit_;
  bool quoted_ = false;
};

// ---------------------------------------------------------------- a formula: expression and parameters
class FormulaEditor : public PieceEditor {
 public:
  explicit FormulaEditor(const front::ValueInfo &info) : info_(info) {
    QGridLayout *g = new QGridLayout(this);
    g->setContentsMargins(0, 0, 0, 0);
    expression_ = new QLineEdit;
    expression_->setPlaceholderText(tr("expression in x, e.g. a0*x^a1"));
    g->addWidget(new QLabel(tr("f(x) =")), 0, 0);
    g->addWidget(expression_, 0, 1, 1, 3);
    params_ = new QWidget;
    paramsLayout_ = new QHBoxLayout(params_);
    paramsLayout_->setContentsMargins(0, 0, 0, 0);
    QToolButton *add = new QToolButton;
    add->setText("+");
    add->setToolTip(tr("Add a parameter of the formula"));
    g->addWidget(new QLabel(tr("with")), 1, 0);
    g->addWidget(params_, 1, 1);
    g->addWidget(add, 1, 2);
    g->setColumnStretch(3, 1);
    connect(expression_, &QLineEdit::textEdited, this, &PieceEditor::changed);
    connect(add, &QToolButton::clicked, this, [this]() {
      addParam(QString("a%1").arg(names_.size()), QString());
      emit changed();
    });
  }
  bool set(const front::ValueParts &p) override {
    expression_->setText(Q(p.expression));
    for (QWidget *w : rows_) delete w;
    rows_.clear();
    names_.clear();
    values_.clear();
    for (const auto &x : p.params) addParam(Q(x.first), Q(x.second));
    return true;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    p.expression = S(expression_->text());
    for (int i = 0; i < names_.size(); i++)
      if (!names_[i]->text().trimmed().isEmpty()) p.params.push_back({S(names_[i]->text()), S(values_[i]->text())});
    return p;
  }

 private:
  void addParam(const QString &name, const QString &value) {
    QWidget *row = new QWidget;
    QHBoxLayout *h = Row(row);
    QLineEdit *n = new QLineEdit(name);
    n->setMaximumWidth(n->fontMetrics().horizontalAdvance("abcdef"));
    QLineEdit *v = NumberLine(tr("value"));
    v->setText(value);
    h->addWidget(n);
    h->addWidget(new QLabel("="));
    h->addWidget(v);
    paramsLayout_->addWidget(row);
    rows_ << row;
    names_ << n;
    values_ << v;
    connect(n, &QLineEdit::textEdited, this, &PieceEditor::changed);
    connect(v, &QLineEdit::textEdited, this, &PieceEditor::changed);
  }
  front::ValueInfo info_;
  QLineEdit *expression_;
  QWidget *params_;
  QHBoxLayout *paramsLayout_;
  QList<QWidget *> rows_;
  QList<QLineEdit *> names_, values_;
};

// ---------------------------------------------------------------- a structured value: alternative and parts
class FormEditor : public PieceEditor {
 public:
  FormEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info) {
    QVBoxLayout *v = new QVBoxLayout(this);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(2);
    alternative_ = new QComboBox;
    alternative_->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    for (const front::AlternativeInfo &a : info.alternatives)
      alternative_->addItem(a.keyword.empty() ? Q(a.doc) : Q(a.keyword), Q(a.keyword));  // an unkeyed one by its doc
    for (int i = 0; i < alternative_->count(); i++) alternative_->setItemData(i, Q(info.alternatives[i].doc), Qt::ToolTipRole);
    doc_ = new QLabel;
    doc_->setWordWrap(true);
    doc_->setStyleSheet("color: palette(mid);");
    QWidget *top = new QWidget;
    QHBoxLayout *h = Row(top);
    h->addWidget(alternative_);
    h->addWidget(doc_, 1);
    v->addWidget(top);
    top->setVisible(info.alternatives.size() > 1 || !info.alternatives.empty());
    alternative_->setVisible(info.alternatives.size() > 1);
    pages_ = new QStackedWidget;
    for (const front::AlternativeInfo &a : info.alternatives) {
      QWidget *page = new QWidget;
      QFormLayout *f = new QFormLayout(page);
      f->setContentsMargins(12, 0, 0, 0);
      f->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
      QList<Part> parts;
      for (const front::PartInfo &p : a.parts) {
        Part part;
        part.info = p;
        EditorContext partContext = c;
        partContext.defaultText = p.optional ? p.defaultText : std::string();
        part.editor = MakePieceEditor(p.value, partContext, p.optional);
        part.editor->setToolTip(Q(p.doc));
        QString title = Q(p.name);
        title.replace('_', ' ');
        if (p.optional) {
          part.use = new QCheckBox(title);
          part.use->setToolTip(Q(p.doc) + (p.defaultText.empty() ? QString() : tr(" (left out: %1)").arg(Q(p.defaultText))));
          connect(part.use, &QCheckBox::toggled, part.editor, &QWidget::setEnabled);
          connect(part.use, &QCheckBox::toggled, this, &PieceEditor::changed);
          part.editor->setEnabled(false);
          f->addRow(part.use, part.editor);
        } else {
          QLabel *label = new QLabel(title);
          label->setToolTip(Q(p.doc));
          f->addRow(label, part.editor);
        }
        connect(part.editor, &PieceEditor::changed, this, &PieceEditor::changed);
        parts << part;
      }
      if (a.parts.empty()) f->addRow(new QLabel(tr("(nothing more to set)")));
      parts_ << parts;
      pages_->addWidget(page);
    }
    v->addWidget(pages_);
    connect(alternative_, QOverload<int>::of(&QComboBox::activated), this, [this](int i) {
      showAlternative(i);
      emit changed();
    });
    showAlternative(0);
  }
  bool set(const front::ValueParts &p) override {
    const int at = alternative_->findData(Q(p.alternative));
    if (at < 0) return false;
    alternative_->setCurrentIndex(at);
    showAlternative(at);
    for (Part &part : parts_[at]) {
      const front::ValueParts *given = p.Part(part.info.name);
      if (part.use) {
        part.use->blockSignals(true);
        part.use->setChecked(given != nullptr);
        part.use->blockSignals(false);
        part.editor->setEnabled(given != nullptr);
      }
      front::ValueParts empty;
      empty.type = part.info.value.type;
      if (!part.editor->set(given ? *given : empty)) return false;
    }
    return true;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    const int at = alternative_->currentIndex();
    if (at < 0) return p;
    p.alternative = S(alternative_->currentData().toString());
    for (const Part &part : parts_[at]) {
      if (part.use && !part.use->isChecked()) continue;
      p.parts.push_back({part.info.name, part.editor->get()});
    }
    return p;
  }
  void refresh() override {
    for (QList<Part> &parts : parts_)
      for (Part &part : parts) part.editor->refresh();
  }
  void prefill() override {
    const int at = alternative_->currentIndex();
    if (at >= 0)
      for (Part &part : parts_[at])
        if (!part.use) part.editor->prefill();
  }

 private:
  struct Part {
    front::PartInfo info;
    PieceEditor *editor = nullptr;
    QCheckBox *use = nullptr;
  };
  void showAlternative(int i) {
    if (i < 0 || i >= (int)info_.alternatives.size()) return;
    pages_->setCurrentIndex(i);
    // the alternative's doc beside the choice explains the choice; a value of one form only has no choice, and its doc
    // is the field's (its box or label says it): the tip only, not a grey line repeated on every entry of a repeatable
    // field (natural_neon's isotopes; qt_gui.md section 19.8)
    doc_->setText(Q(info_.alternatives[i].doc));
    // (an unkeyed alternative is listed in the choice by its doc already: no second copy beside it)
    const bool explains = info_.alternatives.size() > 1 && !info_.alternatives[i].keyword.empty();
    doc_->setVisible(explains);
    setToolTip(explains ? QString() : Q(info_.alternatives[i].doc));
    // the stack takes the height of the page shown, not of the tallest
    for (int k = 0; k < pages_->count(); k++)
      pages_->widget(k)->setSizePolicy(QSizePolicy::Preferred, k == i ? QSizePolicy::Preferred : QSizePolicy::Ignored);
    pages_->setVisible(!info_.alternatives[i].parts.empty());
  }
  front::ValueInfo info_;
  QComboBox *alternative_;
  QLabel *doc_;
  QStackedWidget *pages_;
  QList<QList<Part>> parts_;
};

// ---------------------------------------------------------------- a grid of energies or angles
class GridEditor : public PieceEditor {
 public:
  GridEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info), kinds_(front::GridKinds(info.none)) {
    QGridLayout *g = new QGridLayout(this);
    g->setContentsMargins(0, 0, 0, 0);
    kind_ = new QComboBox;
    for (const front::GridKind &k : kinds_) {
      kind_->addItem(Q(k.word), Q(k.word));
      kind_->setItemData(kind_->count() - 1, Q(k.doc), Qt::ToolTipRole);
    }
    lo_ = NumberLine(tr("from"));
    hi_ = NumberLine(tr("to"));
    list_ = new QLineEdit;
    list_->setPlaceholderText(tr("values separated by blanks"));
    path_ = new PathEditor(front::ValueInfo(), c);
    unit_ = UnitCombo(Words(info.units));
    frame_ = new QComboBox;
    for (const std::string &f : front::Frames()) frame_->addItem(Q(f), Q(f));
    spacing_ = new QComboBox;
    spacing_->addItem(tr("step"));
    spacing_->addItem(tr("number of points"));
    step_ = NumberLine(tr("step"));
    points_ = NumberLine(tr("points"), true);
    log_ = new QCheckBox(tr("evenly in log"));
    tolerance_ = NumberLine(tr("tolerance"));
    maxPoints_ = NumberLine(tr("max points"), true);
    startStep_ = NumberLine(tr("start step"));
    QWidget *span = new QWidget;
    QHBoxLayout *s = Row(span);
    s->addWidget(kind_);
    s->addWidget(lo_);
    s->addWidget(toLabel_ = new QLabel(tr("to")));
    s->addWidget(hi_);
    s->addWidget(list_, 1);
    s->addWidget(path_, 1);
    s->addWidget(unit_);
    s->addWidget(frame_);
    s->addStretch(0);
    QWidget *more = new QWidget;
    QHBoxLayout *m = Row(more);
    m->addWidget(spacing_);
    m->addWidget(step_);
    m->addWidget(points_);
    m->addWidget(log_);
    m->addWidget(tolLabel_ = new QLabel(tr("tolerance")));
    m->addWidget(tolerance_);
    m->addWidget(maxLabel_ = new QLabel(tr("max points")));
    m->addWidget(maxPoints_);
    m->addWidget(startLabel_ = new QLabel(tr("start step")));
    m->addWidget(startStep_);
    m->addStretch(1);
    g->addWidget(span, 0, 0);
    g->addWidget(more, 1, 0);
    for (QLineEdit *e : {lo_, hi_, list_, step_, points_, tolerance_, maxPoints_, startStep_})
      connect(e, &QLineEdit::textEdited, this, &PieceEditor::changed);
    for (QComboBox *b : {unit_, frame_}) connect(b, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
    connect(log_, &QCheckBox::toggled, this, &PieceEditor::changed);
    connect(path_, &PieceEditor::changed, this, &PieceEditor::changed);
    connect(kind_, QOverload<int>::of(&QComboBox::activated), this, [this]() {
      showKind();
      emit changed();
    });
    connect(spacing_, QOverload<int>::of(&QComboBox::activated), this, [this]() {
      showKind();
      emit changed();
    });
    showKind();
  }
  bool set(const front::ValueParts &p) override {
    const int at = kind_->findData(Q(p.gridKind));
    if (at < 0) return false;
    kind_->setCurrentIndex(at);
    ShowNumbers(lo_, Q(p.lo));
    ShowNumbers(hi_, Q(p.hi));
    ShowNumbers(list_, Words(p.numbers).join(' '));
    front::ValueParts path;
    path.text = p.text;
    path.quoted = p.quoted;
    path_->set(path);
    SelectWord(unit_, Q(p.unit));
    SelectWord(frame_, Q(p.frame));
    spacing_->setCurrentIndex(p.byPoints ? 1 : 0);
    ShowNumbers(step_, Q(p.step));
    points_->setText(Q(p.points));
    log_->setChecked(p.log);
    tolerance_->setText(Q(p.tolerance));
    maxPoints_->setText(Q(p.maxPoints));
    startStep_->setText(Q(p.startStep));
    showKind();
    return true;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    p.gridKind = S(Word(kind_));
    p.lo = S(NumbersOf(lo_));
    p.hi = S(NumbersOf(hi_));
    for (const QString &t : NumbersOf(list_).split(' ', Qt::SkipEmptyParts)) p.numbers.push_back(S(t));
    const front::ValueParts path = path_->get();
    p.text = path.text;
    p.quoted = path.quoted;
    p.unit = S(Word(unit_));
    p.frame = S(Word(frame_));
    p.byPoints = spacing_->currentIndex() == 1;
    p.step = S(NumbersOf(step_));
    p.points = S(points_->text());
    p.log = log_->isChecked();
    p.tolerance = S(tolerance_->text());
    p.maxPoints = S(maxPoints_->text());
    p.startStep = S(startStep_->text());
    return p;
  }

 private:
  void showKind() {
    const int i = kind_->currentIndex();
    if (i < 0) return;
    const front::GridKind &k = kinds_[i];
    for (QWidget *w : {(QWidget *)lo_, (QWidget *)toLabel_, (QWidget *)hi_}) w->setVisible(k.span);
    list_->setVisible(k.list);
    path_->setVisible(k.file);
    unit_->setVisible(k.unit);
    frame_->setVisible(k.unit);
    spacing_->setVisible(k.spacing);
    step_->setVisible(k.spacing && spacing_->currentIndex() == 0);
    points_->setVisible(k.spacing && spacing_->currentIndex() == 1);
    log_->setVisible(k.spacing && spacing_->currentIndex() == 1);
    for (QWidget *w : {(QWidget *)tolLabel_, (QWidget *)tolerance_, (QWidget *)maxLabel_, (QWidget *)maxPoints_,
                       (QWidget *)startLabel_, (QWidget *)startStep_})
      w->setVisible(k.adaptive);
  }
  front::ValueInfo info_;
  std::vector<front::GridKind> kinds_;
  QComboBox *kind_, *unit_, *frame_, *spacing_;
  QLineEdit *lo_, *hi_, *list_, *step_, *points_, *tolerance_, *maxPoints_, *startStep_;
  PathEditor *path_;
  QCheckBox *log_;
  QLabel *toLabel_, *tolLabel_, *maxLabel_, *startLabel_;
};

// ---------------------------------------------------------------- the columns of a data file
class ColumnsEditor : public PieceEditor {
 public:
  ColumnsEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info), context_(c), roles_(front::ColumnRoles()) {
    QVBoxLayout *v = new QVBoxLayout(this);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(2);
    fromFile_ = new QCheckBox(tr("Read the roles from the data file's header line"));
    columns_ = new QWidget;
    columnsLayout_ = new QHBoxLayout(columns_);
    columnsLayout_->setContentsMargins(0, 0, 0, 0);
    QToolButton *add = new QToolButton, *remove = new QToolButton;
    add->setText("+");
    add->setToolTip(tr("Add a column"));
    remove->setText("-");
    remove->setToolTip(tr("Remove the last column"));
    QWidget *row = new QWidget;
    QHBoxLayout *r = Row(row);
    r->addWidget(columns_);
    r->addWidget(add);
    r->addWidget(remove);
    r->addStretch(1);
    v->addWidget(row);
    v->addWidget(fromFile_);
    connect(fromFile_, &QCheckBox::toggled, this, [this](bool on) {
      columns_->setEnabled(!on);
      emit changed();
    });
    connect(add, &QToolButton::clicked, this, [this]() {
      addColumn(QString(), QString());
      emit changed();
    });
    connect(remove, &QToolButton::clicked, this, [this]() {
      if (columns.isEmpty()) return;
      delete columns.takeLast().host;
      emit changed();
    });
  }
  bool set(const front::ValueParts &p) override {
    for (const Column &c : columns) delete c.host;
    columns.clear();
    fromFile_->setChecked(p.fromFile);
    for (const auto &c : p.columns) addColumn(Q(c.first), Q(c.second));
    return true;
  }
  void prefill() override {
    if (!columns.isEmpty() || !context_.quantity) return;
    front::ValueParts p;
    for (const auto &c : front::DefaultColumns(S(context_.quantity()))) p.columns.push_back(c);
    set(p);
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    p.fromFile = fromFile_->isChecked();
    for (const Column &c : columns)
      p.columns.push_back({S(Word(c.role)), c.unit->isVisible() || !Word(c.unit).isEmpty() ? S(Word(c.unit)) : std::string()});
    return p;
  }

 private:
  struct Column {
    QWidget *host;
    QComboBox *role, *unit;
  };
  QString DimensionOf(const QString &role) const {
    for (const auto &r : roles_)
      if (Q(r.first) == role) return Q(r.second);
    return QString();
  }
  void addColumn(const QString &role, const QString &unit) {
    Column c;
    c.host = new QWidget;
    QVBoxLayout *v = new QVBoxLayout(c.host);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(1);
    c.role = new QComboBox;
    for (const auto &r : roles_) c.role->addItem(Q(r.first), Q(r.first));
    c.unit = new QComboBox;
    v->addWidget(new QLabel(tr("column %1").arg(columns.size() + 1)));
    v->addWidget(c.role);
    v->addWidget(c.unit);
    columnsLayout_->addWidget(c.host);
    SelectWord(c.role, role.isEmpty() ? c.role->itemData(0).toString() : role);
    auto units = [this, c](const QString &keep) {
      const QString dim = DimensionOf(Word(c.role));
      FillCombo(c.unit, dim.isEmpty() ? QStringList() : Words(front::DimensionUnits(S(dim))), {}, false);
      if (!dim.isEmpty()) SelectWord(c.unit, keep.isEmpty() ? c.unit->itemData(0).toString() : keep);
      c.unit->setVisible(!dim.isEmpty());
    };
    units(unit);
    connect(c.role, QOverload<int>::of(&QComboBox::activated), this, [this, units]() {
      units(QString());
      emit changed();
    });
    connect(c.unit, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
    columns << c;
  }
  front::ValueInfo info_;
  EditorContext context_;
  std::vector<std::pair<std::string, std::string>> roles_;
  QCheckBox *fromFile_;
  QWidget *columns_;
  QHBoxLayout *columnsLayout_;
  QList<Column> columns;
};

// ---------------------------------------------------------------- a reaction: entrance and exit pairs
/// The entrance and exit drop-downs of one reaction (the pairs of <levels>; the exits of the entrance pair's compound
/// nucleus, and its total capture).
class PairPicker : public QWidget {
 public:
  PairPicker(const EditorContext &c, std::function<void()> changed) : context_(c), changed_(std::move(changed)) {
    QHBoxLayout *h = Row(this);
    entrance_ = new QComboBox;
    exit_ = new QComboBox;
    entrance_->setToolTip(tr("The entrance pair (beam and target, or a beta decay), from the Particle Pairs tab"));
    exit_->setToolTip(tr("The exit pair: the reactions the entrance pair makes in the same compound nucleus"));
    entrance_->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    exit_->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    h->addWidget(new QLabel(tr("entrance")));
    h->addWidget(entrance_);
    h->addWidget(new QLabel(tr("exit")));
    h->addWidget(exit_);
    connect(entrance_, QOverload<int>::of(&QComboBox::activated), this, [this]() {
      fillExits(exit_->currentData().toInt());
      changed_();
    });
    connect(exit_, QOverload<int>::of(&QComboBox::activated), this, [this]() { changed_(); });
    refresh();
  }
  void refresh() {
    pairs_ = front::Pairs(context_.holder->project());
    const int keep = entrance_->currentData().toInt(), keepExit = exit_->currentData().toInt();
    entrance_->blockSignals(true);
    entrance_->clear();
    for (const front::PairChoice &p : pairs_)
      if (p.kind != 10) entrance_->addItem(QString("%1: %2").arg(p.key).arg(Q(p.label)), p.key);
    const int at = entrance_->findData(keep);
    entrance_->setCurrentIndex(at >= 0 ? at : (entrance_->count() ? 0 : -1));
    entrance_->blockSignals(false);
    fillExits(keepExit);
  }
  /// Show pair keys `in`, `out` (-1 total capture); false when the project has not them.
  bool setPairs(int in, int out) {
    const int at = entrance_->findData(in);
    if (in <= 0 || out == 0 || at < 0) return false;
    entrance_->setCurrentIndex(at);
    fillExits(out);
    return exit_->currentData().toInt() == out;
  }
  void prefill() {
    refresh();
    if (entrance_->currentIndex() < 0 && entrance_->count()) entrance_->setCurrentIndex(0);
    fillExits(0);
    if (exit_->count()) exit_->setCurrentIndex(0);
  }
  int entrance() const { return entrance_->currentData().toInt(); }
  int exit() const { return exit_->currentData().toInt(); }
  /// The label of the reaction chosen ("" when none is).
  QString label() const {
    const int in = entrance(), out = exit();
    return in > 0 && out != 0 ? Q(front::ReactionLabel(context_.holder->project(), in, out > 0 ? out : 0)) : QString();
  }

 private:
  void fillExits(int keep) {
    const int in = entrance_->currentData().toInt();
    std::string nucleus;
    for (const front::PairChoice &p : pairs_)
      if (p.key == in) nucleus = p.nucleus;
    exit_->blockSignals(true);
    exit_->clear();
    bool photon = false;
    for (const front::PairChoice &p : pairs_) {
      if (p.nucleus != nucleus) continue;
      photon = photon || p.kind == 10;
      exit_->addItem(QString("%1: %2").arg(p.key).arg(Q(front::ReactionLabel(context_.holder->project(), in, p.key))), p.key);
    }
    if (photon) exit_->addItem(tr("total capture: %1").arg(Q(front::ReactionLabel(context_.holder->project(), in, 0))), -1);
    const int at = exit_->findData(keep);
    exit_->setCurrentIndex(at >= 0 ? at : -1);
    exit_->blockSignals(false);
  }
  EditorContext context_;
  std::function<void()> changed_;
  std::vector<front::PairChoice> pairs_;
  QComboBox *entrance_, *exit_;
};

/// A reaction expression (section 3.5): one reaction by its pairs, or a sum or a ratio of terms, each a weight, a
/// reaction by its pairs (or as written, when the project has not its pairs) and an optional fixed angle.
class ReactionEditor : public PieceEditor {
 public:
  ReactionEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info), context_(c) {
    QVBoxLayout *v = new QVBoxLayout(this);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(2);
    QWidget *top = new QWidget;
    QHBoxLayout *h = Row(top);
    mode_ = new QComboBox;
    mode_->setObjectName("reactionMode");
    // the operations of the grammar (front::ReactionExpression::op), in words
    mode_->addItem(tr("one reaction"), "single");
    mode_->addItem(tr("sum of reactions"), "sum");
    mode_->addItem(tr("ratio of two"), "ratio");
    mode_->setToolTip(tr("A data set of one reaction, of a weighted sum of reactions (a total, unresolved states, a "
                         "natural target) or of the ratio of two"));
    single_ = new PairPicker(c, [this]() { emit changed(); });
    h->addWidget(mode_);
    h->addWidget(single_);
    h->addStretch(1);
    v->addWidget(top);
    termsHost_ = new QWidget;
    termsLayout_ = new QVBoxLayout(termsHost_);
    termsLayout_->setContentsMargins(16, 0, 0, 0);
    termsLayout_->setSpacing(2);
    v->addWidget(termsHost_);
    addTerm_ = new QPushButton(tr("Add term"));
    addTerm_->setObjectName("addReactionTerm");
    QHBoxLayout *more = new QHBoxLayout;
    more->setContentsMargins(16, 0, 0, 0);
    more->addWidget(addTerm_);
    more->addStretch(1);
    v->addLayout(more);
    connect(addTerm_, &QPushButton::clicked, this, [this]() {
      addTerm()->picker->prefill();
      showMode();
      emit changed();
    });
    connect(mode_, QOverload<int>::of(&QComboBox::activated), this, [this]() {
      // a sum or a ratio starts from the reaction shown, and a second term
      if (Word(mode_) != "single" && terms_.isEmpty()) {
        Term *first = addTerm();
        first->picker->refresh();
        first->picker->setPairs(single_->entrance(), single_->exit());
        addTerm()->picker->prefill();
      }
      while (Word(mode_) == "ratio" && terms_.size() > 2) removeTerm(terms_.last());
      while (Word(mode_) == "ratio" && terms_.size() < 2) addTerm()->picker->prefill();
      showMode();
      emit changed();
    });
    showMode();
  }
  ~ReactionEditor() override { qDeleteAll(terms_); }
  void refresh() override {
    single_->refresh();
    for (Term *t : terms_)
      if (t->picker) t->picker->refresh();
  }
  bool set(const front::ValueParts &p) override {
    clearTerms();
    front::ReactionExpression e;
    std::string why;
    if (!front::ParseReactionExpression(context_.holder->project(), S(Q(p.text)), e, why)) return false;
    const bool plain = e.op == "single" && e.terms.size() == 1 && e.terms[0].scale.empty() && e.terms[0].angle.empty();
    if (plain) {
      SelectWord(mode_, "single");
      single_->refresh();
      const bool ok = single_->setPairs(e.terms[0].entrance, e.terms[0].exit);
      showMode();
      return ok;
    }
    SelectWord(mode_, Q(e.op));
    for (const front::ReactionTermInfo &t : e.terms) {
      Term *row = addTerm();
      row->original = t;
      row->weight->setText(Q(t.scale));
      if (!row->picker->setPairs(t.entrance, t.exit)) {
        // the project has not its pairs (another compound nucleus): the reaction as written
        row->asText = true;
        row->picker->setVisible(false);
        row->text->setVisible(true);
        row->text->setText(Q(t.reaction));
      }
      row->useAngle->setChecked(!t.angle.empty());
      row->angle->setText(Q(t.angle));
      if (!t.angle.empty()) {
        SelectWord(row->angleUnit, Q(t.angleUnit));
        SelectWord(row->frame, Q(t.angleFrame));
      }
      row->showAngle();
    }
    showMode();
    return true;
  }
  void prefill() override {
    SelectWord(mode_, "single");
    single_->prefill();
    showMode();
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    if (Word(mode_) == "single") {
      p.text = S(single_->label());
      return p;
    }
    front::ReactionExpression e;
    e.op = S(Word(mode_));
    for (const Term *t : terms_) {
      front::ReactionTermInfo term;
      term.scale = S(t->weight->text());
      if (t->asText) {
        term.reaction = S(t->text->text());
      } else if (t->picker->entrance() == t->original.entrance && t->picker->exit() == t->original.exit &&
                 !t->original.reaction.empty()) {
        term.reaction = t->original.reaction;  // the pairs as read: the reaction as written
      } else {
        term.reaction = S(t->picker->label());
      }
      if (t->useAngle->isChecked()) {
        term.angle = S(t->angle->text());
        term.angleUnit = S(Word(t->angleUnit));
        term.angleFrame = S(Word(t->frame));
      }
      if (!term.reaction.empty()) e.terms.push_back(term);
    }
    p.text = front::ReactionExpressionText(e);
    return p;
  }

 private:
  struct Term {
    QWidget *host = nullptr;
    QLineEdit *weight = nullptr, *text = nullptr, *angle = nullptr;
    PairPicker *picker = nullptr;
    QCheckBox *useAngle = nullptr;
    QComboBox *angleUnit = nullptr, *frame = nullptr;
    front::ReactionTermInfo original;
    bool asText = false;  ///< the reaction as written (the project has not its pairs)
    void showAngle() {
      for (QWidget *w : {(QWidget *)angle, (QWidget *)angleUnit, (QWidget *)frame}) w->setVisible(useAngle->isChecked());
    }
  };
  Term *addTerm() {
    Term *t = new Term;
    t->host = new QWidget;
    t->host->setObjectName("reactionTerm");
    QHBoxLayout *h = Row(t->host);
    t->weight = NumberLine(tr("1"));
    t->weight->setToolTip(tr("The term's weight (left empty: 1)"));
    t->weight->setMaximumWidth(t->weight->fontMetrics().horizontalAdvance("0.000000"));
    t->picker = new PairPicker(context_, [this]() { emit changed(); });
    t->text = new QLineEdit;
    t->text->setVisible(false);
    t->text->setToolTip(tr("A reaction of pairs this project has not (written as text)"));
    t->useAngle = new QCheckBox(tr("at a fixed angle"));
    t->useAngle->setToolTip(tr("The term at a fixed angle (a ratio of differential cross sections)"));
    t->angle = NumberLine(tr("angle"));
    t->angle->setMaximumWidth(t->angle->fontMetrics().horizontalAdvance("000.0000"));
    t->angleUnit = UnitCombo(Words(front::DimensionUnits("angle")));
    t->frame = new QComboBox;
    for (const std::string &f : front::Frames()) t->frame->addItem(Q(f), Q(f));
    QToolButton *remove = new QToolButton;
    remove->setText(QString::fromUtf8("−"));
    remove->setToolTip(tr("Remove this term"));
    h->addWidget(t->weight);
    h->addWidget(new QLabel(QString::fromUtf8("×")));
    h->addWidget(t->picker);
    h->addWidget(t->text, 1);
    h->addWidget(t->useAngle);
    h->addWidget(t->angle);
    h->addWidget(t->angleUnit);
    h->addWidget(t->frame);
    h->addWidget(remove);
    h->addStretch(1);
    t->showAngle();
    termsLayout_->addWidget(t->host);
    terms_ << t;
    for (QLineEdit *e : {t->weight, t->text, t->angle}) connect(e, &QLineEdit::textEdited, this, &PieceEditor::changed);
    for (QComboBox *b : {t->angleUnit, t->frame}) connect(b, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
    connect(t->useAngle, &QCheckBox::toggled, this, [this, t]() {
      t->showAngle();
      emit changed();
    });
    connect(remove, &QToolButton::clicked, this, [this, t]() {
      if (terms_.size() <= 2) return;  // a sum or a ratio has two terms at least
      removeTerm(t);
      showMode();
      emit changed();
    });
    return t;
  }
  void removeTerm(Term *t) {
    terms_.removeAll(t);
    t->host->deleteLater();
    delete t;
  }
  void clearTerms() {
    while (!terms_.isEmpty()) removeTerm(terms_.last());
  }
  void showMode() {
    const QString m = Word(mode_);
    single_->setVisible(m == "single");
    termsHost_->setVisible(m != "single");
    addTerm_->setVisible(m == "sum");
  }
  front::ValueInfo info_;
  EditorContext context_;
  QComboBox *mode_;
  PairPicker *single_;
  QWidget *termsHost_;
  QVBoxLayout *termsLayout_;
  QPushButton *addTerm_;
  QList<Term *> terms_;
};

// ---------------------------------------------------------------- a spin, with or without its parity
/// J by a spin box of halves (empty below 0: an optional spin left out), and for a spin_parity a +/- toggle.
class SpinEditor : public PieceEditor {
 public:
  SpinEditor(const front::ValueInfo &info) : info_(info), parity_(Is(info, "spin_parity")) {
    QHBoxLayout *h = Row(this);
    j_ = new QDoubleSpinBox;
    j_->setDecimals(1);
    j_->setSingleStep(0.5);
    j_->setRange(-0.5, 40.0);
    j_->setSpecialValueText(tr("(none)"));  // -0.5: not given
    j_->setValue(-0.5);
    j_->setToolTip(tr("J, in halves"));
    h->addWidget(j_);
    if (parity_) {
      sign_ = new QToolButton;
      sign_->setCheckable(true);
      sign_->setToolTip(tr("The parity: click to switch between + and -"));
      showSign(true);
      connect(sign_, &QToolButton::toggled, this, [this](bool on) {
        showSign(on);
        emit changed();
      });
      h->addWidget(sign_);
    }
    h->addStretch(1);
    connect(j_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PieceEditor::changed);
  }
  bool set(const front::ValueParts &p) override {
    QString t = Q(p.text).trimmed();
    j_->blockSignals(true);
    if (sign_) sign_->blockSignals(true);
    bool ok = true;
    if (t.isEmpty()) {
      j_->setValue(-0.5);
    } else {
      if (sign_ && (t.endsWith('+') || t.endsWith('-'))) {
        showSign(t.endsWith('+'));
        t.chop(1);
      }
      double j = 0;
      if (t.endsWith("/2")) j = t.chopped(2).toInt(&ok) / 2.0;
      else j = t.toInt(&ok);
      j_->setValue(j);
    }
    j_->blockSignals(false);
    if (sign_) sign_->blockSignals(false);
    return ok;
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    if (j_->value() < 0) return p;
    const int twice = int(std::lround(2.0 * j_->value()));
    p.text = S(twice % 2 ? QString("%1/2").arg(twice) : QString::number(twice / 2));
    if (sign_) p.text += sign_->isChecked() ? "+" : "-";
    return p;
  }
  void prefill() override { j_->setValue(parity_ ? 0.5 : 0.0); }

 private:
  void showSign(bool plus) {
    sign_->setChecked(plus);
    sign_->setText(plus ? "+" : QString::fromUtf8("−"));
  }
  front::ValueInfo info_;
  bool parity_;
  QDoubleSpinBox *j_;
  QToolButton *sign_ = nullptr;
};

// ---------------------------------------------------------------- a pair of the project
/// A pair by its label ("p+15N", "g+16O[6.13]"): a drop-down of the project's pairs (front::Pairs), still editable
/// (a pair of a nucleus the project has not yet, a spelling of its own: kept as written until one is chosen).
class PairEditor : public PieceEditor {
 public:
  PairEditor(const front::ValueInfo &info, const EditorContext &c) : info_(info), context_(c) {
    QHBoxLayout *h = Row(this);
    combo_ = new QComboBox;
    combo_->setEditable(true);
    combo_->setInsertPolicy(QComboBox::NoInsert);
    combo_->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    combo_->setMinimumContentsLength(14);
    combo_->setToolTip(tr("A pair of the Particle Pairs tab (or its label typed)"));
    h->addWidget(combo_);
    h->addStretch(1);
    refresh();
    connect(combo_, QOverload<int>::of(&QComboBox::activated), this, &PieceEditor::changed);
    connect(combo_->lineEdit(), &QLineEdit::textEdited, this, &PieceEditor::changed);
  }
  void refresh() override {
    const QString keep = combo_->currentText();
    combo_->blockSignals(true);
    combo_->clear();
    if (context_.holder)
      for (const front::PairChoice &p : front::Pairs(context_.holder->project()))
        combo_->addItem(Q(p.label), Q(p.label));
    combo_->setEditText(keep);
    combo_->blockSignals(false);
  }
  bool set(const front::ValueParts &p) override {
    refresh();
    combo_->blockSignals(true);
    const int at = combo_->findData(Q(p.text));
    if (at >= 0) combo_->setCurrentIndex(at);
    combo_->setEditText(Q(p.text));
    combo_->blockSignals(false);
    return true;
  }
  void prefill() override {
    refresh();
    if (combo_->count()) {
      combo_->setCurrentIndex(0);
      combo_->setEditText(combo_->itemText(0));
    }
  }
  front::ValueParts get() const override {
    front::ValueParts p;
    p.type = info_.type;
    p.text = S(combo_->currentText());
    return p;
  }

 private:
  front::ValueInfo info_;
  EditorContext context_;
  QComboBox *combo_;
};

/// The editors that hold an empty state of their own: an optional field of these types gets a "use" switch.
bool NeedsUseSwitch(const front::ValueInfo &info) {
  return Is(info, "nuisance") || Is(info, "structured") || Is(info, "grid") || Is(info, "columns") || Is(info, "formula") ||
         Is(info, "range") || Is(info, "number_with_unit");
}

}  // namespace

PieceEditor *MakePieceEditor(const front::ValueInfo &info, const EditorContext &c, bool optional, QWidget *parent) {
  PieceEditor *e = nullptr;
  if (Is(info, "enum") || Is(info, "unit") || Is(info, "reference")) {
    front::FieldInfo probe;
    probe.value = info;
    probe.type = info.type;
    probe.values = info.values;
    e = new ChoiceEditor(info, c, optional, IsQuantityField(probe, c.kind));
  } else if (Is(info, "number_with_unit") || Is(info, "range")) {
    e = new NumberUnitEditor(info, c);
  } else if (Is(info, "nuisance")) {
    e = new NuisanceEditor(info, c);
  } else if (Is(info, "structured")) {
    e = new FormEditor(info, c);
  } else if (Is(info, "path")) {
    e = new PathEditor(info, c);
  } else if (Is(info, "formula")) {
    e = new FormulaEditor(info);
  } else if (Is(info, "grid")) {
    e = new GridEditor(info, c);
  } else if (Is(info, "columns")) {
    e = new ColumnsEditor(info, c);
  } else if (Is(info, "reaction_expr") && c.holder) {
    e = new ReactionEditor(info, c);
  } else if (Is(info, "spin") || Is(info, "spin_parity")) {
    e = new SpinEditor(info);
  } else if (Is(info, "pair") && c.holder) {
    e = new PairEditor(info, c);
  } else {
    e = new WordEditor(info, optional ? c.defaultText : std::string());
  }
  if (parent) e->setParent(parent);
  return e;
}

bool IsQuantityField(const front::FieldInfo &field, const std::string &kind) {
  if (field.value.type != "enum" || field.value.values.empty()) return false;
  const std::vector<std::string> names = front::QuantityNames(kind);
  return !names.empty() && field.value.values == names;
}

QStringList RecordNamesOfKinds(const ProjectHolder *holder, const std::vector<std::string> &kinds) {
  QStringList out;
  if (!holder) return out;
  for (const std::string &kind : kinds) {
    const front::KindInfo *k = holder->kind(kind);
    if (!k) continue;
    for (const std::string &section : k->sections)
      for (const front::Record &r : holder->records(Q(section)))
        if (r.kind == kind && !out.contains(Q(r.name))) out << Q(r.name);
  }
  return out;
}

QString OneBlank(const QString &text) {
  QString out;
  bool quote = false, blank = false;
  for (QChar c : text) {
    if (c == '"') quote = !quote;
    if (!quote && (c == ' ' || c == '\t')) {
      blank = true;
      continue;
    }
    if (blank && !out.isEmpty()) out += ' ';
    blank = false;
    out += c;
  }
  return out;
}

// ---------------------------------------------------------------- FieldEditor
FieldEditor::FieldEditor(const front::FieldInfo &field, const EditorContext &context, QWidget *parent)
    : QWidget(parent), field_(field), context_(context) {
  setObjectName("field_" + Q(field.key));
  QVBoxLayout *v = new QVBoxLayout(this);
  v->setContentsMargins(0, 0, 0, 0);
  v->setSpacing(1);
  QWidget *row = new QWidget;
  QHBoxLayout *h = Row(row);
  const bool optional = !field.required;
  if (optional && NeedsUseSwitch(field.value)) {
    use_ = new QCheckBox;
    use_->setToolTip(tr("Give this field (else it is left out%1)")
                         .arg(field.defaultText.empty() ? QString() : tr(": %1").arg(Q(field.defaultText))));
    h->addWidget(use_, 0, Qt::AlignTop);
    connect(use_, &QCheckBox::toggled, this, &FieldEditor::useToggled);
  }
  EditorContext pieceContext = context;
  pieceContext.defaultText = optional ? field.defaultText : std::string();
  piece_ = MakePieceEditor(field.value, pieceContext, optional);
  textEdit_ = new QLineEdit;
  textEdit_->setPlaceholderText(Q(field.grammar));
  textEdit_->setVisible(false);
  textToggle_ = new QToolButton;
  textToggle_->setObjectName("editAsText");
  textToggle_->setIcon(QIcon(":/actions/edit-text.svg"));
  textToggle_->setText(tr("Text"));
  textToggle_->setToolButtonStyle(Qt::ToolButtonIconOnly);
  textToggle_->setCheckable(true);
  textToggle_->setToolTip(expert::Enabled()
                              ? tr("<b>Edit this field as text</b>: the value as the project file writes it (its grammar: "
                                   "<tt>%1</tt>). Click again to return to the widgets.")
                                    .arg(Q(field.grammar).toHtmlEscaped())
                              : tr("<b>Edit this field as text</b>: the value as the project file writes it. Click again "
                                   "to return to the widgets."));
  textToggle_->setAutoRaise(true);
  textToggle_->setFixedSize(24, 22);  // one size on every row (qt_gui.md section 19.9)
  h->addWidget(piece_, 1);
  h->addWidget(textEdit_, 1);
  h->addWidget(textToggle_, 0, Qt::AlignTop);
  message_ = new QLabel;
  message_->setWordWrap(true);
  message_->setVisible(false);
  message_->setObjectName("message_" + Q(field.key));
  v->addWidget(row);
  v->addWidget(message_);
  connect(piece_, &PieceEditor::changed, this, &FieldEditor::pieceChanged);
  connect(textToggle_, &QToolButton::toggled, this, &FieldEditor::setTextMode);
  connect(textEdit_, &QLineEdit::textEdited, this, [this]() {
    edited_ = true;
    check();
    emit changed();
  });
}

void FieldEditor::setText(const QString &text) {
  loading_ = true;
  original_ = text;
  edited_ = false;
  textEdit_->setText(text);
  bool represented = true;
  front::ValueParts parts;
  parts.type = field_.value.type;
  if (!text.trimmed().isEmpty()) {
    std::string why;
    represented = front::ParseFieldValue(context_.kind, field_.key, S(text), context_.quantity ? S(context_.quantity()) : "",
                                         parts, why) &&
                  piece_->set(parts);
  } else {
    piece_->set(parts);
  }
  if (use_) {
    use_->blockSignals(true);
    use_->setChecked(!text.trimmed().isEmpty());
    use_->blockSignals(false);
    piece_->setEnabled(use_->isChecked());
  }
  textMode_ = !represented;
  textToggle_->blockSignals(true);
  textToggle_->setChecked(textMode_);
  textToggle_->blockSignals(false);
  piece_->setVisible(!textMode_);
  textEdit_->setVisible(textMode_);
  loading_ = false;
  check();
}

bool FieldEditor::isUsed() const { return !use_ || use_->isChecked(); }

QString FieldEditor::text() const {
  if (!isUsed()) return QString();
  if (textMode_) return edited_ ? textEdit_->text().trimmed() : original_;
  if (!edited_) return original_;
  return Q(front::ValueText(field_.value, piece_->get()));
}

void FieldEditor::pieceChanged() {
  if (loading_) return;
  edited_ = true;
  textEdit_->setText(text());
  check();
  emit changed();
}

void FieldEditor::useToggled(bool on) {
  piece_->setEnabled(on);
  if (loading_) return;
  edited_ = true;
  check();
  emit changed();
}

void FieldEditor::setTextMode(bool on) {
  if (on == textMode_) return;
  if (on) {
    textEdit_->setText(text());
  } else {
    front::ValueParts parts;
    std::string why;
    const QString t = textEdit_->text().trimmed();
    if (!t.isEmpty() && !(front::ParseFieldValue(context_.kind, field_.key, S(t), context_.quantity ? S(context_.quantity()) : "",
                                                   parts, why) &&
                          piece_->set(parts))) {
      showMessage(tr("Kept as text: %1").arg(why.empty() ? tr("the widgets cannot show this value") : Q(why)));
      textToggle_->blockSignals(true);
      textToggle_->setChecked(true);
      textToggle_->blockSignals(false);
      return;
    }
  }
  textMode_ = on;
  piece_->setVisible(!on);
  textEdit_->setVisible(on);
  textToggle_->blockSignals(true);
  textToggle_->setChecked(on);
  textToggle_->blockSignals(false);
}

void FieldEditor::typeText(const QString &text) {
  setTextMode(true);
  if (use_) use_->setChecked(true);
  textEdit_->setText(text);
  edited_ = true;
  check();
  emit changed();
}

void FieldEditor::setApplies(bool on) {
  if (applies_ == on) return;
  applies_ = on;
  check();
}

void FieldEditor::showMessage(const QString &message, bool error) {
  message_->setText(message);
  message_->setStyleSheet(error ? "color: #c0392b;" : "color: palette(mid);");
  message_->setVisible(!message.isEmpty());
}

QString FieldEditor::message() const { return message_->text(); }

bool FieldEditor::check() {
  const QString t = text();
  if (t.trimmed().isEmpty()) {
    const bool required = field_.required && applies_;  // required where it applies (qt_gui.md section 19.8)
    // the grammar in the expert view only (the label's tip has it there): a plain "Required" otherwise
    showMessage(!required ? QString() : expert::Enabled() ? tr("Required: %1").arg(Q(field_.grammar)) : tr("Required"),
                required);
    return !required;
  }
  front::ValueParts parts;
  std::string why;
  if (front::ParseFieldValue(context_.kind, field_.key, S(t), context_.quantity ? S(context_.quantity()) : "", parts, why)) {
    showMessage(QString());
    return true;
  }
  showMessage(Q(why));
  return false;
}

void FieldEditor::prefill() {
  loading_ = true;
  if (use_) use_->setChecked(true);
  piece_->setEnabled(true);
  piece_->prefill();
  loading_ = false;
  edited_ = true;
  textEdit_->setText(text());
  check();
}

void FieldEditor::setUsed(bool on) {
  if (use_) use_->setChecked(on);
}

void FieldEditor::refresh() {
  const bool was = loading_;
  loading_ = true;
  piece_->refresh();
  loading_ = was;
}
