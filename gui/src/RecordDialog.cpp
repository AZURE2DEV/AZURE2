#include "RecordDialog.h"
#include "ExpertMode.h"

#include <QComboBox>
#include <QDialogButtonBox>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QRegularExpression>
#include <QScrollArea>
#include <QTextStream>
#include <QTimer>
#include <QVBoxLayout>

#include "ProjectHolder.h"
#include "ValueEditor.h"

namespace front = azure2::front;

namespace {

QString Q(const std::string &s) { return QString::fromStdString(s); }

QString Title(const front::FieldInfo &f) { return Q(f.title).replace("J^pi", QString::fromUtf8("Jπ")); }

QString Tip(const front::FieldInfo &f) {
  // the grammar in the expert view only (qt_gui.md section 19.9)
  QString t = "<b>" + Q(f.title).toHtmlEscaped() + "</b> &mdash; " + Q(f.doc).toHtmlEscaped() +
              (expert::Enabled() ? "<br><tt>" + Q(f.grammar).toHtmlEscaped() + "</tt>" : QString());
  if (!f.defaultText.empty()) t += "<br>Left out: <tt>" + Q(f.defaultText).toHtmlEscaped() + "</tt>.";
  return t;
}

/// A message of the validator without its head ("dataset 'x', field 'y': "): the form shows it beside the field.
QString Plain(const std::string &message) {
  const QString m = Q(message);
  const int field = m.indexOf("field '");
  const int at = field >= 0 ? m.indexOf("': ", field) : -1;
  return at > 0 ? m.mid(at + 3) : m;
}

/// The whole record as one text, to tell whether it changed between two presses of OK.
QString TextOf(const front::Record &r) {
  QString t = Q(r.kind) + " " + Q(r.name);
  for (const front::Field &f : r.fields) t += "\n" + Q(f.key) + " " + Q(f.value);
  return t;
}

}  // namespace

RecordDialog::RecordDialog(ProjectHolder *holder, const QString &section, const QStringList &kinds,
                           std::function<bool(const front::Record &, QString &)> writer, QWidget *parent)
    : QDialog(parent), holder_(holder), section_(section), writer_(std::move(writer)) {
  setObjectName("recordDialog");
  setWindowTitle(tr("Add a record"));
  QVBoxLayout *v = new QVBoxLayout(this);
  QFormLayout *top = new QFormLayout;
  kind_ = new QComboBox;
  kind_->setObjectName("recordKind");
  for (const QString &k : kinds) {
    const front::KindInfo *info = holder_->kind(k.toStdString());
    kind_->addItem(info ? Q(info->title) : k, k);
  }
  name_ = new QLineEdit;
  name_->setObjectName("recordName");
  top->addRow(tr("Kind"), kind_);
  top->addRow(tr("Name"), name_);
  kind_->setVisible(kinds.size() > 1);
  top->labelForField(kind_)->setVisible(kinds.size() > 1);
  v->addLayout(top);
  doc_ = new QLabel;
  doc_->setWordWrap(true);
  v->addWidget(doc_);
  QWidget *host = new QWidget;
  fields_ = new QFormLayout(host);
  fields_->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  QScrollArea *scroll = new QScrollArea;
  scroll->setWidgetResizable(true);
  scroll->setWidget(host);
  v->addWidget(scroll, 1);
  preview_ = new QPlainTextEdit;
  preview_->setObjectName("filePreview");
  preview_->setReadOnly(true);
  preview_->setLineWrapMode(QPlainTextEdit::NoWrap);
  preview_->setMaximumHeight(preview_->fontMetrics().height() * 7);
  fileNote_ = new QLabel;
  fileNote_->setWordWrap(true);
  fileNote_->setStyleSheet("color: palette(mid);");
  status_ = new QLabel;
  status_->setObjectName("recordDialogStatus");
  status_->setWordWrap(true);
  v->addWidget(status_);
  QDialogButtonBox *box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  another_ = box->addButton(tr("Add another"), QDialogButtonBox::ActionRole);
  another_->setObjectName("addAnother");
  another_->setToolTip(tr("Write this record and keep the dialog open, with the same choices, for the next one"));
  connect(box, &QDialogButtonBox::accepted, this, &RecordDialog::accept);
  connect(box, &QDialogButtonBox::rejected, this, &QDialog::reject);
  connect(another_, &QPushButton::clicked, this, &RecordDialog::addAnother);
  v->addWidget(box);
  connect(kind_, QOverload<int>::of(&QComboBox::activated), this, [this]() {
    name_->setText(uniqueName(kind_->currentData().toString()));
    build(QString(), {});
  });
  name_->setText(uniqueName(kind_->currentData().toString()));
  build(QString(), {});
  resize(860, 720);
}

QString RecordDialog::uniqueName(const QString &stem) const {
  QStringList names;
  for (const front::Record &r : holder_->records(section_)) names << Q(r.name);
  for (const QString &w : written_) names << w;
  QString base = stem, name = stem;
  // a name that ends in a number goes on counting ("s3" -> "s4")
  int i = 2;
  const int digits = [&]() {
    int n = 0;
    while (n < stem.size() && stem.at(stem.size() - 1 - n).isDigit()) n++;
    return n;
  }();
  if (digits > 0 && digits < stem.size()) {
    base = stem.left(stem.size() - digits);
    i = stem.right(digits).toInt() + 1;
    if (!names.contains(stem)) return stem;
    for (;; i++) {
      name = base + QString::number(i);
      if (!names.contains(name)) return name;
    }
  }
  for (; names.contains(name); i++) name = QString("%1_%2").arg(stem).arg(i);
  return name;
}

QString RecordDialog::quantityText() const {
  for (FieldEditor *e : editors_)
    if (Q(e->field().key) == quantityKey_) return e->text().trimmed();
  return QString();
}

FieldEditor *RecordDialog::fieldEditor(const QString &key) const {
  for (FieldEditor *e : editors_)
    if (Q(e->field().key) == key) return e;
  return nullptr;
}

void RecordDialog::build(const QString &chosenQuantity, const QMap<QString, QString> &keep) {
  rebuilding_ = true;
  // the preview's widgets are reused: take them out of the form before its rows go
  preview_->setParent(nullptr);
  fileNote_->setParent(nullptr);
  for (FieldEditor *e : editors_) delete e;
  editors_.clear();
  while (fields_->rowCount()) fields_->removeRow(0);
  const std::string k = kind_->currentData().toString().toStdString();
  const front::KindInfo *info = holder_->kind(k);
  doc_->setText(info ? Q(info->doc) : QString());
  quantityKey_.clear();
  fileKey_.clear();
  columnsKey_.clear();
  if (!info) {
    rebuilding_ = false;
    return;
  }
  for (const front::FieldInfo &f : info->fields) {
    if (IsQuantityField(f, k)) quantityKey_ = Q(f.key);
    if (f.value.type == "columns") columnsKey_ = Q(f.key);
    if (f.value.type == "path" && f.required && fileKey_.isEmpty()) fileKey_ = Q(f.key);
  }
  EditorContext context;
  context.holder = holder_;
  context.kind = k;
  context.quantity = [this]() { return quantityText(); };
  QString quantity = chosenQuantity;
  if (quantity.isEmpty() && editing_) quantity = Q(original_.Value(quantityKey_.toStdString()));
  std::vector<front::QuantityFieldRule> rules;
  const bool known = front::QuantityFieldRules(quantity.toStdString(), rules);
  auto given = [this](const std::string &key) {
    if (!editing_) return false;
    for (const front::Field &f : original_.fields)
      if (f.key == key) return true;
    return false;
  };
  struct Fresh {
    FieldEditor *editor;
    bool kept, required;
  };
  QList<Fresh> fresh;
  for (const front::FieldInfo &f : info->fields) {
    if (f.repeatable && !given(f.key)) continue;  // an epoch, an isotope: added in the form below the table
    if (f.repeatable) continue;                   // kept as written (record() takes them over)
    const front::QuantityFieldRule *rule = nullptr;
    for (const front::QuantityFieldRule &r : rules)
      if (r.field == f.key) rule = &r;
    bool wanted = (f.required && f.appliesWhen.empty()) || f.summary || given(f.key);
    if (f.byQuantity) {
      wanted = given(f.key) || (rule && (rule->required || f.summary)) || (!known && f.summary);
      if (known && !rule && !given(f.key)) wanted = false;
    }
    if (Q(f.key) == Q(info->includeField) && !given(f.key)) wanted = false;  // the table's "Use" box
    if (!wanted) continue;
    FieldEditor *e = new FieldEditor(f, context);
    e->setToolTip(Tip(f));
    editors_ << e;
    const bool required = (f.required && f.appliesWhen.empty()) || (rule && rule->required && !rule->conditional);
    fields_->addRow(Title(f) + (required ? " *" : ""), e);
    const QString key = Q(f.key);
    QString text;
    bool kept = false;
    if (key == quantityKey_ && !chosenQuantity.isEmpty()) {
      text = chosenQuantity;
      kept = true;
    } else if (keep.contains(key) && !keep.value(key).isEmpty() && !f.unitOfQuantity && !f.value.unitOfQuantity &&
               f.value.type != "columns") {
      text = keep.value(key);
      kept = true;
    } else if (editing_ && given(f.key) && chosenQuantity.isEmpty()) {
      text = Q(original_.Value(f.key));
      kept = true;
    }
    e->setText(text);
    fresh << Fresh{e, kept, required};
    if (key == fileKey_) {
      fields_->addRow(QString(), preview_);
      fields_->addRow(QString(), fileNote_);
    }
  }
  for (const Fresh &x : fresh) {
    // a required field starts from the project's first choice; an optional one is left out until given
    if (!x.kept && x.required && (!editing_ || !chosenQuantity.isEmpty())) {
      x.editor->prefill();
      x.editor->showMessage(QString());  // a value still to be typed (a grid's numbers) is not an error yet
    }
    const QString key = Q(x.editor->field().key);
    if (key == quantityKey_)
      connect(x.editor, &FieldEditor::changed, this, [this]() {
        if (!rebuilding_) QTimer::singleShot(0, this, [this]() { quantityChanged(); });
      });
    if (key == fileKey_) connect(x.editor, &FieldEditor::changed, this, [this]() { showFile(); });
  }
  rebuilding_ = false;
  showFile();
  // a quantity prefilled that requires fields: its rows too
  if (chosenQuantity.isEmpty() && !editing_ && !quantityKey_.isEmpty()) {
    std::vector<front::QuantityFieldRule> first;
    bool more = false;
    if (front::QuantityFieldRules(quantityText().toStdString(), first))
      for (const front::QuantityFieldRule &r : first)
        more = more || (r.required && !r.conditional && !fieldEditor(Q(r.field)));
    if (more) quantityChanged();
  }
}

void RecordDialog::quantityChanged() {
  QMap<QString, QString> texts;
  for (FieldEditor *e : editors_) texts[Q(e->field().key)] = e->text();
  build(quantityText(), texts);
}

void RecordDialog::showFile() {
  FieldEditor *file = fileKey_.isEmpty() ? nullptr : fieldEditor(fileKey_);
  preview_->setVisible(file != nullptr);
  fileNote_->setVisible(file != nullptr);
  if (!file) return;
  QString path = file->text().trimmed();
  if (path.startsWith('"') && path.endsWith('"') && path.size() > 1) path = path.mid(1, path.size() - 2);
  const QString full = QDir(holder_->directory()).absoluteFilePath(path);
  QFile f(full);
  if (path.isEmpty() || !f.open(QIODevice::ReadOnly | QIODevice::Text)) {
    preview_->setPlainText(QString());
    fileNote_->setText(path.isEmpty() ? tr("Choose the data file: its first lines are shown here.")
                                      : tr("%1 cannot be read (relative to the project's directory).").arg(path));
    return;
  }
  QTextStream in(&f);
  QStringList lines;
  int columns = -1, rows = 0;
  while (!in.atEnd()) {
    const QString line = in.readLine();
    if (lines.size() < 6) lines << line;
    const QStringList tokens = line.split(QRegularExpression("[\\s,]+"), Qt::SkipEmptyParts);
    bool numbers = !tokens.isEmpty();
    for (const QString &t : tokens) {
      bool ok = false;
      t.toDouble(&ok);
      numbers = numbers && ok;
    }
    if (!numbers) continue;
    rows++;
    if (columns < 0) columns = tokens.size();
  }
  preview_->setPlainText(lines.join("\n"));
  if (columns < 0) {
    fileNote_->setText(tr("No line of numbers in %1.").arg(path));
    return;
  }
  fileNote_->setText(tr("%n row(s) of numbers", nullptr, rows) + tr(", %n column(s) each", nullptr, columns));
  // the columns as many as the file has: the quantity's first ones, then columns the engine skips
  FieldEditor *cols = columnsKey_.isEmpty() ? nullptr : fieldEditor(columnsKey_);
  if (!cols) return;
  front::ValueParts parts;
  std::string why;
  const bool read = !cols->text().isEmpty() &&
                    front::ParseFieldValue(kind_->currentData().toString().toStdString(), columnsKey_.toStdString(),
                                           cols->text().toStdString(), quantityText().toStdString(), parts, why);
  if (read && (parts.fromFile || (int)parts.columns.size() == columns)) return;
  front::ValueParts made;
  made.type = cols->field().value.type;
  made.columns = front::ColumnsForFile(quantityText().toStdString(), columns);
  cols->setText(Q(front::ValueText(cols->field().value, made)));
  fileNote_->setText(fileNote_->text() + tr(": the columns below made %1; check their roles.").arg(columns));
}

front::Record RecordDialog::record() const {
  front::Record r;
  r.kind = kind_->currentData().toString().toStdString();
  r.name = name_->text().trimmed().toStdString();
  QMap<QString, QString> texts;
  for (FieldEditor *e : editors_) texts[Q(e->field().key)] = e->text().trimmed();
  QStringList done;
  if (editing_) {
    // the record's own order and comments; a field the dialog does not show (a repeatable one) as written
    for (const front::Field &f : original_.fields) {
      const QString key = Q(f.key);
      if (!texts.contains(key)) {
        r.fields.push_back(f);
        continue;
      }
      if (done.contains(key)) continue;
      done << key;
      if (!texts.value(key).isEmpty()) r.fields.push_back({f.key, texts.value(key).toStdString(), f.comment});
    }
  }
  for (FieldEditor *e : editors_) {
    const QString key = Q(e->field().key);
    if (done.contains(key)) continue;
    done << key;
    const QString t = texts.value(key);
    const bool required = e->field().required && e->field().appliesWhen.empty();
    if (!t.isEmpty()) r.fields.push_back({e->field().key, t.toStdString(), ""});
    else if (required && !editing_) r.fields.push_back({e->field().key, "?", ""});
  }
  return r;
}

bool RecordDialog::check() {
  const front::Record r = record();
  front::Project copy = front::Project::FromText(holder_->project().Text(), holder_->project().Path());
  std::string why;
  front::Record probe = r;
  if (editing_) probe.name = original_.name;  // the copy's record of that name is replaced
  for (FieldEditor *e : editors_) e->check();
  if (!copy.SetRecord(section_.toStdString(), probe, why)) {
    status_->setText("<span style=\"color:#c0392b\">" + Q(why).toHtmlEscaped() + "</span>");
    return false;
  }
  QStringList general;
  QMap<QString, QStringList> byField;
  for (const front::Diagnostic &d : copy.Validate().diagnostics)
    if (d.severity == "error" && d.record == probe.name) {
      if (fieldEditor(Q(d.field))) byField[Q(d.field)] << Plain(d.message);
      else general << (d.field.empty() ? QString() : Q(d.field) + ": ") + Plain(d.message);
    }
  for (auto it = byField.begin(); it != byField.end(); ++it) fieldEditor(it.key())->showMessage(it.value().join(" "));
  if (byField.isEmpty() && general.isEmpty()) {
    status_->clear();
    return true;
  }
  QStringList titles;
  for (auto it = byField.begin(); it != byField.end(); ++it) titles << Title(fieldEditor(it.key())->field());
  status_->setText("<span style=\"color:#c0392b\">" +
                   tr("The engine refuses the record%1%2. Press again to write it as it is (the diagnostics name it "
                      "until it is fixed).")
                       .arg(titles.isEmpty() ? QString() : tr(": see %1").arg(titles.join(", ")).toHtmlEscaped())
                       .arg(general.isEmpty() ? QString() : ("; " + general.join("; ")).toHtmlEscaped()) +
                   "</span>");
  return false;
}

bool RecordDialog::commit(QString &why) {
  const front::Record r = record();
  if (r.name.empty()) {
    why = tr("the record needs a name");
    status_->setText(why);
    return false;
  }
  if (!check() && TextOf(r) != lastChecked_) {
    lastChecked_ = TextOf(r);
    why = tr("the engine refuses the record");
    return false;
  }
  if (!writer_(r, why)) {
    status_->setText("<span style=\"color:#c0392b\">" + why.toHtmlEscaped() + "</span>");
    return false;
  }
  written_ << Q(r.name);
  lastChecked_.clear();
  return true;
}

void RecordDialog::accept() {
  if (editing_ && TextOf(record()) == TextOf(original_)) {
    QDialog::accept();  // nothing changed: not a byte moves
    return;
  }
  QString why;
  if (commit(why)) QDialog::accept();
}

void RecordDialog::addAnother() {
  QString why;
  if (!commit(why)) return;
  const QString was = written_.last();
  name_->setText(uniqueName(was));
  status_->setText(tr("%1 added. The next one, with the same choices: change what differs (the name, the data file, "
                      "...) and add it.")
                       .arg(was));
}

void RecordDialog::setRecord(const front::Record &record) {
  editing_ = true;
  original_ = record;
  setWindowTitle(tr("Edit %1").arg(Q(record.name)));
  another_->setVisible(false);
  const int at = kind_->findData(Q(record.kind));
  if (at >= 0) kind_->setCurrentIndex(at);
  kind_->setEnabled(false);
  name_->setText(Q(record.name));
  build(QString(), {});
}
