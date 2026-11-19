#include "RecordsPanel.h"

#include <QCheckBox>
#include <QDir>
#include <QComboBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSplitter>
#include <QTableWidget>
#include <QToolButton>
#include <QVBoxLayout>

namespace front = azure2::front;

namespace {

QString Q(const std::string &s) { return QString::fromStdString(s); }

/// What a record is, in one line: the values of its first fields (a data set: reaction, quantity, file).
QString Summary(const front::Record &r) {
  QStringList parts;
  for (const front::Field &f : r.fields) {
    if (parts.size() == 3) break;
    if (f.key == "kind") continue;
    parts << Q(f.key) + " " + Q(f.value);
  }
  return parts.join("   ");
}

QString ToolTip(const front::FieldInfo &f) {
  QString t = "<b>" + Q(f.key).toHtmlEscaped() + "</b> &mdash; " + Q(f.doc).toHtmlEscaped() + "<br><tt>" +
              Q(f.grammar).toHtmlEscaped() + "</tt>";
  if (f.required) t += "<br>Required.";
  if (!f.requiredWhen.empty()) t += "<br>Required when " + Q(f.requiredWhen).toHtmlEscaped() + ".";
  if (!f.appliesWhen.empty()) t += "<br>Applies when " + Q(f.appliesWhen).toHtmlEscaped() + ".";
  if (f.byQuantity) t += "<br>The record's quantity decides whether it applies.";
  if (!f.defaultText.empty()) t += "<br>Default: <tt>" + Q(f.defaultText).toHtmlEscaped() + "</tt>.";
  return t;
}

}  // namespace

RecordsPanel::RecordsPanel(ProjectHolder *holder, const QString &section, const QString &help, QWidget *parent)
    : QWidget(parent), holder_(holder), section_(section) {
  setObjectName("records_" + section);
  list_ = new QTableWidget(0, 3);
  list_->setObjectName("recordList_" + section);
  list_->setHorizontalHeaderLabels({tr("Name"), tr("Kind"), tr("Record")});
  list_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
  list_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
  list_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
  list_->verticalHeader()->setVisible(false);
  list_->verticalHeader()->setDefaultSectionSize(fontMetrics().height() + 6);
  list_->setSelectionBehavior(QAbstractItemView::SelectRows);
  list_->setSelectionMode(QAbstractItemView::SingleSelection);
  list_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  list_->setShowGrid(false);
  connect(list_, &QTableWidget::itemSelectionChanged, this, &RecordsPanel::selectionChanged);

  addButton_ = new QPushButton(tr("+"));
  addButton_->setMaximumSize(28, 28);
  addButton_->setToolTip(tr("Add a record"));
  removeButton_ = new QPushButton(tr("-"));
  removeButton_->setMaximumSize(28, 28);
  removeButton_->setToolTip(tr("Remove the selected record"));
  duplicateButton_ = new QPushButton(tr("Duplicate"));
  duplicateButton_->setToolTip(tr("A copy of the selected record under a new name"));
  connect(addButton_, &QPushButton::clicked, this, &RecordsPanel::addRecord);
  connect(removeButton_, &QPushButton::clicked, this, &RecordsPanel::removeRecord);
  connect(duplicateButton_, &QPushButton::clicked, this, &RecordsPanel::duplicateRecord);
  QHBoxLayout *buttons = new QHBoxLayout;
  buttons->addWidget(addButton_);
  buttons->addWidget(removeButton_);
  buttons->addWidget(duplicateButton_);
  buttons->addStretch(1);

  QWidget *top = new QWidget;
  QVBoxLayout *topLayout = new QVBoxLayout(top);
  topLayout->setContentsMargins(0, 0, 0, 0);
  if (!help.isEmpty()) {
    QLabel *helpLabel = new QLabel(help);
    helpLabel->setWordWrap(true);
    topLayout->addWidget(helpLabel);
  }
  topLayout->addWidget(list_);
  topLayout->addLayout(buttons);

  kindLabel_ = new QLabel;
  kindLabel_->setWordWrap(true);
  everyField_ = new QCheckBox(tr("Every field"));
  everyField_->setToolTip(tr("Show the fields set rarely too (each is shown anyway when the record sets it)."));
  connect(everyField_, &QCheckBox::toggled, this, &RecordsPanel::showEveryField);
  applyButton_ = new QPushButton(tr("&Apply"));
  applyButton_->setObjectName("applyRecord_" + section);
  applyButton_->setToolTip(tr("Write the record into the project (the engine's writer: only the changed lines move)."));
  connect(applyButton_, &QPushButton::clicked, this, &RecordsPanel::apply);
  QPushButton *revert = new QPushButton(tr("Revert"));
  connect(revert, &QPushButton::clicked, this, &RecordsPanel::selectionChanged);
  QHBoxLayout *formButtons = new QHBoxLayout;
  formButtons->addWidget(everyField_);
  formButtons->addStretch(1);
  formButtons->addWidget(revert);
  formButtons->addWidget(applyButton_);

  formHost_ = new QWidget;
  QScrollArea *scroll = new QScrollArea;
  scroll->setWidgetResizable(true);
  scroll->setWidget(formHost_);
  QWidget *bottom = new QWidget;
  QVBoxLayout *bottomLayout = new QVBoxLayout(bottom);
  bottomLayout->setContentsMargins(0, 0, 0, 0);
  bottomLayout->addWidget(kindLabel_);
  bottomLayout->addWidget(scroll, 1);
  bottomLayout->addLayout(formButtons);

  QSplitter *splitter = new QSplitter(Qt::Vertical);
  splitter->addWidget(top);
  splitter->addWidget(bottom);
  splitter->setStretchFactor(0, 1);
  splitter->setStretchFactor(1, 2);
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  layout->addWidget(splitter);

  connect(holder_, &ProjectHolder::opened, this, &RecordsPanel::reload);
  connect(holder_, &ProjectHolder::sectionEdited, this, [this](const QString &s) {
    if (s == section_ && !reloading_) reload();
  });
  reload();
}

int RecordsPanel::recordCount() const { return list_->rowCount(); }

void RecordsPanel::reload() {
  const QString selected = list_->currentRow() >= 0 && list_->item(list_->currentRow(), 0)
                               ? list_->item(list_->currentRow(), 0)->text()
                               : QString();
  const std::vector<front::Record> records = holder_->records(section_);
  list_->blockSignals(true);
  list_->setRowCount((int)records.size());
  int keep = -1;
  for (int i = 0; i < (int)records.size(); i++) {
    const front::Record &r = records[i];
    list_->setItem(i, 0, new QTableWidgetItem(Q(r.name)));
    list_->setItem(i, 1, new QTableWidgetItem(Q(r.kind)));
    list_->setItem(i, 2, new QTableWidgetItem(Summary(r)));
    const QString active = Q(r.Value("active"));
    if (active == "no")
      for (int c = 0; c < 3; c++) list_->item(i, c)->setForeground(palette().color(QPalette::Disabled, QPalette::Text));
    if (Q(r.name) == selected) keep = i;
  }
  list_->blockSignals(false);
  if (keep < 0 && !records.empty()) keep = 0;
  if (keep >= 0)
    list_->selectRow(keep);
  else
    clearForm();
  selectionChanged();
  addButton_->setEnabled(!holder_->kindsOfSection(section_).isEmpty());
}

bool RecordsPanel::selectRecord(const QString &name) {
  for (int i = 0; i < list_->rowCount(); i++)
    if (list_->item(i, 0)->text() == name) {
      list_->selectRow(i);
      return true;
    }
  return false;
}

void RecordsPanel::selectionChanged() {
  const int row = list_->currentRow();
  removeButton_->setEnabled(row >= 0);
  duplicateButton_->setEnabled(row >= 0);
  applyButton_->setEnabled(row >= 0);
  if (row < 0 || !list_->item(row, 0)) {
    clearForm();
    return;
  }
  const QString name = list_->item(row, 0)->text();
  for (const front::Record &r : holder_->records(section_))
    if (Q(r.name) == name) {
      buildForm(r);
      return;
    }
  clearForm();
}

void RecordsPanel::clearForm() {
  widgets_.clear();
  nameEdit_ = nullptr;
  delete formHost_->layout();
  for (QWidget *w : formHost_->findChildren<QWidget *>(QString(), Qt::FindDirectChildrenOnly)) delete w;
  form_ = nullptr;
  kindLabel_->clear();
}

void RecordsPanel::buildForm(const front::Record &record) {
  clearForm();
  shown_ = record;
  form_ = new QFormLayout(formHost_);
  form_->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  const front::KindInfo *kind = holder_->kind(record.kind);
  kindLabel_->setText(kind ? "<b>" + Q(kind->title).toHtmlEscaped() + "</b> (" + Q(record.kind) + "): " +
                                 Q(kind->doc).toHtmlEscaped()
                           : tr("<b>%1</b>: a kind the schema does not know").arg(Q(record.kind)));
  nameEdit_ = new QLineEdit(Q(record.name));
  nameEdit_->setObjectName("recordName");
  nameEdit_->setToolTip(tr("The record's name: other records refer to it by this name."));
  form_->addRow(tr("Name"), nameEdit_);
  if (!kind) return;
  const QString quantity = Q(record.Value("quantity"));
  for (const front::FieldInfo &f : kind->fields) {
    const QString key = Q(f.key);
    QStringList values;
    for (const front::Field &x : record.fields)
      if (x.key == f.key) values << Q(x.value);
    const bool set = !values.isEmpty();
    QWidget *editor = nullptr;
    if (f.repeatable) {
      QPlainTextEdit *edit = new QPlainTextEdit(values.join("\n"));
      edit->setPlaceholderText(tr("one value per line: %1").arg(Q(f.grammar)));
      edit->setMaximumHeight(70);
      editor = edit;
    } else if (f.type == "enum" && !f.values.empty()) {
      QComboBox *combo = new QComboBox;
      if (!f.required || !set) combo->addItem(QString());  // nothing chosen: the field is not written
      for (const std::string &v : f.values) combo->addItem(Q(v));
      if (set && combo->findText(values.first()) < 0) combo->addItem(values.first());  // an alias as written
      combo->setCurrentText(set ? values.first() : QString());
      if (f.key == "quantity") connect(combo, &QComboBox::currentTextChanged, this, &RecordsPanel::quantityChanged);
      editor = combo;
    } else if (f.type == "unit") {
      QComboBox *combo = new QComboBox;
      combo->setEditable(true);
      const std::vector<std::string> units = f.unitOfQuantity ? front::QuantityUnits(quantity.toStdString()) : f.units;
      for (const std::string &u : units) combo->addItem(Q(u));
      combo->setCurrentText(set ? values.first() : QString());
      editor = combo;
    } else {
      QLineEdit *edit = new QLineEdit(set ? values.first() : QString());
      edit->setPlaceholderText(Q(f.grammar));
      editor = edit;
      if (f.type == "path") {
        QWidget *host = new QWidget;
        QHBoxLayout *h = new QHBoxLayout(host);
        h->setContentsMargins(0, 0, 0, 0);
        h->addWidget(edit, 1);
        QToolButton *browse = new QToolButton;
        browse->setText("...");
        connect(browse, &QToolButton::clicked, this, [this, edit]() {
          const QString file = QFileDialog::getOpenFileName(this, tr("Choose a file"), holder_->directory());
          if (!file.isEmpty()) edit->setText(QDir(holder_->directory()).relativeFilePath(file));
        });
        h->addWidget(browse);
        widgets_.insert(key, edit);
        edit->setObjectName("field_" + key);
        edit->setToolTip(ToolTip(f));
        QLabel *label = new QLabel(Q(f.title) + (f.required ? " *" : ""));
        label->setToolTip(ToolTip(f));
        form_->addRow(label, host);
        const bool show = set || !f.rare || everyField_->isChecked();
        label->setVisible(show);
        host->setVisible(show);
        continue;
      }
    }
    editor->setObjectName("field_" + key);
    editor->setToolTip(ToolTip(f));
    widgets_.insert(key, editor);
    QLabel *label = new QLabel(Q(f.title) + (f.required ? " *" : ""));
    label->setToolTip(ToolTip(f));
    form_->addRow(label, editor);
    const bool show = set || !f.rare || everyField_->isChecked();
    label->setVisible(show);
    editor->setVisible(show);
  }
  // fields the schema does not list for the kind: shown as written, kept on Apply (the engine names them)
  for (const front::Field &x : record.fields) {
    bool known = false;
    for (const front::FieldInfo &f : kind->fields) known = known || f.key == x.key;
    if (known || widgets_.contains(Q(x.key))) continue;
    QLineEdit *edit = new QLineEdit(Q(x.value));
    edit->setObjectName("field_" + Q(x.key));
    widgets_.insert(Q(x.key), edit);
    form_->addRow(new QLabel(Q(x.key) + tr(" (not a field of %1)").arg(Q(record.kind))), edit);
  }
}

void RecordsPanel::showEveryField(bool) {
  if (list_->currentRow() >= 0) buildForm(shown_);
}

void RecordsPanel::quantityChanged(const QString &quantity) {
  QComboBox *unit = qobject_cast<QComboBox *>(widgets_.value("unit", nullptr));
  if (!unit) return;
  const QString was = unit->currentText();
  unit->clear();
  for (const std::string &u : front::QuantityUnits(quantity.toStdString())) unit->addItem(Q(u));
  unit->setCurrentText(was);
}

QString RecordsPanel::valueOf(QWidget *w) const {
  if (QComboBox *c = qobject_cast<QComboBox *>(w)) return c->currentText().trimmed();
  if (QLineEdit *e = qobject_cast<QLineEdit *>(w)) return e->text().trimmed();
  if (QPlainTextEdit *p = qobject_cast<QPlainTextEdit *>(w)) return p->toPlainText();
  return QString();
}

bool RecordsPanel::applyForm(QString &why) {
  if (!nameEdit_) {
    why = tr("no record is shown");
    return false;
  }
  front::Record r;
  r.kind = shown_.kind;
  r.name = shown_.name;
  // the record's fields in its own order, then the new ones in the schema's
  QStringList done;
  auto add = [&r](const QString &key, const QString &value, const std::string &comment) {
    if (value.isEmpty()) return;
    r.fields.push_back({key.toStdString(), value.toStdString(), comment});
  };
  const front::KindInfo *kind = holder_->kind(shown_.kind);
  auto isRepeatable = [kind](const std::string &key) {
    if (kind)
      for (const front::FieldInfo &f : kind->fields)
        if (f.key == key) return f.repeatable;
    return false;
  };
  for (const front::Field &x : shown_.fields) {
    const QString key = Q(x.key);
    if (done.contains(key)) continue;
    done << key;
    QWidget *w = widgets_.value(key, nullptr);
    if (!w) {
      r.fields.push_back(x);
      continue;
    }
    if (isRepeatable(x.key)) {
      for (const QString &line : valueOf(w).split('\n')) add(key, line.trimmed(), "");
      continue;
    }
    add(key, valueOf(w), x.comment);
  }
  if (kind)
    for (const front::FieldInfo &f : kind->fields) {
      const QString key = Q(f.key);
      if (done.contains(key) || !widgets_.contains(key)) continue;
      if (f.repeatable) {
        for (const QString &line : valueOf(widgets_.value(key)).split('\n')) add(key, line.trimmed(), "");
        continue;
      }
      add(key, valueOf(widgets_.value(key)), "");
    }
  reloading_ = true;
  bool ok = holder_->setRecord(section_, r, why);
  const QString newName = nameEdit_->text().trimmed();
  if (ok && newName != Q(shown_.name)) ok = holder_->renameRecord(section_, Q(shown_.name), newName, why);
  reloading_ = false;
  reload();
  if (ok) {
    selectRecord(newName);
    emit edited();
  }
  return ok;
}

void RecordsPanel::apply() {
  QString why;
  if (!applyForm(why)) QMessageBox::warning(this, tr("Record not written"), why);
}

QString RecordsPanel::uniqueName(const QString &stem) const {
  QStringList names;
  for (const front::Record &r : holder_->records(section_)) names << Q(r.name);
  QString name = stem;
  for (int i = 2; names.contains(name); i++) name = QString("%1_%2").arg(stem).arg(i);
  return name;
}

void RecordsPanel::addRecord() {
  const QStringList kinds = holder_->kindsOfSection(section_);
  if (kinds.isEmpty()) return;
  QString kind = kinds.first();
  if (kinds.size() > 1) {
    bool ok = false;
    kind = QInputDialog::getItem(this, tr("Add a record"), tr("Kind:"), kinds, 0, false, &ok);
    if (!ok) return;
  }
  bool ok = false;
  const QString name = QInputDialog::getText(this, tr("Add a record"), tr("Name:"), QLineEdit::Normal, uniqueName(kind), &ok);
  if (!ok || name.trimmed().isEmpty()) return;
  front::Record r;
  r.kind = kind.toStdString();
  r.name = name.trimmed().toStdString();
  // the required fields as '?' placeholders: the engine names each one until it is given (E_PLACEHOLDER)
  if (const front::KindInfo *k = holder_->kind(r.kind))
    for (const front::FieldInfo &f : k->fields)
      if (f.required) r.fields.push_back({f.key, "?", ""});
  QString why;
  if (!holder_->setRecord(section_, r, why)) {
    QMessageBox::warning(this, tr("Record not added"), why);
    return;
  }
  selectRecord(name.trimmed());
  emit edited();
}

void RecordsPanel::duplicateRecord() {
  if (list_->currentRow() < 0) return;
  front::Record r = shown_;
  r.name = uniqueName(Q(shown_.name) + "_copy").toStdString();
  QString why;
  if (!holder_->setRecord(section_, r, why)) {
    QMessageBox::warning(this, tr("Record not added"), why);
    return;
  }
  selectRecord(Q(r.name));
  emit edited();
}

void RecordsPanel::removeRecord() {
  if (list_->currentRow() < 0) return;
  const QString name = list_->item(list_->currentRow(), 0)->text();
  QString why;
  if (!holder_->removeRecord(section_, name, why)) QMessageBox::warning(this, tr("Record not removed"), why);
  else emit edited();
}
