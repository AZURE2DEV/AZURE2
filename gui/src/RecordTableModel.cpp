#include "RecordTableModel.h"

#include <QApplication>
#include <QComboBox>
#include <QLineEdit>
#include <QPalette>

#include "ProjectHolder.h"
#include "ValueEditor.h"

namespace front = azure2::front;

namespace {

QString Q(const std::string &s) { return QString::fromStdString(s); }

const int ChoicesRole = Qt::UserRole + 1;
const int EditTextRole = Qt::UserRole + 2;

/// The record as it stands in the file, for the row's tooltip.
QString RecordText(const front::Record &r) {
  QString t = "<pre>" + Q(r.kind).toHtmlEscaped() + " " + Q(r.name).toHtmlEscaped();
  for (const front::Field &f : r.fields) t += "\n    " + Q(f.key).leftJustified(16).toHtmlEscaped() + Q(f.value).toHtmlEscaped();
  return t + "</pre>";
}

bool InPlace(const front::FieldInfo &f) {
  if (f.repeatable) return false;
  const std::string &t = f.value.type;
  return t == "enum" || t == "unit" || t == "reference" || t == "number" || t == "integer" || t == "nuisance" || t == "path";
}

}  // namespace

RecordTableModel::RecordTableModel(ProjectHolder *holder, const QString &section, const QStringList &kinds, QObject *parent)
    : QAbstractTableModel(parent), holder_(holder), section_(section), kinds_(kinds) {
  bool include = false;
  for (const QString &k : kinds_)
    if (const front::KindInfo *info = holder_->kind(k.toStdString())) include = include || !info->includeField.empty();
  if (include) {
    Column c;
    c.type = Column::Include;
    c.title = tr("Use");
    columns_.push_back(c);
  }
  Column name;
  name.type = Column::Name;
  name.title = tr("Name");
  columns_.push_back(name);
  if (kinds_.size() > 1) {
    Column c;
    c.type = Column::Kind;
    c.title = tr("Kind");
    columns_.push_back(c);
  }
  for (const QString &k : kinds_) {
    const front::KindInfo *info = holder_->kind(k.toStdString());
    if (!info) continue;
    for (const std::string &key : info->summaryFields) {
      if (columnOf(Column::Field, key) >= 0) continue;
      for (const front::FieldInfo &f : info->fields) {
        if (f.key != key) continue;
        Column c;
        c.type = Column::Field;
        c.key = key;
        c.field = f;
        c.title = Q(f.title);
        columns_.push_back(c);
        if (f.value.type == "nuisance") {
          Column fit = c;
          fit.type = Column::Fit;
          fit.title = tr("Fit?");
          columns_.push_back(fit);
        }
      }
    }
  }
  reload();
}

int RecordTableModel::columnOf(Column::Type type, const std::string &key) const {
  for (std::size_t i = 0; i < columns_.size(); i++)
    if (columns_[i].type == type && (key.empty() || columns_[i].key == key)) return (int)i;
  return -1;
}

int RecordTableModel::rowOf(const QString &name) const {
  for (std::size_t i = 0; i < records_.size(); i++)
    if (Q(records_[i].name) == name) return (int)i;
  return -1;
}

void RecordTableModel::reload() {
  beginResetModel();
  records_.clear();
  cells_.clear();
  for (const front::Record &r : holder_->records(section_)) {
    if (!kinds_.contains(Q(r.kind))) continue;
    records_.push_back(r);
    std::vector<QString> row;
    for (const Column &c : columns_) {
      QString text;
      if (c.type == Column::Name) text = Q(r.name);
      else if (c.type == Column::Kind) text = holder_->kind(r.kind) ? Q(holder_->kind(r.kind)->title) : Q(r.kind);
      else if (c.type == Column::Field) {
        QStringList values;
        for (const front::Field &f : r.fields)
          if (f.key == c.key) values << Q(front::ShortValueText(r.kind, f.key, f.value));
        text = values.join("; ");
      }
      row.push_back(text);
    }
    cells_.push_back(row);
  }
  endResetModel();
}

int RecordTableModel::rowCount(const QModelIndex &parent) const { return parent.isValid() ? 0 : (int)records_.size(); }
int RecordTableModel::columnCount(const QModelIndex &parent) const { return parent.isValid() ? 0 : (int)columns_.size(); }

const front::KindInfo *RecordTableModel::kindOf(int row) const { return holder_->kind(records_[row].kind); }

std::string RecordTableModel::quantityOf(int row) const {
  const front::KindInfo *k = kindOf(row);
  if (!k) return std::string();
  for (const front::FieldInfo &f : k->fields)
    if (IsQuantityField(f, k->kind)) return records_[row].Value(f.key);
  return std::string();
}

bool RecordTableModel::included(int row) const {
  const front::KindInfo *k = kindOf(row);
  if (!k || k->includeField.empty()) return true;
  return records_[row].Value(k->includeField) != k->includeOff;
}

QVariant RecordTableModel::data(const QModelIndex &index, int role) const {
  if (!index.isValid() || index.row() >= (int)records_.size()) return QVariant();
  const int row = index.row();
  const Column &c = columns_[index.column()];
  const front::Record &r = records_[row];
  const bool has = std::any_of(r.fields.begin(), r.fields.end(), [&c](const front::Field &f) { return f.key == c.key; });
  switch (role) {
    case Qt::DisplayRole:
      if (c.type == Column::Include || c.type == Column::Fit) return QVariant();
      return cells_[row][index.column()];
    case Qt::EditRole:
      return editText(index);
    case Qt::CheckStateRole:
      if (c.type == Column::Include) {
        const front::KindInfo *k = kindOf(row);
        if (!k || k->includeField.empty()) return QVariant();
        return included(row) ? Qt::Checked : Qt::Unchecked;
      }
      if (c.type == Column::Fit && has) {
        front::ValueParts p;
        std::string why;
        if (!front::ParseFieldValue(r.kind, c.key, r.Value(c.key), quantityOf(row), p, why) || p.reference || p.profiled)
          return QVariant();
        return p.fit ? Qt::Checked : Qt::Unchecked;
      }
      return QVariant();
    case Qt::ToolTipRole:
      if (c.type == Column::Include) return tr("Use this record (unchecked: it is left out of the run)");
      if (c.type == Column::Fit) return tr("Vary %1 in a fit").arg(Q(c.field.title));
      return RecordText(r);
    case Qt::ForegroundRole:
      if (!included(row)) return QApplication::palette().color(QPalette::Disabled, QPalette::Text);
      return QVariant();
    case ChoicesRole:
      return choices(index);
    case EditTextRole:
      return editText(index);
  }
  return QVariant();
}

QVariant RecordTableModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if (orientation != Qt::Horizontal || section < 0 || section >= (int)columns_.size()) return QVariant();
  const Column &c = columns_[section];
  if (role == Qt::DisplayRole) return c.title;
  if (role == Qt::ToolTipRole && c.type == Column::Field) return Q(c.field.doc);
  return QVariant();
}

Qt::ItemFlags RecordTableModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::NoItemFlags;
  Qt::ItemFlags f = Qt::ItemIsEnabled | Qt::ItemIsSelectable;
  const Column &c = columns_[index.column()];
  const front::KindInfo *k = kindOf(index.row());
  if (c.type == Column::Include && k && !k->includeField.empty()) f |= Qt::ItemIsUserCheckable;
  if (c.type == Column::Fit && data(index, Qt::CheckStateRole).isValid()) f |= Qt::ItemIsUserCheckable;
  if (c.type == Column::Name) f |= Qt::ItemIsEditable;
  if (c.type == Column::Field && InPlace(c.field)) {
    // a field of this record's kind only (a column may be another kind's)
    bool ofKind = false;
    if (k)
      for (const front::FieldInfo &x : k->fields) ofKind = ofKind || x.key == c.key;
    if (ofKind) f |= Qt::ItemIsEditable;
  }
  return f;
}

QStringList RecordTableModel::choices(const QModelIndex &index) const {
  const Column &c = columns_[index.column()];
  if (c.type != Column::Field) return QStringList();
  const front::ValueInfo &v = c.field.value;
  QStringList words;
  if (v.type == "enum") {
    for (const std::string &w : v.values) words << Q(w);
  } else if (v.type == "unit") {
    for (const std::string &u : v.unitOfQuantity ? front::QuantityUnits(quantityOf(index.row())) : v.units) words << Q(u);
  } else if (v.type == "reference") {
    words = RecordNamesOfKinds(holder_, v.values);
  } else {
    return QStringList();
  }
  if (!c.field.required) words.prepend(QString());
  return words;
}

QString RecordTableModel::editText(const QModelIndex &index) const {
  const Column &c = columns_[index.column()];
  const front::Record &r = records_[index.row()];
  if (c.type == Column::Name) return Q(r.name);
  if (c.type != Column::Field) return QString();
  const QString value = Q(r.Value(c.key));
  if (c.field.value.type == "nuisance") {
    front::ValueParts p;
    std::string why;
    if (front::ParseFieldValue(r.kind, c.key, value.toStdString(), quantityOf(index.row()), p, why) && !p.reference)
      return Q(p.number);
  }
  if (c.field.value.type == "reference") {
    front::ValueParts p;
    std::string why;
    if (front::ParseFieldValue(r.kind, c.key, value.toStdString(), "", p, why)) return Q(p.text);
  }
  return value;
}

bool RecordTableModel::writeField(int row, const std::string &key, const std::string &value, QString &why) {
  front::Record r = records_[row];
  bool found = false;
  std::vector<front::Field> fields;
  for (const front::Field &f : r.fields) {
    if (f.key == key) {
      if (found) continue;
      found = true;
      if (!value.empty()) fields.push_back({key, value, f.comment});
      continue;
    }
    fields.push_back(f);
  }
  if (!found && !value.empty()) fields.push_back({key, value, ""});
  r.fields = fields;
  writing_ = true;
  const bool ok = holder_->setRecord(section_, r, why);
  writing_ = false;
  if (!ok) return false;
  // this row read again (the others did not move)
  for (const front::Record &x : holder_->records(section_))
    if (x.name == r.name) records_[row] = x;
  for (std::size_t c = 0; c < columns_.size(); c++)
    if (columns_[c].type == Column::Field) {
      QStringList values;
      for (const front::Field &f : records_[row].fields)
        if (f.key == columns_[c].key) values << Q(front::ShortValueText(r.kind, f.key, f.value));
      cells_[row][c] = values.join("; ");
    }
  emit dataChanged(index(row, 0), index(row, (int)columns_.size() - 1));
  emit recordWritten(Q(r.name));
  return true;
}

bool RecordTableModel::setData(const QModelIndex &idx, const QVariant &value, int role) {
  if (!idx.isValid()) return false;
  const int row = idx.row();
  const Column &c = columns_[idx.column()];
  const front::Record &r = records_[row];
  QString why;
  if (c.type == Column::Include && role == Qt::CheckStateRole)
    return setIncluded({Q(r.name)}, value.toInt() == Qt::Checked, why) || (emit refused(why), false);
  if (c.type == Column::Fit && role == Qt::CheckStateRole) {
    front::ValueParts p;
    std::string error;
    if (!front::ParseFieldValue(r.kind, c.key, r.Value(c.key), quantityOf(row), p, error)) return false;
    p.fit = value.toInt() == Qt::Checked;
    return writeField(row, c.key, front::FieldValueText(r.kind, c.key, p), why) || (emit refused(why), false);
  }
  if (role != Qt::EditRole) return false;
  if (c.type == Column::Name) {
    const QString to = value.toString().trimmed();
    if (to.isEmpty() || to == Q(r.name)) return false;
    writing_ = true;
    const bool ok = holder_->renameRecord(section_, Q(r.name), to, why);
    writing_ = false;
    if (!ok) {
      emit refused(why);
      return false;
    }
    records_[row].name = to.toStdString();
    cells_[row][idx.column()] = to;
    emit dataChanged(idx, idx);
    emit recordWritten(to);
    return true;
  }
  if (c.type != Column::Field) return false;
  std::string text = value.toString().trimmed().toStdString();
  if (c.field.value.type == "nuisance") {
    front::ValueParts p;
    std::string error;
    if (!front::ParseFieldValue(r.kind, c.key, r.Value(c.key), quantityOf(row), p, error) || p.reference) return false;
    if (p.number == text) return false;
    p.number = text;
    text = front::FieldValueText(r.kind, c.key, p);
  }
  if (text == r.Value(c.key)) return false;
  if (!text.empty()) {
    front::ValueParts p;
    std::string error;
    if (!front::ParseFieldValue(r.kind, c.key, text, quantityOf(row), p, error)) {
      emit refused(Q(c.field.title) + ": " + Q(error));
      return false;
    }
  } else if (c.field.required) {
    emit refused(tr("%1 is required").arg(Q(c.field.title)));
    return false;
  }
  return writeField(row, c.key, text, why) || (emit refused(why), false);
}

bool RecordTableModel::setIncluded(const QStringList &names, bool on, QString &why) {
  for (const QString &name : names) {
    const int row = rowOf(name);
    if (row < 0) continue;
    const front::KindInfo *k = kindOf(row);
    if (!k || k->includeField.empty() || included(row) == on) continue;
    // on: the field left out (its default is on), so that a record left out and used again is as it was
    if (!writeField(row, k->includeField, on ? std::string() : k->includeOff, why)) return false;
  }
  return true;
}

QWidget *RecordCellDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const {
  const QStringList choices = index.data(ChoicesRole).toStringList();
  if (choices.isEmpty()) {
    QLineEdit *e = new QLineEdit(parent);
    e->setFrame(false);
    return e;
  }
  QComboBox *c = new QComboBox(parent);
  for (const QString &w : choices) c->addItem(w.isEmpty() ? tr("(none)") : w, w);
  Q_UNUSED(option);
  return c;
}

void RecordCellDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const {
  const QString text = index.data(EditTextRole).toString();
  if (QComboBox *c = qobject_cast<QComboBox *>(editor)) {
    int at = c->findData(text);
    if (at < 0) {
      c->addItem(text, text);
      at = c->count() - 1;
    }
    c->setCurrentIndex(at);
    c->showPopup();
  } else if (QLineEdit *e = qobject_cast<QLineEdit *>(editor)) {
    e->setText(text);
  }
}

void RecordCellDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const {
  if (QComboBox *c = qobject_cast<QComboBox *>(editor))
    model->setData(index, c->currentData().toString(), Qt::EditRole);
  else if (QLineEdit *e = qobject_cast<QLineEdit *>(editor))
    model->setData(index, e->text(), Qt::EditRole);
}
