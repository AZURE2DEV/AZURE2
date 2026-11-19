#include "ConfigDialogs.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QScrollArea>
#include <QSpinBox>
#include <QVBoxLayout>

ConfigDialog::ConfigDialog(ProjectHolder *holder, const QString &title, const QStringList &groups,
                           const QStringList &exclude, QWidget *parent)
    : QDialog(parent), holder_(holder) {
  setWindowTitle(title);
  const QJsonArray table = QJsonDocument::fromJson(QByteArray::fromStdString(azure2::front::ConfigTableJson())).array();
  QMap<QString, QFormLayout *> forms;
  QStringList order;
  QWidget *host = new QWidget;
  QVBoxLayout *hostLayout = new QVBoxLayout(host);
  for (const QJsonValue &v : table) {
    const QJsonObject row = v.toObject();
    const QString name = row["name"].toString(), group = row["group"].toString(), kind = row["kind"].toString();
    if (kind == "version") continue;  // the file's format: written by the converter, not edited
    if (!groups.isEmpty() && !groups.contains(group)) continue;
    if (groups.isEmpty() && exclude.contains(group)) continue;
    if (!forms.contains(group)) {
      QGroupBox *box = new QGroupBox(group);
      forms.insert(group, new QFormLayout(box));
      hostLayout->addWidget(box);
    }
    const bool option = row["option"].toBool();
    const bool absent = holder_->configValue(name, "\x01") == "\x01";
    QString value = holder_->configValue(name);
    QWidget *editor = nullptr;
    if (kind == "bool") {
      QCheckBox *box = new QCheckBox;
      bool on = option ? row["default"].toBool() : false;
      if (!absent) on = value == "1" || value == "true";
      box->setChecked(on);
      value = on ? "1" : "0";
      editor = box;
    } else if (kind == "word") {
      QComboBox *combo = new QComboBox;
      for (const QJsonValue &w : row["values"].toArray()) combo->addItem(w.toString());
      combo->setCurrentText(value);
      editor = combo;
    } else if (kind == "int") {
      QSpinBox *spin = new QSpinBox;
      spin->setRange(row["minimum"].toInt(), 1000000);
      spin->setValue(absent ? row["minimum"].toInt() : value.toInt());
      value = QString::number(spin->value());
      editor = spin;
    } else {
      QLineEdit *edit = new QLineEdit(value);
      editor = edit;
    }
    QString tip = "<b>" + name + "</b>: " + row["doc"].toString().toHtmlEscaped();
    if (!row["absent"].toString().isEmpty()) tip += "<br>Without the entry: " + row["absent"].toString().toHtmlEscaped();
    if (option) tip += QString("<br>Also the program's flag %1 / %2 and pyazr's keyword %3.")
                           .arg(row["flag"].toString(), row["negation"].toString(), name);
    editor->setToolTip(tip);
    editor->setObjectName("config_" + name);
    QLabel *label = new QLabel(row["title"].toString().isEmpty() ? name : row["title"].toString());
    label->setToolTip(tip);
    forms.value(group)->addRow(label, editor);
    editors_.insert(name, editor);
    initial_.insert(name, absent ? QString("\x01") + value : value);
    rows_.insert(name, row);
  }
  hostLayout->addStretch(1);
  QScrollArea *scroll = new QScrollArea;
  scroll->setWidgetResizable(true);
  scroll->setWidget(host);
  QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(buttons, &QDialogButtonBox::accepted, this, &ConfigDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, this, &ConfigDialog::reject);
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->addWidget(scroll, 1);
  machine_ = new QFormLayout;
  layout->addLayout(machine_);
  layout->addWidget(buttons);
  resize(640, 520);
}

void ConfigDialog::addMachineRow(const QString &label, QWidget *editor) { machine_->addRow(label, editor); }

bool ConfigDialog::write(QString &why) {
  for (auto it = editors_.constBegin(); it != editors_.constEnd(); ++it) {
    QString value;
    if (QCheckBox *b = qobject_cast<QCheckBox *>(it.value())) value = b->isChecked() ? "1" : "0";
    else if (QComboBox *c = qobject_cast<QComboBox *>(it.value())) value = c->currentText();
    else if (QSpinBox *s = qobject_cast<QSpinBox *>(it.value())) value = QString::number(s->value());
    else if (QLineEdit *e = qobject_cast<QLineEdit *>(it.value())) value = e->text().trimmed();
    QString was = initial_.value(it.key());
    const bool absent = was.startsWith(QChar(1));
    if (absent) was = was.mid(1);
    if (value == was) continue;  // unchanged: the file keeps its line, or keeps having none
    if (absent && value.isEmpty()) continue;
    if (!holder_->setConfig(it.key(), value, why, value.isEmpty() && rows_.value(it.key())["kind"].toString() == "path")) {
      why = it.key() + ": " + why;
      return false;
    }
  }
  return true;
}

void ConfigDialog::accept() {
  QString why;
  if (!write(why)) {
    QMessageBox::warning(this, tr("Not written"), why);
    return;
  }
  QDialog::accept();
}
