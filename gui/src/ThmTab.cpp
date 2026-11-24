#include "ThmTab.h"

#include <QComboBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QScrollArea>
#include <QSplitter>
#include <QTabWidget>
#include <QVBoxLayout>

#include "ExpertMode.h"
#include "SectionTextPanel.h"

namespace {

QString EditorText(QWidget *w) {
  if (QComboBox *c = qobject_cast<QComboBox *>(w)) return c->currentData().toString();
  if (QLineEdit *e = qobject_cast<QLineEdit *>(w)) return e->text().trimmed();
  return QString();
}

void SetEditorText(QWidget *w, const QString &value) {
  if (QComboBox *c = qobject_cast<QComboBox *>(w)) {
    int at = c->findData(value);
    if (at < 0 && !value.isEmpty()) {  // a value the table does not list: shown, the engine names it
      c->addItem(value, value);
      at = c->count() - 1;
    }
    c->setCurrentIndex(at < 0 ? 0 : at);
  } else if (QLineEdit *e = qobject_cast<QLineEdit *>(w)) {
    e->setText(value);
  }
}

}  // namespace

ThmTab::ThmTab(ProjectHolder *holder, QWidget *parent) : QWidget(parent), holder_(holder) {
  kinds_ = azure2::front::ThmKinds();
  message_ = new QLabel;
  message_->setObjectName("thmMessage");
  message_->setWordWrap(true);

  expert_ = expert::Enabled();
  QFormLayout *optionsForm = new QFormLayout;
  QFormLayout *experimentForm = new QFormLayout;
  optionsForm_ = optionsForm;
  experimentForm_ = experimentForm;
  for (int k = 0; k < 2 && k < (int)kinds_.size(); k++)
    for (const azure2::front::FieldInfo &f : kinds_[k].fields) {
      QWidget *editor = editorFor(f, k == 1);
      QLabel *label = new QLabel(QString::fromStdString(f.title) + ":");
      label->setToolTip(editor->toolTip());
      (k == 0 ? optionsForm : experimentForm)->addRow(label, editor);
      (k == 0 ? optionEditors_ : experimentEditors_).insert(QString::fromStdString(f.key), editor);
      if (k == 1 && !f.appliesWhen.empty())
        presentWith_.insert(QString::fromStdString(f.key),
                            QString::fromStdString(f.appliesWhen.substr(0, f.appliesWhen.find(' '))));
    }
  QGroupBox *optionsBox = new QGroupBox(kinds_.empty() ? tr("Options") : QString::fromStdString(kinds_[0].title));
  optionsBox->setToolTip(kinds_.empty() ? QString() : QString::fromStdString(kinds_[0].doc));
  optionsBox->setLayout(optionsForm);

  experimentList_ = new QListWidget;
  experimentList_->setObjectName("thmExperiments");
  connect(experimentList_, &QListWidget::currentRowChanged, this, &ThmTab::fill);
  QPushButton *add = new QPushButton(tr("Add..."));
  QPushButton *remove = new QPushButton(tr("Remove"));
  connect(add, &QPushButton::clicked, this, [this]() {
    const QString name = QInputDialog::getText(this, tr("New THM experiment"), tr("Name (letters, digits, _):"));
    QString why;
    if (!name.isEmpty() && !addExperiment(name, why)) QMessageBox::warning(this, tr("THM"), why);
  });
  connect(remove, &QPushButton::clicked, this, [this]() {
    QString why;
    if (!currentExperiment().isEmpty() && !removeExperiment(currentExperiment(), why))
      QMessageBox::warning(this, tr("THM"), why);
  });
  QWidget *listPane = new QWidget;
  QVBoxLayout *lp = new QVBoxLayout(listPane);
  lp->setContentsMargins(0, 0, 0, 0);
  lp->addWidget(new QLabel(tr("Experiments:")));
  lp->addWidget(experimentList_, 1);
  QHBoxLayout *lb = new QHBoxLayout;
  lb->addWidget(add);
  lb->addWidget(remove);
  lp->addLayout(lb);
  QWidget *experimentFormHost = new QWidget;
  experimentFormHost->setLayout(experimentForm);
  QScrollArea *experimentScroll = new QScrollArea;
  experimentScroll->setWidget(experimentFormHost);
  experimentScroll->setWidgetResizable(true);
  QSplitter *experimentSplit = new QSplitter;
  experimentSplit->addWidget(listPane);
  experimentSplit->addWidget(experimentScroll);
  experimentSplit->setStretchFactor(1, 3);
  QGroupBox *experimentBox =
      new QGroupBox(kinds_.size() > 1 ? QString::fromStdString(kinds_[1].title) : tr("Experiments"));
  experimentBox->setToolTip(kinds_.size() > 1 ? QString::fromStdString(kinds_[1].doc) : QString());
  (new QVBoxLayout(experimentBox))->addWidget(experimentSplit);

  QWidget *formPage = new QWidget;
  QVBoxLayout *fp = new QVBoxLayout(formPage);
  fp->addWidget(optionsBox);
  fp->addWidget(experimentBox, 1);
  fp->addWidget(message_);
  text_ = new SectionTextPanel(holder_, "thm",
                               tr("The whole <thm> section as text (the Trojan horse experiments: one line per "
                                  "experiment, key=value).  The form and this text edit the same lines."));
  QTabWidget *pages = new QTabWidget;
  pages->setObjectName("thmPages");
  pages->addTab(formPage, tr("Form"));
  pages->addTab(text_, tr("Text"));
  (new QVBoxLayout(this))->addWidget(pages);
  connect(holder_, &ProjectHolder::opened, this, &ThmTab::reload);
  connect(holder_, &ProjectHolder::sectionEdited, this, [this](const QString &s) {
    if (s == "thm" && !filling_) reload();
  });
  reload();
}

QWidget *ThmTab::editorFor(const azure2::front::FieldInfo &f, bool experiment) {
  const QString key = QString::fromStdString(f.key);
  QString tip = key + ": " + QString::fromStdString(f.doc);
  if (!f.grammar.empty()) tip += "\n" + QString::fromStdString(f.grammar);
  if (!f.appliesWhen.empty()) tip += "\n" + tr("Applies when %1.").arg(QString::fromStdString(f.appliesWhen));
  QWidget *editor;
  if (f.type == "enum") {
    QComboBox *c = new QComboBox;
    c->addItem(f.defaultText.empty() ? tr("(not given)") : tr("(not given: %1)").arg(QString::fromStdString(f.defaultText)),
               QString());
    for (const std::string &v : f.values) c->addItem(QString::fromStdString(v), QString::fromStdString(v));
    connect(c, &QComboBox::activated, this, [this, key, experiment, c]() {
      QString why;
      const bool ok = experiment ? setExperimentKey(currentExperiment(), key, EditorText(c), why)
                                 : setOption(key, EditorText(c), why);
      if (!ok) message_->setText(why);
    });
    editor = c;
  } else {
    QLineEdit *e = new QLineEdit;
    e->setPlaceholderText(!f.defaultText.empty() ? tr("default %1").arg(QString::fromStdString(f.defaultText))
                                                  : QString::fromStdString(f.grammar));
    connect(e, &QLineEdit::editingFinished, this, [this, key, experiment, e]() {
      if (filling_) return;
      QString why;
      const bool ok = experiment ? setExperimentKey(currentExperiment(), key, EditorText(e), why)
                                 : setOption(key, EditorText(e), why);
      if (!ok) message_->setText(why);
    });
    editor = e;
  }
  editor->setObjectName(QString("thm_%1_%2").arg(experiment ? "experiment" : "option", key));
  editor->setToolTip(tip);
  return editor;
}

QStringList ThmTab::experiments() const {
  QStringList out;
  for (const azure2::front::ThmEntry &e : azure2::front::ThmEntries(holder_->project())) {
    const QString name = QString::fromStdString(e.experiment);
    if (!name.isEmpty() && !out.contains(name)) out << name;
  }
  return out;
}

QString ThmTab::currentExperiment() const {
  return experimentList_->currentItem() ? experimentList_->currentItem()->text() : QString();
}

void ThmTab::reload() {
  filling_ = true;
  const QString selected = currentExperiment();
  experimentList_->clear();
  experimentList_->addItems(experiments());
  const QList<QListWidgetItem *> again = experimentList_->findItems(selected, Qt::MatchExactly);
  experimentList_->setCurrentRow(again.isEmpty() ? (experimentList_->count() ? 0 : -1) : experimentList_->row(again[0]));
  QMap<QString, QString> options;
  for (const azure2::front::ThmEntry &e : azure2::front::ThmEntries(holder_->project()))
    if (e.experiment.empty()) options[QString::fromStdString(e.key)] = QString::fromStdString(e.value);
  for (auto it = optionEditors_.constBegin(); it != optionEditors_.constEnd(); ++it) {
    SetEditorText(it.value(), options.value(it.key()));
    // an option the project gives is shown in either view; the others are expert options (qt_gui.md section 19.4)
    optionsForm_->setRowVisible(it.value(), expert_ || options.contains(it.key()));
  }
  filling_ = false;
  fill();
}

void ThmTab::setExpertMode(bool on) {
  expert_ = on;
  reload();
}

bool ThmTab::rowShown(const QString &key, bool experiment) const {
  QWidget *e = (experiment ? experimentEditors_ : optionEditors_).value(key, nullptr);
  return e && !e->isHidden();
}

void ThmTab::fill() {
  filling_ = true;
  const std::string name = currentExperiment().toStdString();
  QMap<QString, QString> values;
  for (const azure2::front::ThmEntry &e : azure2::front::ThmEntries(holder_->project()))
    if (!name.empty() && e.experiment == name) values[QString::fromStdString(e.key)] = QString::fromStdString(e.value);
  for (auto it = experimentEditors_.constBegin(); it != experimentEditors_.constEnd(); ++it) {
    SetEditorText(it.value(), values.value(it.key()));
    // the keys an occasional THM user sets, and any key the experiment gives, in either view; the rest in expert view
    const bool standard = QStringList{"beam", "target", "spectator", "Ebeam", "ps", "background"}.contains(it.key());
    experimentForm_->setRowVisible(it.value(), expert_ || standard || values.contains(it.key()));
    const QString needs = presentWith_.value(it.key());
    it.value()->setEnabled(!name.empty() && (needs.isEmpty() || values.contains(needs)));
  }
  filling_ = false;
}

bool ThmTab::setOption(const QString &key, const QString &value, QString &why) {
  std::string reason;
  filling_ = true;
  const bool ok = azure2::front::SetThmOption(holder_->project(), key.toStdString(), value.toStdString(), reason,
                                              value.isEmpty());
  filling_ = false;
  if (!ok) {
    why = QString::fromStdString(reason);
    return false;
  }
  if (QWidget *e = optionEditors_.value(key)) SetEditorText(e, value);
  message_->clear();
  emit holder_->sectionEdited("thm");
  return true;
}

bool ThmTab::setExperimentKey(const QString &experiment, const QString &key, const QString &value, QString &why) {
  if (experiment.isEmpty()) {
    why = tr("no experiment is selected");
    return false;
  }
  std::string reason;
  filling_ = true;
  const bool ok = azure2::front::SetThmExperimentKey(holder_->project(), experiment.toStdString(), key.toStdString(),
                                                     value.toStdString(), reason, value.isEmpty());
  filling_ = false;
  if (!ok) {
    why = tr("experiment[%1]: %2").arg(experiment, QString::fromStdString(reason));
    return false;
  }
  message_->clear();
  emit holder_->sectionEdited("thm");
  selectExperiment(experiment);
  return true;
}

bool ThmTab::addExperiment(const QString &name, QString &why) {
  // an experiment is its line: made with its first key (the beam's), the rest set in the form
  if (experiments().contains(name)) {
    why = tr("experiment[%1] exists").arg(name);
    return false;
  }
  if (!setExperimentKey(name, "Ebeam", "1", why)) return false;
  return selectExperiment(name);
}

bool ThmTab::removeExperiment(const QString &name, QString &why) {
  std::string reason;
  if (!azure2::front::SetThmExperimentKey(holder_->project(), name.toStdString(), "", "", reason, true)) {
    why = QString::fromStdString(reason);
    return false;
  }
  emit holder_->sectionEdited("thm");
  reload();
  return true;
}

bool ThmTab::selectExperiment(const QString &name) {
  reload();
  const QList<QListWidgetItem *> found = experimentList_->findItems(name, Qt::MatchExactly);
  if (found.isEmpty()) return false;
  experimentList_->setCurrentItem(found[0]);
  fill();
  return true;
}

void ThmTab::applyText() { text_->apply(); }
