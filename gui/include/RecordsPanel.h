#ifndef RECORDSPANEL_H
#define RECORDSPANEL_H

#include <QMap>
#include <QWidget>

#include "ProjectHolder.h"

QT_BEGIN_NAMESPACE
class QCheckBox;
class QFormLayout;
class QLabel;
class QLineEdit;
class QPushButton;
class QTableWidget;
QT_END_NAMESPACE

/*!
 * The records of one format-3 section (docs/upgrade/qt_gui.md section 2): a list (name, kind, what the record is) and
 * the form of the chosen record, built from the schema's row of its kind (azure2::front::RecordKinds): a closed choice
 * is a drop-down, a unit one of its dimension's (of the record's quantity for a data set), a path has a file chooser;
 * every label carries the field's doc and grammar.  Rare fields are one click away ("Every field") unless the record
 * sets them.  Apply writes the record through the engine's writer (ProjectHolder::setRecord): only its changed lines
 * move.  Nothing is checked here: the engine's diagnostics (Validate) are shown in the window's panel.
 */
class RecordsPanel : public QWidget {
  Q_OBJECT

 public:
  RecordsPanel(ProjectHolder *holder, const QString &section, const QString &help, QWidget *parent = nullptr);
  QString section() const { return section_; }
  /// Select the record `name` in the list (from a diagnostic).
  bool selectRecord(const QString &name);
  int recordCount() const;
  /// Apply the form to the record it shows (the tests and Save call it); false with the engine's reason.
  bool applyForm(QString &why);
  /// The form's widget of a field, for the tests (nullptr when the form has none).
  QWidget *fieldWidget(const QString &key) const { return widgets_.value(key, nullptr); }

 signals:
  void edited();

 public slots:
  void reload();

 private slots:
  void selectionChanged();
  void addRecord();
  void duplicateRecord();
  void removeRecord();
  void apply();
  void showEveryField(bool on);
  void quantityChanged(const QString &quantity);

 private:
  void buildForm(const azure2::front::Record &record);
  void clearForm();
  QString valueOf(QWidget *w) const;
  QString uniqueName(const QString &stem) const;

  ProjectHolder *holder_;
  QString section_;
  QTableWidget *list_;
  QPushButton *addButton_, *removeButton_, *duplicateButton_, *applyButton_;
  QCheckBox *everyField_;
  QWidget *formHost_;
  QFormLayout *form_ = nullptr;
  QLineEdit *nameEdit_ = nullptr;
  QLabel *kindLabel_;
  QMap<QString, QWidget *> widgets_;  ///< field key -> its editor
  azure2::front::Record shown_;        ///< the record as last read
  bool reloading_ = false;
};

#endif
