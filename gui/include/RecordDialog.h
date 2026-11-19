#ifndef RECORDDIALOG_H
#define RECORDDIALOG_H

#include <QDialog>
#include <QList>
#include <QMap>
#include <QStringList>
#include <functional>

#include "azure2/front/Engine.h"

class FieldEditor;
class ProjectHolder;

QT_BEGIN_NAMESPACE
class QComboBox;
class QFormLayout;
class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QPushButton;
QT_END_NAMESPACE

/*!
 * The pop-up of a record, to add one or to edit one (docs/upgrade/qt_gui.md section 15.10), as dev's AddSegDataDialog
 * and AddSegTestDialog were for the segments: every field a record of the kind is usually given (the kind's required
 * fields and its summary fields: for a data set the reaction, quantity, unit, data file, columns, windows, the
 * normalization, target, beam and detector; for an extrapolation the reaction, quantity, the energy and angle grids
 * and the target effects), the fields the quantity requires, and in edit mode every field the record gives.  Each
 * value by its own widgets (FieldEditor), every choice from the facade.  A data file's first lines are shown under
 * it and the columns made as many as the file has.  OK checks the record with the engine: a refusal is shown under
 * its field and the dialog stays (OK again writes it as it is); "Add another" writes it and keeps the dialog open
 * with the same choices for the next record.
 */
class RecordDialog : public QDialog {
  Q_OBJECT

 public:
  /// `kinds`: the kinds the dialog may make (a drop-down when several); `writer` writes a record (Add: the panel's
  /// add; Edit: its replace), returning false with the reason.
  RecordDialog(ProjectHolder *holder, const QString &section, const QStringList &kinds,
               std::function<bool(const azure2::front::Record &, QString &)> writer, QWidget *parent = nullptr);
  /// Edit mode: the record shown, its name and kind fixed until changed; every field it gives has a row.
  void setRecord(const azure2::front::Record &record);
  /// The record as the dialog holds it (a field left empty: "?" in add mode, left out in edit mode).
  azure2::front::Record record() const;
  /// The engine's errors of the project with the record written, by field ("" the record's own); shown under the
  /// fields.  True when there is none.
  bool check();
  /// Write the record (check first): false with the reason when it is not written.
  bool commit(QString &why);
  FieldEditor *fieldEditor(const QString &key) const;
  QPlainTextEdit *filePreview() const { return preview_; }
  QLineEdit *nameEdit() const { return name_; }
  QPushButton *addAnotherButton() const { return another_; }
  QLabel *statusLabel() const { return status_; }
  /// The names written by this dialog, in order.
  QStringList written() const { return written_; }

 public slots:
  void accept() override;
  void addAnother();

 private:
  void build(const QString &quantity, const QMap<QString, QString> &keep);
  QString quantityText() const;
  void quantityChanged();
  void showFile();
  QString uniqueName(const QString &stem) const;

  ProjectHolder *holder_;
  QString section_;
  std::function<bool(const azure2::front::Record &, QString &)> writer_;
  QComboBox *kind_;
  QLineEdit *name_;
  QLabel *doc_, *status_, *fileNote_;
  QFormLayout *fields_;
  QPlainTextEdit *preview_;
  QPushButton *another_;
  QList<FieldEditor *> editors_;
  QString quantityKey_, fileKey_, columnsKey_;
  bool editing_ = false;
  azure2::front::Record original_;
  QString lastChecked_;  ///< the record's text when OK last found errors (OK again writes it)
  QStringList written_;
  bool rebuilding_ = false;
};

#endif
