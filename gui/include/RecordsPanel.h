#ifndef RECORDSPANEL_H
#define RECORDSPANEL_H

#include <QList>
#include <QMap>
#include <QWidget>

#include "ProjectHolder.h"

class FieldEditor;
class RecordTableModel;

QT_BEGIN_NAMESPACE
class QAction;
class QCheckBox;
class QFormLayout;
class QLabel;
class QLineEdit;
class QPushButton;
class QSortFilterProxyModel;
class QTableView;
class QVBoxLayout;
QT_END_NAMESPACE

/*!
 * The records of one format-3 section, or of some of its kinds (docs/upgrade/qt_gui.md section 11): a table with real
 * columns (RecordTableModel: the include switch, the name, the summary fields of the schema; sortable, filtered,
 * several rows selected; cells edited in place) and the form of the chosen record, built from the schema's row of its
 * kind: every value by its own widgets (FieldEditor: drop-downs of units, words, record names and pairs; numbers;
 * "Fit"; the prior), a field that belongs with another in its box, a box closed unless the record sets it.  Apply
 * writes the record through the engine's writer (only its changed lines move) after the engine has checked it; an
 * unedited record is written back byte for byte.
 */
class RecordsPanel : public QWidget {
  Q_OBJECT

 public:
  /// `kinds`: the kinds of the section this panel shows (empty: all of them).
  RecordsPanel(ProjectHolder *holder, const QString &section, const QString &help, QWidget *parent = nullptr,
               const QStringList &kinds = QStringList());
  QString section() const { return section_; }
  QStringList kinds() const { return kinds_; }
  /// Select the record `name` in the table (from a diagnostic); false when this panel does not show it.
  bool selectRecord(const QString &name);
  int recordCount() const;
  /// Apply the form to the record it shows (the tests and Save call it); false with the engine's reason.  With
  /// `refuseErrors`, a record the engine would refuse is not written (`why` names the field).
  bool applyForm(QString &why, bool refuseErrors = false);
  /// The editor of a field (the first, for a repeatable field), for the tests; nullptr when the form has none.
  QWidget *fieldWidget(const QString &key) const;
  FieldEditor *fieldEditor(const QString &key, int occurrence = 0) const;
  /// Set a field's text in the form as the user would by its text line (an edit).
  bool setFieldText(const QString &key, const QString &text);
  /// What Add, Duplicate, Remove, Move and the form's name do, without their dialogs (the tests; false with the reason).
  bool addNamedRecord(const QString &kind, const QString &name, QString &why);
  bool duplicateShown(QString &newName, QString &why);
  bool removeShown(QString &why);
  bool renameShown(const QString &newName, QString &why);
  bool moveShown(int delta, QString &why);
  /// Switch the selected records on or off.
  bool includeSelected(bool on, QString &why);
  /// Select several records (the table's selection).
  void selectRecords(const QStringList &names);
  QStringList selectedNames() const;
  QStringList recordNames() const;
  /// The records of the project that name the shown record (a data set's target, ...), as "section/name".
  QStringList usedBy() const;
  /// Make the records `names` name the shown record (and only them), in their field that may name its kind.
  bool assignShownTo(const QStringList &names, QString &why);
  RecordTableModel *model() const { return model_; }
  QTableView *table() const { return table_; }
  /// The text of the filter line.
  void setFilter(const QString &text);
  /// The boxes of the form, by their field's key (for the tests and the screenshots).
  QWidget *box(const QString &key) const { return boxes_.value(key, nullptr); }
  QString shownName() const;
  /// The preview of the record's kernels (nullptr or hidden when it has none).
  QWidget *preview() const { return preview_; }
  /// The kernels the preview draws, as the facade gave them last (empty when none).
  const azure2::front::KernelPreview &kernelPreview() const { return kernels_; }
  /// Draw the kernels for user `user` ("data: <name>", "" the first) at lab energy `labEnergy` (MeV; NaN the default).
  void setPreviewPoint(const QString &user, double labEnergy);
  /// Why the kernels are not drawn ("" when they are).
  QString previewWhy() const { return previewWhy_; }
  /// The fields of the shown record that its quantity (as the form holds it) hides, and of those the ones that hold a
  /// value Apply will drop (the quantity was changed in the form).
  QStringList hiddenByQuantity() const { return hiddenByQuantity_; }
  /// Whether a field's row is shown on the form (the tests).
  bool fieldRowShown(const QString &key) const;

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
  void revert();
  void showEveryField(bool on);
  void assignDialog();
  void contextMenu(const QPoint &at);

 private:
  void buildForm(const azure2::front::Record &record);
  void clearForm();
  void addFieldRows(const azure2::front::KindInfo &kind, const azure2::front::Record &record, QFormLayout *form,
                    const azure2::front::FieldInfo &field, bool inBox);
  FieldEditor *makeEditor(const azure2::front::FieldInfo &field, const QString &text);
  QString uniqueName(const QString &stem) const;
  QString quantityText() const;
  void quantityChanged();
  void updateButtons();
  void showUsedBy();
  /// Draw the record's kernels (Doppler, beam spread, energy resolution, beam profile) as the form holds it, from the
  /// engine's own (front::KernelPreviewOf).
  void updatePreview();
  /// Show the fields the record's quantity uses and hide the others (front::QuantityFieldRules).
  void applyQuantityRules();
  bool quantityEdited() const;
  /// The record as the form holds it.
  azure2::front::Record formRecord() const;
  bool formEdited() const;
  /// The engine's errors of the project with `record` written, by field ("" for the record).
  QMap<QString, QStringList> errorsWith(const azure2::front::Record &record) const;

  ProjectHolder *holder_;
  QString section_;
  QStringList kinds_;
  RecordTableModel *model_;
  QSortFilterProxyModel *proxy_;
  QTableView *table_;
  QLineEdit *filter_;
  QPushButton *addButton_, *removeButton_, *duplicateButton_, *upButton_, *downButton_, *applyButton_, *revertButton_;
  QPushButton *includeButton_, *excludeButton_;
  QAction *addAction_, *removeAction_, *duplicateAction_, *upAction_, *downAction_;
  QCheckBox *everyField_;
  QWidget *formHost_;
  QVBoxLayout *formLayout_ = nullptr;
  QLineEdit *nameEdit_ = nullptr;
  QLabel *kindLabel_;
  QLabel *usedBy_;
  QPushButton *assignButton_;
  QLabel *status_;
  QWidget *preview_ = nullptr;
  azure2::front::KernelPreview kernels_;
  QString previewUser_;
  QString previewWhy_;
  double previewEnergy_ = 0;  ///< lab MeV; 0: the default
  QWidget *formScroll_ = nullptr;
  QWidget *bottom_ = nullptr;     ///< the form's half of the splitter (hidden while the section has no record)
  QLabel *emptyHint_ = nullptr;   ///< "No record yet: Add..." under the empty table
  /// The rows of the fields the quantity decides: key -> (form, label, field widget).
  struct QuantityRow {
    QFormLayout *form = nullptr;
    QLabel *label = nullptr;
    QWidget *row = nullptr;
    azure2::front::FieldInfo field;
  };
  QMap<QString, QuantityRow> quantityRows_;
  QStringList hiddenByQuantity_;
  QMap<QString, QList<FieldEditor *>> editors_;  ///< field key -> its editors (several for a repeatable field)
  QMap<QString, QWidget *> boxes_;               ///< field key -> its box (a structured field and its group)
  QList<QWidget *> rareRows_;                    ///< rows of rare fields the record does not set
  azure2::front::Record shown_;                   ///< the record as last read
  bool haveShown_ = false;
  bool reloading_ = false;
  bool previewPending_ = false;
  QString quantityKey_;
};

#endif
