#ifndef THMEXPERIMENTSPAGE_H
#define THMEXPERIMENTSPAGE_H

#include <QList>
#include <QStringList>
#include <QWidget>

#include "ThmSettings.h"

QT_BEGIN_NAMESPACE
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QListWidgetItem;
class QPushButton;
class QTableWidget;
QT_END_NAMESPACE

class PairsModel;
class SegmentsDataModel;

/*!
 * The Experiments page of the THM workspace: the experiment[<name>] lines of
 * the <thm> block (docs/source/theory/thm_implementation.rst, "THM
 * experiments").  A table of the experiments and an editor of the selected
 * one: its name, its segments (THM data segments with a free norm that no
 * other experiment has), the background, and the optional three-body
 * reaction (beam, target, spectator, lab beam energy -- all four or none),
 * with the binding and quasi-free energies AZURE2 prints for it.  Keys the
 * page does not show are kept as written.
 */
class ThmExperimentsPage : public QWidget {
  Q_OBJECT

 public:
  ThmExperimentsPage(const QStringList &experimentLines, SegmentsDataModel *segments, PairsModel *pairs,
                     QWidget *parent = 0);

  /// The experiment lines to write (ThmExperimentRecord::compose): the
  /// lines read, verbatim, if nothing changed.
  QStringList experimentLines() const;
  const QList<ThmExperimentRecord> &records() const { return records_; }
  /*! What AZURE2 would refuse, with its message: the engine's parser and
      check on the lines, then the startup checks against the project
      (segment exists, is THM, has a free norm; the reaction gives the
      entrance pair of the segments).  "" if none. */
  QString check() const;

  /// The nuclide names of the engine's built-in table, in (Z, A) order.
  static QStringList nuclideNames();
  /*! The binding and quasi-free energies AZURE2 prints for a record with
      complete kinematics, or "" (and the reason in *error when the engine
      would refuse the reaction, e.g. it does not give the entrance pair). */
  QString derivedInfo(const ThmExperimentRecord &record, QString *error = nullptr) const;

  /// Editing, as the buttons and the editor do it (for the tests too).
  void addExperiment();
  void removeCurrent();
  void selectExperiment(int row);
  int currentRow() const { return current_; }
  /// Sets the segments of the selected experiment as given, eligible or not.
  void setSegmentsOfCurrent(const QList<int> &segments);

  QTableWidget *experimentTable;
  QPushButton *addButton;
  QPushButton *removeButton;
  QGroupBox *editorBox;
  QLineEdit *nameEdit;
  QListWidget *segmentList;  ///< checkable; item data Qt::UserRole = segment number
  QComboBox *backgroundCombo;
  QGroupBox *kinematicsBox;  ///< checkable: the four keys as a unit
  QComboBox *beamCombo;
  QComboBox *targetCombo;
  QComboBox *spectatorCombo;
  QLineEdit *beamEnergyEdit;
  QLabel *derivedLabel;

 private slots:
  void tableSelectionChanged();
  void nameEdited(const QString &text);
  void segmentItemChanged(QListWidgetItem *item);
  void backgroundChanged(int index);
  void kinematicsEdited();

 private:
  void loadEditor();
  void fillSegmentList();
  void storeSegments(const QList<int> &segments);
  void refreshRow(int row);
  void refreshTable();
  QString segmentLabel(int key) const;

  SegmentsDataModel *segments_;
  PairsModel *pairs_;
  QStringList oldLines_;
  QList<ThmExperimentRecord> oldRecords_;
  QList<ThmExperimentRecord> records_;
  int current_ = -1;
  bool loading_ = false;
};

#endif
