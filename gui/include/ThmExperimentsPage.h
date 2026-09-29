#ifndef THMEXPERIMENTSPAGE_H
#define THMEXPERIMENTSPAGE_H

#include <QList>
#include <QStringList>
#include <QWidget>
#include <functional>

#include "ThmExperiment.h"
#include "ThmSettings.h"

QT_BEGIN_NAMESPACE
class QCheckBox;
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QListWidgetItem;
class QPushButton;
class QSpinBox;
class QTableWidget;
QT_END_NAMESPACE

class PairsModel;
class SegmentsDataModel;
struct ThmSpectatorWindow;

/*!
 * The Experiments page of the THM workspace: the experiment[<name>] lines of
 * the <thm> block (docs/source/theory/thm_implementation.rst, "THM
 * experiments").  A table of the experiments and an editor of the selected
 * one: its name, its segments (THM data segments with a free norm that no
 * other experiment has), the background, and the optional three-body
 * reaction (beam, target, spectator, lab beam energy -- all four or none),
 * with the binding and quasi-free energies AZURE2 prints for it, and the
 * Coulomb line shape of the spectator (lineshape=on) with zeta at the ends of
 * the data, and the spectator-momentum window (ps=, psNodes=) with the mean
 * spectator energy <T_s> the engine's ThmSpectatorWindow gives.  Keys the
 * page does not show are kept as written.
 */
class ThmExperimentsPage : public QWidget {
  Q_OBJECT

 public:
  /// `projectDir`: where relative data file names resolve; `brune`: the
  /// project uses the Brune parameterization (lineshape=on needs it).
  ThmExperimentsPage(const QStringList &experimentLines, SegmentsDataModel *segments, PairsModel *pairs,
                     const QString &projectDir = QString(), bool brune = true, QWidget *parent = 0);

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
  /*! For a record with lineshape=on and a complete reaction: zeta of every
      exit pair at the lowest and highest point energy of its segments (the
      engine's ThmLineshape on the engine's data reader), or "" and the
      reason in *error when AZURE2 would refuse it (no Brune, E_sF <= 0). */
  QString lineshapeInfo(const ThmExperimentRecord &record, QString *error = nullptr) const;
  /*! For a record with a ps window and a complete reaction: the window as
      the engine builds it (BuildThmSpectatorWindow on the engine's parse of
      the line, the table read by ReadThmPsTable, mu_sx from the reaction),
      described with <T_s>; or "" and the reason in *error when AZURE2 would
      refuse it (a bad table, a spectator energy for the same pair). */
  QString windowInfo(const ThmExperimentRecord &record, QString *error = nullptr,
                     ThmSpectatorWindow *window = nullptr) const;
  /// The ps= value the spectator-momentum controls describe ("" for a point).
  QString psText() const;
  /// A file chosen for ps=table: relative to the project directory when inside it.
  static QString projectRelative(const QString &file, const QString &projectDir);
  /// The spectator energy per entrance pair (the Model page's values); a ps
  /// window is refused together with a non-zero one for its pair.
  void setSpectatorEnergy(std::function<double(int pairKey)> energy) { spectatorEnergy_ = energy; }
  /// Recomputes the derived text of the selected experiment (other pages changed).
  void refreshDerived();
  /// The lowest and highest c.m. point energy of the given data segments, as
  /// ESegment::FillData reads them; false if a file cannot be read.
  bool pointRange(const QList<int> &segments, double &lo, double &hi) const;

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
  QCheckBox *lineshapeCheck;  ///< lineshape=on; enabled with a complete reaction
  QGroupBox *psBox;          ///< spectator momentum; enabled with a complete reaction
  QComboBox *psKindCombo;    ///< item data: delta | hulthen | gauss | table
  QLineEdit *psMinEdit, *psMaxEdit;  ///< MeV/c
  QCheckBox *psCustomCheck;  ///< Hulthen a, b other than the deuteron's
  QLineEdit *psAEdit, *psBEdit;      ///< fm^-1
  QLineEdit *psFwhmEdit;     ///< MeV/c
  QLineEdit *psTableEdit;
  QPushButton *psTableButton;
  QSpinBox *psNodesSpin;     ///< psNodes=, 1-64, default 16
  QLabel *derivedLabel;

 private slots:
  void tableSelectionChanged();
  void nameEdited(const QString &text);
  void segmentItemChanged(QListWidgetItem *item);
  void backgroundChanged(int index);
  void kinematicsEdited();
  void lineshapeToggled(bool on);
  void psEdited();
  void psNodesChanged(int n);
  void chooseTable();

 private:
  struct Reaction {
    ThmNuclide beam, target, spectator, horse;
    int pairKey = 0;
    double beamEnergy = 0.0, bind = 0.0, exa = 0.0;
    double mX = 0.0;  ///< mass of x = Trojan horse - spectator (u), as the engine takes it
  };
  /// The reaction of a record with all four keys, as EData::SetupThmExperiments checks it.
  bool reaction(const ThmExperimentRecord &x, Reaction &out, QString *error) const;
  void showDerived(const ThmExperimentRecord &r);
  void loadEditor();
  void loadPs(const ThmExperimentRecord &r);
  void showPsRows();
  void fillSegmentList();
  void storeSegments(const QList<int> &segments);
  void refreshRow(int row);
  void refreshTable();
  QString segmentLabel(int key) const;

  SegmentsDataModel *segments_;
  PairsModel *pairs_;
  QString projectDir_;
  bool brune_;
  QStringList oldLines_;
  QList<ThmExperimentRecord> oldRecords_;
  QList<ThmExperimentRecord> records_;
  QList<QWidget *> psWindowRow_, psHulthenRow_, psGaussRow_, psTableRow_, psNodesRow_;
  std::function<double(int)> spectatorEnergy_ = [](int) { return 0.0; };
  int current_ = -1;
  bool loading_ = false;
};

#endif
