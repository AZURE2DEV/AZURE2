#ifndef THMEXPERIMENTSPAGE_H
#define THMEXPERIMENTSPAGE_H

#include <QDialog>
#include <QList>
#include <QStringList>
#include <QWidget>
#include <QVector>
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
class ThmNumberSpin;
class ThmDistortion;
struct ThmSpectatorWindow;

/*!
 * The Experiments page of the THM workspace: the experiment[<name>] lines of
 * the <thm> block (docs/source/theory/thm_implementation.rst, "THM
 * experiments").  A compact list of the experiments on the left and an
 * editor of the selected one on the right, in sections: Experiment (name,
 * segments -- THM data segments with a free norm that no other experiment
 * has --, background, the exit angle: angle-integrated or a theta= window),
 * Three-body reaction (beam, target, spectator, lab
 * beam energy -- all four or none -- with the binding and quasi-free
 * energies AZURE2 prints for it, and the line shape, lineshape=on, with zeta
 * at the ends of the data) and Spectator momentum window (ps=, psNodes=, with the
 * mean spectator energy <T_s> the engine's ThmSpectatorWindow gives) and
 * Distortion (distortion= and its keys, with R(E) at the ends of the data
 * from the engine's ThmDistortion).  When B(x+s) from the masses and the
 * entrance pair's B (field 32) disagree, the reaction section says so.
 * Numbers are spin boxes with their unit that give back the text they were
 * read with until changed (ThmNumberSpin).  Keys the page does not show are
 * kept as written.
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
  QString lineshapeInfo(const ThmExperimentRecord &record, QString *error = nullptr,
                        QString *compact = nullptr) const;
  /// Everything derivedInfo gives for the selected experiment, and the reason
  /// AZURE2 would refuse it, as the tooltips of the page's values show it.
  QString derivedText() const { return derivedText_; }
  /*! For a record with a ps window and a complete reaction: the window as
      the engine builds it (BuildThmSpectatorWindow on the engine's parse of
      the line, the table read by ReadThmPsTable, mu_sx from the reaction),
      described with <T_s>; or "" and the reason in *error when AZURE2 would
      refuse it (a bad table, a spectator energy for the same pair). */
  QString windowInfo(const ThmExperimentRecord &record, QString *error = nullptr,
                     ThmSpectatorWindow *window = nullptr) const;
  /*! For a record with a distortion: R(E) (coulomb, optical: the engine's
      ThmDistortion set up as EData::BuildThmGroups does, evaluated directly)
      or w(E) (a table, read by ThmWeightTable) at the lowest and highest point
      energy of its segments, and the whole text; "" and the reason in
      *error, in the engine's words, when AZURE2 would refuse it. */
  QString distortionInfo(const ThmExperimentRecord &record, QString *error = nullptr, double *lo = nullptr,
                         double *hi = nullptr, double *rLo = nullptr, double *rHi = nullptr) const;
  /*! When B(x+s) from the masses of a record's reaction differs from its
      entrance pair's B (field 32) by more than 1 keV, as the engine warns:
      the one-line note (*detail: the engine's whole warning); else "". */
  QString bindingMismatch(const ThmExperimentRecord &record, QString *detail = nullptr) const;
  /// The same for every experiment whose entrance pair is `pairKey`, with
  /// that pair's B taken as `pairB` (the Channels page).  "" if none.
  QString bindingMismatchOfPair(int pairKey, double pairB) const;
  /// The entrance pair's B (field 32) per pair key; the Channels page's
  /// values in the workspace (default: the project's).
  void setPairBinding(std::function<double(int pairKey)> binding) { pairBinding_ = binding; }
  /// Sets opticalAA= (channel 0) or opticalSF= (1) of the selected
  /// experiment to ten numbers, as the Woods-Saxon dialog does on OK.
  void setOpticalText(int channel, const QString &tenNumbers);
  /// Sets opticalAA= / opticalSF= of the selected experiment to a global
  /// optical potential, `name` or `name:extrapolate` (ThmOptical.h).
  void setGlobalText(int channel, const QString &value);
  /*! For a record whose channel (0 a + A, 1 s + F) has a global optical
      potential and a complete reaction: its ten numbers at the lowest and
      highest point energy (the projectile lab energies in elab; for a + A
      both ends are E_aA), as the engine evaluates them; false if they cannot
      be evaluated (no reaction, no data, a projectile the model does not
      describe). */
  bool globalEnds(const ThmExperimentRecord &record, int channel, double elab[2], double p[2][10],
                  double *lo = nullptr, double *hi = nullptr) const;
  /// The ps= value the spectator-momentum controls describe ("" for a point).
  QString psText() const;
  /// A file chosen for ps=table: relative to the project directory when inside it.
  static QString projectRelative(const QString &file, const QString &projectDir);
  /// The spectator energy per entrance pair (the Model page's values); a ps
  /// window is refused together with a non-zero one for its pair.
  void setSpectatorEnergy(std::function<double(int pairKey)> energy) { spectatorEnergy_ = energy; }
  /// entranceL of the <thm> block (the Model page's value); a theta window is
  /// refused together with entranceL=coherent.
  void setEntranceL(std::function<QString()> entranceL) { entranceL_ = entranceL; }
  /// coulombIntegral of the <thm> block (the Model page's value); refused with
  /// vertexModel=dw and with a computed R(E) whose a + A wave is distorted.
  void setCoulombIntegral(std::function<bool()> coulombIntegral) { coulombIntegral_ = coulombIntegral; }
  /// The engine's refusal of a theta window with entranceL=coherent (EData::BuildThmGroups).
  static QString coherentRefusal();
  /// The theta= value the exit-angle controls describe ("all" for angle-integrated).
  QString thetaText() const;
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
  QWidget *editorBox;
  QLineEdit *nameEdit;
  QListWidget *segmentList;  ///< checkable; item data Qt::UserRole = segment number
  QComboBox *backgroundCombo;
  /// theta=: the exit angle, item data all | window; the window's ends in degrees.
  QComboBox *thetaCombo;
  ThmNumberSpin *thetaMinEdit, *thetaMaxEdit;
  QGroupBox *kinematicsBox;  ///< checkable: the four keys as a unit
  QComboBox *beamCombo;
  QComboBox *targetCombo;
  QComboBox *spectatorCombo;
  ThmNumberSpin *beamEnergyEdit;  ///< MeV; the minimum (shown as a dash) = not given
  QCheckBox *lineshapeCheck;  ///< lineshape=on; enabled with a complete reaction
  QGroupBox *psBox;          ///< spectator momentum; enabled with a complete reaction
  QComboBox *psKindCombo;    ///< item data: delta | hulthen | gauss | table
  ThmNumberSpin *psMinEdit, *psMaxEdit;  ///< MeV/c
  QCheckBox *psCustomCheck;  ///< Hulthen a, b other than the deuteron's
  ThmNumberSpin *psAEdit, *psBEdit;  ///< fm^-1
  ThmNumberSpin *psFwhmEdit;         ///< MeV/c
  QLineEdit *psTableEdit;
  QPushButton *psTableButton;
  QSpinBox *psNodesSpin;     ///< psNodes=, 1-64, default 16
  /// Distortion (distortion= ...): the kind, and for coulomb/optical the
  /// spectator angle, E_ref, the ratio, the bound state and (optical) the
  /// two channels; for a table its file.
  QGroupBox *distortionBox;
  QComboBox *distortionCombo;   ///< item data: none | coulomb | optical | table
  QComboBox *angleKindCombo;    ///< item data: qf | lab | cm
  ThmNumberSpin *angleEdit;     ///< deg
  ThmNumberSpin *distortionRefEdit;  ///< MeV; the minimum (shown as "auto") = not given
  QComboBox *ratioCombo;        ///< item data: dwpw | dw
  QComboBox *boundCombo;        ///< item data: whittaker | yukawa
  ThmNumberSpin *rminEdit;      ///< fm; the minimum (a dash) = not given
  QComboBox *opticalCombo[2];   ///< a + A, s + F; item data: plane | coulomb | global | ws
  QComboBox *globalCombo[2];    ///< the global potential (item data: its name); shown for global
  QPushButton *opticalButton[2];  ///< "Edit..." the ten Woods-Saxon numbers, or the global potential's
  QLineEdit *distortionTableEdit;
  QPushButton *distortionTableButton;
  /// Derived values (the reaction, zeta, the window, R at the ends), and why AZURE2 would refuse the experiment.
  QLabel *bindingValue, *qfValue, *zetaValue, *meanTsValue, *distortionValue;
  /// B(x+s) from the masses differs from the pair's B: icon and one line in the reaction section.
  QLabel *bindingWarningIcon, *bindingWarningLabel;
  QLabel *messageLabel;
  QLabel *messageIcon;

 private slots:
  void tableSelectionChanged();
  void nameEdited(const QString &text);
  void segmentItemChanged(QListWidgetItem *item);
  void backgroundChanged(int index);
  void thetaEdited();
  void kinematicsEdited();
  void lineshapeToggled(bool on);
  void psEdited();
  void psNodesChanged(int n);
  void chooseTable();
  void distortionKindChanged();
  void distortionEdited();
  void opticalKindChanged();
  void globalNameChanged();
  void chooseDistortionTable();
  void editOptical(int channel);

 private:
  struct Reaction {
    ThmNuclide beam, target, spectator, horse;
    int pairKey = 0;
    double beamEnergy = 0.0, bind = 0.0, exa = 0.0;
    double mX = 0.0;  ///< mass of x = Trojan horse - spectator (u), as the engine takes it
    bool horseIsBeam = true;
    ThmNuclide other;  ///< A, the nucleus that is not the Trojan horse
  };
  /// The reaction of a record with all four keys, as EData::SetupThmExperiments checks it.
  bool reaction(const ThmExperimentRecord &x, Reaction &out, QString *error) const;
  void showDerived(const ThmExperimentRecord &r);
  QString derivedText_;
  void loadEditor();
  void loadPs(const ThmExperimentRecord &r);
  void loadDistortion(const ThmExperimentRecord &r);
  void showDistortionRows();
  void updateDistortionItems(bool complete);
  /// The c.m. energies of the points of the segments, in the engine's order.
  bool pointEnergies(const QList<int> &segments, QVector<double> &energies) const;
  void showPsRows();
  void loadTheta(const ThmExperimentRecord &r);
  void showThetaRows();
  /// The window item is not offered with entranceL=coherent (unless it is the current one).
  void updateThetaItems();
  QList<QWidget *> thetaWindowRow_;
  std::function<QString()> entranceL_ = []() { return QString("incoherent"); };
  std::function<bool()> coulombIntegral_ = []() { return false; };
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
  QList<QWidget *> distortionComputedRows_, distortionOpticalRow_, distortionTableRow_, distortionValueRow_;
  QString lastOptical_[2];  ///< the ten numbers last shown per channel (kept across plane/coulomb)
  QString lastGlobal_[2];   ///< the global potential last shown per channel (name[:extrapolate])
  /// d.kin (and d.eAA) of a record with a complete reaction, as EData::BuildThmGroups sets them.
  bool fillKinematics(const ThmExperimentRecord &record, ThmDistortion &d) const;
  /// The Edit button's tooltip for a global potential: its numbers at the data ends.
  QString globalSummary(const ThmExperimentRecord &record, int channel) const;
  std::function<double(int)> pairBinding_;
  /// distortionInfo of the last line asked for (the setup takes milliseconds to a second).
  mutable QString cacheKey_, cacheText_, cacheError_;
  mutable double cache_[4] = {0.0, 0.0, 0.0, 0.0};
  std::function<double(int)> spectatorEnergy_ = [](int) { return 0.0; };
  int current_ = -1;
  bool loading_ = false;
};

/*!
 * The ten numbers of a Woods-Saxon optical potential (opticalAA= /
 * opticalSF= V,R,a,W,RW,aW,WD,RD,aD,RC) in a compact form: real volume,
 * imaginary volume, imaginary surface (depth, radius, diffuseness) and the
 * Coulomb radius, each a spin box with its unit.  The text gives back each
 * number as it was written until it is changed.
 */
class ThmOpticalDialog : public QDialog {
 public:
  explicit ThmOpticalDialog(const QString &title, const QString &tenNumbers, QWidget *parent = nullptr);
  /// V,R,a,W,RW,aW,WD,RD,aD,RC as the fields give them.
  QString text() const;
  ThmNumberSpin *fields[10];
};

/*!
 * A global optical potential of a channel (opticalAA= / opticalSF=
 * <name>[:extrapolate]): its reference and validity range, the ten numbers
 * it gives at the two ends of the data (read only), and whether it may be
 * used outside its validity range.
 */
class ThmGlobalOpticalDialog : public QDialog {
 public:
  /// `value` is name[:extrapolate]; `elab`/`p` the ends (have = false: not available).
  ThmGlobalOpticalDialog(const QString &title, const QString &value, bool have, const double elab[2],
                         const double p[2][10], double lo, double hi, QWidget *parent = nullptr);
  /// name or name:extrapolate.
  QString text() const;
  QCheckBox *extrapolateCheck;
  QLabel *valueLabels[2][10];

 private:
  QString name_;
};

#endif
