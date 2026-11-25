#ifndef POTENTIALPANEL_H
#define POTENTIALPANEL_H

#include <functional>
#include <vector>

#include <QList>
#include <QString>
#include <QWidget>

#include "azure2/front/Engine.h"

class PhysicsWorker;
class ProjectHolder;

QT_BEGIN_NAMESPACE
class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QTableWidget;
class QTimer;
QT_END_NAMESPACE

/*!
 * The Potential panel of the Particle Pairs tab (docs/upgrade/qt_gui.md section 19): the `potential` record of the pair
 * selected in the tab, wave by wave, with what the potential creates and what it gives by itself.
 *
 *  - on/off: the record switched off is kept in the file as comment lines (Project::SetRecordOn), so that switching it
 *    on restores it byte for byte;
 *  - a table of the partial waves: every `wave` line of the record and a row per l without one, each with its mode
 *    (central potential, exterior only, exterior + interior = `background interior`, none = `only bound_state`), its
 *    depth, width, B_bg, `subtract n`, `capture interior`, and whether the wave is in the calculation at all (an active
 *    level of its J^pi, or includeAllJPi: front::PartialWaveStates);
 *  - the poles: the box states of the selected wave at the channel radius (front::PotentialWaveOf: E_p, E_x, gamma_p^2,
 *    theta^2, nodes, the subtracted one marked) beside the project's levels of the same J^pi;
 *  - tuning: the wave's own a and r and bound states (read out as it is edited), the depth for a chosen a (the branch of
 *    the current number of bound states), the depth and Gaussian width for a chosen a and r (the same branch), the depth
 *    that binds a level below threshold with the potential's ANC beside the level's.
 * The engine's messages of the edited record (Validate of a copy, the set-up's checks of the wave) are shown as they
 * come.  Apply writes the record with Project::SetRecord (an unchanged wave keeps its text).  Every computation runs
 * in a worker thread (PhysicsWorker); the window never waits for it.
 */
class PotentialPanel : public QWidget {
  Q_OBJECT

 public:
  /// What the panel needs of the window: the Levels tab's unsaved state and how to write the Pairs and Levels tabs'
  /// rows into the project before the engine reads it.
  struct Hooks {
    std::function<bool(QString &why)> commitLevels;
    std::function<std::vector<azure2::front::LevelRow>()> levelRows;
    std::function<bool()> includeAllJPi;
    std::function<int()> maxOrbitalMomentum;  ///< the key's l as the engine will read it (-1: the rows')
    std::function<int()> maxL;                ///< the Levels tab's maximum l (the waves listed without a line)
  };
  /// The modes of a wave (the Mode column's items, in order).
  enum Mode { Central = 0, Exterior = 1, ExteriorInterior = 2, PointCoulomb = 3 };

  PotentialPanel(ProjectHolder *holder, Hooks hooks, QWidget *parent = nullptr);
  ~PotentialPanel() override;

  /// Show the potential of pair `pairKey` (<levels> key; 0 none), `label` in nuclear notation.
  void setPair(int pairKey, const QString &label);
  int pair() const { return pair_; }
  /// Read the record again (after the project changed elsewhere); the edits not applied are dropped.
  void reload();

  /// The record as the panel holds it (the waves as edited).
  azure2::front::Record editedRecord() const;
  bool hasRecord() const { return hasRecord_; }
  bool recordOn() const { return on_; }
  /// Write the edited record into the project; false with the reason.
  bool apply(QString &why);
  /// A computation is running.
  bool busy() const;

  // ---- the widgets (the tests drive them) ----
  QCheckBox *onCheck() const { return onCheck_; }
  QTableWidget *wavesTable() const { return waves_; }
  QTableWidget *polesTable() const { return poles_; }
  QTableWidget *levelsTable() const { return levels_; }
  QLabel *messagesLabel() const { return messages_; }
  QLabel *readOutLabel() const { return readOut_; }
  QLabel *waveTitle() const { return waveTitle_; }
  QLineEdit *lengthEdit() const { return lengthEdit_; }
  QPushButton *lengthButton() const { return lengthButton_; }
  QLineEdit *rangeAEdit() const { return rangeA_; }
  QLineEdit *rangeREdit() const { return rangeR_; }
  QPushButton *rangeButton() const { return rangeButton_; }
  QComboBox *bindCombo() const { return bindCombo_; }
  QPushButton *bindButton() const { return bindButton_; }
  QLabel *tuningLabel() const { return tuning_; }
  QPushButton *useDepthButton() const { return useDepth_; }
  QPushButton *applyButton() const { return apply_; }
  QPushButton *revertButton() const { return revert_; }
  QPushButton *addWaveButton() const { return addWave_; }
  QPushButton *createButton() const { return create_; }
  /// The last report of the selected wave (valid when `reportOk()`), the last depth found.
  const azure2::front::PotentialWaveReport &report() const { return report_; }
  bool reportOk() const { return reportOk_; }
  const azure2::front::PotentialDepth &depthFound() const { return depth_; }
  /// The partial waves of the table's row `row` as the engine will take them.
  std::vector<azure2::front::PartialWaveState> waveStates(int row) const;

  /// Set a row's mode (as the Mode column's drop-down does).
  void setMode(int row, Mode mode);
  /// Add a wave line (l, J; "" every J) as "Add wave..." does; returns its row.
  int addWave(int l, const QString &J);

 signals:
  /// A computation's answer is on screen (the selected wave's report, a depth, the record's messages).
  void computed();

 public slots:
  /// Compute the selected wave again (after an edit, a change of the levels).
  void refresh();
  /// The Levels tab changed (a level, the limits, includeAllJPi): the waves' states and the selected wave again.
  void levelsChanged();

 private:
  struct Row {
    bool line = false;         ///< a `wave` field of the record (else: the central potential in this l)
    int field = -1;            ///< its place among the record's fields (-1: a new line)
    QString original;          ///< the field's text as read ("" for a new line)
    bool edited = false;
    azure2::front::PotentialWave wave;
  };
  void buildRows();
  void fillTable();
  void fillRow(int row);
  int currentRow() const;
  int requestTwoJ(const Row &row) const;
  QString rowText(const Row &row) const;
  void cellEdited(int row, int column);
  void modeChosen(int row, int index);
  void captureToggled(int row, bool on);
  void switchOn(bool on);
  void showReport();
  void showMessages();
  void findLength();
  void findBinding();
  void findRange();
  void useDepth();
  void createPotential();
  void editRecord();
  void addWaveAsked();
  void removeWaveLine();
  void revert();
  void schedule();
  void updateStates();
  QString jpiText(int twoJ, int parity) const;
  /// The panel's help (InTabDocs.cpp), shown by its "?" button.
  static const QString helpText;

  ProjectHolder *holder_;
  Hooks hooks_;
  PhysicsWorker *worker_;
  QTimer *debounce_;
  int pair_ = 0;
  QString pairLabel_;
  bool hasRecord_ = false, on_ = true;
  azure2::front::Record record_;  ///< as read
  QList<Row> rows_;
  bool filling_ = false;
  int reportRow_ = -1;
  azure2::front::PotentialWaveReport report_;
  bool reportOk_ = false;
  QString reportWhy_;
  QStringList recordMessages_;
  azure2::front::PotentialDepth depth_;
  int depthRow_ = -1;
  bool depthWithWidth_ = false;  ///< the depth found comes with a width (depth and range)
  int ticket_ = 0;

  QLabel *title_, *help_, *messages_, *readOut_, *tuning_, *waveTitle_, *busy_;
  QCheckBox *onCheck_;
  QPushButton *create_, *edit_, *addWave_, *removeWave_, *apply_, *revert_;
  QTableWidget *waves_, *poles_, *levels_;
  QDoubleSpinBox *eMax_;
  QLineEdit *lengthEdit_, *rangeA_, *rangeR_;
  QPushButton *lengthButton_, *bindButton_, *useDepth_, *rangeButton_;
  QComboBox *bindCombo_;
};

#endif
