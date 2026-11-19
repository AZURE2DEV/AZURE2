#ifndef LEVELSTAB_H
#define LEVELSTAB_H

#include <QWidget>
#include <QMap>
#include <QSet>
#include <QStringList>
#include "azure2/front/Engine.h"
#include <QTabWidget>
#include <QItemSelection>
#include <QLineEdit>
#include <QSpinBox>
#include <QCheckBox>
#include <QLabel>
#include <QPushButton>
#include <QTableView>
#include <QSortFilterProxyModel>
#include <QSignalMapper>
#include <QPointer>
#include <QFileDialog>
#include <QClipboard>
#include <QApplication>
#include <QVector>
#include "PairsModel.h"
#include "LevelsModel.h"
#include "ChannelsModel.h"
#include "ChannelDetails.h"
#include "AddLevelDialog.h"
#include "LevelsUndo.h"
#include "LevelsViews.h"

class InfoDialog;
class ProjectHolder;
class QUndoStack;
class QComboBox;

/*!
 * The Levels tab: the level scheme, and the channels each J^pi group carries. Channels are derived from the levels and the pair list rather than entered directly.
 */
class LevelsTab : public QWidget {
  Q_OBJECT

 public:
  LevelsTab(QWidget *parent = 0);
  void setPairsModel(PairsModel *);
  void updateChannelsLevelAdded(int levelIndex);
  void updateChannelsLevelDeleted(int levelIndex);
  void updateChannelsLevelEdited(int levelIndex);
  QList<ChannelsData> calculateChannels(int levelIndex);
  QList<ChannelsData> calculateChannels(const LevelsData &level, int levelIndex);
  /*! The J^pi that give at least one channel under the current pairs and limits (every pair). */
  QList<SpinParityOption> spinParityOptions();
  /*! The same for the pairs of one compound nucleus ((0, 0) = every pair). */
  QList<SpinParityOption> spinParityOptions(const Nucleus &nucleus);
  /*! The <levels> rows as the engine read them (azure2::front::Project::Levels) into the pairs, levels and channels;
   *  false with `why` for a row the GUI cannot place. */
  bool readLevels(const std::vector<azure2::front::LevelRow> &rows, QString &why);
  /*! The <levels> rows of the pairs, levels and channels, in writtenLevelOrder, a blank row after each level; a
   *  value the GUI did not change keeps the token it was read with (LevelsFile.cpp). */
  std::vector<azure2::front::LevelRow> levelRows();
  /*! A panel shown below the levels (gamma decays, isospin; AZURESetup). */
  void addExtraPanel(const QString &title, QWidget *panel);
  /*! Model rows of the levels in the order writeNuclearFile writes them: compound nuclei in order of first
   *  appearance, then J, parity and energy. */
  QList<int> writtenLevelOrder();
  /*!
   * How AZURE2 numbers the parameters of these levels when it reads a file that lists them in `order` (model
   * rows), following CNuc::Fill and CNuc::FillMnParams: only active levels exist; J-groups (nucleus, J, parity)
   * are numbered in order of first appearance and an active level's energy is energy_N, N counting the levels
   * group by group; a width is width_N_M, M the channel's position among the channels (pair, l, s) of the
   * level's J-group in order of first appearance.  level[row] = N (0: inactive), channel[channel row] = M
   * (0: the level is inactive).  The Fitting tab names its "Level N ..." parameters, and so the
   * <parameterSettings> block, in this numbering, which is the one the engine reads them in.
   */
  struct EngineNumbering {
    QVector<int> level;
    QVector<int> channel;
  };
  EngineNumbering engineNumbering(const QList<int> &order);
  void reset();
  /// Select level `row` (of the levels model) and remove it as the Remove button does (the tests).
  void removeLevelAt(int row);
  LevelsModel *getLevelsModel() { return levelsModel; }
  /*! Hard-sphere scattering in every J^pi up to the maximum orbital momentum, levels or not
   *  (includeAllJPi in the .azr <config> block). On for a new project, off for a file that lacks the key. */
  bool includeAllJPi() const { return includeAllJPiCheck->isChecked(); }
  void setIncludeAllJPi(bool on) { includeAllJPiCheck->setChecked(on); }
  int maxL() const { return maxLSpin->value(); }
  void setMaxL(int l) { maxLSpin->setValue(l); }
  /*! The J^pi the pairs allow up to the maximum orbital momentum that no active level carries. */
  QStringList missingSpinParities(QStringList *links = nullptr);
  ChannelsModel *getChannelsModel() { return channelsModel; }
  /*! The <levels> row the tab writes for channel `ch` of `level` (LevelsFile.cpp), with level ID `id`. */
  azure2::front::LevelRow channelRow(const ChannelsData &ch, const LevelsData &level, const PairsData &pair,
                                     const PairsData &entrance, int id) const;

  // ---- on the go (docs/upgrade/qt_gui.md section 11; LevelsTabEdits.cpp, LevelsTabNavigation.cpp) ----
  /*! The undo stack of this tab and the Particle Pairs tab. */
  QUndoStack *undoStack() const { return undo_; }
  QTableView *getLevelsView() const { return levelsView; }
  QTableView *getChannelsView() const { return channelsView; }
  LevelsFilterProxy *getLevelsProxy() const { return levelsModelProxy; }
  /*! The selected levels (model rows, in the order shown) and the current one (-1: none). */
  QList<int> selectedLevelRows() const;
  int currentLevelRow() const;
  /*! Select these levels (model rows; the first becomes current), show the first. */
  void selectLevelRows(const QList<int> &rows);
  /*! The add bar's fields and its button (the tests drive them). */
  HalfSpinBox *addJSpin() const { return addJ_; }
  ParityToggle *addParityToggle() const { return addPi_; }
  QLineEdit *addEnergyEdit() const { return addE_; }
  QPushButton *addBarButton() const { return addBarButton_; }
  QLineEdit *searchEdit() const { return search_; }
  QComboBox *jpiFilterCombo() const { return jpiFilter_; }
  QCheckBox *onlyIncludedCheck() const { return onlyIncluded_; }
  QLabel *statusLabel() const { return status_; }
  QLabel *missingLabel() const { return missingJPiLabel; }
  /*! Edits the tab checks (in place, the add bar): false with the reason, the models unchanged. */
  bool addLevelFromBar(QString &why);
  bool setLevelEnergy(int row, const QString &text, QString &why);
  bool setLevelJpi(int row, double j, int parity, QString &why);
  bool setChannelValue(int row, const QString &text, QString &why);
  /*! Remove these levels (model rows) with their channels; `ask`: after one confirmation naming them. */
  bool removeLevels(QList<int> rows, bool ask);
  /*! Fix or free the energies (or the widths of every channel) of these levels. */
  void setFixed(const QList<int> &rows, bool energies, bool fixed);
  /*! Record the state as an undo entry now if it changed (the event loop does it after every change). */
  void recordUndo();
  /*! The quantities of a channel (model row) at its level's energy, from the facade (azure2::front::ChannelAtLevel). */
  azure2::front::ChannelQuantities channelQuantities(int channelRow) const;
  bool eventFilter(QObject *object, QEvent *event) override;

 public slots:
  void addLevel();
  void addLevel(LevelsData level, bool fromFile);
  void removeLevel();
  void editLevel();
  void addLevelsFromIAEA();
  void updateButtons(const QItemSelection &selection);
  void updateFilter(const QItemSelection &selection);
  void updateChannelsPairAddedEdited();
  void updateChannelsPairRemoved(int pairIndex);
  void updateDetails(const QItemSelection &selection);
  void updateMissingSpinParities();
  void updateReducedWidth(const QString &string);
  void calculateWignerLimit();
  void showInfo(int which = 0, QString title = "");
  void fixAllWidths();
  void fixAllEnergies();
  void exportLatexTable();
  /*! Show the nucleus column and filter only when the pairs form more than one compound nucleus. */
  void updateNucleusColumn();
  /*! Show the levels of one nucleus (0 = all). */
  void filterByNucleus(int index);
  void removeSelectedLevels();
  void duplicateSelectedLevels();
  void addLevelFromBarClicked();
  /*! Put the cursor in the add bar's energy (Insert). */
  void focusAddBar();
  /*! The add bar set to a J^pi (a link of the missing-J^pi text: "jpi:<2J>:<pi>"). */
  void fillAddBar(const QString &link);
  void addBackgroundPole();
  void showLevelsMenu(const QPoint &pos);
  void showChannelsMenu(const QPoint &pos);
  void undoLater();
  void currentLevelChanged();
  void refreshChannelInfo();
  void updateJpiFilter();
  void jpiFilterChosen(int index);
  void updateCount();

 signals:
  void readNewPair(PairsData, int, bool);
  void readExistingPair(PairsData, int, bool);

 private:
  QTabWidget *extras = nullptr;
  /// The state of the tabs as it was after the last recorded action, and the stack (LevelsUndo.h).
  LevelsState currentState() const;
  void applyState(const LevelsState &state);
  void finishLoading();
  QLayout *buildFilterBar();
  QLayout *buildAddBar();
  void setUpShortcuts();
  /// The window's project (the facade's functions that read the project: channel radii, the relativistic option).
  ProjectHolder *holder() const;
  void showRefusal(const QModelIndex &index, const QString &why);
  void setStatus(const QString &text, bool error = false);
  QString levelName(const LevelsData &level) const;
  /// Copy the widths, fix flags and file tokens of `from`'s channels onto the channels of level `to` (same pair, s, l).
  void copyChannelValues(int from, int to);
  QUndoStack *undo_ = nullptr;
  LevelsState lastState_;
  bool loading_ = false, restoring_ = false, undoPending_ = false;
  QString undoText_;
  HalfSpinBox *addJ_ = nullptr;
  ParityToggle *addPi_ = nullptr;
  QLineEdit *addE_ = nullptr;
  QComboBox *addNucleus_ = nullptr;
  QPushButton *addBarButton_ = nullptr;
  QPushButton *duplicateButton_ = nullptr;
  QLineEdit *search_ = nullptr;
  QComboBox *jpiFilter_ = nullptr;
  QCheckBox *onlyIncluded_ = nullptr;
  QLabel *count_ = nullptr;
  QLabel *status_ = nullptr;
  QSpinBox *maxLSpin;
  QSpinBox *maxMultSpin;
  QSpinBox *maxNumMultSpin;
  QCheckBox *includeAllJPiCheck;
  QLabel *missingJPiLabel;
  QPushButton *addLevelButton;
  QPushButton *removeLevelButton;
  QPushButton *iaeaQueryButton;
  QPushButton *fixAllWidthsButton;
  QPushButton *fixAllEnergiesButton;
  QPushButton *exportLatexButton;
  QComboBox *nucleusFilterCombo;
  PairsModel *pairsModel;
  LevelsModel *levelsModel;
  ChannelsModel *channelsModel;
  QTableView *levelsView;
  QTableView *channelsView;
  LevelsFilterProxy *levelsModelProxy;
  QSortFilterProxyModel *proxyModel;
  ChannelDetails *channelDetails;
  // Inputs for the physical Wigner-limit calculation of the selected channel.
  bool wignerApplicable_ = false;
  double wignerZ1_ = 0., wignerZ2_ = 0., wignerRedMass_ = 0.;
  double wignerRadius_ = 0., wignerEcm_ = 0.;
  int wignerL_ = 0;
  QSignalMapper *mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};


#endif
