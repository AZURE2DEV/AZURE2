#ifndef LEVELSTAB_H
#define LEVELSTAB_H

#include <QWidget>
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

class InfoDialog;

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
  bool writeNuclearFile(QTextStream &outStream);
  bool readNuclearFile(QTextStream &inStream);
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
  LevelsModel *getLevelsModel() { return levelsModel; }
  /*! Hard-sphere scattering in every J^pi up to the maximum orbital momentum, levels or not
   *  (includeAllJPi in the .azr <config> block). On for a new project, off for a file that lacks the key. */
  bool includeAllJPi() const { return includeAllJPiCheck->isChecked(); }
  void setIncludeAllJPi(bool on) { includeAllJPiCheck->setChecked(on); }
  int maxL() const { return maxLSpin->value(); }
  void setMaxL(int l) { maxLSpin->setValue(l); }
  /*! The J^pi the pairs allow up to the maximum orbital momentum that no active level carries. */
  QStringList missingSpinParities();
  ChannelsModel *getChannelsModel() { return channelsModel; }

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

 signals:
  void readNewPair(PairsData, int, bool);
  void readExistingPair(PairsData, int, bool);

 private:

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
  QSortFilterProxyModel *levelsModelProxy;
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
