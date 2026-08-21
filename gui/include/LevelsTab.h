#ifndef LEVELSTAB_H
#define LEVELSTAB_H

#include <QWidget>
#include <QItemSelection>
#include <QLineEdit>
#include <QSpinBox>
#include <QPushButton>
#include <QTableView>
#include <QSortFilterProxyModel>
#include <QSignalMapper>
#include <QPointer>
#include <QFileDialog>
#include <QClipboard>
#include <QApplication>
#include "PairsModel.h"
#include "LevelsModel.h"
#include "ChannelsModel.h"
#include "ChannelDetails.h"

class InfoDialog;

class LevelsTab : public QWidget {
  Q_OBJECT

 public:
  LevelsTab(QWidget *parent = 0);
  void setPairsModel(PairsModel *);
  void updateChannelsLevelAdded(int levelIndex);
  void updateChannelsLevelDeleted(int levelIndex);
  void updateChannelsLevelEdited(int levelIndex);
  QList<ChannelsData> calculateChannels(int levelIndex);
  bool writeNuclearFile(QTextStream &outStream);
  bool readNuclearFile(QTextStream &inStream);
  void reset();
  LevelsModel *getLevelsModel() { return levelsModel; }
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
  void updateReducedWidth(const QString &string);
  void calculateWignerLimit();
  void showInfo(int which = 0, QString title = "");
  void fixAllWidths();
  void fixAllEnergies();
  void exportLatexTable();

 signals:
  void readNewPair(PairsData, int, bool);
  void readExistingPair(PairsData, int, bool);

 private:
  int calculateCompoundNucleus(int &massNumber, int &atomicNumber);

  QSpinBox *maxLSpin;
  QSpinBox *maxMultSpin;
  QSpinBox *maxNumMultSpin;
  QPushButton *addLevelButton;
  QPushButton *removeLevelButton;
  QPushButton *iaeaQueryButton;
  QPushButton *fixAllWidthsButton;
  QPushButton *fixAllEnergiesButton;
  QPushButton *exportLatexButton;
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
