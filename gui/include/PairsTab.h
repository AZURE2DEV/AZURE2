#ifndef PAIRSTAB_H
#define PAIRSTAB_H

#include <QWidget>
#include <QItemSelection>
#include <QSignalMapper>
#include <QTableView>
#include <QPushButton>
#include <QSignalMapper>
#include <QPointer>

#include "PairsModel.h"
#include "AbundancesModel.h"
#include "ChannelRadiiModel.h"
#include <QGroupBox>
#include <QLabel>
#include <QShowEvent>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QTextStream>
#include "AddPairDialog.h"

QT_BEGIN_NAMESPACE

class QPushButton;

QT_END_NAMESPACE

class InfoDialog;

/*!
 * The Particle Pairs tab: the pairs of the compound nucleus, and the source every other tab takes its pair list from.
 */
class PairsTab : public QWidget {
  Q_OBJECT

 public:
  PairsTab(QWidget *parent = 0);
  PairsModel *getPairsModel();
  AbundancesModel *getAbundancesModel();
  /// The per-channel radii (<channelRadii> block, 2026-11-10).
  ChannelRadiiModel *getChannelRadiiModel();
  /// Read the <channelRadii> block body (after its opening tag) up to </channelRadii>.
  bool readChannelRadii(QTextStream &);
  /// Write the block body; the caller writes the tags, and only when hasChannelRadii().
  bool writeChannelRadii(QTextStream &);
  bool hasChannelRadii() const;
  /// The list of every particle channel of the calculation (ChannelRadiiListModel), and where it comes from: the
  /// Levels tab's channels of its active levels (AZURESetup).  Refreshed when the tab is shown and after a file read.
  ChannelRadiiListModel *getChannelRadiiListModel() { return channelRadiiListModel; }
  void setChannelEnumerator(std::function<QList<ChannelKey>()> enumerate) { channelEnumerator = std::move(enumerate); }
  /// Configure > Runtime Options "Per-channel radii": the box is shown when this is on, and always for a project
  /// whose file carries <channelRadii> rows (so model input is never hidden).
  void setPerChannelRadiiOption(bool on);
  bool perChannelRadiiOption() const { return perChannelRadiiOn; }
  bool parseOldECSection(QTextStream &);
  /*! Read the <abundances> block ("[<sample>:] Z A abundance vary sigma" and "[<sample>:] sum sigma" lines) up to its closing tag. */
  bool readAbundances(QTextStream &);
  /*! Write the block's lines (no tags): every sample's nuclei, then its own sum-to-one width if any. */
  bool writeAbundances(QTextStream &);
  /*! Should the file carry an <abundances> block: more than one nucleus, or a non-default entry. */
  bool hasAbundanceBlock() const;
  double abundanceSumSigma() const;
  void setAbundanceSumSigma(double sigma);
  /*! A new file: default abundances, no sum-to-one prior. */
  void reset();

 public slots:
  void addPair();
  void addPair(PairsData pair, int pairIndex, bool fromFile);
  void editPair();
  void editPair(PairsData pair, int pairIndex, bool fromFile);
  void removePair();
  void updateButtons(const QItemSelection &selection);
  void showInfo(int which = 0, QString title = "");
  void updateNucleusColumn();
  void updateAbundancesBox();
  /// Rebuild the channel list from the enumerator and update the box (visibility, the note on broader rows).
  void refreshChannelRadii();
  /// Remove the selected channel's own radius (its pair's, or a broader row's, applies again).
  void resetChannelRadius();
  /*! The sample selector follows the model's samples; the sum width shown is the chosen sample's. */
  void updateSampleCombo();
  void sampleChosen(int index);
  void addSample();
  void removeSample();
  void sampleSumSigmaEdited(double sigma);
  /*! "Masses from AME2020...": compare every pair with the mass table and fill the ticked ones (MassCheckDialog). */
  void checkMasses();

 signals:
  void pairAdded(int);
  void pairRemoved(int);
  void pairEdited(int);

 private:
  PairsModel *pairsModel;
  QTableView *pairsView;
  AbundancesModel *abundancesModel;
  QTableView *abundancesView;
  QGroupBox *abundancesBox;
  QDoubleSpinBox *sumSigmaSpin;
  QComboBox *sampleCombo;
  QPushButton *addSampleButton;
  QPushButton *removeSampleButton;
  QDoubleSpinBox *sampleSumSigmaSpin;
  ChannelRadiiModel *channelRadiiModel;
  ChannelRadiiListModel *channelRadiiListModel;
  QTableView *channelRadiiView;
  QGroupBox *channelRadiiBox;
  QLabel *channelRadiiNote;
  QPushButton *resetRadiusButton;
  std::function<QList<ChannelKey>()> channelEnumerator;
  bool perChannelRadiiOn = false;
  void updateChannelRadiiBox();

 protected:
  void showEvent(QShowEvent *event) override;

 private:
  QPushButton *addButton;
  QPushButton *deleteButton;
  QPushButton *massButton;
  QSignalMapper *mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};

#endif
