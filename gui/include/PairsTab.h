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
#include <QGroupBox>
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
  bool parseOldECSection(QTextStream &);
  /*! Read the <abundances> block ("[<sample>:] Z A abundance vary sigma" and "[<sample>:] sum sigma" lines) up to its closing tag. */
  bool readAbundances(QTextStream &);
  /*! Write the block's lines (no tags): every sample's nuclei, then its own sum-to-one width if any. */
  bool writeAbundances(QTextStream &);
  /*! Should the file carry an <abundances> block: more than one nucleus, or a non-default entry. */
  bool hasAbundanceBlock() const;
  double abundanceSumSigma() const;
  void setAbundanceSumSigma(double sigma);
  /*! Relativistic two-body kinematics ("relativistic=1" in <config>). */
  bool relativistic() const;
  void setRelativistic(bool on);
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
  /*! The sample selector follows the model's samples; the sum width shown is the chosen sample's. */
  void updateSampleCombo();
  void sampleChosen(int index);
  void addSample();
  void removeSample();
  void sampleSumSigmaEdited(double sigma);

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
  QCheckBox *relativisticCheck;
  QComboBox *sampleCombo;
  QPushButton *addSampleButton;
  QPushButton *removeSampleButton;
  QDoubleSpinBox *sampleSumSigmaSpin;
  QPushButton *addButton;
  QPushButton *deleteButton;
  QSignalMapper *mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};

#endif
