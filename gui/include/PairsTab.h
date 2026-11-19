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
#include <QGroupBox>
#include <QLabel>
#include <QShowEvent>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QTextStream>
#include "AddPairDialog.h"

QT_BEGIN_NAMESPACE

class QPushButton;
class QTabWidget;

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
  /// A panel shown below the pairs (the per-channel radii of <channelRadii>, AZURESetup).
  void addExtraPanel(const QString &title, QWidget *panel);
  /*! A new file: no pairs. */
  void reset();
  /// Select pair `row` and remove it as the Remove button does (the tests).
  void removePairAt(int row);

 public slots:
  void addPair();
  void addPair(PairsData pair, int pairIndex, bool fromFile);
  void editPair();
  void editPair(PairsData pair, int pairIndex, bool fromFile);
  void removePair();
  void updateButtons(const QItemSelection &selection);
  void showInfo(int which = 0, QString title = "");
  void updateNucleusColumn();
  /*! "Masses from AME2020...": compare every pair with the mass table and fill the ticked ones (MassCheckDialog). */
  void checkMasses();

 signals:
  void pairAdded(int);
  void pairRemoved(int);
  void pairEdited(int);

 private:
  PairsModel *pairsModel;
  QTableView *pairsView;
  QTabWidget *extras;

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
