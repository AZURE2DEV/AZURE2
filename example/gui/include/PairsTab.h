#ifndef PAIRSTAB_H
#define PAIRSTAB_H

#include <QWidget>
#include <QItemSelection>
#include <QSignalMapper>
#include "PairsModel.h"
#include "AddPairDialog.h"

QT_BEGIN_NAMESPACE

class QPushButton;

QT_END_NAMESPACE

class InfoDialog;

class PairsTab : public QWidget {
  Q_OBJECT

 public:
  PairsTab(QWidget *parent = 0);
  PairsModel *getPairsModel();
  bool parseOldECSection(QTextStream&);

 public slots:
  void addPair();
  void addPair(PairsData pair,int pairIndex,bool fromFile);
  void editPair();
  void editPair(PairsData pair,int pairIndex,bool fromFile);
  void removePair();
  void updateButtons(const QItemSelection &selection);
  void showInfo(int which=0,QString title="");

 signals:
  void pairAdded(int);
  void pairRemoved(int);
  void pairEdited(int);

 private:
  PairsModel *pairsModel;
  QTableView *pairsView;
  QPushButton *addButton;
  QPushButton *deleteButton;
  QSignalMapper* mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};

#endif
