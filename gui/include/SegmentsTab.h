#ifndef SEGMENTSTAB_H
#define SEGMENTSTAB_H

#include <QWidget>
#include <QLineEdit>
#include <QItemSelection>
#include <QTableView>
#include <QPushButton>
#include <QSignalMapper>
#include <QPointer>
#include <QComboBox>
#include <QSortFilterProxyModel>
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"
#include "AddSegDataDialog.h"
#include "AddSegTestDialog.h"

class InfoDialog;
class PairsModel;

class SegmentsTab : public QWidget {
  Q_OBJECT

 public:
  SegmentsTab(QWidget *parent = 0);
  SegmentsTestModel* getSegmentsTestModel();
  SegmentsDataModel* getSegmentsDataModel();
  void reset();
  /*QLineEdit *getSegDataFileText() const {return segDataFileText;};
    QLineEdit *getSegTestFileText() const {return segTestFileText;};*/

 public slots:
  void addSegDataLine();
  void addSegDataLine(SegmentsDataData line);
  void addSegTestLine();
  void addSegTestLine(SegmentsTestData line);
  void editSegDataLine();
  void editSegTestLine();
  void deleteSegDataLine();
  void deleteSegTestLine();
  void moveSegDataLineUp();
  void moveSegDataLineDown();
  void moveSegTestLineUp();
  void moveSegTestLineDown();
  void updateSegDataButtons(const QItemSelection &selection);
  void updateSegTestButtons(const QItemSelection &selection);
  void checkAllSegData();
  void uncheckAllSegData();
  void getDataFromExfor();
  void filterSegDataByPairs();
  /*void openSegDataFile();
  void openSegDataFile(QString filename);
  void saveSegDataFile();
  void saveAsSegDataFile();*/
  /*bool readSegDataFile(QString filename);*/
  bool readSegDataFile(QTextStream& inStream);
  /*bool writeSegDataFile(QString filename);*/
  bool writeSegDataFile(QTextStream& outStream);
  /*void openSegTestFile();
  void openSegTestFile(QString filename);
  void saveSegTestFile();
  void saveAsSegTestFile();*/
  /*bool readSegTestFile(QString filename);*/
  bool readSegTestFile(QTextStream& inStream);
  /*bool writeSegTestFile(QString filename);*/
  bool writeSegTestFile(QTextStream& outStream);
  void setPairsModel(PairsModel* model);
  void updateFilterComboboxes(PairsModel* model);
  void updateFilterComboboxes(); // Update filters using stored model
  void showInfo(int which=0,QString title="");

 private:
  void moveSegDataLine(unsigned int upDown);
  void moveSegTestLine(unsigned int upDown);

  /*QLineEdit *segDataFileText;*/
  PairsModel *pairsModel;
  SegmentsDataModel *segmentsDataModel;
  QTableView *segmentsDataView;
  QPushButton *segDataAddButton;
  //QPushButton *segDataEditButton;
  QPushButton *segDataDeleteButton;
  QPushButton *segDataUpButton;
  QPushButton *segDataDownButton;
  QPushButton *segDataCheckAllButton;
  QPushButton *segDataUncheckAllButton;
  QPushButton *segDataExforButton;
  QComboBox *segDataEntranceFilter;
  QComboBox *segDataExitFilter;
  /*QLineEdit *segTestFileText;*/
  SegmentsTestModel *segmentsTestModel;
  QTableView *segmentsTestView;
  QPushButton *segTestAddButton;
  //QPushButton *segTestEditButton;
  QPushButton *segTestDeleteButton;
  QPushButton *segTestUpButton;
  QPushButton *segTestDownButton;
  QSignalMapper* mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};

#endif
