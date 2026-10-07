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
#include <QVector>
#include <QSortFilterProxyModel>
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"
#include "AddSegDataDialog.h"
#include "AddSegTestDialog.h"

class InfoDialog;
class PairsModel;

/*!
 * The Segments tab: the data segments the fit runs against and the test segments an extrapolation is evaluated on, held in two separate models.
 */
class SegmentsTab : public QWidget {
  Q_OBJECT

 public:
  SegmentsTab(QWidget *parent = 0);
  SegmentsTestModel *getSegmentsTestModel();
  SegmentsDataModel *getSegmentsDataModel();
  void reset();
  /// THM on (AZURESetup::setThmEnabled): the segment dialogs show their THM tick.
  void setThmEnabled(bool on) { thmEnabled_ = on; }

  /*! Moves data (test) segment line `from` to `to` (0-based rows) or removes
      line `row`, as the up/down and delete buttons do, and emits
      dataSegmentsRenumbered (testSegmentsRenumbered) so that what refers to
      the lines by number can follow them.  False if the move is refused. */
  bool moveDataSegment(int from, int to);
  void deleteDataSegment(int row);
  bool moveTestSegment(int from, int to);
  void deleteTestSegment(int row);

 signals:
  /*! The lines were renumbered: line k (1-based) before is line newNumber[k - 1]
      now, 0 if it was deleted.  A line added at the end (newNumber covers the
      lines before it) is signalled too. */
  void dataSegmentsRenumbered(const QVector<int> &newNumber);
  void testSegmentsRenumbered(const QVector<int> &newNumber);
  /*QLineEdit *getSegDataFileText() const {return segDataFileText;};
    QLineEdit *getSegTestFileText() const {return segTestFileText;};*/

 public slots:
  void addSegDataLine();
  // fromFile: a line read from the project is kept as it is, as the engine
  // reads it; only one added here is refused as a duplicate.
  void addSegDataLine(SegmentsDataData line, bool fromFile = false);
  void addSegTestLine();
  void addSegTestLine(SegmentsTestData line, bool fromFile = false);
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
  bool readSegDataFile(QTextStream &inStream);
  /*bool writeSegDataFile(QString filename);*/
  bool writeSegDataFile(QTextStream &outStream);
  /*void openSegTestFile();
  void openSegTestFile(QString filename);
  void saveSegTestFile();
  void saveAsSegTestFile();*/
  /*bool readSegTestFile(QString filename);*/
  bool readSegTestFile(QTextStream &inStream);
  /*bool writeSegTestFile(QString filename);*/
  bool writeSegTestFile(QTextStream &outStream);
  void setPairsModel(PairsModel *model);
  void updateFilterComboboxes(PairsModel *model);
  void updateFilterComboboxes();  // Update filters using stored model
  void showInfo(int which = 0, QString title = "");

 private:
  void moveSegDataLine(unsigned int upDown);
  void moveSegTestLine(unsigned int upDown);

  /*QLineEdit *segDataFileText;*/
  PairsModel *pairsModel;
  SegmentsDataModel *segmentsDataModel;
  QTableView *segmentsDataView;
  QPushButton *segDataAddButton;
  // QPushButton *segDataEditButton;
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
  // QPushButton *segTestEditButton;
  QPushButton *segTestDeleteButton;
  QPushButton *segTestUpButton;
  QPushButton *segTestDownButton;
  QSignalMapper *mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
  bool thmEnabled_ = false;
};

#endif
