#ifndef PLOTTAB_H
#define PLOTTAB_H

#include <QWidget>
#include <QSortFilterProxyModel>
#include <QSignalMapper>
#include <QPointer>
#include <QDialog>

#include <qwt_symbol.h>

class Config;
class AZUREPlot;
class PlotEntry;
class SegmentsDataModel;
class SegmentsTestModel;
class LevelsModel;
class InfoDialog;

QT_BEGIN_NAMESPACE

class QRadioButton;
class QListView;
class QListWidget;
class QCheckBox;
class QPushButton;
class QComboBox;
class QLineEdit;
class QSpinBox;
class QLabel;

QT_END_NAMESPACE

class ChannelFilterProxy : public QSortFilterProxyModel {
  Q_OBJECT
 public:
  ChannelFilterProxy(QWidget* parent = 0) : QSortFilterProxyModel(parent),
    entranceFilter_(-1), exitFilter_(-1) {};
  void setEntranceFilter(int pairIndex);
  void setExitFilter(int pairIndex);

 protected:
  bool entranceExitAccepts(int source_row, const QModelIndex& source_parent) const;
  int entranceFilter_;
  int exitFilter_;
};

class SegTestProxyModel : public ChannelFilterProxy {
 public:
  SegTestProxyModel(QWidget* parent = 0) : ChannelFilterProxy(parent) {};
  QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const;
  bool filterAcceptsRow(int source_row, const QModelIndex &source_parent) const;
};

class SegDataProxyModel : public ChannelFilterProxy {
 public:
  SegDataProxyModel(QWidget* parent = 0) : ChannelFilterProxy(parent) {};
  QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const;
  bool filterAcceptsRow(int source_row, const QModelIndex &source_parent) const;
};

class CurveStyleDialog : public QDialog {
  Q_OBJECT
 public:
  CurveStyleDialog(PlotEntry* entry, QWidget* parent = 0);

 private slots:
  void chooseColor();
  void chooseFitColor();
  void accept() override;

 private:
  void updateColorButton();
  void updateFitColorButton();

  PlotEntry* entry_;
  QLineEdit* labelEdit_;
  QPushButton* colorButton_;
  QPushButton* fitColorButton_;
  QColor color_;
  QColor fitColor_;
  QComboBox* symbolCombo_;
  QSpinBox* symbolSizeSpin_;
  QSpinBox* lineWidthSpin_;
};

class PlotTab : public QWidget {

  Q_OBJECT

 public:
  PlotTab(Config& config, SegmentsDataModel* dataModel, SegmentsTestModel* testModel, LevelsModel* levelsModel = 0, QWidget* parent = 0);
  QList<PlotEntry*> getDataSegments();
  QList<PlotEntry*> getTestSegments();
  void reset();

 public slots:
  void draw();
  void xAxisTypeChanged();
  void yAxisTypeChanged();
  void xAxisLogScaleChanged(bool);
  void yAxisLogScaleChanged(bool);
  void gridToggled(bool);
  void legendToggled(bool);
  void levelsToggled(bool);
  void bandToggled(bool);
  void dataChannelFilterChanged();
  void testChannelFilterChanged();
  void rebuildChannelFilterOptions();
  void rebuildCurveList();
  void editSelectedCurveStyle();
  void showInfo(int which=0,QString title="");

 public:
  friend class AZUREPlot;

 private:
  void populateChannelFilter(QComboBox* combo, QAbstractItemModel* model, int column);

  Config& configure;
  AZUREPlot* azurePlot;
  QListView* dataSegmentSelectorList;
  QListView* testSegmentSelectorList;
  QComboBox* dataInChannelCombo;
  QComboBox* dataOutChannelCombo;
  QComboBox* testInChannelCombo;
  QComboBox* testOutChannelCombo;
  QListWidget* curveList;
  QRadioButton* yAxisXSButton;
  QRadioButton* yAxisSFButton;
  QComboBox* xAxisTypeCombo;
  QCheckBox* xAxisIsLogCheck;
  QCheckBox* yAxisIsLogCheck;
  QCheckBox* gridCheck;
  QCheckBox* legendCheck;
  QCheckBox* levelsCheck;
  QCheckBox* bandCheck;
  SegTestProxyModel* segTestProxyModel;
  SegDataProxyModel* segDataProxyModel;
  QPushButton* refreshButton;
  QPushButton* clearButton;
  QPushButton* exportButton;
  QPushButton* printButton;
  QSignalMapper* mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};

#endif
