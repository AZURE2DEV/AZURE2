#ifndef AZUREPLOT_H
#define AZUREPLOT_H

#include <qwt_plot.h>
#include <qwt_text.h>
#include <qwt_symbol.h>
#include <qwt_plot_zoomer.h>
#include <QColor>

class QwtPlotCurve;
class QwtPlotIntervalCurve;
class QwtPlotGrid;
class QwtPlotMarker;
class QwtLegend;
class PlotTab;
class LevelsModel;

struct PlotPoint {
  double energy;
  double excitationEnergy;
  double angle;
  double fitCrossSection;
  double fitSFactor;
  double dataCrossSection;
  double dataErrorCrossSection;
  double dataSFactor;
  double dataErrorSFactor;
  // 1-sigma analytic uncertainty band on the calculation (0 if no band file).
  double fitCrossSectionError;
  double fitSFactorError;
};

class AZUREZoomer : public QwtPlotZoomer {
 public:
  AZUREZoomer(QWidget *canvas) : QwtPlotZoomer(canvas) {};
 protected:
  QwtText trackerTextF( const QPointF &pos ) const;

};

class PlotEntry {
 public:
  PlotEntry(int type, int entranceKey, int exitKey, int index, QString filename);
  ~PlotEntry();

  int type() const {return type_;};
  int entranceKey() const {return entranceKey_;};
  int exitKey() const {return exitKey_;};
  QString filename() const {return filename_;};

  QString label() const {return label_;};
  void setLabel(const QString& label) {label_ = label;};

  QColor color() const {return color_;};
  void setColor(const QColor& c) {color_ = c;};

  QColor fitColor() const {return fitColor_;};
  void setFitColor(const QColor& c) {fitColor_ = c;};

  QwtSymbol::Style symbolStyle() const {return symbolStyle_;};
  void setSymbolStyle(QwtSymbol::Style s) {symbolStyle_ = s;};

  int symbolSize() const {return symbolSize_;};
  void setSymbolSize(int s) {symbolSize_ = s;};

  int lineWidth() const {return lineWidth_;};
  void setLineWidth(int w) {lineWidth_ = w;};

  bool readData();
  bool hasBand() const {return hasBand_;};
  void attach(QwtPlot*, int xAxisType, int yAxisType, bool showBand=false);
  void detach();

  static QString labelFromFilename(const QString& filename);

 public:
  friend class AZUREPlot;

 private:
  void sortPointsByXAxis(int xAxisType);
  bool readBandData(QVector<double>& xsErr, QVector<double>& sErr);
  int type_;
  int entranceKey_;
  int exitKey_;
  int index_;
  QString filename_;
  QString label_;
  QColor color_;
  QColor fitColor_;
  QwtSymbol::Style symbolStyle_;
  int symbolSize_;
  int lineWidth_;
  bool hasBand_;
  QwtPlotCurve* dataCurve_;
  QwtPlotIntervalCurve* dataErrorCurve_;
  QwtPlotCurve* fitCurve_;
  QwtPlotIntervalCurve* bandCurve_;
  QVector<PlotPoint> points_;
};

class AZUREPlot : public QwtPlot {

  Q_OBJECT

 public:
  AZUREPlot(PlotTab* plotTab, QWidget* parent = 0);
  void setXAxisLog(bool set);
  void setYAxisLog(bool set);
  void setXAxisType(unsigned int type);
  void setYAxisType(unsigned int type);
  void setGridVisible(bool visible);
  void setLegendVisible(bool visible);
  void setBandVisible(bool visible);
  bool isBandVisible() const {return bandVisible;};
  void setLevelsModel(LevelsModel* model);
  void setLevelsVisible(bool visible);
  void refreshLevelMarkers();

  const QList<PlotEntry*>& getEntries() const {return entries;};
  void redrawEntries();

 public slots:
  void clearEntries();
  void draw(QList<PlotEntry*> newEntries);
  void update();
  void exportPlot();
  void print();

 private:
  void clearLevelMarkers();

  unsigned int xAxisType;
  unsigned int yAxisType;
  QList<PlotEntry*> entries;
  AZUREZoomer* zoomer;
  PlotTab* containingTab;
  QwtPlotGrid* grid;
  QwtLegend* legend;
  LevelsModel* levelsModel;
  bool levelsVisible;
  bool bandVisible;
  QList<QwtPlotMarker*> levelMarkers;
};

#endif
