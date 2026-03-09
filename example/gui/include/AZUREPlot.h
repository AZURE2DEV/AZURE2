#ifndef AZUREPLOT_H
#define AZUREPLOT_H

#include <qwt_plot.h>
#include <qwt_symbol.h>
#include <qwt_plot_zoomer.h>

class QwtPlotCurve;
class QwtPlotIntervalCurve;
class PlotTab;

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

  bool readData();
  void attach(QwtPlot*,int,int,QwtSymbol::Style);
  void detach();

 public:
  friend class AZUREPlot;

 private:
  bool hasNegative_;
  int type_;
  int entranceKey_;
  int exitKey_;
  int index_;
  QString filename_;
  QwtPlotCurve* dataCurve_;
  QwtPlotIntervalCurve* dataErrorCurve_;
  QwtPlotCurve* fitCurve_;
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
  void clearEntries();

 public slots:
  void draw(QList<PlotEntry*> newEntries);
  void update();
  void exportPlot();
  void print();

 private:
  unsigned int xAxisType;
  unsigned int yAxisType;
  QList<PlotEntry*> entries;
  AZUREZoomer* zoomer;
  PlotTab* containingTab;
};

#endif
