#include <QBrush>
#include <QCheckBox>
#include <QDialog>
#include <QFileDialog>
#include <QListView>
#include <QPushButton>
#include <QMessageBox>
#include <QImageWriter>
#include <QtPrintSupport/QPrinter>
#include <QtPrintSupport/QPrintDialog>
#include <QPageLayout>
#include <QTextStream>
#include <QFileInfo>

#include "AZUREPlot.h"
#include "PlotTab.h"
#include "LevelsModel.h"
#include <qwt_text.h>
#include "qwt_plot_curve.h"
#include "qwt_plot_grid.h"
#include "qwt_plot_marker.h"
#include "qwt_legend.h"
#include "qwt_plot_intervalcurve.h"
#include "qwt_interval_symbol.h"
#include "qwt_scale_engine.h"
#include "qwt_plot_panner.h"
#include "qwt_plot_renderer.h"
#include <iostream>

#include <math.h>
#include <algorithm>

QwtText AZUREZoomer::trackerTextF(const QPointF &pos) const {
  QString text;
  switch (rubberBand()) {
    case HLineRubberBand:
      text = QString::asprintf("%.4g", pos.y());
      break;
    case VLineRubberBand:
      text = QString::asprintf("%.4g", pos.x());
      break;
    default:
      text = QString::asprintf("%.4g, %.4g", pos.x(), pos.y());
  }
  return QwtText(text);
}

static const QColor kDefaultPalette[] = {
    QColor(0x1f, 0x77, 0xb4),  // blue
    QColor(0xd6, 0x27, 0x28),  // red
    QColor(0x2c, 0xa0, 0x2c),  // green
    QColor(0xff, 0x7f, 0x0e),  // orange
    QColor(0x94, 0x67, 0xbd),  // purple
    QColor(0x8c, 0x56, 0x4b),  // brown
    QColor(0xe3, 0x77, 0xc2),  // pink
    QColor(0x7f, 0x7f, 0x7f),  // gray
    QColor(0xbc, 0xbd, 0x22),  // olive
    QColor(0x17, 0xbe, 0xcf)   // cyan
};
static const int kDefaultPaletteSize = sizeof(kDefaultPalette) / sizeof(QColor);

static const QwtSymbol::Style kSymbolCycle[] = {
    QwtSymbol::Ellipse,
    QwtSymbol::Rect,
    QwtSymbol::Diamond,
    QwtSymbol::Triangle,
    QwtSymbol::DTriangle,
    QwtSymbol::UTriangle,
    QwtSymbol::Cross,
    QwtSymbol::XCross,
    QwtSymbol::Hexagon,
    QwtSymbol::Star1};
static const int kSymbolCycleSize = sizeof(kSymbolCycle) / sizeof(QwtSymbol::Style);

QString PlotEntry::labelFromFilename(const QString &filename) {
  int slash = filename.lastIndexOf('/');
  if (slash < 0) return filename;
  return filename.mid(slash + 1);
}

PlotEntry::PlotEntry(int type, int entranceKey, int exitKey, int index, QString filename) :
  type_(type),
  entranceKey_(entranceKey),
  exitKey_(exitKey),
  index_(index),
  filename_(filename),
  color_(Qt::black),
  // For type==0 (data + fit) default the fit/calculation line to red so it
  // stands out from the data points. For type==1 (fit-only) the fitColor_
  // gets initialized to match color_ later in AZUREPlot::draw().
  fitColor_(type == 0 ? QColor(0xd6, 0x27, 0x28) : QColor()),
  symbolStyle_(QwtSymbol::Ellipse),
  symbolSize_(6),
  lineWidth_(2),
  hasBand_(false),
  dataCurve_(NULL),
  dataErrorCurve_(NULL),
  fitCurve_(NULL),
  bandCurve_(NULL) {
  label_ = labelFromFilename(filename);
}

PlotEntry::~PlotEntry() {
  if (dataCurve_) delete dataCurve_;
  if (dataErrorCurve_) delete dataErrorCurve_;
  if (fitCurve_) delete fitCurve_;
  if (bandCurve_) delete bandCurve_;
}

static bool isFiniteAndPositive(double v) {
  if (v != v) return false;             // NaN
  if (!std::isfinite(v)) return false;  // +/-inf
  if (v <= 0.) return false;
  return true;
}

bool PlotEntry::readData() {
  QFile file(filename_);
  if (!file.open(QIODevice::ReadOnly)) {
    return false;
  }
  QTextStream inStream(&file);
  QString line("");
  bool previousLineBreak = false;
  bool foundBlock = false;
  int blockNumber = 0;

  // Read the sibling ".band" file (analytic uncertainty band) for this block, if
  // present.  bandXsErr[n]/bandSErr[n] correspond to raw line n of this block,
  // so they are assigned before the (order-preserving) validity filtering below.
  QVector<double> bandXsErr, bandSErr;
  hasBand_ = readBandData(bandXsErr, bandSErr);
  int rawLine = 0;

  while (!inStream.atEnd() && !foundBlock) {
    line = inStream.readLine();
    if (line.trimmed().isEmpty()) {
      if (!previousLineBreak) {
        previousLineBreak = true;
        continue;
      }
      if (blockNumber == index_) {
        foundBlock = true;
        break;
      } else {
        blockNumber++;
        previousLineBreak = false;
        points_.clear();
        rawLine = 0;
        continue;
      }
    }
    if (!inStream.atEnd() && !foundBlock) {
      QTextStream in(&line);
      PlotPoint newPoint = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      if (type_ == 0) {
        in >> newPoint.energy >> newPoint.excitationEnergy >> newPoint.angle >> newPoint.fitCrossSection >> newPoint.fitSFactor >> newPoint.dataCrossSection >> newPoint.dataErrorCrossSection >> newPoint.dataSFactor >> newPoint.dataErrorSFactor;
      } else {
        in >> newPoint.energy >> newPoint.excitationEnergy >> newPoint.angle >> newPoint.fitCrossSection >> newPoint.fitSFactor;
      }
      // Attach this raw line's band uncertainty (in file order, before filtering).
      if (blockNumber == index_ && rawLine < bandXsErr.size()) {
        newPoint.fitCrossSectionError = bandXsErr[rawLine];
        newPoint.fitSFactorError = bandSErr[rawLine];
      }
      rawLine++;
      // Discard points that would break log scale: NaN, inf, zero, negative,
      // or (for data points) error bars that drop the lower bound to <= 0.
      // A quantity that is legitimately negative -- an analyzing power -- is
      // exempt from everything but the finiteness test, since for it the
      // positivity guard would delete the negative half of the distribution.
      bool valid;
      if (allowNonPositive_) {
        valid = std::isfinite(newPoint.fitCrossSection) &&
            std::isfinite(newPoint.fitSFactor);
        if (valid && type_ == 0)
          valid = std::isfinite(newPoint.dataCrossSection) &&
              std::isfinite(newPoint.dataSFactor) &&
              std::isfinite(newPoint.dataErrorCrossSection) &&
              std::isfinite(newPoint.dataErrorSFactor);
      } else {
        valid = isFiniteAndPositive(newPoint.fitCrossSection) &&
            isFiniteAndPositive(newPoint.fitSFactor);
        if (valid && type_ == 0) {
          valid = isFiniteAndPositive(newPoint.dataCrossSection) &&
              isFiniteAndPositive(newPoint.dataSFactor) &&
              std::isfinite(newPoint.dataErrorCrossSection) &&
              std::isfinite(newPoint.dataErrorSFactor) &&
              (newPoint.dataCrossSection - newPoint.dataErrorCrossSection) > 0. &&
              (newPoint.dataSFactor - newPoint.dataErrorSFactor) > 0.;
        }
      }
      if (valid) points_.push_back(newPoint);
    }
  }
  inStream.flush();
  file.close();
  if (!foundBlock) {
    return false;
  }
  return !points_.isEmpty();
}

// Read the sibling "<primary>.band" file for this entry's block index.  The band
// file has the same block/line structure as the primary output file, with
// columns: energy, excitation, angle, xs, dXS, S, dS.  Fills xsErr/sErr (the dXS
// and dS columns) for every raw line of the block, in file order.  Returns true
// if a band file existed and yielded a non-empty block.
bool PlotEntry::readBandData(QVector<double> &xsErr, QVector<double> &sErr) {
  xsErr.clear();
  sErr.clear();
  QFile file(filename_ + ".band");
  if (!file.open(QIODevice::ReadOnly)) return false;
  QTextStream inStream(&file);
  QString line("");
  bool previousLineBreak = false;
  bool foundBlock = false;
  int blockNumber = 0;
  while (!inStream.atEnd() && !foundBlock) {
    line = inStream.readLine();
    if (line.trimmed().isEmpty()) {
      if (!previousLineBreak) {
        previousLineBreak = true;
        continue;
      }
      if (blockNumber == index_) {
        foundBlock = true;
        break;
      }
      blockNumber++;
      previousLineBreak = false;
      xsErr.clear();
      sErr.clear();
      continue;
    }
    if (blockNumber == index_) {
      QTextStream in(&line);
      double e, ex, ang, xs, dxs, s, ds;
      in >> e >> ex >> ang >> xs >> dxs >> s >> ds;
      xsErr.push_back(dxs);
      sErr.push_back(ds);
    }
  }
  file.close();
  return (foundBlock || blockNumber == index_) && !xsErr.isEmpty();
}

void PlotEntry::sortPointsByXAxis(int xAxisType) {
  std::sort(points_.begin(), points_.end(), [xAxisType](const PlotPoint &a, const PlotPoint &b) {
    switch (xAxisType) {
      case 0: return a.energy < b.energy;
      case 1: return a.excitationEnergy < b.excitationEnergy;
      case 2: return a.angle < b.angle;
      default: return a.energy < b.energy;
    }
  });
}

void PlotEntry::attach(QwtPlot *plot, int xAxisType, int yAxisType, bool showBand) {
  sortPointsByXAxis(xAxisType);
  QVector<QPointF> fit(points_.size());
  if (type_ == 0) {
    QVector<QPointF> data(points_.size());
    QVector<QwtIntervalSample> error(points_.size());

    if (xAxisType == 0 && yAxisType == 0) {
      for (int i = 0; i < points_.size(); i++) {
        data[i] = QPointF(points_[i].energy, points_[i].dataCrossSection);
        fit[i] = QPointF(points_[i].energy, points_[i].fitCrossSection);
        error[i] = QwtIntervalSample(points_[i].energy,
                                     QwtInterval(points_[i].dataCrossSection - points_[i].dataErrorCrossSection,
                                                 points_[i].dataCrossSection + points_[i].dataErrorCrossSection));
      }
    } else if (xAxisType == 0 && yAxisType == 1) {
      for (int i = 0; i < points_.size(); i++) {
        data[i] = QPointF(points_[i].energy, points_[i].dataSFactor);
        fit[i] = QPointF(points_[i].energy, points_[i].fitSFactor);
        error[i] = QwtIntervalSample(points_[i].energy,
                                     QwtInterval(points_[i].dataSFactor - points_[i].dataErrorSFactor,
                                                 points_[i].dataSFactor + points_[i].dataErrorSFactor));
      }
    } else if (xAxisType == 1 && yAxisType == 0) {
      for (int i = 0; i < points_.size(); i++) {
        data[i] = QPointF(points_[i].excitationEnergy, points_[i].dataCrossSection);
        fit[i] = QPointF(points_[i].excitationEnergy, points_[i].fitCrossSection);
        error[i] = QwtIntervalSample(points_[i].excitationEnergy,
                                     QwtInterval(points_[i].dataCrossSection - points_[i].dataErrorCrossSection,
                                                 points_[i].dataCrossSection + points_[i].dataErrorCrossSection));
      }
    } else if (xAxisType == 1 && yAxisType == 1) {
      for (int i = 0; i < points_.size(); i++) {
        data[i] = QPointF(points_[i].excitationEnergy, points_[i].dataSFactor);
        fit[i] = QPointF(points_[i].excitationEnergy, points_[i].fitSFactor);
        error[i] = QwtIntervalSample(points_[i].excitationEnergy,
                                     QwtInterval(points_[i].dataSFactor - points_[i].dataErrorSFactor,
                                                 points_[i].dataSFactor + points_[i].dataErrorSFactor));
      }
    } else if (xAxisType == 2 && yAxisType == 0) {
      for (int i = 0; i < points_.size(); i++) {
        data[i] = QPointF(points_[i].angle, points_[i].dataCrossSection);
        fit[i] = QPointF(points_[i].angle, points_[i].fitCrossSection);
        error[i] = QwtIntervalSample(points_[i].angle,
                                     QwtInterval(points_[i].dataCrossSection - points_[i].dataErrorCrossSection,
                                                 points_[i].dataCrossSection + points_[i].dataErrorCrossSection));
      }
    } else if (xAxisType == 2 && yAxisType == 1) {
      for (int i = 0; i < points_.size(); i++) {
        data[i] = QPointF(points_[i].angle, points_[i].dataSFactor);
        fit[i] = QPointF(points_[i].angle, points_[i].fitSFactor);
        error[i] = QwtIntervalSample(points_[i].angle,
                                     QwtInterval(points_[i].dataSFactor - points_[i].dataErrorSFactor,
                                                 points_[i].dataSFactor + points_[i].dataErrorSFactor));
      }
    }

    dataCurve_ = new QwtPlotCurve;
    dataCurve_->setTitle(label_);
    dataCurve_->setRenderHint(QwtPlotItem::RenderAntialiased);
    dataCurve_->setStyle(QwtPlotCurve::NoCurve);
    QwtSymbol *symbol = new QwtSymbol(symbolStyle_);
    symbol->setSize(symbolSize_);
    symbol->setPen(QPen(color_));
    symbol->setColor(color_);

    dataErrorCurve_ = new QwtPlotIntervalCurve;
    dataErrorCurve_->setRenderHint(QwtPlotItem::RenderAntialiased);
    dataErrorCurve_->setStyle(QwtPlotIntervalCurve::NoCurve);
    dataErrorCurve_->setItemAttribute(QwtPlotItem::Legend, false);
    QwtIntervalSymbol *errorBar =
        new QwtIntervalSymbol(QwtIntervalSymbol::Bar);
    errorBar->setWidth(8);
    errorBar->setPen(QPen(color_));

    dataCurve_->setSymbol(symbol);
    dataCurve_->setSamples(data);

    dataErrorCurve_->setSymbol(errorBar);
    dataErrorCurve_->setSamples(error);
  } else {
    if (xAxisType == 0 && yAxisType == 0) {
      for (int i = 0; i < points_.size(); i++)
        fit[i] = QPointF(points_[i].energy, points_[i].fitCrossSection);
    } else if (xAxisType == 0 && yAxisType == 1) {
      for (int i = 0; i < points_.size(); i++)
        fit[i] = QPointF(points_[i].energy, points_[i].fitSFactor);
    } else if (xAxisType == 1 && yAxisType == 0) {
      for (int i = 0; i < points_.size(); i++)
        fit[i] = QPointF(points_[i].excitationEnergy, points_[i].fitCrossSection);
    } else if (xAxisType == 1 && yAxisType == 1) {
      for (int i = 0; i < points_.size(); i++)
        fit[i] = QPointF(points_[i].excitationEnergy, points_[i].fitSFactor);
    } else if (xAxisType == 2 && yAxisType == 0) {
      for (int i = 0; i < points_.size(); i++)
        fit[i] = QPointF(points_[i].angle, points_[i].fitCrossSection);
    } else if (xAxisType == 2 && yAxisType == 1) {
      for (int i = 0; i < points_.size(); i++)
        fit[i] = QPointF(points_[i].angle, points_[i].fitSFactor);
    }
  }

  // Analytic 1-sigma uncertainty band: a translucent tube around the
  // calculation curve, drawn behind everything else.  The band half-width is the
  // cross-section or S-factor uncertainty depending on the current y axis.
  if (showBand && hasBand_) {
    QVector<QwtIntervalSample> band(points_.size());
    bool any = false;
    for (int i = 0; i < points_.size(); i++) {
      double err = (yAxisType == 0) ? points_[i].fitCrossSectionError
                                    : points_[i].fitSFactorError;
      double yc = fit[i].y();
      double lo = yc - err, hi = yc + err;
      if (err > 0.) any = true;
      if (lo <= 0.) lo = yc * 1.e-6;  // keep the tube positive on log-scale axes
      band[i] = QwtIntervalSample(fit[i].x(), QwtInterval(lo, hi));
    }
    if (any) {
      bandCurve_ = new QwtPlotIntervalCurve;
      bandCurve_->setRenderHint(QwtPlotItem::RenderAntialiased);
      bandCurve_->setStyle(QwtPlotIntervalCurve::Tube);
      bandCurve_->setItemAttribute(QwtPlotItem::Legend, false);
      QColor bandColor = (type_ == 0 && fitColor_.isValid()) ? fitColor_ : color_;
      bandColor.setAlpha(60);
      bandCurve_->setBrush(QBrush(bandColor));
      bandCurve_->setPen(QPen(Qt::NoPen));
      bandCurve_->setSamples(band);
      bandCurve_->setZ(5);  // behind data (10/20) and calculation (30)
      bandCurve_->attach(plot);
    }
  }

  fitCurve_ = new QwtPlotCurve;
  fitCurve_->setRenderHint(QwtPlotItem::RenderAntialiased);
  fitCurve_->setStyle(QwtPlotCurve::Lines);
  // For type==0 (data + calculation), use the dedicated fit color so the
  // calculation is visually separated from the data points. For type==1
  // (calculation-only) the entry has a single color_, which is the fit color.
  QColor lineColor = (type_ == 0 && fitColor_.isValid()) ? fitColor_ : color_;
  fitCurve_->setPen(QPen(lineColor, lineWidth_));

  if (type_ == 0) {
    // Avoid duplicate legend entries: keep only the data curve entry.
    fitCurve_->setItemAttribute(QwtPlotItem::Legend, false);
  } else {
    fitCurve_->setTitle(label_);
  }

  fitCurve_->setSymbol(new QwtSymbol(QwtSymbol::NoSymbol));
  fitCurve_->setSamples(fit);

  // Draw the calculation line on top of the data points and error bars.
  if (dataErrorCurve_) dataErrorCurve_->setZ(10);
  if (dataCurve_) dataCurve_->setZ(20);
  fitCurve_->setZ(30);

  if (dataCurve_) dataCurve_->attach(plot);
  if (dataErrorCurve_) dataErrorCurve_->attach(plot);
  if (fitCurve_) fitCurve_->attach(plot);
}

void PlotEntry::detach() {
  if (dataCurve_) {
    dataCurve_->detach();
    delete dataCurve_;
    dataCurve_ = NULL;
  }
  if (dataErrorCurve_) {
    dataErrorCurve_->detach();
    delete dataErrorCurve_;
    dataErrorCurve_ = NULL;
  }
  if (fitCurve_) {
    fitCurve_->detach();
    delete fitCurve_;
    fitCurve_ = NULL;
  }
  if (bandCurve_) {
    bandCurve_->detach();
    delete bandCurve_;
    bandCurve_ = NULL;
  }
}

AZUREPlot::AZUREPlot(PlotTab *plotTab, QWidget *parent) :
  QwtPlot(parent),
  xAxisType(0),
  yAxisType(0),
  containingTab(plotTab),
  levelsModel(NULL),
  levelsVisible(false),
  bandVisible(false) {
  setCanvasBackground(QColor(Qt::white));
  setAutoReplot(true);

  grid = new QwtPlotGrid;
  grid->setMajorPen(QPen(QColor(220, 220, 220), 0, Qt::SolidLine));
  grid->setMinorPen(QPen(QColor(240, 240, 240), 0, Qt::DotLine));
  grid->enableXMin(true);
  grid->enableYMin(true);
  grid->attach(this);
  grid->setVisible(false);

  legend = new QwtLegend;
  legend->setDefaultItemMode(QwtLegendData::ReadOnly);
  insertLegend(legend, QwtPlot::RightLegend);

  zoomer = new AZUREZoomer(canvas());
  zoomer->setRubberBandPen(QColor(Qt::black));
  zoomer->setTrackerPen(QColor(Qt::black));
  zoomer->setMousePattern(QwtEventPattern::MouseSelect2,
                          Qt::RightButton, Qt::ControlModifier);
  zoomer->setMousePattern(QwtEventPattern::MouseSelect3,
                          Qt::RightButton);

  QwtPlotPanner *panner = new QwtPlotPanner(canvas());
  panner->setMouseButton(Qt::MiddleButton);
}


void AZUREPlot::setXAxisLog(bool set) {
  setAxisAutoScale(QwtPlot::xBottom, true);
  setAxisAutoScale(QwtPlot::yLeft, true);
  if (set)
    setAxisScaleEngine(QwtPlot::xBottom, new QwtLogScaleEngine);
  else
    setAxisScaleEngine(QwtPlot::xBottom, new QwtLinearScaleEngine);
  zoomer->setZoomBase(false);
};

void AZUREPlot::setYAxisLog(bool set) {
  setAxisAutoScale(QwtPlot::xBottom, true);
  setAxisAutoScale(QwtPlot::yLeft, true);
  if (set) {
    QwtLogScaleEngine *scaleEngine = new QwtLogScaleEngine;
    scaleEngine->setMargins(0.5, 0.5);
    setAxisScaleEngine(QwtPlot::yLeft, scaleEngine);
  } else
    setAxisScaleEngine(QwtPlot::yLeft, new QwtLinearScaleEngine);
  zoomer->setZoomBase(false);
};

void AZUREPlot::setXAxisType(unsigned int type) {
  QwtText text;
  if (type == 0)
    text = QwtText(QString("Center of Mass Energy [MeV]"));
  else if (type == 1)
    text = QwtText(QString("Excitation Energy [MeV]"));
  else
    text = QwtText(QString("Center of Mass Angle [degrees]"));
  setAxisTitle(QwtPlot::xBottom, text);

  xAxisType = type;
  update();
}

void AZUREPlot::setYAxisType(unsigned int type) {
  QwtText text = (type == 0) ? QwtText(QString("Cross Section [b]")) : QwtText(QString("S-Factor [MeV b]"));
  setAxisTitle(QwtPlot::yLeft, text);
  yAxisType = type;
  update();
}

void AZUREPlot::setGridVisible(bool visible) {
  grid->setVisible(visible);
  replot();
}

void AZUREPlot::setLevelsModel(LevelsModel *model) {
  levelsModel = model;
}

void AZUREPlot::setLevelsVisible(bool visible) {
  levelsVisible = visible;
  refreshLevelMarkers();
}

void AZUREPlot::clearLevelMarkers() {
  for (int i = 0; i < levelMarkers.size(); i++) {
    levelMarkers[i]->detach();
    delete levelMarkers[i];
  }
  levelMarkers.clear();
}

void AZUREPlot::refreshLevelMarkers() {
  clearLevelMarkers();
  if (!levelsVisible || !levelsModel) {
    replot();
    return;
  }
  // Levels are not meaningful when plotting against an angle.
  if (xAxisType == 2) {
    replot();
    return;
  }

  // Collect (entranceKey -> threshold) for CoM-energy axis.
  // threshold = excitationEnergy - CoM energy, taken from the first plot point
  // of any loaded entry with that entrance pair.
  QMap<int, double> thresholdByEntrance;
  if (xAxisType == 0) {
    for (int i = 0; i < entries.size(); i++) {
      PlotEntry *e = entries[i];
      if (e->points_.isEmpty()) continue;
      const PlotPoint &p = e->points_.first();
      double thr = p.excitationEnergy - p.energy;
      if (!thresholdByEntrance.contains(e->entranceKey())) {
        thresholdByEntrance.insert(e->entranceKey(), thr);
      }
    }
    // If no entries are loaded yet, fall back to threshold=0 so the user still
    // sees something at the excitation-energy values.
    if (thresholdByEntrance.isEmpty()) thresholdByEntrance.insert(-1, 0.);
  }

  QList<LevelsData> levels = levelsModel->getLevels();
  for (int i = 0; i < levels.size(); i++) {
    const LevelsData &lvl = levels[i];
    QString jpi = levelsModel->getSpinLabel(lvl);

    QList<double> xPositions;
    if (xAxisType == 1) {
      xPositions.append(lvl.energy);
    } else {  // xAxisType==0
      QMap<int, double>::const_iterator it = thresholdByEntrance.constBegin();
      for (; it != thresholdByEntrance.constEnd(); ++it) {
        xPositions.append(lvl.energy - it.value());
      }
    }

    for (int k = 0; k < xPositions.size(); k++) {
      QwtPlotMarker *m = new QwtPlotMarker;
      m->setLineStyle(QwtPlotMarker::VLine);
      QPen pen(QColor(80, 80, 80));
      pen.setStyle(Qt::DashLine);
      m->setLinePen(pen);
      m->setXValue(xPositions[k]);
      QwtText t(jpi);
      t.setColor(QColor(40, 40, 40));
      QFont f = t.font();
      f.setPointSize(f.pointSize() > 0 ? f.pointSize() : 9);
      f.setBold(true);
      t.setFont(f);
      m->setLabel(t);
      m->setLabelAlignment(Qt::AlignRight | Qt::AlignTop);
      m->setLabelOrientation(Qt::Vertical);
      m->setItemAttribute(QwtPlotItem::Legend, false);
      m->attach(this);
      levelMarkers.append(m);
    }
  }
  replot();
}

void AZUREPlot::setLegendVisible(bool visible) {
  if (visible) {
    if (!legend) {
      legend = new QwtLegend;
      legend->setDefaultItemMode(QwtLegendData::ReadOnly);
    }
    insertLegend(legend, QwtPlot::RightLegend);
  } else {
    insertLegend(NULL);
    legend = NULL;
  }
  replot();
}

void AZUREPlot::setBandVisible(bool visible) {
  bandVisible = visible;
  redrawEntries();
}

void AZUREPlot::draw(QList<PlotEntry *> newEntries) {
  clearEntries();

  int paletteIndex = 0;
  int symbolIndex = 0;
  for (int i = 0; i < newEntries.size(); i++) {
    if (newEntries[i]->readData()) {
      // Assign default color/symbol if not already customized.
      if (!newEntries[i]->color().isValid() || newEntries[i]->color() == Qt::black) {
        newEntries[i]->setColor(kDefaultPalette[paletteIndex % kDefaultPaletteSize]);
        paletteIndex++;
      }
      if (newEntries[i]->type() == 0) {
        newEntries[i]->setSymbolStyle(kSymbolCycle[symbolIndex % kSymbolCycleSize]);
        symbolIndex++;
      } else {
        // For calculation-only entries, the fit takes color_ (no separate
        // data color exists). Mirror it into fitColor_ so attach() works
        // uniformly and the dialog can edit it via the data-color field.
        if (!newEntries[i]->fitColor().isValid())
          newEntries[i]->setFitColor(newEntries[i]->color());
      }
      newEntries[i]->attach(this, xAxisType, yAxisType, bandVisible);
      entries.push_back(newEntries[i]);
    } else
      delete newEntries[i];
  }
  setAxisAutoScale(QwtPlot::xBottom, true);
  setAxisAutoScale(QwtPlot::yLeft, true);
  refreshLevelMarkers();
  replot();
  zoomer->setZoomBase(false);
}

void AZUREPlot::redrawEntries() {
  for (int i = 0; i < entries.size(); i++) {
    entries[i]->detach();
    entries[i]->attach(this, xAxisType, yAxisType, bandVisible);
  }
  refreshLevelMarkers();
  replot();
  zoomer->setZoomBase(false);
}

void AZUREPlot::update() {
  setAxisAutoScale(QwtPlot::xBottom, true);
  setAxisAutoScale(QwtPlot::yLeft, true);
  redrawEntries();
}

void AZUREPlot::exportPlot() {
#ifndef QT_NO_PRINTER
  QString fileName = "AZUREPlot.pdf";
#else
  QString fileName = "AZUREPlot.png";
#endif

#ifndef QT_NO_FILEDIALOG
  const QList<QByteArray> imageFormats =
      QImageWriter::supportedImageFormats();

  QStringList filter;
  filter += "PDF Documents (*.pdf)";
#ifndef QWT_NO_SVG
  filter += "SVG Documents (*.svg)";
#endif
  filter += "Postscript Documents (*.ps)";

  if (imageFormats.size() > 0) {
    QString imageFilter("Images (");
    for (int i = 0; i < imageFormats.size(); i++) {
      if (i > 0) imageFilter += " ";
      imageFilter += "*.";
      imageFilter += imageFormats[i];
    }
    imageFilter += ")";

    filter += imageFilter;
  }
  fileName = QFileDialog::getSaveFileName(
      this, "Export File Name", fileName,
      filter.join(";;"), NULL, QFileDialog::DontConfirmOverwrite);
#endif
  if (!fileName.isEmpty()) {
    QwtPlotRenderer renderer;
    renderer.setDiscardFlag(QwtPlotRenderer::DiscardBackground, true);
    renderer.renderDocument(this, fileName, QSizeF(300, 200), 85);
  }
}

void AZUREPlot::print() {
  QPrinter printer(QPrinter::HighResolution);
  QString docName("AZUREPlot");
  printer.setDocName(docName);

  printer.setCreator("AZURE2");
  printer.setPageOrientation(QPageLayout::Landscape);

  QPrintDialog dialog(&printer);
  if (dialog.exec()) {
    QwtPlotRenderer renderer;
    renderer.renderTo(this, printer);
  }
}

void AZUREPlot::clearEntries() {
  for (int i = 0; i < entries.size(); i++) {
    entries[i]->detach();
    delete entries[i];
  }
  entries.clear();
  refreshLevelMarkers();
  replot();
}
