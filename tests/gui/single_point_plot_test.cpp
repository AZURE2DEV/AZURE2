// The plotter with a one-point segment (2026-11-15): a calculation of one point has no line to draw, so it is shown
// as an open marker (its uncertainty band as a bar); with two points or more it stays a line.  Headless
// (QT_QPA_PLATFORM=offscreen).  AZURE2_SINGLE_POINT_PNG=<file> saves the one-point plot for a look.
#include <QApplication>
#include <QFile>
#include <QTemporaryDir>
#include <QTextStream>
#include <iostream>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_symbol.h>
#include <vector>
#include "AZUREPlot.h"
#include "Config.h"

// The GUI library calls into the command-line front end; a test links without it (as the other GUI tests).
struct SegPairs {int firstPair; int secondPair;};
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

static int failures = 0;
static void ok(const char *what, bool pass) {
  std::cout << (pass ? "  ok    " : "  FAIL  ") << what << std::endl;
  if (!pass) failures++;
}

// The calculation curve of a plot: the curve drawn on top (z = 30, AZUREPlot.cpp).
static const QwtPlotCurve *fitCurve(QwtPlot &plot) {
  for (QwtPlotItem *item : plot.itemList(QwtPlotItem::Rtti_PlotCurve))
    if (item->z() == 30) return static_cast<const QwtPlotCurve *>(item);
  return nullptr;
}

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QTemporaryDir tmp;
  // AZUREOut columns: E, Ex, angle, fit sigma, fit S, data sigma, data error, data S, data S error; the block ends
  // with blank lines.
  const QString one = tmp.filePath("one.out"), two = tmp.filePath("two.out");
  {
    QFile f(one);
    f.open(QIODevice::WriteOnly);
    QTextStream(&f) << "1.0 3.0 90 2.0e-3 1.0 2.2e-3 2.0e-4 1.1 0.1\n\n\n";
  }
  {
    QFile f(two);
    f.open(QIODevice::WriteOnly);
    QTextStream(&f) << "1.0 3.0 90 2.0e-3 1.0 2.2e-3 2.0e-4 1.1 0.1\n1.2 3.2 90 3.0e-3 1.4 2.9e-3 2.0e-4 1.3 0.1\n\n\n";
  }

  PlotEntry single(0, 1, 2, 0, one);
  ok("one-point file read", single.readData() && single.numPoints() == 1);
  QwtPlot plot1;
  single.attach(&plot1, 0, 0, false);
  const QwtPlotCurve *c1 = fitCurve(plot1);
  ok("one point: the calculation is a marker, not a line",
     c1 && c1->style() == QwtPlotCurve::NoCurve && c1->symbol() && c1->symbol()->style() == QwtSymbol::Diamond &&
         c1->dataSize() == 1);
  ok("... an open one, larger than the data symbol", c1 && c1->symbol() && c1->symbol()->brush().style() == Qt::NoBrush &&
                                                         c1->symbol()->size().width() > single.symbolSize());

  PlotEntry pair(0, 1, 2, 0, two);
  ok("two-point file read", pair.readData() && pair.numPoints() == 2);
  QwtPlot plot2;
  pair.attach(&plot2, 0, 0, false);
  const QwtPlotCurve *c2 = fitCurve(plot2);
  ok("two points: a line, no marker",
     c2 && c2->style() == QwtPlotCurve::Lines && (!c2->symbol() || c2->symbol()->style() == QwtSymbol::NoSymbol));

  if (qEnvironmentVariableIsSet("AZURE2_SINGLE_POINT_PNG")) {
    plot1.resize(500, 380);
    plot1.replot();
    plot1.show();
    app.processEvents();
    plot1.grab().save(qEnvironmentVariable("AZURE2_SINGLE_POINT_PNG"));
  }
  single.detach();
  pair.detach();
  std::cout << (failures ? "FAILED" : "PASSED") << std::endl;
  return failures ? 1 : 0;
}
