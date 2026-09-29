#ifndef THMDIAGNOSTICSPAGE_H
#define THMDIAGNOSTICSPAGE_H

#include <QWidget>
#include <functional>

#include "ThmDiagnostics.h"

QT_BEGIN_NAMESPACE
class QComboBox;
class QGridLayout;
class QLabel;
class QProgressBar;
class QPushButton;
QT_END_NAMESPACE

class ThmPlotWidget;

/*!
 * The Diagnostics page of the THM workspace: read-only plots, computed by the
 * engine (ThmDiagnostics) for one THM data segment when Compute is pressed --
 * never on an edit.  The engine runs on the project as the workspace would
 * leave it (its pages' values included), in a temporary copy, off the GUI
 * thread; the busy bar runs meanwhile.  Panels:
 *   entrance vertex  |M_l(E)|^2 per entrance l of the chosen J^pi, nodes marked;
 *   HOES / on-shell  the two cross sections of the segment's channel, log scale;
 *   line shape       |N_C|^2 of the levels near the data, and zeta(E)
 *                    (only for an experiment with lineshape=on);
 *   weight           w(E) (only with a weight table);
 *   p_s window       w(p) = |phi(p)|^2 p^2 over the window with its quadrature
 *                    nodes (only for an experiment with ps=), and the vertex
 *                    panel adds <|M_l|^2> over the window, dashed.
 * Every panel carries the reaction as a bold label inside the axes.  Later
 * stages add panels here (distortion), not controls.
 */
class ThmDiagnosticsPage : public QWidget {
  Q_OBJECT

 public:
  /// A THM data segment the page offers: its line and how it is listed.
  struct Target {
    int segment = 0;
    QString text;      ///< combo text
    QString reaction;  ///< rich text, the bold label of the plots
  };
  /// `prepare` fills a request for a segment (the project text and the
  /// engine flags as the workspace would leave them), or returns the reason
  /// it cannot; `targets` lists the segments to offer.
  ThmDiagnosticsPage(std::function<QString(int segment, ThmDiagnosticsRequest &)> prepare,
                     std::function<QList<Target>()> targets, QWidget *parent = 0);
  ~ThmDiagnosticsPage();

  /// Refills the segment list (the workspace calls it when the page is shown).
  void refreshTargets();
  /// Starts a computation for the selected segment (what Compute does).
  void compute();
  /// Computes in the calling thread and shows the result (for the tests).
  bool computeNow();
  bool busy() const { return thread_ != nullptr; }
  const ThmDiagnosticsResult &result() const { return result_; }
  /// Shows a result (the thread's, or computeNow's).
  void showResult(const ThmDiagnosticsResult &result);

  QComboBox *segmentCombo;
  QPushButton *computeButton;
  QProgressBar *busyBar;
  QLabel *statusLabel;
  QComboBox *vertexGroupCombo;  ///< J^pi of the vertex panel
  ThmPlotWidget *vertexPlot;
  ThmPlotWidget *hoesPlot;
  ThmPlotWidget *lineshapePlot;
  ThmPlotWidget *zetaPlot;
  ThmPlotWidget *weightPlot;
  ThmPlotWidget *windowPlot;
  QWidget *vertexPanel;  ///< vertexGroupCombo + vertexPlot

 signals:
  /// Emitted when a computation has finished and its result is shown.
  void computed();

 private slots:
  void threadFinished();
  void drawVertex();

 private:
  void clearPlots(const QString &message);
  /// The optional panels that are shown, in reading order after the first row.
  void layoutPanels();
  QString reactionOf(int segment) const;

  std::function<QString(int, ThmDiagnosticsRequest &)> prepare_;
  std::function<QList<Target>()> targets_;
  QList<Target> list_;
  ThmDiagnosticsThread *thread_ = nullptr;
  ThmDiagnosticsResult result_;
  QString computedText_;  ///< the project text of the shown result
  int computedSegment_ = 0;
  QGridLayout *grid_ = nullptr;
};

#endif
