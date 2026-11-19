#ifndef MCMCTAB_H
#define MCMCTAB_H

#include <QWidget>

#include "ProjectHolder.h"

class PlotWidget;

QT_BEGIN_NAMESPACE
class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTableWidget;
QT_END_NAMESPACE

/*!
 * The MCMC tab (docs/upgrade/qt_gui.md section 7), shown only when the engine's command table has `mcmc`
 * (azure2::front::HasCommand): a command removed from the engine takes the tab with it.  The priors table is the
 * <parameters> block of <mcmc> (the free parameters of a session of the project, azure2::front::ReadParameters,
 * merged with the rows the project has; written by azure2::front::SetSamplerPriors in the grammar of the sampler's
 * reader); the sampling settings are the command's options (--walkers, --steps, --parameter-spread, --energy-spread,
 * --sampler-threads, --basis, --overwrite-samples) from the engine's table of options; Run starts the program's mcmc
 * command in the window's worker thread, its steps arrive as progress records (stage "mcmc") through the facade's
 * observer.  Burn-in is not the engine's: samples.mcmc holds every step, and a reader (pyazr) drops the burn-in.
 *
 * Results (docs/upgrade/qt_gui.md section 11): the chain of the last run (azure2::front::ReadSamplerChain, the
 * engine's reader of samples.mcmc), the burn-in steps to drop chosen here: per parameter the mean, the standard
 * deviation and the 16, 50 and 84 % quantiles over the kept steps and walkers; the acceptance fraction (the share of
 * walker moves after the burn-in that changed the position); a plot of a parameter's trace (every walker against the
 * step), of its histogram, or of the walkers' mean log posterior against the step.
 */
class MCMCTab : public QWidget {
  Q_OBJECT

 public:
  explicit MCMCTab(ProjectHolder *holder, QWidget *parent = nullptr);
  QTableWidget *priorTable() const { return table_; }
  /// Fill the priors table from a session of the saved project; false with `why`.
  bool readParameters(const QString &projectFile, QString &why);
  /// The <mcmc> <parameters> rows of the table (checked rows) into the project, when the table was edited.
  void apply();
  /// The command line options of the sampler's run.
  QStringList options() const;
  void setWalkers(int walkers);
  void setSteps(int steps);
  void setSpread(double percent);
  void setThreads(int threads);
  bool setBasis(const QString &basis);
  void setFreshStart(bool fresh);
  void setRunning(bool running);
  QString acceptance() const;
  int lastStep() const { return lastStep_; }
  /// Read the chain of the last run and fill the Results page; false with `why`.
  bool readResults(QString &why);
  /// The burn-in: the first steps left out of the statistics.
  void setBurnIn(int steps);
  QTableWidget *resultsTable() const { return results_; }
  PlotWidget *resultsPlot() const { return resultsPlot_; }
  /// The share of moves after the burn-in that changed a walker's position (NaN before a chain is read).
  double acceptanceFraction() const { return accepted_; }
  /// What the plot shows: "trace:<parameter index>", "histogram:<parameter index>" or "logprob".
  bool showResultPlot(const QString &what);

 public slots:
  void reload();
  void showProgress(const QString &stage, qlonglong step, qlonglong total, double logPosterior, double logLikelihood);
  void appendLine(const QString &line, bool transient);

 signals:
  void parametersRequested();
  void runRequested();
  void stopRequested();

 private:
  ProjectHolder *holder_;
  QTableWidget *table_;
  QLabel *note_;
  QSpinBox *walkers_, *steps_, *threads_;
  QDoubleSpinBox *spread_, *energySpread_;
  QComboBox *basis_;
  QCheckBox *fresh_;
  QPushButton *run_, *stop_;
  QProgressBar *bar_;
  QLabel *status_, *acceptance_;
  QPlainTextEdit *log_;
  bool dirty_ = false, filling_ = false;
  int lastStep_ = 0;
  // the Results page
  void fillResults();
  QSpinBox *burnIn_ = nullptr;
  QComboBox *plotChoice_ = nullptr;
  QLabel *resultsNote_ = nullptr;
  QTableWidget *results_ = nullptr;
  PlotWidget *resultsPlot_ = nullptr;
  azure2::front::SamplerChain chain_;
  QStringList chainNames_;
  double accepted_;
};

#endif
