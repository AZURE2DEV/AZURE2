#ifndef PHYSICSTAB_H
#define PHYSICSTAB_H

#include <atomic>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <QList>
#include <QString>
#include <QWidget>

#include "PhysicsForms.h"
#include "azure2/front/Engine.h"

class HeatmapWidget;
class PhysicsWorker;
class ProjectHolder;
class QCheckBox;
class QComboBox;
class QLabel;
class QListWidget;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QLineEdit;
class QStackedWidget;
class QTableWidget;

/*!
 * The project's session as the Physics tab's project pages share it (docs/upgrade/qt_gui.md section 15.3): which
 * saved file and which parameters, the session's pairs, J^pi groups, poles and THM experiments once it is open, and
 * the run of a page's job on it in the session's worker thread (the session opened or re-opened first when the file,
 * its time stamp or the parameters changed).
 */
class PhysicsContext : public QObject {
  Q_OBJECT

 public:
  enum Parameters { StartValues, LastFit, ChosenFile };
  struct Info {
    bool ok = false;
    QString why;
    std::string key;
    std::string parametrization;
    bool extrapolations = false;
    std::vector<azure2::front::PhysicsPair> pairs;
    std::vector<azure2::front::PhysicsGroup> groups;
    std::vector<azure2::front::PhysicsPole> poles;
    std::vector<std::string> experiments;
  };
  using Job = std::function<bool(azure2::front::PhysicsSession &, const std::atomic<bool> &cancel, std::string &why)>;

  PhysicsContext(ProjectHolder *holder, PhysicsWorker *worker, QObject *parent = nullptr);
  /// The saved project file ("" when the project has not been saved: the pages need one).
  QString projectFile() const;
  /// What the project pages need, in words, when there is no saved project; "" when there is one.
  QString needs() const;
  void setParameters(Parameters choice, const QString &file = QString());
  Parameters parameters() const { return choice_; }
  /// The parameter file in use ("" the start values).
  QString parameterFile() const;
  /// The last fit's param.sav (the facade's LastParameterFile), "" when there is none.
  QString lastFit() const;
  /// What the session is opened on: the file, its time stamp, the parameter file.
  std::string key() const;
  const Info &info() const { return info_; }
  /// The info is for the present key.
  bool current() const { return info_.key == key() && !info_.key.empty(); }
  /// Read the session's pairs, groups, poles and experiments again (after opening it as needed).
  void refresh();
  /// refresh() unless the info is current.
  void ensure();
  /// Run `job` on the session (opened first as needed) in the worker; `done(ok, why)` in the GUI thread.
  void run(const QString &lane, Job job, std::function<void(bool ok, const QString &why)> done);
  /// The pair's label in nuclear notation ("p+12C") from the project's <levels>, else "pair k".
  QString pairLabel(int key) const;
  ProjectHolder *holder() const { return holder_; }

 signals:
  void infoChanged();

 private:
  ProjectHolder *holder_;
  PhysicsWorker *worker_;
  Parameters choice_ = StartValues;
  QString chosenFile_;
  Info info_;
};

namespace physicsqt {

/// "1/2+" for J = 0.5, parity +1.
QString JPi(double J, int parity);

/// A page that needs the project: the message what it needs, or its content.
class ProjectPage : public QWidget {
  Q_OBJECT
 public:
  ProjectPage(PhysicsContext *context, QWidget *parent = nullptr);
  /// The message shown in place of the content ("" shows the content).
  void setNeeds(const QString &text);
  QString needsText() const;
  bool contentShown() const;
  /// The page became visible or the project changed: show the message or the content, read the session's info.
  virtual void activate();

 protected:
  void setContent(QWidget *content);
  virtual void infoChanged() {}
  PhysicsContext *context_;
  QStackedWidget *stack_;
  QLabel *needs_;
  std::string shownKey_;  ///< the session key the page last filled itself from
};

/// Page 1 (Functions) and page 6 (Stopping): a row of the physics table or of the stopping table, swept.
class FunctionsPage : public QWidget {
  Q_OBJECT
 public:
  FunctionsPage(bool stopping, PhysicsWorker *worker, QWidget *parent = nullptr);
  QComboBox *functionChoice() const { return function_; }
  QComboBox *returnChoice() const { return return_; }
  ArgumentForm *form() const { return form_; }
  PlotPanel *panel() const { return panel_; }
  /// Evaluate the form and replace the curves (add: overlay).
  void draw(bool add);
  void clear();
  /// The last evaluation's grid and values (for the tests).
  QVector<double> lastGrid() const { return curves_.isEmpty() ? QVector<double>() : curves_.back().grid; }

 private:
  struct Curve {
    QStringList parts;  ///< the fixed arguments in words: the label is the parts that differ between curves
    QVector<double> grid;
    azure2::front::PhysicsValues values;
  };
  void functionChanged();
  void redraw();
  bool stopping_;
  PhysicsWorker *worker_;
  std::vector<azure2::front::PhysicsFunction> rows_;
  QComboBox *function_, *return_;
  QLabel *doc_;
  ArgumentForm *form_;
  PlotPanel *panel_;
  QList<Curve> curves_;
};

/// Page 2: the session's channel functions of a pair.
class ChannelsPage : public ProjectPage {
  Q_OBJECT
 public:
  explicit ChannelsPage(PhysicsContext *context, QWidget *parent = nullptr);
  QComboBox *pair() const { return pair_; }
  QSpinBox *l() const { return l_; }
  QComboBox *quantity() const { return quantity_; }
  SweepFields *energies() const { return energies_; }
  PlotPanel *panel() const { return panel_; }
  void draw(bool add);

 protected:
  void infoChanged() override;

 private:
  struct Curve {
    QString label;
    azure2::front::PhysicsCurves answer;
  };
  void redraw();
  QComboBox *pair_, *quantity_;
  QSpinBox *l_, *twoJ_;
  QLineEdit *radius_;
  SweepFields *energies_;
  PlotPanel *panel_;
  QList<Curve> curves_;
};

/// Page 3: one J^pi group of the R-matrix with the poles.
class RMatrixPage : public ProjectPage {
  Q_OBJECT
 public:
  explicit RMatrixPage(PhysicsContext *context, QWidget *parent = nullptr);
  QComboBox *group() const { return group_; }
  QComboBox *quantity() const { return quantity_; }
  SweepFields *energies() const { return energies_; }
  QListWidget *seriesList() const { return series_; }
  QTableWidget *poleTable() const { return poles_; }
  PlotPanel *panel() const { return panel_; }
  const azure2::front::PhysicsCurves &answer() const { return answer_; }
  int answerGroup() const { return answerGroup_; }
  void compute();
  /// Zoom the plot on the pole of row `row` of the table.
  void centreOnPole(int row);

 protected:
  void infoChanged() override;

 private:
  void groupChanged();
  void fillSeries();
  void redraw();
  void fillPoles();
  QComboBox *group_, *quantity_;
  SweepFields *energies_;
  QListWidget *series_;
  QTableWidget *poles_;
  QLabel *facts_;
  PlotPanel *panel_;
  azure2::front::PhysicsCurves answer_;
  int answerGroup_ = 0;
};

/// Page 4: the real-energy pole map.
class PoleMapPage : public ProjectPage {
  Q_OBJECT
 public:
  explicit PoleMapPage(PhysicsContext *context, QWidget *parent = nullptr);
  QComboBox *kind() const { return kind_; }
  QComboBox *group() const { return group_; }
  QCheckBox *logColour() const { return log_; }
  SweepFields *energies() const { return energies_; }
  HeatmapWidget *map() const { return map_; }
  PlotPanel *cut() const { return cut_; }
  QLabel *readout() const { return readout_; }
  QLabel *explanation() const { return explanation_; }
  void compute();
  void showCut(int row);

 protected:
  void infoChanged() override;

 private:
  void kindChanged();
  QComboBox *kind_, *group_;
  QCheckBox *log_;
  SweepFields *energies_;
  HeatmapWidget *map_;
  PlotPanel *cut_;
  QLabel *readout_, *explanation_;
  QVector<double> cutMarkers_;
};

/// Page 5: the THM experiments' vertex, line shape and distortion.
class ThmPage : public ProjectPage {
  Q_OBJECT
 public:
  explicit ThmPage(PhysicsContext *context, QWidget *parent = nullptr);
  QComboBox *experiment() const { return experiment_; }
  QComboBox *part() const { return part_; }
  QComboBox *quantity() const { return quantity_; }
  SweepFields *energies() const { return energies_; }
  PlotPanel *panel() const { return panel_; }
  void draw(bool add);

 protected:
  void infoChanged() override;

 private:
  struct Curve {
    QString label;
    azure2::front::PhysicsCurves answer;
  };
  void partChanged();
  void redraw();
  QComboBox *experiment_, *part_, *quantity_;
  SweepFields *energies_;
  QLabel *facts_;
  PlotPanel *panel_;
  QList<Curve> curves_;
};

/// A phase shift in degrees, defined modulo 180, made continuous for drawing (a step of more than 90 degrees removed
/// as the wrap; the curve starts in (-90, 90]).  Display only (qt_gui.md section 19.8).
QVector<double> UnwrappedPhaseDegrees(const QVector<double> &deg);

}  // namespace physicsqt

/*!
 * The Physics tab (docs/upgrade/qt_gui.md section 15; the engine side docs/upgrade/physics_tab.md): six pages in a
 * list -- Functions, Channels of the project, R-matrix of the project, Pole map, THM, Stopping -- each drawn from
 * the facade's physics calls.  Pages 1 and 6 need no project.  The computations run in two worker threads with a busy
 * bar and Cancel; the window never waits.
 */
class PhysicsTab : public QWidget {
  Q_OBJECT

 public:
  enum Page { Functions, Channels, RMatrix, PoleMap, Thm, Stopping };
  explicit PhysicsTab(ProjectHolder *holder, QWidget *parent = nullptr);
  ~PhysicsTab() override;
  void showPage(int page);
  /// Standard or expert view (qt_gui.md section 19.4): the pole map, and the THM page of a project without <thm>,
  /// are not listed in standard view.
  void setExpertMode(bool on);
  QListWidget *pages() const { return pages_; }
  bool busy() const;
  void cancel();
  PhysicsContext *context() const { return context_; }
  QComboBox *parameterChoice() const { return parameters_; }
  physicsqt::FunctionsPage *functionsPage() const { return functions_; }
  physicsqt::FunctionsPage *stoppingPage() const { return stopping_; }
  physicsqt::ChannelsPage *channelsPage() const { return channels_; }
  physicsqt::RMatrixPage *rmatrixPage() const { return rmatrix_; }
  physicsqt::PoleMapPage *poleMapPage() const { return poleMap_; }
  physicsqt::ThmPage *thmPage() const { return thm_; }

 public slots:
  /// The project was opened, closed or saved: the project pages read the session again when shown.
  void projectChanged();

 private:
  void busyChanged();
  void parametersChosen();
  ProjectHolder *holder_;
  PhysicsWorker *sessionWorker_, *freeWorker_;
  PhysicsContext *context_;
  QListWidget *pages_;
  QStackedWidget *stack_;
  QWidget *sessionBar_;
  QComboBox *parameters_;
  QLabel *sessionStatus_;
  QWidget *busyStrip_;
  QProgressBar *busyBar_;
  physicsqt::FunctionsPage *functions_, *stopping_;
  physicsqt::ChannelsPage *channels_;
  physicsqt::RMatrixPage *rmatrix_;
  physicsqt::PoleMapPage *poleMap_;
  physicsqt::ThmPage *thm_;
};

#endif
