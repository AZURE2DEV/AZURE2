#ifndef STARTPAGE_H
#define STARTPAGE_H

#include <QStringList>
#include <QWidget>

QT_BEGIN_NAMESPACE
class QListWidget;
class QCheckBox;
class QListWidgetItem;
class QPushButton;
QT_END_NAMESPACE

/*!
 * The window's page while no project is open (docs/upgrade/qt_gui.md section 11): New and Open, the recent projects,
 * the example projects found on this computer (examplePlaces()), and one line per tab saying what it is for.  A
 * click on a project opens it (openFile); the window shows the tabs from then on.
 */
class StartPage : public QWidget {
  Q_OBJECT

 public:
  explicit StartPage(QWidget *parent = nullptr);
  /// The recent projects (newest first); files that are gone are left out.
  void setRecent(const QStringList &files);
  /// The tab titles and what each is for, shown under "The tabs".
  void setTabNotes(const QList<QPair<QString, QString>> &notes);
  QListWidget *recentList() const { return recent_; }
  QListWidget *exampleList() const { return examples_; }
  QPushButton *newButton() const { return new_; }
  QPushButton *openButton() const { return open_; }
  /// "Show expert options" (qt_gui.md section 19.4): the same switch as View > Show Expert Options.
  QCheckBox *expertBox() const { return expert_; }
  void setExpertChecked(bool on);

  /// One line of the examples' manifest (examples/index.json, qt_gui.md section 19.9): a group's heading (project
  /// empty), an entry (one reaction) or a variant of the entry before it.
  struct Example {
    QString group, reaction, source, what, project;  ///< project: absolute path ("" for a heading)
    bool variant = false;
  };
  /// The manifest of the first place (examplePlaces()) that has one, in its order; empty when none has one.  `why`:
  /// what could not be read.
  static QList<Example> readManifest(QString *why = nullptr);
  static QList<Example> readManifest(const QString &root, QString *why = nullptr);
  /// The example projects: the manifest's entries and variants in its order, or, where there is no manifest, the
  /// *.azr of examples/*/ (and of $AZURE2_EXAMPLES/*/), sorted by folder.
  static QStringList findExamples();
  /// Where examples are looked for: $AZURE2_EXAMPLES; a folder holding examples/ beside or above the program (an
  /// installed share/azure2 or a build directory inside a checkout); the checkout the program was built from.
  static QStringList examplePlaces();

 signals:
  void newRequested();
  void openRequested();
  void openFile(const QString &file);
  void expertToggled(bool on);

 private:
  QListWidget *recent_;
  QListWidget *examples_;
  QWidget *tabNotes_;
  QPushButton *new_, *open_;
  QCheckBox *expert_ = nullptr;
};

#endif
