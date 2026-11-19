#ifndef STARTPAGE_H
#define STARTPAGE_H

#include <QStringList>
#include <QWidget>

QT_BEGIN_NAMESPACE
class QListWidget;
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

  /// The example projects (*.azr) of the directories examplePlaces() finds, sorted by folder.
  static QStringList findExamples();
  /// Where examples are looked for: $AZURE2_EXAMPLES; examples/ and tests/upgrade/ beside or above the program (an
  /// installed share/azure2/examples or a build directory inside a checkout); the checkout the program was built from.
  static QStringList examplePlaces();

 signals:
  void newRequested();
  void openRequested();
  void openFile(const QString &file);

 private:
  QListWidget *recent_;
  QListWidget *examples_;
  QWidget *tabNotes_;
  QPushButton *new_, *open_;
};

#endif
