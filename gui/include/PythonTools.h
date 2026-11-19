#ifndef PYTHONTOOLS_H
#define PYTHONTOOLS_H

#include <QObject>
#include <QPointer>
#include <QProcess>
#include <QString>
#include <QStringList>

QT_BEGIN_NAMESPACE
class QMenu;
class QWidget;
QT_END_NAMESPACE

/*!
 * The Tools menu (docs/upgrade/qt_gui.md section 9): what lives in pyazr and not in the engine (ENDF-6 export and
 * import, the Coulomb table's build, EXFOR and ENSDF imports, the web page of the project) started as
 * `<python> -m pyazr ...`, a child process whose output goes to the log pane.  Nothing of it is re-implemented in
 * C++.  The interpreter is the one chosen in Configure > Directories (a setting of this computer, not of the
 * project; default python3); every item is disabled, with the reason in its tool tip, when `import pyazr` fails
 * with it.  The web page is `python -m pyazr ui FILE --json`, which binds to 127.0.0.1; its process is stopped with
 * the window.
 */
class PythonTools : public QObject {
  Q_OBJECT

 public:
  explicit PythonTools(QWidget *window);
  ~PythonTools() override;
  /// The interpreter of this computer's settings (QSettings "pythonInterpreter"), default python3.
  static QString interpreter();
  static void setInterpreter(const QString &python);
  /// Whether `import pyazr` works with the interpreter; `why` says what failed (the interpreter's own message).
  bool available(QString &why);
  /// Fill the menu: one item per tool, enabled or disabled with its reason.
  void fillMenu(QMenu *menu, const QString &projectFile);
  /// Start `python -m pyazr <arguments>` in `directory`; false with `why` when one runs already or it cannot start.
  bool start(const QStringList &arguments, const QString &directory, QString &why);
  bool running() const { return process_ && process_->state() != QProcess::NotRunning; }
  /// The web page of the project: `pyazr ui FILE --json`; the URL is opened in the browser when the server says it.
  bool openWebPage(const QString &projectFile, QString &why);
  /// The URL the web page's server gave (empty before).
  QString webUrl() const { return webUrl_; }
  void stopWebPage();

 signals:
  void line(const QString &text, bool transient);
  void finished(int exitCode);
  void webPageReady(const QString &url);

 private:
  QWidget *window_;
  QPointer<QProcess> process_, web_;
  QString checkedInterpreter_, unavailable_;
  bool checked_ = false;
  QString webUrl_, webOutput_;
};

#endif
