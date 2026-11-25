#ifndef NOTESBAR_H
#define NOTESBAR_H

#include <QStringList>
#include <QToolBar>

QT_BEGIN_NAMESPACE
class QAction;
class QDialog;
class QLabel;
class QPlainTextEdit;
class QPushButton;
QT_END_NAMESPACE

/*!
 * The notes about a project just opened, in one bar under the tool bar (docs/upgrade/qt_gui.md section 19.2): what the
 * engine's check found, what the converter assumed, what the tabs could not read.  Never a modal box: the bar says
 * how many notes there are and the first one, "Details" lists them all in a window that does not block the main one,
 * "Dismiss" hides the bar until the next project is opened.  The window clears it when a project opens.
 */
class NotesBar : public QToolBar {
  Q_OBJECT

 public:
  explicit NotesBar(QWidget *parent = nullptr);
  /// Forget every note and hide the bar (a project opened or closed).
  void clearNotes();
  /// Add a note: `summary` one line for the bar, `details` the full text for the Details window (empty: the summary).
  /// `diagnostics`: the note is about the engine's check, so Details also offers the Diagnostics.
  /// `warning`: something the user must know before saving (a section the tabs could not read): put first, the bar
  /// shows the warning icon.
  void addNote(const QString &summary, const QString &details = QString(), bool diagnostics = false, bool warning = false);
  bool hasWarning() const { return warning_; }
  QStringList notes() const { return summaries_; }
  /// The Details window (null before the first click), for the tests.
  QDialog *detailsWindow() const { return details_; }
  QAction *detailsAction() const { return detailsAction_; }
  QAction *dismissAction() const { return dismissAction_; }
  /// "Locate data files...": shown while a note says that data files are not found (the window answers it).
  QAction *locateAction() const { return locateAction_; }
  void showLocate(bool on);

 signals:
  void diagnosticsRequested();
  void locateRequested();

 private:
  void showDetails();
  void updateText();

  QStringList summaries_;
  QStringList texts_;
  bool aboutDiagnostics_ = false;
  bool warning_ = false;
  QLabel *icon_ = nullptr;
  QLabel *text_ = nullptr;
  QAction *detailsAction_ = nullptr;
  QAction *dismissAction_ = nullptr;
  QAction *locateAction_ = nullptr;
  QDialog *details_ = nullptr;
  QPlainTextEdit *detailsText_ = nullptr;
  QPushButton *diagnosticsButton_ = nullptr;
};

#endif
