// The GUI's file dialogs on macOS are Qt's own, not the native NSOpenPanel (gui/include/FileDialogPolicy.h).
//
// Report of 2026-11-05 from a Mac user: the "Choose..." buttons (Run tab parameter file, directory chooser, ...)
// failed with the native macOS dialog and worked once a call passed QFileDialog::DontUseNativeDialog -- one button
// at a time, so the others stayed broken.  The fix sets Qt::AA_DontUseNativeDialogs for the whole application on
// macOS, which every QFileDialog::get...() call honours.
//
// Checks:
//   (a) the platform default: Qt dialogs on macOS, native ones elsewhere (the attribute after applyFileDialogPolicy());
//   (b) with the attribute set, a QFileDialog builds Qt's own widgets (its list view exists).  On macOS this test runs
//       on the cocoa platform (tests/gui/CMakeLists.txt), where a native helper exists, so (b) fails there on a build
//       without the policy;
//   (c) AZURE2_NATIVE_DIALOGS=1 / 0 force either choice on any platform; on macOS (cocoa) =1 really gives the native
//       dialog (no widgets), showing that (b) discriminates;
//   (d) an unrecognised value falls back to the platform default;
//   (e) the static QFileDialog::getOpenFileName -- the call every "Choose..." button makes -- opens a Qt dialog
//       under the policy (found and cancelled by a timer).
#include <QApplication>
#include <QFileDialog>
#include <QListView>
#include <QTimer>
#include <cstdio>
#include <cstdlib>

#include "FileDialogPolicy.h"

static int failures = 0;
static void check(bool ok, const char *what) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) failures++;
}

static bool qtWidgetsInDialog() {
  QFileDialog dialog;
  return dialog.findChild<QListView *>("listView") != nullptr;
}

static bool attributeSet() { return QCoreApplication::testAttribute(Qt::AA_DontUseNativeDialogs); }

int main(int argc, char **argv) {
#ifdef Q_OS_MACOS
  const bool mac = true;
#else
  const bool mac = false;
#endif
  qunsetenv("AZURE2_NATIVE_DIALOGS");
  applyFileDialogPolicy();  // as main.cpp: before the application exists
  QApplication app(argc, argv);
  std::printf("platform plugin: %s\n", qPrintable(QGuiApplication::platformName()));

  // (a)
  check(attributeSet() == mac, mac ? "(a) macOS default: native dialogs off" : "(a) default here: native dialogs on");
  check(useNativeFileDialogs() == !mac, "(a) useNativeFileDialogs() agrees with the platform default");

  // (b)
  if (mac) check(qtWidgetsInDialog(), "(b) macOS default: a QFileDialog is Qt's own (list view present)");
  qputenv("AZURE2_NATIVE_DIALOGS", "0");
  applyFileDialogPolicy();
  check(attributeSet(), "(c) AZURE2_NATIVE_DIALOGS=0 sets AA_DontUseNativeDialogs");
  check(qtWidgetsInDialog(), "(b) with the attribute set, a QFileDialog is Qt's own (list view present)");

  // (c)
  qputenv("AZURE2_NATIVE_DIALOGS", "1");
  applyFileDialogPolicy();
  check(!attributeSet(), "(c) AZURE2_NATIVE_DIALOGS=1 clears AA_DontUseNativeDialogs");
  if (mac && QGuiApplication::platformName() == "cocoa")
    check(!qtWidgetsInDialog(), "(c) macOS, AZURE2_NATIVE_DIALOGS=1: the native panel is used (no Qt list view)");

  // (d)
  qputenv("AZURE2_NATIVE_DIALOGS", "maybe");
  applyFileDialogPolicy();
  check(attributeSet() == mac, "(d) an unrecognised AZURE2_NATIVE_DIALOGS value gives the platform default");

  // (e) the static call path, under the macOS policy forced on every platform
  qputenv("AZURE2_NATIVE_DIALOGS", "0");
  applyFileDialogPolicy();
  bool sawQtDialog = false;
  QTimer::singleShot(200, [&sawQtDialog]() {
    for (QWidget *w : QApplication::topLevelWidgets())
      if (QFileDialog *d = qobject_cast<QFileDialog *>(w))
        if (d->isVisible()) {
          sawQtDialog = d->findChild<QListView *>("listView") != nullptr;
          d->reject();
        }
  });
  QTimer::singleShot(5000, []() {  // safety net: never hang the test
    for (QWidget *w : QApplication::topLevelWidgets())
      if (QDialog *d = qobject_cast<QDialog *>(w)) d->reject();
  });
  const QString chosen = QFileDialog::getOpenFileName(nullptr, "file_dialog_policy_test");
  check(chosen.isEmpty() && sawQtDialog, "(e) QFileDialog::getOpenFileName opens Qt's dialog under the policy");

  qunsetenv("AZURE2_NATIVE_DIALOGS");
  if (failures) {
    std::printf("FAIL: %d check(s)\n", failures);
    return 1;
  }
  std::printf("PASS: file dialog policy\n");
  return 0;
}
