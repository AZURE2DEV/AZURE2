#include <QMainWindow>
#include <QApplication>
#include <QResource>
#include <QOperatingSystemVersion>
#include <clocale>
#include <cstdlib>
#include <cstring>
#include "AZURESetup.h"
#include "AzureTheme.h"
#include "FileDialogPolicy.h"
#include <QTimer>
#include <QIcon>
#include <QDir>

inline void initResource() { Q_INIT_RESOURCE(AZURESetup); }

// The setup program AZURESetup: `AZURESetup [project.azr]`.  AZURE2_SCREENSHOT_DIR=<dir> grabs every tab of the
// opened project into <dir> and quits (QT_QPA_PLATFORM=offscreen needs no display).
int main(int argc, char *argv[]) {
#ifdef Q_OS_MACOS
  if (QOperatingSystemVersion::current() > QOperatingSystemVersion(QOperatingSystemVersion::MacOS, 10, 8)) {
    // fix Mac OS X 10.9 (mavericks) font issue
    // https://bugreports.qt-project.org/browse/QTBUG-32789
    QFont::insertSubstitution(".Lucida Grande UI", "Lucida Grande");
  }
  // Fix for macOS font warnings with .AppleSystemUIFont
  QFont::insertSubstitution(".AppleSystemUIFont", "Helvetica Neue");
  QFont::insertSubstitution(".SF NS Text", "Helvetica Neue");
#endif

  // Qt's own file dialogs on macOS, where the native ones failed (FileDialogPolicy.h). Set before the
  // application exists, independent of the theme (AZURE2_THEME=system keeps it).
  applyFileDialogPolicy();

  QApplication app(argc, argv);
  // numbers with a decimal point whatever the user's language (see AZURESetup's constructor)
  std::setlocale(LC_NUMERIC, "C");

#ifdef WIN_SPACING
  QFont font = app.font();
  font.setPointSizeF(10.5);
  app.setFont(font);
#endif

  initResource();
  applyAzureTheme(app);
  app.setWindowIcon(QIcon(":/azure-icon.png"));   // title bar, task bar, dock
  QCoreApplication::setOrganizationName("jina");
  QCoreApplication::setApplicationName("azure2");

  AZURESetup azureSetup;
  azureSetup.show();

  QString filename = "";
  for (int i = 1; i < argc; i++)
    if (strncmp(argv[i], "--", 2) != 0) {
      filename = argv[i];
      break;
    }
  if (!filename.trimmed().isEmpty()) azureSetup.open(filename);

  // AZURE2_SCREENSHOT_DIR=<dir>: save every tab of the opened project as <dir>/tab<n>_<name>.png and exit.
  // For documentation figures and layout reviews without a display (run under xvfb-run).
  if (const char *shotDir = std::getenv("AZURE2_SCREENSHOT_DIR")) {
    QString dir = QString::fromLocal8Bit(shotDir);
    QTimer::singleShot(1500, [&azureSetup, &app, dir]() {
      azureSetup.saveTabScreenshots(dir);
      app.quit();
    });
  }

  return app.exec();
}
