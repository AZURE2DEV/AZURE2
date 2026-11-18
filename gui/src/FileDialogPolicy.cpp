#include "FileDialogPolicy.h"

#include <QCoreApplication>
#include <QString>
#include <cstdlib>

bool useNativeFileDialogs() {
  if (const char *choice = std::getenv("AZURE2_NATIVE_DIALOGS")) {
    const QString c = QString::fromLocal8Bit(choice).trimmed().toLower();
    if (c == "1" || c == "yes" || c == "on" || c == "true") return true;
    if (c == "0" || c == "no" || c == "off" || c == "false") return false;
  }
#ifdef Q_OS_MACOS
  return false;
#else
  return true;
#endif
}

void applyFileDialogPolicy() {
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, !useNativeFileDialogs());
}
