#include <QtGui>
#include <iostream>
#include "AZURESetup.h"

inline void initResource() { Q_INIT_RESOURCE(AZURESetup); }

int start_gui(int argc, char *argv[]) {

#ifdef Q_OS_MACX
    if ( QSysInfo::MacintoshVersion > QSysInfo::MV_10_8 )
    {
        // fix Mac OS X 10.9 (mavericks) font issue
        // https://bugreports.qt-project.org/browse/QTBUG-32789
        QFont::insertSubstitution(".Lucida Grande UI", "Lucida Grande");
    }
#endif

  QApplication app(argc, argv);
 
#ifdef WIN_SPACING
  QFont font = app.font();
  font.setPointSizeF(10.5);
  app.setFont(font);
#endif
 
  initResource();
  QCoreApplication::setOrganizationName("jina");
  QCoreApplication::setApplicationName("azure2");

  AZURESetup azureSetup;
  azureSetup.show();
  
  QString filename="";
  for(int i=1;i<argc;i++) 
    if(strncmp(argv[i],"--",2)!=0) {
      filename=argv[i];
      break;
    }
  if(!filename.trimmed().isEmpty()) azureSetup.open(filename);
  
  return app.exec();
}
