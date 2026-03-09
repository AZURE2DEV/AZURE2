#ifndef INFODIALOG_H
#define INFODIALOG_H

#include <QDialog>

class InfoDialog : public QDialog {
Q_OBJECT

 public:
  InfoDialog(const QString&, QWidget* parent=0, QString title="");
  
};

#endif
