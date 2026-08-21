#ifndef INFODIALOG_H
#define INFODIALOG_H

#include <QDialog>

/*!
 * Per-tab help shown by the info button.
 */
class InfoDialog : public QDialog {
  Q_OBJECT

 public:
  InfoDialog(const QString &, QWidget *parent = 0, QString title = "");
};

#endif
