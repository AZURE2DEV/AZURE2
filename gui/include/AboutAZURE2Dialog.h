#ifndef ABOUTAZURE2DIALOG_H
#define ABOUTAZURE2DIALOG_H

#include <QDialog>

QT_BEGIN_NAMESPACE
QT_END_NAMESPACE

/*!
 * The About box.
 */
class AboutAZURE2Dialog : public QDialog {
  Q_OBJECT

 public:
  /// logo: the resource of the logo to show, the one the window shows for its project (AZURESetup::logoFor).
  AboutAZURE2Dialog(QWidget *parent = 0, const QString &logo = ":/azure-icon.png");
};

#endif
