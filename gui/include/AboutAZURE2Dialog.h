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
  /// polarized: show the logo with the polarization arrow (the project holds a polarization segment).
  AboutAZURE2Dialog(QWidget *parent = 0, bool polarized = false);
};

#endif
