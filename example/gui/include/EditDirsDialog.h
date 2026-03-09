#ifndef EDITDIRSDIALOG_H
#define EDITDIRSDIALOG_H

#include <QDialog>

QT_BEGIN_NAMESPACE

class QPushButton;
class QGroupBox;
class QComboBox;
class QLineEdit;

QT_END_NAMESPACE

class EditDirsDialog : public QDialog {

Q_OBJECT

 public:
  EditDirsDialog(QWidget *parent = 0);
  QLineEdit *outputDirectoryText;
  QLineEdit *checksDirectoryText;

 private slots:
  void setChooseDirectory(QLineEdit*);

 private:
  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
