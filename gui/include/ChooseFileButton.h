#ifndef CHOOSEFILEBUTTON_H
#define CHOOSEFILEBUTTON_H

#include <QPushButton>
#include <QWidget>
#include <QLineEdit>
#include <QString>

/*!
 * Push button that opens a file dialog and writes the chosen path into an associated line edit.
 */
class ChooseFileButton : public QPushButton {
  Q_OBJECT;

 public:
  ChooseFileButton(const QString &text, QWidget *parent = 0);
  void setLineEdit(QLineEdit *lineEdit);

 public slots:
  void click();

 signals:
  void clicked(QLineEdit *lineEdit);

 private:
  QLineEdit *thisLineEdit;
};

#endif
