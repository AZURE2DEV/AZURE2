#ifndef CHOOSEFILEBUTTON_H
#define CHOOSEFILEBUTTON_H

#include <QPushButton>
#include <QWidget>
#include <QLineEdit>
#include <QString>

class ChooseFileButton : public QPushButton {

  Q_OBJECT;

 public:
  ChooseFileButton(const QString& text, QWidget *parent = 0);
  void setLineEdit(QLineEdit* lineEdit);

 public slots:
  void click();

 signals:
  void clicked(QLineEdit *lineEdit); 

 private: 
  QLineEdit *thisLineEdit;

};

#endif
