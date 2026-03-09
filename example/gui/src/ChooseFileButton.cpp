#include "ChooseFileButton.h"

ChooseFileButton::ChooseFileButton(const QString& text, QWidget *parent) :
  QPushButton(text, parent) {
    connect(this,SIGNAL(clicked()),this,SLOT(click()));
};

void ChooseFileButton::setLineEdit(QLineEdit* lineEdit) {
  thisLineEdit=lineEdit;
};

void ChooseFileButton::click() {
  emit(clicked(thisLineEdit));
};

