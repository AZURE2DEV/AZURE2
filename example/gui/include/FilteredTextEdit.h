#ifndef FILTEREDTEXTEDIT_H
#define FILTEREDTEXTEDIT_H

#include <iostream>
#include <QtGui>

class FilteredTextEdit : public QTextEdit {
  Q_OBJECT

 public:
  FilteredTextEdit(QWidget *parent = 0) : 
  QTextEdit(parent), filtered_(false) {
    QFont font("Courier");
    font.setStyleHint(QFont::TypeWriter);
    setCurrentFont(font);
  };
  void SetMouseFiltered(bool filtered) { filtered_=filtered;};
  bool IsMouseFiltered() const {return filtered_;};
 public slots:
  void write(QString string) {
    if(string[0]=='\r') {
      QTextCursor cursor = textCursor();
      cursor.select(QTextCursor::LineUnderCursor);
      setTextCursor(cursor);
      string.remove(0,1);
    }  
    insertPlainText(string);
  };
  void mousePressEvent(QMouseEvent *event) {
    if(!filtered_) QTextEdit::mousePressEvent(event);
  };
  void mouseDoubleClickEvent(QMouseEvent *event) {
    if(!filtered_) QTextEdit::mouseDoubleClickEvent(event);
  }
 private:
  bool filtered_;
};

#endif
