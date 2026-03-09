#ifndef TEXTEDITBUFFER_H
#define TEXTEDITBUFFER_H

#include <QWidget>
#include <QString>
#include <string>
#include <streambuf>
#include <assert.h>
#include <iostream>

class QTextEdit;

class TextEditBuffer : public QWidget, public std::streambuf {

Q_OBJECT

 public:
  TextEditBuffer(std::size_t buff_size = 256, QWidget* parent=0) : QWidget(parent), buffer_(buff_size+1) {
    char* base = &buffer_.front();
    setp(base, base + buffer_.size() -1);
  };
 signals:
  void updateLog(QString);

 protected:
  virtual int_type overflow(int_type ch) {
    if(ch !=traits_type::eof()){
      assert(std::less_equal<char*>()(pptr(),epptr()));
      *pptr()=ch;
      pbump(1);
      if(writeToTextEdit()) return traits_type::to_int_type(ch);
    }
    return traits_type::eof();
  };
  virtual int sync() {
    return writeToTextEdit() ? 0 : -1;
  };
 private:
  bool writeToTextEdit() {
    std::string tempString(pbase(),pptr());
    pbump(pbase()-pptr());
    if(!tempString.empty()) emit updateLog(QString::fromStdString(tempString));
    return true;
  };
 private:
  std::vector<char> buffer_;
};

#endif
