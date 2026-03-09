#ifndef CHANNELDETAILS_H
#define CHANNELDETAILS_H

#include <QWidget>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;

QT_END_NAMESPACE

class ChannelDetails : public QWidget {
  Q_OBJECT

 public:
  ChannelDetails(QWidget *parent = 0);
  void setNormParam(int which);
  QLineEdit *reducedWidthText;
  QLabel *details;
 private:
  QLabel *normParam;
  QLabel *normUnits;
};

#endif
