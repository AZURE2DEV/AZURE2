#ifndef CHANNELDETAILS_H
#define CHANNELDETAILS_H

#include <QWidget>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;
class QPushButton;

QT_END_NAMESPACE

/*!
 * Per-channel detail panel beside the channels table.
 */
class ChannelDetails : public QWidget {
  Q_OBJECT

 public:
  ChannelDetails(QWidget *parent = 0);
  void setNormParam(int which);
  QLineEdit *reducedWidthText;
  QLabel *details;
  QPushButton *wignerButton;
  QLineEdit *wignerLimitText;

 private:
  QLabel *normParam;
  QLabel *normUnits;
};

#endif
