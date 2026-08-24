#ifndef CHANNELDETAILS_H
#define CHANNELDETAILS_H

#include <QWidget>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;
class QRadioButton;
class QButtonGroup;
class QPushButton;

QT_END_NAMESPACE

class ChannelDetails : public QWidget {
  Q_OBJECT

 public:
  ChannelDetails(QWidget *parent = 0);
  void setNormParam(int which);
  void setConventionChoice(bool isParticle, bool isRWA);
  QLineEdit *reducedWidthText;
  QLabel *details;
  QRadioButton *physicalButton;
  QRadioButton *rwaButton;
  QPushButton *wignerButton;
  QLineEdit *wignerLimitText;

 private:
  QLabel *normParam;
  QLabel *normUnits;
  QButtonGroup *conventionGroup;
  int normParamWhich_;
};

#endif
