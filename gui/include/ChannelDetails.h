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
  /// A particle channel whose width column is a reduced width amplitude
  /// (field 33, set in the THM workspace): the value is labelled MeV^(1/2).
  void setWidthIsAmplitude(bool amplitude);
  bool widthIsAmplitude() const { return widthIsAmplitude_; }
  QLineEdit *reducedWidthText;
  QLabel *details;
  QPushButton *wignerButton;
  QLineEdit *wignerLimitText;

 private:
  QLabel *normParam;
  QLabel *normUnits;
  int normParamWhich_;
  bool widthIsAmplitude_ = false;
};

#endif
