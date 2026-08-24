#ifndef ADDPAIRDIALOG_H
#define ADDPAIRDIALOG_H

#include <QDialog>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;
class QCheckBox;
class QGroupBox;
class QComboBox;
class QPushButton;

QT_END_NAMESPACE

class AddPairDialog : public QDialog {
  Q_OBJECT

 public:
  AddPairDialog(QWidget *parent = 0);
  QLineEdit *lightJText;
  QComboBox *lightPiCombo;
  QLineEdit *lightZText;
  QLineEdit *lightMText;
  // QLineEdit *lightGText;
  QLineEdit *heavyJText;
  QComboBox *heavyPiCombo;
  QLineEdit *heavyZText;
  QLineEdit *heavyMText;
  // QLineEdit *heavyGText;
  QLineEdit *excitationEnergyText;
  QLineEdit *seperationEnergyText;
  QLineEdit *channelRadiusText;
  QLineEdit *bindingEnergyText;
  QComboBox *pairTypeCombo;
  QCheckBox *e1Check;
  // QCheckBox *m1Check;
  QCheckBox *e2Check;
  QGroupBox *multBox;

 public slots:
  void updateLightParticle(int index);

 private:
  QLabel *lightJLabel;
  QLabel *lightPiLabel;
  QLabel *lightZLabel;
  QLabel *lightMLabel;
  // QLabel *lightGLabel;
  QLabel *heavyJLabel;
  QLabel *heavyPiLabel;
  QLabel *heavyZLabel;
  QLabel *heavyMLabel;
  // QLabel *heavyGLabel;
  QLabel *excitationEnergyLabel;
  QLabel *seperationEnergyLabel;
  QLabel *channelRadiusLabel;
  QLabel *bindingEnergyLabel;
  QLabel *pairTypeLabel;

  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
