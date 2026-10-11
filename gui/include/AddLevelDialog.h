#ifndef ADDLEVELDIALOG_H
#define ADDLEVELDIALOG_H

#include <QDialog>
#include <QtWidgets>

QT_BEGIN_NAMESPACE

class QLineEdit;
class QLabel;

QT_END_NAMESPACE

/*!
 * One J^pi the level dialog offers, with the number of channels the current pairs and limits give it.
 */
struct SpinParityOption {
  double jValue;
  int piValue;
  int numChannels;
};

/*!
 * Dialog for adding or editing one level.
 *
 * J^pi is picked from a list rather than typed: the Levels tab passes the values that give the level at least one
 * channel under the current pairs and limits, so a J^pi with the wrong integer/half-integer character cannot be entered.
 */
class AddLevelDialog : public QDialog {
  Q_OBJECT

 public:
  AddLevelDialog(QWidget *parent = 0);
  /*!
   * Fills the J^pi list. An empty list (no pairs yet) falls back to J = 0 to 10 in steps of 1/2, both parities. When
   * haveCurrent is set, the current J^pi is selected, and added (marked as having no channels) if it is not in the list.
   */
  void setSpinParityOptions(const QList<SpinParityOption> &options, double currentJ = 0., int currentPi = 1,
                            bool haveCurrent = false);
  double jValue() const;
  int piValue() const;
  static QString spinParityLabel(double jValue, int piValue);

  QLineEdit *energyText;

 private slots:
  void updateOkButton();

 private:
  void addOption(double jValue, int piValue, const QString &suffix);

  QLabel *spinParityLabelWidget;
  QLabel *energyLabel;
  QComboBox *spinParityCombo;
  QLabel *noteLabel;

  QPushButton *okButton;
  QPushButton *cancelButton;
};

#endif
