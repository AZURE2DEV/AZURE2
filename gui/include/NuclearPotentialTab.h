#ifndef NUCLEARPOTENTIALTAB_H
#define NUCLEARPOTENTIALTAB_H

#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QGroupBox>
#include <QPushButton>
#include <QLabel>

class Config;
class QTextStream;

/*!
 * @brief GUI tab for configuring nuclear potential parameters
 *
 * This tab allows users to:
 * - Select the type of nuclear potential (Woods-Saxon, Gaussian)
 * - Set parameters for the selected potential
 * - Preview the potential shape
 * - Apply settings to all CoulFunc calculations
 */
class NuclearPotentialTab : public QWidget {
  Q_OBJECT

 public:
  /*!
   * @brief Constructor
   */
  NuclearPotentialTab(QWidget *parent = nullptr);

  /*!
   * @brief Read nuclear potential settings from text stream
   * Used for loading configuration from file
   */
  bool readPotentialSettings(QTextStream &inStream, Config &config);

  /*!
   * @brief Write nuclear potential settings to text stream
   * Used for saving configuration to file
   */
  bool writePotentialSettings(QTextStream &outStream);

 private slots:
  /*!
   * @brief Handle potential type selection change
   */
  void onPotentialTypeChanged(int index);

  /*!
   * @brief Apply current settings to the global nuclear potential manager
   */
  void onApplySettings();

  /*!
   * @brief Reset to default Woods-Saxon potential
   */
  void onResetToDefault();

  /*!
   * @brief Update GUI when values change
   */
  void onParameterChanged();

 private:
  /*!
   * @brief Create UI components
   */
  void createUI();

  /*!
   * @brief Update parameter labels and descriptions based on potential type
   */
  void updateParameterLabels();

  /*!
   * @brief Load current settings from nuclear potential manager into GUI
   */
  void loadCurrentSettings();

  /*!
   * @brief Validate parameter values
   */
  bool validateParameters();

  // Potential type selection
  QLabel *potentialTypeLabel_;
  QComboBox *potentialTypeCombo_;

  // Woods-Saxon parameters
  QGroupBox *woodsSaxonGroup_;
  QLabel *v0Label_;
  QLineEdit *v0Input_;
  QLabel *rLabel_;
  QLineEdit *rInput_;
  QLabel *aLabel_;
  QLineEdit *aInput_;
  QLabel *v0UnitLabel_;
  QLabel *rUnitLabel_;
  QLabel *aUnitLabel_;

  // Gaussian parameters
  QGroupBox *gaussianGroup_;
  QLabel *gaussV0Label_;
  QLineEdit *gaussV0Input_;
  QLabel *gaussR0Label_;
  QLineEdit *gaussR0Input_;
  QLabel *gaussV0UnitLabel_;
  QLabel *gaussR0UnitLabel_;

  // Control buttons
  QPushButton *applyButton_;
  QPushButton *resetButton_;
};

#endif  // NUCLEARPOTENTIALTAB_H
