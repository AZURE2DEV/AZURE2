#ifndef NUCLEARPOTENTIALTAB_H
#define NUCLEARPOTENTIALTAB_H

#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QGroupBox>
#include <QPushButton>
#include <QLabel>
#include <QCheckBox>

#include "NuclearPotentialManager.h"

class Config;
class PairsModel;
class QTextStream;

/*!
 * @brief GUI tab for configuring nuclear potential parameters
 *
 * A nuclear potential belongs to a particle pair -- it bends the radial wave
 * functions of that channel and no other -- so the tab edits one pair at a
 * time.  The selector at the top chooses which: "Default" is the setting every
 * pair inherits when it has none of its own, and each pair below it can be
 * switched on, off, or given a different shape independently.
 *
 * This tab allows users to:
 * - Choose which particle pair (or the default) is being edited
 * - Switch the hybrid model on or off for that pair alone
 * - Select the type of nuclear potential (Woods-Saxon, Gaussian)
 * - Set parameters for the selected potential
 * - Apply settings to the CoulFunc calculations of that pair
 */
class NuclearPotentialTab : public QWidget {
  Q_OBJECT

public:
  /*!
   * @brief Constructor
   */
  NuclearPotentialTab(QWidget* parent = nullptr);

  /*!
   * @brief Give the tab the pair list it offers in the selector.
   * A pair's key is its 1-based position, which is what the .azr stores and
   * what PPair::GetPairKey() reports.
   */
  void setPairsModel(PairsModel* model);

  /*!
   * @brief Read nuclear potential settings from text stream
   * Used for loading configuration from file
   */
  bool readPotentialSettings(QTextStream& inStream, Config& config);

  /*!
   * @brief Write nuclear potential settings to text stream
   * Used for saving configuration to file
   */
  bool writePotentialSettings(QTextStream& outStream);

private slots:
  /*!
   * @brief Handle potential type selection change
   */
  void onPotentialTypeChanged(int index);

  /*!
   * @brief Switch to editing another pair.  What is on screen is committed to
   * the pair being left, so switching never discards typing.
   */
  void onPairSelectionChanged(int index);

  /*!
   * @brief Turn the hybrid model on or off for the pair on screen
   */
  void onEnabledToggled(bool enabled);

  /*!
   * @brief Drop this pair's own setting so it follows the default again
   */
  void onUseDefaultForPair();

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
   * @brief Rebuild the pair selector from the pairs model
   */
  void refreshPairCombo();

  /*!
   * @brief Load the setting for one pair key (0 = the default) into the widgets
   */
  void loadSettingsFor(int pairKey);

  /*!
   * @brief Write what is on screen into the manager for the pair on screen.
   * Returns false (and warns) if the values do not validate.
   */
  bool commitCurrent();

  /*!
   * @brief Restate, under the buttons, which pairs are switched on
   */
  void refreshSummary();

  /*!
   * @brief The setting currently on screen
   */
  NuclearPotentialSetting settingFromWidgets() const;

  /*!
   * @brief Validate parameter values
   */
  bool validateParameters();

  // Which pair is being edited
  PairsModel* pairsModel_ = nullptr;
  QLabel* pairLabel_;
  QComboBox* pairCombo_;
  QCheckBox* enabledCheck_;
  QLabel* summaryLabel_;
  QPushButton* useDefaultButton_;
  ///The pair key on screen; 0 is the default every unnamed pair inherits.
  int currentPairKey_ = 0;
  ///Set while the widgets are being filled, so the change slots stay quiet.
  bool loading_ = false;

  // Potential type selection
  QLabel* potentialTypeLabel_;
  QComboBox* potentialTypeCombo_;

  // Woods-Saxon parameters
  QGroupBox* woodsSaxonGroup_;
  QLabel* v0Label_;
  QLineEdit* v0Input_;
  QLabel* rLabel_;
  QLineEdit* rInput_;
  QLabel* aLabel_;
  QLineEdit* aInput_;
  QLabel* v0UnitLabel_;
  QLabel* rUnitLabel_;
  QLabel* aUnitLabel_;

  // Gaussian parameters
  QGroupBox* gaussianGroup_;
  QLabel* gaussV0Label_;
  QLineEdit* gaussV0Input_;
  QLabel* gaussR0Label_;
  QLineEdit* gaussR0Input_;
  QLabel* gaussV0UnitLabel_;
  QLabel* gaussR0UnitLabel_;

  // Control buttons
  QPushButton* applyButton_;
  QPushButton* resetButton_;
};

#endif // NUCLEARPOTENTIALTAB_H
