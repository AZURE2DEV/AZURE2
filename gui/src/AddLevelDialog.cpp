#include "AddLevelDialog.h"

AddLevelDialog::AddLevelDialog(QWidget *parent) :
  QDialog(parent) {
  energyLabel = new QLabel(tr("Excitation Energy [MeV]:"));
  energyText = new QLineEdit;
  energyText->setMinimumWidth(energyText->fontMetrics().horizontalAdvance(QString(14, QChar('0'))));

  // Rich text for the superscript pi; QLabel renders it, QComboBox items would not.
  spinParityLabelWidget = new QLabel(tr("Spin-Parity (J<sup>%1</sup>):").arg(QChar(0x03C0)));
  spinParityLabelWidget->setTextFormat(Qt::RichText);
  spinParityCombo = new QComboBox;
  spinParityCombo->setSizeAdjustPolicy(QComboBox::AdjustToContents);
  spinParityCombo->setMaxVisibleItems(20);
  spinParityCombo->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  spinParityCombo->setObjectName("spinParityCombo");

  // Both boxes carry object names; tests look the J^pi box up by its name.
  nucleusLabelWidget = new QLabel(tr("Compound Nucleus:"));
  nucleusCombo = new QComboBox;
  nucleusCombo->setObjectName("nucleusCombo");
  nucleusCombo->setSizeAdjustPolicy(QComboBox::AdjustToContents);
  nucleusCombo->setToolTip(tr("The compound nucleus this level belongs to; its channels come from that nucleus's pairs only."));
  nucleusLabelWidget->setVisible(false);
  nucleusCombo->setVisible(false);
  connect(nucleusCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(nucleusChanged(int)));

  noteLabel = new QLabel;
  noteLabel->setWordWrap(true);
  noteLabel->setVisible(false);

  cancelButton = new QPushButton(tr("Cancel"));
  okButton = new QPushButton(tr("Accept"));
  okButton->setDefault(true);

  QFormLayout *formLayout = new QFormLayout;
  formLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  formLayout->addRow(nucleusLabelWidget, nucleusCombo);
  formLayout->addRow(energyLabel, energyText);
  formLayout->addRow(spinParityLabelWidget, spinParityCombo);

  QHBoxLayout *buttonBox = new QHBoxLayout;
  buttonBox->addStretch(1);
  buttonBox->addWidget(cancelButton);
  buttonBox->addWidget(okButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(formLayout);
  mainLayout->addWidget(noteLabel);
  mainLayout->addLayout(buttonBox);
  mainLayout->setSizeConstraint(QLayout::SetFixedSize);

  setLayout(mainLayout);

  connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
  connect(energyText, SIGNAL(textChanged(const QString &)), this, SLOT(updateOkButton()));

  setSpinParityOptions(QList<SpinParityOption>());
  updateOkButton();

  setWindowTitle(tr("Add a Level"));
}

QString AddLevelDialog::spinParityLabel(double jValue, int piValue) {
  // Same notation as the Levels table (LevelsModel::getSpinLabel).
  QString spin;
  if (((int)(jValue * 2)) % 2 != 0 && jValue != 0.)
    spin = QString("%1/2").arg((int)(jValue * 2));
  else
    spin = QString("%1").arg(jValue);
  // Superscript plus/minus (U+207A/B): plain text, so it renders in the combo box without a rich-text delegate.
  return spin + (piValue == -1 ? QString(QChar(0x207B)) : QString(QChar(0x207A)));
}

void AddLevelDialog::addOption(double jValue, int piValue, const QString &suffix) {
  QString text = spinParityLabel(jValue, piValue);
  if (!suffix.isEmpty()) text += "  " + suffix;
  spinParityCombo->addItem(text);
  int row = spinParityCombo->count() - 1;
  spinParityCombo->setItemData(row, jValue, Qt::UserRole);
  spinParityCombo->setItemData(row, piValue, Qt::UserRole + 1);
}

void AddLevelDialog::setSpinParityOptions(const QList<SpinParityOption> &options, double currentJ, int currentPi,
                                          bool haveCurrent) {
  spinParityCombo->clear();
  QList<SpinParityOption> list = options;
  bool fallback = list.isEmpty();
  if (fallback) {
    for (int twoJ = 0; twoJ <= 20; twoJ++)
      for (int pi : {1, -1}) list.append({twoJ / 2.0, pi, 0});
    noteLabel->setTextFormat(Qt::RichText);
    noteLabel->setText(tr("No particle pairs are defined yet, so every J<sup>%1</sup> up to 10 is listed. "
                          "Define the pairs first to list only the values that give channels.")
                         .arg(QChar(0x03C0)));
  }
  noteLabel->setVisible(fallback);

  int current = -1;
  for (const SpinParityOption &option : list) {
    if (haveCurrent && option.jValue == currentJ && option.piValue == currentPi) current = spinParityCombo->count();
    addOption(option.jValue, option.piValue, QString());
  }
  if (haveCurrent && current == -1) {
    // A level read from a file, or one whose limits have since been lowered: keep its J^pi selectable.
    current = 0;
    for (int i = 0; i < spinParityCombo->count(); i++) {
      double j = spinParityCombo->itemData(i, Qt::UserRole).toDouble();
      int pi = spinParityCombo->itemData(i, Qt::UserRole + 1).toInt();
      if (j < currentJ || (j == currentJ && pi > currentPi)) current = i + 1;
    }
    spinParityCombo->insertItem(current, spinParityLabel(currentJ, currentPi) + "  " + tr("(no channels)"));
    spinParityCombo->setItemData(current, currentJ, Qt::UserRole);
    spinParityCombo->setItemData(current, currentPi, Qt::UserRole + 1);
  }
  spinParityCombo->setCurrentIndex(current == -1 ? 0 : current);
}

double AddLevelDialog::jValue() const { return spinParityCombo->currentData(Qt::UserRole).toDouble(); }

int AddLevelDialog::piValue() const { return spinParityCombo->currentData(Qt::UserRole + 1).toInt(); }

void AddLevelDialog::updateOkButton() {
  bool ok = false;
  energyText->text().trimmed().toDouble(&ok);
  okButton->setEnabled(ok && spinParityCombo->count() > 0);
}

void AddLevelDialog::setNucleusOptions(const QStringList &labels, const QList<QList<SpinParityOption> > &optionsPerNucleus,
                                       int current, bool editable) {
  nucleusCombo->blockSignals(true);
  nucleusCombo->clear();
  nucleusOptions = optionsPerNucleus;
  for (const QString &label : labels) nucleusCombo->addItem(label);
  if (current < 0 || current >= labels.size()) current = 0;
  if (!labels.isEmpty()) nucleusCombo->setCurrentIndex(current);
  nucleusCombo->blockSignals(false);
  bool several = labels.size() > 1;
  nucleusLabelWidget->setVisible(several);
  nucleusCombo->setVisible(several);
  nucleusCombo->setEnabled(editable);
  if (!several || !editable) return;
  // The J^pi list of the chosen nucleus; an edit keeps the list the caller set with its current J^pi.
  nucleusChanged(current);
}

int AddLevelDialog::nucleusIndex() const { return nucleusCombo->count() > 0 ? nucleusCombo->currentIndex() : 0; }

void AddLevelDialog::nucleusChanged(int index) {
  if (index < 0 || index >= nucleusOptions.size()) return;
  setSpinParityOptions(nucleusOptions.at(index));
  updateOkButton();
}
