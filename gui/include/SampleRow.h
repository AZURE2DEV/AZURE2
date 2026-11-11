#ifndef SAMPLEROW_H
#define SAMPLEROW_H

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleValidator>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QWidget>

#include "SegmentKeywords.h"

/*!
 * The "Sample:" row of the segment dialogs: the transmission / self-shielded-yield block of the line
 * ("transmission n [vary sigma]", "yield n [vary sigma]", SegLine.h) and, for a data segment, its thickness group
 * ("thicknessgroup <label>").  A combo None / Transmission / Self-shielded yield, the thickness n in atoms/barn and, with
 * \p withFit (data segments), a "Fit" check box (the vary flag of the fit form, n the parameter
 * segment_<key>_thickness), the prior width sigma (atoms/barn, empty = no prior) and an editable thickness-group combo.
 * A test segment has no fit parameter (EData::MakePoints ignores a fit form there), so its row shows kind and n only
 * and keeps whatever fit form the line had.
 *
 * Built without Q_OBJECT (lambda connections), so it needs no moc of its own.
 */
class SampleRow : public QWidget {
 public:
  explicit SampleRow(bool withFit, QWidget *parent = nullptr) : QWidget(parent), withFit_(withFit) {
    kindCombo = new QComboBox;
    kindCombo->setObjectName("sampleKindCombo");
    kindCombo->addItem(tr("None"), QString());
    kindCombo->addItem(tr("Transmission"), QString("transmission"));
    kindCombo->addItem(tr("Self-shielded yield"), QString("yield"));
    kindCombo->setToolTip(tr("A thick-sample neutron observable of the segment's cross section sigma:\n"
                             "Transmission: the data are T = exp(-n sigma) (\"transmission n\" on the line; the\n"
                             "segment is the entrance pair's total cross section, \"total 1\").\n"
                             "Self-shielded yield: the data are Y0 = (sigma/sigma_t)(1 - exp(-n sigma_t)) of a\n"
                             "capture or reaction segment (\"yield n\"), sigma_t the total cross section.\n"
                             "None: the data are the cross section itself."));
    thicknessText = new QLineEdit;
    thicknessText->setObjectName("sampleThicknessText");
    QDoubleValidator *positive = new QDoubleValidator(this);
    positive->setBottom(0.0);
    thicknessText->setValidator(positive);
    thicknessText->setMaximumWidth(90);
    thicknessText->setToolTip(tr("The sample thickness n in atoms/barn (> 0)."));
    QHBoxLayout *layout = new QHBoxLayout;
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(kindCombo);
    layout->addWidget(new QLabel(tr("n [atoms/b]:")));
    layout->addWidget(thicknessText);
    if (withFit_) {
      fitCheck = new QCheckBox(tr("Fit"));
      fitCheck->setObjectName("sampleFitCheck");
      fitCheck->setToolTip(tr("Make n the fit parameter segment_<key>_thickness (the fit form \"transmission n vary\n"
                              "sigma\" with vary 1).  An engine older than November 2026 ignores the fit form."));
      sigmaText = new QLineEdit;
      sigmaText->setObjectName("sampleSigmaText");
      QDoubleValidator *nonNegative = new QDoubleValidator(this);
      nonNegative->setBottom(0.0);
      sigmaText->setValidator(nonNegative);
      sigmaText->setMaximumWidth(80);
      sigmaText->setPlaceholderText(tr("none"));
      sigmaText->setToolTip(tr("Gaussian prior width of the fitted thickness in atoms/barn: the fit adds\n"
                               "((n - n0)/sigma)^2, n0 the value above.  Empty or 0 = no prior."));
      groupCombo = new QComboBox;
      groupCombo->setObjectName("thicknessGroupCombo");
      groupCombo->setEditable(true);
      groupCombo->setInsertPolicy(QComboBox::NoInsert);
      groupCombo->setMinimumContentsLength(6);
      groupCombo->setToolTip(tr("Thickness group (\"thicknessgroup <label>\" on the segment line).\n"
                                "Every active segment with this label uses ONE sample thickness: that of the group's\n"
                                "first active segment in file order (its n, Fit and prior), e.g. a transmission and a\n"
                                "capture-yield measurement on the same sample.  Each member needs its own\n"
                                "transmission or yield.  Empty = no group."));
      layout->addWidget(fitCheck);
      layout->addWidget(new QLabel(tr("Prior sigma:")));
      layout->addWidget(sigmaText);
      layout->addWidget(new QLabel(tr("Group:")));
      layout->addWidget(groupCombo, 1);
    } else {
      layout->addStretch(1);
    }
    setLayout(layout);
    connect(kindCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) { updateEnabled(); });
    updateEnabled();
  }

  /*! Show the block read from the line. */
  void setBlock(const SegmentKeywords::SampleBlock &b) {
    loaded_ = b;
    const int i = kindCombo->findData(b.kind);
    kindCombo->setCurrentIndex(i >= 0 ? i : 0);
    thicknessText->setText(b.n);
    if (withFit_) {
      fitCheck->setChecked(b.fit && b.vary);
      sigmaText->setText(b.fit ? b.sigma : QString());
    }
    updateEnabled();
  }

  /*! The block the row describes.  The fit form is written when Fit is checked, a prior is given or the line had
   *  it (a fixed fitted thickness, vary 0, stays one); a test segment's row keeps the line's fit form as read. */
  SegmentKeywords::SampleBlock block() const {
    SegmentKeywords::SampleBlock b;
    b.kind = kindCombo->currentData().toString();
    if (b.kind.isEmpty()) return b;
    b.n = thicknessText->text().trimmed();
    if (withFit_) {
      const QString sg = sigmaText->text().trimmed();
      b.vary = fitCheck->isChecked();
      b.fit = b.vary || (!sg.isEmpty() && sg.toDouble() != 0.0) || loaded_.fit;
      b.sigma = sg.isEmpty() ? QString("0") : sg;
      if (loaded_.fit && b.sigma == "0" && sg.isEmpty() && !loaded_.sigma.isEmpty() && loaded_.sigma.toDouble() == 0.0)
        b.sigma = loaded_.sigma;  // "0.0" read stays "0.0"
    } else {
      b.fit = loaded_.fit;
      b.vary = loaded_.vary;
      b.sigma = loaded_.sigma;
    }
    return b;
  }

  /*! True when the line the row was filled from had a transmission / yield block. */
  bool hasLoadedBlock() const { return !loaded_.kind.isEmpty(); }
  /*! The thickness-group label (data segments), "" for none. */
  QString groupLabel() const { return (withFit_ && !block().kind.isEmpty()) ? groupCombo->currentText().trimmed() : QString(); }
  /*! Offer the project's thickness-group labels and show the line's own. */
  void setGroupChoices(const QStringList &labels) {
    if (!withFit_) return;
    const QString keep = groupCombo->currentText();
    groupCombo->blockSignals(true);
    groupCombo->clear();
    groupCombo->addItem(QString());
    for (const QString &l : labels)
      if (!l.isEmpty()) groupCombo->addItem(l);
    groupCombo->setCurrentIndex(0);
    groupCombo->setEditText(keep);
    groupCombo->blockSignals(false);
  }
  void setGroupLabel(const QString &label) {
    if (!withFit_) return;
    const int i = groupCombo->findText(label);
    if (i >= 0)
      groupCombo->setCurrentIndex(i);
    else
      groupCombo->setEditText(label);
  }

  QComboBox *kindCombo;
  QLineEdit *thicknessText;
  QCheckBox *fitCheck = nullptr;
  QLineEdit *sigmaText = nullptr;
  QComboBox *groupCombo = nullptr;

 private:
  void updateEnabled() {
    const bool on = !kindCombo->currentData().toString().isEmpty();
    thicknessText->setEnabled(on);
    if (withFit_) {
      fitCheck->setEnabled(on);
      sigmaText->setEnabled(on);
      groupCombo->setEnabled(on);
    }
  }
  bool withFit_;
  SegmentKeywords::SampleBlock loaded_;
};

#endif
