#ifndef TARGETINTTAB_H
#define TARGETINTTAB_H

#include <QWidget>
#include <QItemSelection>
#include <QTableView>
#include <QSignalMapper>
#include <QPointer>
#include "TargetIntModel.h"
#include "SegmentsDataModel.h"
#include "PairsModel.h"
#include "AddTargetIntDialog.h"

QT_BEGIN_NAMESPACE

class QPushButton;
class QLineEdit;

QT_END_NAMESPACE

class InfoDialog;

/*!
 * The Target Integration tab: beam convolution, target thickness and attenuation coefficients, written to the target-effects file.
 */
class TargetIntTab : public QWidget {
  Q_OBJECT

 public:
  TargetIntTab(QWidget *parent = 0);
  /*! The segments and pairs, so an effect's entrance pair (and its c.m./lab energy ratio) can be looked up. */
  void setModels(SegmentsDataModel *segments, PairsModel *pairs);
  /*! c.m./lab energy ratio m_t/(m_p+m_t) of the entrance pair shared by every segment in a key list
   *  ("1,3-5"); 0 when the pairs differ or a key is unknown. */
  double conversionFactorFor(const QString &segmentsList) const;
  /// The adiabatic-limit screening potential (eV) of the segments' common entrance pair, or -1 when the
  /// list is empty, names a missing segment, mixes pairs, or the pair is not charged on charged.
  double adiabaticScreeningFor(const QString &segmentsList) const;
  TargetIntModel *getTargetIntModel();
  bool writeFile(QTextStream &outStream);
  bool readFile(QTextStream &inStream);
  void reset();

 public slots:
  void addLine();
  void addLine(TargetIntData line);
  void editLine();
  void deleteLine();
  void updateButtons(const QItemSelection &selection);
  void showInfo(int which = 0, QString title = "");

 private:
  TargetIntModel *targetIntModel;
  SegmentsDataModel *segmentsModel_ = nullptr;
  PairsModel *pairsModel_ = nullptr;
  QTableView *targetIntView;
  QPushButton *addButton;
  QPushButton *deleteButton;
  QSignalMapper *mapper;
  QPushButton *infoButton[5];
  static const std::vector<QString> infoText;
  QPointer<InfoDialog> infoDialog[5];
};

#endif
