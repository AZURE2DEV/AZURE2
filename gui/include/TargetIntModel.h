#ifndef TARGETINTMODEL_H
#define TARGETINTMODEL_H

#include <QAbstractTableModel>
#include <QList>

Q_DECLARE_METATYPE(QList<double>);

struct TargetIntData {
  static const int SIZE = 34;  // column 33: derived, display-only summary of the effect types
  int isActive;
  QString segmentsList;
  int numPoints;
  bool isConvolution;
  double sigma;
  bool isTargetIntegration;
  double density;
  QString stoppingPowerEq;
  int numParameters;
  QList<double> parameters;
  bool isQCoefficients;
  QList<double> qCoefficients;
  bool isConvCoefficients;
  QList<double> convCoefficients;
  QString convolutionEq;
  bool isStraggling;
  double stragglingCoefficient;
  double resonanceWidthMultiplier;
  double pointsPerWidth;
  // Optional restriction of the effect to lab-energy windows.
  QString applyRanges;        // "lo1-hi1,lo2-hi2"; empty = whole segment
  double transitionWidth;     // MeV; 0 = hard edges
  double autoTolerance;       // relative; 0 = always apply
  // Optional beam-profile kernel: an absolute (not point-centred) beam energy
  // profile of skewed Gaussians, a detector-resolution window per point, and
  // the detailed-balance weight of an inverse photodissociation measurement.
  bool isBeamProfile = false;
  QList<double> beamProfile;         // flattened (xi, omega, alpha, weight) quadruples
  double beamTpcSigma = 0.;          // detector energy resolution, MeV (lab)
  double beamTruncation = 0.;        // zero each component beyond this many s.d.; 0 = none
  bool beamPhotodissociation = false;  // weight the average with the detailed-balance factor
  // Optional user-defined numerical resolution function (SAMMY UDR): a file
  // of time-of-flight delay distributions, the nominal flight path, and the
  // burst / channel widths convolved on top.
  bool isUdr = false;
  QString udrFile;
  double udrFlightPath = 0.;         // metres
  double udrBurstFwhm = 0.;          // ns, Gaussian FWHM; 0 = none
  double udrChannelWidth = 0.;       // ns, rectangular; 0 = none
  bool udrCentred = true;            // re-centre the tabulated function on its centroid
};

/*!
 * Table model behind the Target Integration tab.
 */
class TargetIntModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  TargetIntModel(QObject *parent = 0);
  int rowCount(const QModelIndex &parent) const;
  int columnCount(const QModelIndex &parent) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());
  Qt::ItemFlags flags(const QModelIndex &index) const;
  QList<TargetIntData> getLines() const { return targetIntList; };

 private:
  QList<TargetIntData> targetIntList;
};

#endif
