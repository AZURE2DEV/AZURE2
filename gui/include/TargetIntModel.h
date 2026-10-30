#ifndef TARGETINTMODEL_H
#define TARGETINTMODEL_H

#include <QAbstractTableModel>
#include <QList>

Q_DECLARE_METATYPE(QList<double>);

struct TargetIntData {
  static const int SIZE = 39;  // column 33: derived, display-only summary of the effect types; 34: sigmaLab;
                               // 35-38: electron screening (on, U_e eV, vary, prior sigma eV)
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
  // A sigma(E) expression (isConvCoefficients) is lab in, lab out when this is set (the `sigmalab`
  // keyword), as the constant sigma is; false is the pre-2026-10 convention, c.m. in, c.m. out.
  bool sigmaLab = false;
  // Electron screening (the `screening Ue vary sigma` block): the bare cross section of the listed
  // segments times exp(pi eta U_e / E); U_e a fit parameter when screeningVary, with a Gaussian prior of
  // width screeningSigma (eV, 0 = none) about screeningUe.
  bool isScreening = false;
  double screeningUe = 0.;     // eV
  bool screeningVary = false;
  double screeningSigma = 0.;  // eV
  // A fitted target density (`densityfit vary sigma [max]`, `densitygroup <label>`), kept as the line's text
  // and written back unchanged: the tab does not edit it yet (the engine and pyazr do).
  QString densityBlock;
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
  /// The densityfit / densitygroup text of a row (no table column: it is carried, not edited, by the tab).
  void setDensityBlock(int row, const QString &block) {
    if (row >= 0 && row < targetIntList.size()) targetIntList[row].densityBlock = block;
  }

 private:
  QList<TargetIntData> targetIntList;
};

#endif
