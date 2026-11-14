// The four data check files (lines 6-9 of <config>): the points (data.chk), their energy-dependent values
// (lomatrixandpene.chk), Legendre polynomials (legendre.chk) and Coulomb amplitudes (coulombamplitudes.chk), each
// to checks/ or to the screen.
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "CompoundSystem.h"
#include "Config.h"
#include "WorkingCopy.h"

/*!
 * Prints the data point after the object is filled or points are created.
 */

namespace {
// (position, point) of every point of the segments without components, 0-based, in data order: the rows of the check
// files.
std::vector<std::pair<int, int>> PlainPointRows(const DataSet &data) {
  std::vector<std::pair<int, int>> rows;
  for (int s = 0; s < data.NumSegments(); s++)
    if (!data.GetSegments()[s].HasComponents())
      for (int p = 0; p < data.GetSegments()[s].NumPoints(); p++) rows.push_back({s, p});
  return rows;
}
}  // namespace

void WorkingCopy::PrintData(const Config &configure) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_DATA) {
    std::string outfile = configure.checkdir + "data.chk";
    fbuffer.open(configure.OpenPath(outfile).c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_DATA)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_DATA) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_DATA)) {
    out << std::endl
        << "************************************" << std::endl
        << "*            Segments              *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(11) << "Segment #"
        << std::setw(17) << "Segment Key #"
        << std::setw(17) << "Entrance Key #"
        << std::setw(13) << "Exit Key #"
        << std::setw(12) << "Min Energy"
        << std::setw(12) << "Max Energy"
        << std::setw(11) << "Min Angle"
        << std::setw(11) << "Max Angle"
        << std::setw(25) << "Data File"
        << std::endl;
    for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
      out << std::setw(11) << segment - Data().GetSegments().begin() + 1
          << std::setw(17) << segment->GetSegmentKey()
          << std::setw(17) << segment->GetEntranceKey()
          << std::setw(13) << segment->GetExitKey()
          << std::setw(12) << segment->GetMinEnergy()
          << std::setw(12) << segment->GetMaxEnergy()
          << std::setw(11) << segment->GetMinAngle()
          << std::setw(11) << segment->GetMaxAngle()
          << std::setw(25) << segment->GetDataFile()
          << std::endl;
    }
    out << std::endl
        << "************************************" << std::endl
        << "*               Data               *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(11) << "Segment #"
        << std::setw(14) << "Data Point #"
        << std::setw(15) << "Lab Energy"
        << std::setw(15) << "CM Energy"
        << std::setw(15) << "Angle"
        << std::setw(20) << "Cross Section"
        << std::setw(22) << "Cross Section Error"
        << std::setw(12) << "Map Point"
        << std::setw(18) << "# of Subpoints"
        << std::setw(18) << "Low Sub Energy"
        << std::setw(18) << "High Sub Energy"
        << std::endl;
    for (const std::pair<int, int> &row : PlainPointRows(Data())) {
      const int s = row.first, p = row.second;
      const Segment &seg = Data().GetSegments()[s];
      Unit *pt = seg.GetPoint(*this, p + 1);
      out << std::setw(11) << s + 1
          << std::setw(14) << p + 1
          << std::setw(15) << pt->GetLabEnergy()
          << std::setw(15) << pt->GetCMEnergy()
          << std::setw(15) << pt->GetCMAngle()
          << std::setw(20) << seg.CMValue(State(seg), p)
          << std::setw(22) << seg.CMError(State(seg), p);
      // The mapping column of the file's layout: points are no longer mapped onto one another.
      out << std::setw(12) << "Not Mapped"
          << std::setw(18) << pt->NumSubPoints();
      if (pt->IsTargetEffect() && pt->NumSubPoints() > 0 &&
          Data().GetTargetEffect(pt->GetTargetEffectNum())->IsSubPointEffect()) {
        out << std::setw(18) << pt->GetSubPoint(*this, pt->NumSubPoints())->GetCMEnergy()
            << std::setw(18) << pt->GetSubPoint(*this, 1)->GetCMEnergy();
      }
      out << std::endl;
      if (p == seg.NumPoints() - 1) out << std::endl;
    }
  } else
    configure.outStream << "Could not write data check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

/*!
 * Prints the Legendre polynomials for each point in the WorkingCopy object.
 */

void WorkingCopy::PrintLegendreP(const Config &configure) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_LEGENDRE) {
    std::string outfile = configure.checkdir + "legendre.chk";
    fbuffer.open(configure.OpenPath(outfile).c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_LEGENDRE)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_LEGENDRE) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_LEGENDRE)) {
    out << std::endl
        << "************************************" << std::endl
        << "*       Legendre Polynomials       *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(10) << "Segment #"
        << std::setw(10) << "Point #"
        << std::setw(15) << "CM Energy"
        << std::setw(15) << "Angle"
        << std::setw(5) << "L"
        << std::setw(15) << "Leg. Poly." << std::endl;
    for (const std::pair<int, int> &row : PlainPointRows(Data())) {
      const int s = row.first, p = row.second;
      const Segment &seg = Data().GetSegments()[s];
      Unit *pt = seg.GetPoint(*this, p + 1);
      for (int lOrder = 0; lOrder <= pt->GetMaxLOrder(); lOrder++) {
        out << std::setw(10) << s + 1
            << std::setw(10) << p + 1
            << std::setw(15) << pt->GetCMEnergy()
            << std::setw(15) << pt->GetCMAngle()
            << std::setw(5) << lOrder
            << std::setw(15) << pt->GetLegendreP(lOrder) << std::endl;
      }
    }
  } else
    configure.outStream << "Could not write legendre polynomials check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

/*!
 * Prints the values calculated by Unit::CalcEDependentValues for each point in the entire WorkingCopy object.
 */

void WorkingCopy::PrintEDependentValues(const Config &configure, const CompoundSystem *compound) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_ENERGY_DEP) {
    std::string outfile = configure.checkdir + "lomatrixandpene.chk";
    fbuffer.open(configure.OpenPath(outfile).c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_ENERGY_DEP)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_ENERGY_DEP) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_ENERGY_DEP)) {
    out << std::endl
        << "************************************" << std::endl
        << "*  Lo Matrix and Penetrabilities   *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(10) << "Seg_#"
        << std::setw(10) << "Point_#"
        << std::setw(5) << "j"
        << std::setw(5) << "ch"
        << std::setw(5) << "l"
        << std::setw(15) << "E_chan"
        << std::setw(15) << "sqrt_pene"
        << std::setw(25) << "Lo"
        << std::setw(25) << "expHSP" << std::endl;
    for (const std::pair<int, int> &row : PlainPointRows(Data())) {
      const int s = row.first, p = row.second;
      const Segment &seg = Data().GetSegments()[s];
      Unit *pt = seg.GetPoint(*this, p + 1);
      const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(seg.GetEntranceKey()));
      double inEnergy = pt->GetCMEnergy() + entrancePair->GetSepE();
      for (int j = 1; j <= compound->NumJGroups(); j++) {
        if (compound->GetJGroup(j)->IsInRMatrix() && compound->GetJGroup(j)->GetNucleus() == entrancePair->GetNucleus()) {
          const JGroup *theJGroup = compound->GetJGroup(j);
          for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
            const AChannel *theChannel = theJGroup->GetChannel(ch);
            const PPair *thePair = compound->GetPair(theChannel->GetPairNum());
            int lValue = theChannel->GetL();
            double localEnergy = inEnergy - thePair->GetSepE() - thePair->GetExE();
            complex loElement = pt->GetLoElement(j, ch);
            complex expHSP = pt->GetExpHardSpherePhase(j, ch);
            out << std::setw(10) << s + 1
                << std::setw(10) << p + 1
                << std::setw(5) << j
                << std::setw(5) << ch
                << std::setw(5) << lValue
                << std::setw(15) << localEnergy
                << std::setw(15) << pt->GetSqrtPenetrability(j, ch)
                << std::setw(25) << loElement.real() << " " << loElement.imag()
                << std::setw(25) << expHSP.real() << " " << expHSP.imag() << std::endl;
          }
        }
      }
    }
  } else
    configure.outStream << "Could not write lo-matrix and penetrabilities check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

/*!
 * Prints the values calculated by Unit::CalcCoulombAmplitude for each point in the entire WorkingCopy object.
 */

void WorkingCopy::PrintCoulombAmplitude(const Config &configure, const CompoundSystem *compound) {
  std::streambuf *sbuffer;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & Config::CHECK_COUL_AMPLITUDES) {
    std::string outfile = configure.checkdir + "coulombamplitudes.chk";
    fbuffer.open(configure.OpenPath(outfile).c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & Config::CHECK_COUL_AMPLITUDES)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & Config::CHECK_COUL_AMPLITUDES) && fbuffer.is_open()) ||
      (configure.screenCheckMask & Config::CHECK_COUL_AMPLITUDES)) {
    out << std::endl
        << "************************************" << std::endl
        << "*        Coulomb Amplitudes        *" << std::endl
        << "************************************" << std::endl;
    out << std::setw(10) << "segment #"
        << std::setw(10) << "point #"
        << std::setw(10) << "aa"
        << std::setw(15) << "cmenergy"
        << std::setw(15) << "angle"
        << std::setw(25) << "coulomb amplitude"
        << std::endl;
    for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
      if (segment->GetEntranceKey() == segment->GetExitKey()) {
        for (UnitIterator point = segment->GetPoints(*this).begin(); point < segment->GetPoints(*this).end(); point++) {
          out << std::setw(10) << segment - Data().GetSegments().begin() + 1
              << std::setw(10) << point - segment->GetPoints(*this).begin() + 1
              << std::setw(10) << compound->GetPairNumFromKey(segment->GetEntranceKey())
              << std::setw(15) << point->GetCMEnergy()
              << std::setw(15) << point->GetCMAngle()
              << std::setw(25) << point->GetCoulombAmplitude()
              << std::endl;
        }
      }
    }
  } else
    configure.outStream << "Could not write coulomb amplitudes check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}
