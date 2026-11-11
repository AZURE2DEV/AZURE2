#include "SegmentsDataModel.h"
#include "PairsModel.h"
#include "SegmentKeywords.h"
#include <QColor>

SegmentsDataModel::SegmentsDataModel(QObject *parent) :
  QAbstractTableModel(parent) {
}

int SegmentsDataModel::rowCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return segDataLineList.size();
}

int SegmentsDataModel::columnCount(const QModelIndex &parent) const {
  Q_UNUSED(parent);
  return SegmentsDataData::SIZE;
}

QVariant SegmentsDataModel::data(const QModelIndex &index, int role) const {
  if (!index.isValid()) return QVariant();

  if (index.row() >= segDataLineList.size() || index.row() < 0) return QVariant();

  if (role == Qt::DisplayRole) {
    SegmentsDataData line = segDataLineList.at(index.row());
    // A sample segment (keyword tail "sample <1|name>") shows its sample in front of its reaction.
    const QString natTok = SegmentKeywords::naturalSample(line.extraTokens);
    const QString nat = natTok.isEmpty() ? QString() : (natTok == "1" ? QString("<i>[default sample]</i> ") : QString("<i>[%1]</i> ").arg(natTok));
    if (index.column() == 1) {
      if (line.isAdvanced == 1) {
        // Handle advanced segments (sum/ratio)
        QString result;

        // For RATIO: main segment is numerator, component is denominator
        // For SUM: show all components being added

        if (line.operationType == 1) {  // RATIO
          // Show main segment as numerator
          if (pairsModel->getPairs().size() >= line.entrancePairIndex && pairsModel->getPairs().size() >= line.exitPairIndex) {
            PairsData firstPair = pairsModel->getPairs().at(line.entrancePairIndex - 1);
            PairsData secondPair = pairsModel->getPairs().at(line.exitPairIndex - 1);
            result = pairsModel->getReactionLabel(firstPair, secondPair);
          } else {
            result = "UNDEFINED";
          }

          // Then add the denominator component
          QStringList components = line.componentsList.split(";", Qt::SkipEmptyParts);

          if (components.isEmpty()) {
            result += "<br><font style='color:red;font-weight:bold;'>NO DENOMINATOR</font>";
          } else {
            // For ratio, only use the first component as denominator
            QString component = components[0];
            QStringList parts = component.split(", ");
            if (parts.size() >= 2) {
              int entranceKey = parts[0].split(": ")[1].toInt();
              int exitKey = parts[1].split(": ")[1].toInt();

              // Check if there's an angle specification
              QString angleInfo = "";
              if (parts.size() >= 3 && parts[2].contains("Angle: ")) {
                double angle = parts[2].split(": ")[1].toDouble();
                angleInfo = QString(" (%1°)").arg(angle, 0, 'f', 1);
              }

              if (pairsModel->getPairs().size() >= entranceKey && pairsModel->getPairs().size() >= exitKey) {
                PairsData entrancePair = pairsModel->getPairs().at(entranceKey - 1);
                PairsData exitPair = pairsModel->getPairs().at(exitKey - 1);
                QString reactionLabel = pairsModel->getReactionLabel(entrancePair, exitPair);

                result += " /<br>";
                result += reactionLabel + angleInfo;
              } else {
                result += " /<br>UNDEFINED";
              }
            }
          }
        } else {  // SUM
          // For sum, show the original segment first
          if (pairsModel->getPairs().size() >= line.entrancePairIndex && pairsModel->getPairs().size() >= line.exitPairIndex) {
            PairsData firstPair = pairsModel->getPairs().at(line.entrancePairIndex - 1);
            PairsData secondPair = pairsModel->getPairs().at(line.exitPairIndex - 1);
            result = pairsModel->getReactionLabel(firstPair, secondPair);
          } else {
            result = "UNDEFINED";
          }

          QStringList components = line.componentsList.split(";", Qt::SkipEmptyParts);

          if (components.isEmpty()) {
            result += "<br><font style='color:red;font-weight:bold;'>NO COMPONENTS</font>";
          } else {
            for (int i = 0; i < components.size(); i++) {
              QString component = components[i];
              QStringList parts = component.split(", ");
              if (parts.size() >= 2) {
                int entranceKey = parts[0].split(": ")[1].toInt();
                int exitKey = parts[1].split(": ")[1].toInt();

                // Pick up optional scaling factor so it shows in the table
                QString scalingInfo;
                for (int p = 2; p < parts.size(); p++) {
                  QString tok = parts[p].trimmed();
                  if (tok.startsWith("Scaling:")) {
                    double s = tok.section(":", 1).trimmed().toDouble();
                    scalingInfo = QString("%1 × ").arg(s, 0, 'g', 3);
                    break;
                  }
                }

                if (pairsModel->getPairs().size() >= entranceKey && pairsModel->getPairs().size() >= exitKey) {
                  PairsData entrancePair = pairsModel->getPairs().at(entranceKey - 1);
                  PairsData exitPair = pairsModel->getPairs().at(exitKey - 1);
                  QString reactionLabel = pairsModel->getReactionLabel(entrancePair, exitPair);

                  result += " +<br>";
                  result += scalingInfo + reactionLabel;
                } else {
                  result += " +<br>" + scalingInfo + "UNDEFINED";
                }
              }
            }
          }
        }
        return QString("<div align='left'>%1</div>").arg(nat + result);
      } else if (line.dataType == 3) {
        int i = 0;
        QList<PairsData> pairsList = pairsModel->getPairs();
        for (i = 0; i < pairsList.size(); i++)
          if (pairsList[i].pairType == 10) break;
        if (pairsList.size() >= line.entrancePairIndex && i < pairsList.size()) {
          PairsData firstPair = pairsModel->getPairs().at(line.entrancePairIndex - 1);
          return QString("<div align='left'>%1</div>").arg(nat + pairsModel->getReactionLabelTotalCapture(firstPair));
        } else
          return QString("<div align='left'><font style='color:red;font-weight:bold;'>UNDEFINED</font></div>");
      } else {
        if (pairsModel->getPairs().size() >= line.entrancePairIndex && pairsModel->getPairs().size() >= line.exitPairIndex) {
          PairsData firstPair = pairsModel->getPairs().at(line.entrancePairIndex - 1);
          PairsData secondPair = pairsModel->getPairs().at(line.exitPairIndex - 1);
          return QString("<div align='left'>%1</div>").arg(nat + pairsModel->getReactionLabel(firstPair, secondPair));
        } else
          return QString("<div align='left'><font style='color:red;font-weight:bold;'>UNDEFINED</font></div>");
      }
    } else if (index.column() == 2)
      return QVariant();
    else if (index.column() == 3) {
      if (line.lowEnergy == line.highEnergy)
        return line.lowEnergy;
      else
        return QString("%1-%2").arg(line.lowEnergy).arg(line.highEnergy);
    } else if (index.column() == 4)
      return QVariant();
    else if (index.column() == 5) {
      if (line.lowAngle == line.highAngle)
        return line.lowAngle;
      else
        return QString("%1-%2").arg(line.lowAngle).arg(line.highAngle);
    } else if (index.column() == 6)
      return QVariant();
    else if (index.column() == 7) {
      // a transmission / self-shielded-yield segment: the observable, its thickness and whether it is fitted
      const SegmentKeywords::SampleBlock sb = SegmentKeywords::sampleBlock(line.extraTokens);
      if (!sb.kind.isEmpty()) {
        const QString fit = !sb.fit ? QString() : sb.vary ? tr(", fitted") : tr(", fixed parameter");
        return QString("<center>%1<br><small>n = %2 atoms/b%3</small></center>")
            .arg(sb.kind == "transmission" ? tr("Transmission") : tr("Self-Shielded Yield"), sb.n, fit);
      }
      if (line.dataType == 2) {
        QChar orbital;
        switch (line.phaseL) {
          case 0:
            orbital = 's';
            break;
          case 1:
            orbital = 'p';
            break;
          case 2:
            orbital = 'd';
            break;
          case 3:
            orbital = 'f';
            break;
          case 4:
            orbital = 'g';
            break;
          case 5:
            orbital = 'h';
            break;
          case 6:
            orbital = 'i';
            break;
          default:
            orbital = '?';
        }
        QString tempSpin;
        if (((int)(line.phaseJ * 2)) % 2 != 0 && line.phaseJ != 0.)
          tempSpin = QString("%1/2").arg((int)(line.phaseJ * 2));
        else
          tempSpin = QString("%1").arg(line.phaseJ);
        return QString("<center>Phase Shift [%1<sub>%2</sub>]</center>").arg(orbital).arg(tempSpin);
      } else if (line.dataType == 1)
        return QString(tr("<center>Differential</center>"));
      else if (line.dataType == 3)
        return QString(tr("<center>Total Capture</center>"));
      else if (line.dataType == 4)
        return QString(tr("<center>C.M. Differential</center>"));
      else if (line.dataType == 5)
        return QString(tr("<center>Angle Integrated (E1)</center>"));
      else if (line.dataType == 6)
        return QString(tr("<center>Angle Integrated (E2)</center>"));
      else if (line.dataType == 7)
        return QString(tr("<center>Analyzing Power</center>"));
      else if (line.dataType == 9)
        return QString(tr("<center>Angular Coefficients</center>"));
      else if (line.dataType == 10)
        return QString(tr("<center>Spin: %1</center>")).arg(line.spinObservable);
      else
        return QString(tr("<center>Angle Integrated</center>"));
    } else if (index.column() == 8)
      return line.dataFile;
    else if (index.column() == 9) {
      // a normalization-group member's own norm is not used: grey (its leader's is, see the tooltip)
      if (groupLeader(index.row(), "normgroup") >= 0)
        return QString("<center><font style='color:gray;'>%1</font></center>").arg(line.dataNorm, 0, 'g', 2);
      if (line.varyNorm == 1)
        return QString("<center><font style='color:red;font-weight:bold;'>%1</font></center>").arg(line.dataNorm, 0, 'g', 2);
      else
        return QString("<center>%1</center>").arg(line.dataNorm, 0, 'g', 2);
    } else if (index.column() == 10)
      return QVariant();
    else if (index.column() == 11) {
      if (line.varyNorm == 1)
        return QString(tr("YES"));
      else
        return QString(tr("NO"));
    } else if (index.column() == 12) {
      return QVariant();
    } else if (index.column() == 13) {
      return QVariant();
    } else if (index.column() == 14) {
      if (groupLeader(index.row(), "shiftgroup") >= 0)
        return QString("<center><font style='color:gray;'>%1</font></center>").arg(line.energyShift, 0, 'g', 2);
      if (line.varyEnergyShift == 1)
        return QString("<center><font style='color:red;font-weight:bold;'>%1</font></center>").arg(line.energyShift, 0, 'g', 2);
      else
        return QString("<center>%1</center>").arg(line.energyShift, 0, 'g', 2);
    } else if (index.column() == 15) {
      return line.energyShiftError;
    } else if (index.column() == 16) {
      return line.varyEnergyShift;
    } else if (index.column() == GroupsColumn) {
      return groupsText(index.row());
    }
  } else if (role == Qt::ToolTipRole) {
    if (index.column() == GroupsColumn) {
      const QString tip = groupsToolTip(index.row());
      return tip.isEmpty() ? QVariant() : QVariant(tip);
    }
    const QString keyword = index.column() == 9 ? QString("normgroup") : index.column() == 14 ? QString("shiftgroup") : QString();
    const int leader = keyword.isEmpty() ? -1 : groupLeader(index.row(), keyword);
    if (leader >= 0)
      return tr("Not used: segment %1 takes the %2 of segment %3 (%4 %5).")
          .arg(index.row() + 1)
          .arg(keyword == "normgroup" ? tr("normalization") : tr("energy shift"))
          .arg(leader + 1)
          .arg(keyword)
          .arg(SegmentKeywords::groupLabel(segDataLineList.at(index.row()).extraTokens, keyword));
  } else if (role == Qt::EditRole) {
    SegmentsDataData line = segDataLineList.at(index.row());
    if (index.column() == 1)
      return line.entrancePairIndex;
    else if (index.column() == 2)
      return line.exitPairIndex;
    else if (index.column() == 3)
      return line.lowEnergy;
    else if (index.column() == 4)
      return line.highEnergy;
    else if (index.column() == 5)
      return line.lowAngle;
    else if (index.column() == 6)
      return line.highAngle;
    else if (index.column() == 7)
      return line.dataType;
    else if (index.column() == 8)
      return line.dataFile;
    else if (index.column() == 9)
      return line.dataNorm;
    else if (index.column() == 10)
      return line.dataNormError;
    else if (index.column() == 11)
      return line.varyNorm;
    else if (index.column() == 12)
      return line.phaseJ;
    else if (index.column() == 13)
      return line.phaseL;
    else if (index.column() == 14)
      return line.energyShift;
    else if (index.column() == 15)
      return line.energyShiftError;
    else if (index.column() == 16)
      return line.varyEnergyShift;
    else if (index.column() == 17)
      return line.isAdvanced;
    else if (index.column() == 18)
      return line.operationType;
    else if (index.column() == 19)
      return line.componentsList;
    else if (index.column() == 20)
      return line.isUPOS;
    else if (index.column() == 21)
      return line.secondaryDecayL;
    else if (index.column() == 22)
      return line.finalJ;
    else if (index.column() == 23)
      return line.delta;
  } else if (role == Qt::CheckStateRole && index.column() == 0) {
    SegmentsDataData line = segDataLineList.at(index.row());
    if (line.isActive == 1)
      return Qt::Checked;
    else
      return Qt::Unchecked;
  } else if (role == Qt::TextAlignmentRole)
    return Qt::AlignCenter;

  return QVariant();
}

QVariant SegmentsDataModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if (role != Qt::DisplayRole) return QVariant();
  if (orientation == Qt::Horizontal) {
    switch (section) {
      case 0:
        return tr("");
      case 1:
        return tr("Reaction");
      case 2:
        return QVariant();
      case 3:
        return tr("Energy\nRange");
      case 4:
        return QVariant();
      case 5:
        return tr("Angle\nRange");
      case 6:
        return QVariant();
      case 7:
        return tr("Data Type");
      case 8:
        return tr("Data File");
      case 9:
        return tr("Data\nNorm.");
      case 10:
        return tr("Data\nNorm. Err.");
      case 11:
        return tr("Vary\nNorm.?");
      case 12:
        return QVariant();
      case 13:
        return QVariant();
      case 14:
        return tr("Energy\nShift [MeV]");
      case 15:
        return tr("Energy Shift\nError [MeV]");
      case 16:
        return tr("Vary\nEnergy Shift?");
      case GroupsColumn:
        return tr("Groups");
      default:
        return QVariant();
    }
  } else if (orientation == Qt::Vertical) {
    return section + 1;
  }
  return QVariant();
}

bool SegmentsDataModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  if (index.isValid() && role == Qt::EditRole) {
    int row = index.row();
    SegmentsDataData tempData = segDataLineList.value(row);
    if (index.column() == 0)
      tempData.isActive = value.toInt();
    else if (index.column() == 1)
      tempData.entrancePairIndex = value.toInt();
    else if (index.column() == 2)
      tempData.exitPairIndex = value.toInt();
    else if (index.column() == 3)
      tempData.lowEnergy = value.toDouble();
    else if (index.column() == 4)
      tempData.highEnergy = value.toDouble();
    else if (index.column() == 5)
      tempData.lowAngle = value.toDouble();
    else if (index.column() == 6)
      tempData.highAngle = value.toDouble();
    else if (index.column() == 7)
      tempData.dataType = value.toInt();
    else if (index.column() == 8)
      tempData.dataFile = value.toString();
    else if (index.column() == 9) {
      tempData.dataNorm = value.toDouble();
      // Emit signal for normalization change
      emit normalizationChanged(row, value.toDouble());
    } else if (index.column() == 10) {
      tempData.dataNormError = value.toDouble();
      // Emit signal for normalization error change
      emit normalizationErrorChanged(row, value.toDouble());
    } else if (index.column() == 11) {
      tempData.varyNorm = value.toInt();
      // Emit signal for norm vary flag change
      emit normalizationVaryChanged(row, value.toInt() == 1);
    } else if (index.column() == 12)
      tempData.phaseJ = value.toDouble();
    else if (index.column() == 13)
      tempData.phaseL = value.toInt();
    else if (index.column() == 14) {
      tempData.energyShift = value.toDouble();
      // Emit signal for energy shift change
      emit energyShiftChanged(row, value.toDouble());
    } else if (index.column() == 15) {
      tempData.energyShiftError = value.toDouble();
      // Emit signal for energy shift error change
      emit energyShiftErrorChanged(row, value.toDouble());
    } else if (index.column() == 16) {
      tempData.varyEnergyShift = value.toInt();
      // Emit signal for energy shift vary flag change
      emit energyShiftVaryChanged(row, value.toInt() == 1);
    } else if (index.column() == 17)
      tempData.isAdvanced = value.toInt();
    else if (index.column() == 18)
      tempData.operationType = value.toInt();
    else if (index.column() == 19)
      tempData.componentsList = value.toString();
    else if (index.column() == 20)
      tempData.isUPOS = value.toInt();
    else if (index.column() == 21)
      tempData.secondaryDecayL = value.toInt();
    else if (index.column() == 22)
      tempData.finalJ = value.toDouble();
    else if (index.column() == 23)
      tempData.delta = value.toDouble();
    else
      return false;

    segDataLineList.replace(row, tempData);
    emit(dataChanged(index, index));
    if (index.column() == 0) emitGroupCellsChanged();  // an inactive segment leaves its groups
    return true;
  } else if (role == Qt::CheckStateRole) {
    int row = index.row();
    SegmentsDataData tempData = segDataLineList.value(row);
    if (index.column() == 0) {
      if (value == Qt::Checked)
        tempData.isActive = 1;
      else
        tempData.isActive = 0;
    } else
      return false;
    segDataLineList.replace(row, tempData);
    emit(dataChanged(index, index));
    emitGroupCellsChanged();
    return true;
  }
  return false;
}

void SegmentsDataModel::setLineSpin(int row, const QString &token) {
  if (row < 0 || row >= segDataLineList.size()) return;
  segDataLineList[row].spinObservable = token;
}

void SegmentsDataModel::emitGroupCellsChanged() {
  if (segDataLineList.isEmpty()) return;
  bool any = false;
  for (const SegmentsDataData &line : segDataLineList)
    for (const QString &keyword : SegmentKeywords::groupKeywords())
      if (!SegmentKeywords::groupLabel(line.extraTokens, keyword).isEmpty()) any = true;
  if (any) emit dataChanged(index(0, 0), index(segDataLineList.size() - 1, SegmentsDataData::SIZE - 1));
}

void SegmentsDataModel::setLineExtras(int row, const QString &tokens) {
  if (row < 0 || row >= segDataLineList.size()) return;
  segDataLineList[row].extraTokens = tokens;
  // The tail shows in the reaction (sample) and Groups columns of this row, and a group block changes which row
  // leads its group, so every row's norm, shift and Groups cells may change.
  emit dataChanged(index(0, 0), index(segDataLineList.size() - 1, SegmentsDataData::SIZE - 1));
}

int SegmentsDataModel::groupLeader(int row, const QString &keyword) const {
  if (row < 0 || row >= segDataLineList.size() || segDataLineList.at(row).isActive != 1) return -1;
  const QString label = SegmentKeywords::groupLabel(segDataLineList.at(row).extraTokens, keyword);
  if (label.isEmpty()) return -1;
  for (int i = 0; i < row; i++)
    if (segDataLineList.at(i).isActive == 1 && SegmentKeywords::groupLabel(segDataLineList.at(i).extraTokens, keyword) == label)
      return i;
  return -1;
}

QStringList SegmentsDataModel::groupLabels(const QString &keyword) const {
  QStringList labels;
  for (const SegmentsDataData &line : segDataLineList) {
    const QString label = SegmentKeywords::groupLabel(line.extraTokens, keyword);
    if (!label.isEmpty() && !labels.contains(label)) labels << label;
  }
  labels.sort();
  return labels;
}

static QString groupShortName(const QString &keyword) {
  if (keyword == "normgroup") return QObject::tr("norm");
  if (keyword == "shiftgroup") return QObject::tr("shift");
  if (keyword == "thicknessgroup") return QObject::tr("thickness");
  return QObject::tr("shift fn");
}

QString SegmentsDataModel::groupsText(int row) const {
  if (row < 0 || row >= segDataLineList.size()) return QString();
  QStringList parts;
  for (const QString &keyword : SegmentKeywords::groupKeywords()) {
    const QString label = SegmentKeywords::groupLabel(segDataLineList.at(row).extraTokens, keyword);
    if (!label.isEmpty()) parts << groupShortName(keyword) + ": " + label;
  }
  return parts.join(", ");
}

QString SegmentsDataModel::groupsToolTip(int row) const {
  if (row < 0 || row >= segDataLineList.size()) return QString();
  QStringList lines;
  for (const QString &keyword : SegmentKeywords::groupKeywords()) {
    const QString label = SegmentKeywords::groupLabel(segDataLineList.at(row).extraTokens, keyword);
    if (label.isEmpty()) continue;
    const QString what = keyword == "normgroup"        ? tr("normalization")
                         : keyword == "shiftgroup"     ? tr("constant energy shift")
                         : keyword == "thicknessgroup" ? tr("sample thickness")
                                                       : tr("shift function (each member keeps its own constant shift)");
    if (segDataLineList.at(row).isActive != 1) {
      lines << tr("%1 %2: inactive, not in the group").arg(keyword, label);
      continue;
    }
    const int leader = groupLeader(row, keyword);
    if (leader >= 0) {
      lines << tr("%1 %2: uses the %3 of segment %4").arg(keyword, label, what).arg(leader + 1);
      continue;
    }
    QStringList members;
    for (int i = row + 1; i < segDataLineList.size(); i++)
      if (groupLeader(i, keyword) == row) members << QString::number(i + 1);
    if (members.isEmpty())
      lines << tr("%1 %2: the only active segment of the group, nothing shared").arg(keyword, label);
    else
      lines << tr("%1 %2: first segment of the group; its %3 is used by segment%4 %5")
                   .arg(keyword, label, what, members.size() > 1 ? "s" : "", members.join(", "));
  }
  return lines.join("\n");
}

bool SegmentsDataModel::insertRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if (rows > 0) {
    beginInsertRows(QModelIndex(), position, position + rows - 1);
    for (int row = 0; row < rows; row++) {
      SegmentsDataData tempData;
      segDataLineList.insert(position, tempData);
    }
    endInsertRows();
  }
  return true;
}

bool SegmentsDataModel::removeRows(int position, int rows, const QModelIndex &index) {
  Q_UNUSED(index);
  if (rows > 0) {
    beginRemoveRows(QModelIndex(), position, position + rows - 1);
    for (int row = 0; row < rows; ++row) {
      segDataLineList.removeAt(position);
    }
    endRemoveRows();
  }
  return true;
}

Qt::ItemFlags SegmentsDataModel::flags(const QModelIndex &index) const {
  if (!index.isValid()) return Qt::ItemIsEnabled;
  if (index.column() == 0) return QAbstractTableModel::flags(index) | Qt::ItemIsUserCheckable;
  return QAbstractTableModel::flags(index);
}

// Row of a line equal to this one in every field (components, UPOS and keyword blocks included; not the active
// flag), or -1.  Before 2026-10-31 the components, UPOS fields and keyword blocks were not compared.
int SegmentsDataModel::isSegDataLine(const SegmentsDataData &line) const {
  int foundLine = -1;
  for (int i = 0; i < segDataLineList.size(); i++) {
    SegmentsDataData tempLine = segDataLineList.value(i);
    if (tempLine.entrancePairIndex == line.entrancePairIndex &&
        tempLine.exitPairIndex == line.exitPairIndex &&
        tempLine.lowEnergy == line.lowEnergy &&
        tempLine.highEnergy == line.highEnergy &&
        tempLine.lowAngle == line.lowAngle &&
        tempLine.highAngle == line.highAngle &&
        tempLine.dataType == line.dataType &&
        tempLine.dataFile == line.dataFile &&
        tempLine.dataNorm == line.dataNorm &&
        tempLine.dataNormError == line.dataNormError &&
        tempLine.varyNorm == line.varyNorm &&
        tempLine.phaseJ == line.phaseJ &&
        tempLine.phaseL == line.phaseL &&
        tempLine.spinObservable == line.spinObservable &&
        tempLine.energyShift == line.energyShift &&
        tempLine.energyShiftError == line.energyShiftError &&
        tempLine.varyEnergyShift == line.varyEnergyShift &&
        tempLine.isAdvanced == line.isAdvanced &&
        tempLine.operationType == line.operationType &&
        tempLine.componentsList == line.componentsList &&
        tempLine.isUPOS == line.isUPOS &&
        tempLine.secondaryDecayL == line.secondaryDecayL &&
        tempLine.finalJ == line.finalJ &&
        tempLine.delta == line.delta &&
        tempLine.extraTokens.simplified() == line.extraTokens.simplified()) {
      foundLine = i;
      break;
    }
  }
  return foundLine;
}

void SegmentsDataModel::setPairsModel(PairsModel *model) {
  pairsModel = model;
}

QString SegmentsDataModel::getReactionLabel(const QModelIndex &index) {
  SegmentsDataData line = segDataLineList.at(index.row());
  if (line.dataType == 3) {
    int i = 0;
    QList<PairsData> pairsList = pairsModel->getPairs();
    for (i = 0; i < pairsList.size(); i++)
      if (pairsList[i].pairType == 10) break;
    if (pairsList.size() >= line.entrancePairIndex && i < pairsList.size()) {
      PairsData firstPair = pairsModel->getPairs().at(line.entrancePairIndex - 1);
      return pairsModel->getReactionLabelTotalCapture(firstPair);
    }
    return QString("<font style='color:red;font-weight:bold;'>UNDEFINED</font>");
  } else {
    int numPairs = pairsModel->getPairs().size();
    if (line.entrancePairIndex - 1 >= numPairs ||
        line.exitPairIndex - 1 >= numPairs) return QString("<font style='color:red;font-weight:bold'>UNDEFINED</font>");
    PairsData firstPair = pairsModel->getPairs().at(line.entrancePairIndex - 1);
    PairsData secondPair = pairsModel->getPairs().at(line.exitPairIndex - 1);
    return pairsModel->getReactionLabel(firstPair, secondPair);
  }
}


/*!
 * Moves the line at \p from so that it ends up at \p to, carrying the whole
 * underlying struct.  The up/down buttons used to rebuild the row column by
 * column and silently dropped every field the copy list had fallen behind on
 * (the UPOS block, most recently); moving the struct itself cannot lose
 * anything.
 */

bool SegmentsDataModel::moveLine(int from, int to) {
  if (from < 0 || to < 0 || from >= segDataLineList.size() || to >= segDataLineList.size() || from == to) return false;
  int destination = (to > from) ? to + 1 : to;
  if (!beginMoveRows(QModelIndex(), from, from, QModelIndex(), destination)) return false;
  segDataLineList.move(from, to);
  endMoveRows();
  return true;
}
