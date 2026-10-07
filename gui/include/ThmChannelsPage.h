#ifndef THMCHANNELSPAGE_H
#define THMCHANNELSPAGE_H

#include <QList>
#include <QSet>
#include <QWidget>
#include <functional>

QT_BEGIN_NAMESPACE
class QAction;
class QTableWidget;
QT_END_NAMESPACE

class ChannelsModel;
class LevelsModel;
class PairsModel;
class SegmentsDataModel;
class SegmentsTestModel;

/*!
 * The Channels page of the THM workspace: the THM columns of the <levels>
 * lines.  For each particle pair that is the entrance of a THM segment, the
 * binding energy B of the transferred particle in the Trojan horse (field 32,
 * written on every channel line of the pair); for each particle channel of
 * those pairs, the width input flag (field 33: the width column is a reduced
 * width amplitude in MeV^(1/2) rather than a partial width or ANC).  A pair
 * with B != 0 or a channel with the flag set is listed even when no THM
 * segment uses its pair, so that nothing in the file is out of reach.  A
 * warning icon beside B marks a pair whose B differs from B(x+s) of an
 * experiment's reaction (the masses), as the engine warns at startup.
 */
class ThmChannelsPage : public QWidget {
  Q_OBJECT

 public:
  ThmChannelsPage(PairsModel *pairs, LevelsModel *levels, ChannelsModel *channels, SegmentsDataModel *segmentsData,
                  SegmentsTestModel *segmentsTest, QWidget *parent = 0);

  /// Pair keys (1-based) that are the entrance of a THM data or test segment.
  static QSet<int> thmEntrancePairs(SegmentsDataModel *segmentsData, SegmentsTestModel *segmentsTest);

  /// A binding energy that does not read as a number; "" if none.
  QString check() const;
  /// Writes the edited values into the models; values left as shown are not touched.
  void apply();

  /// Row of a pair key / of a channel (index in ChannelsModel) in the tables, or -1.
  int pairRow(int pairKey) const { return pairRows_.indexOf(pairKey); }
  int channelRow(int channelIndex) const { return channelRows_.indexOf(channelIndex); }
  /// For the tests: set a binding energy / a flag as the widgets do.
  void setBindingText(int pairKey, const QString &text);
  void setReducedWidthFlag(int channelIndex, bool rwa);
  /// B of a pair as the page shows it (the edited value if it reads as a
  /// number, else the project's).
  double bindingOf(int pairKey) const;
  /// The warning for a pair and B ("" none): the Experiments page's
  /// B(x+s) of the reactions on that pair against B.
  void setBindingWarning(std::function<QString(int pairKey, double pairB)> warning);
  /// Recomputes the warning icons (B or the experiments changed).
  void refreshWarnings();
  /// The warning shown beside a pair's B ("" if none), as its tooltip says it.
  QString bindingWarning(int pairKey) const;
  /// The page's fields as text (the B fields and the flags): changes when an edit does.
  QString editState() const;

  QTableWidget *pairTable;     ///< Pair | Nuclei | THM segments | B (MeV) (a line edit)
  QTableWidget *channelTable;  ///< Level | Pair | l | s | Width (as entered) | Amplitude (check: field 33)

 private:
  PairsModel *pairs_;
  ChannelsModel *channels_;
  QList<int> pairRows_;     // pair key per row
  QList<int> channelRows_;  // ChannelsModel row per row
  QList<QAction *> warningActions_;  // per pair row: the icon inside the B field
  std::function<QString(int, double)> warning_;
};

#endif
