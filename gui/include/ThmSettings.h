#ifndef THMSETTINGS_H
#define THMSETTINGS_H

#include <QList>
#include <QMap>
#include <QPair>
#include <QString>
#include <QStringList>

/*!
 * The options of the optional <thm> block (Config::ThmOptions), as the GUI
 * edits them.  Values are canonical: vertex=real is read as perlevel, flags as
 * bool.  Defaults are those of the engine, so a default-constructed object is
 * a project without the block.
 */
struct ThmSettings {
  QString entranceL = "incoherent";  // incoherent | coherent
  QString vertex = "constant";       // constant | perlevel | onshell
  QString kinematics = "lacognata";  // lacognata | triple | kf3body | lambda32
  bool coulombIntegral = false;
  double spectatorEnergy = 0.0;       // MeV, every THM entrance pair
  QMap<int, double> spectatorByPair;  // spectatorEnergy[<pair key>]
  QMap<int, QString> weight;          // weight[<segmentsData line>]
  QMap<int, QString> weightTest;      // weightTest[<segmentsTest line>]
  /// experiment[<name>] ... lines (THM experiments), verbatim, in block
  /// order.  Edited in the THM workspace's Experiments page
  /// (ThmExperimentRecord); a block that has one is never removed as
  /// "all default".
  QStringList experimentLines;

  bool operator==(const ThmSettings &o) const;
  bool operator!=(const ThmSettings &o) const { return !(*this == o); }
  bool isDefault() const { return *this == ThmSettings(); }
  /// The spectator energy of one entrance pair (Config::ThmOptions::SpectatorEnergy).
  double spectatorEnergyOf(int pairKey) const { return spectatorByPair.value(pairKey, spectatorEnergy); }

  /// The non-default keys, canonical key -> canonical value, in the order
  /// they are written (experiment lines are not keys).
  QList<QPair<QString, QString>> keyValues() const;

  /*! Parses the lines between <thm> and </thm> with the rules of
      Config::ReadThmBlock ('#' comments, key=value, same keys and values;
      experiment lines through the engine's own ParseThmExperimentLine and
      CheckThmExperiments).  Returns false and a message on anything the
      engine would refuse.  Weight files are not opened here (see validate). */
  static bool parse(const QStringList &lines, ThmSettings &out, QString *error = nullptr);

  /*! Parses one line; returns false on a line the engine refuses.  A blank or
      comment-only line gives an empty key, and so does an experiment line
      (stored raw in experimentLines, checked by parse). */
  static bool parseLine(const QString &line, QString &key, QString &value, ThmSettings &into);

  /*! The engine's verdict on a set of experiment lines (comments allowed):
      ParseThmExperimentLine on each, then CheckThmExperiments.  An empty
      string, or "<thm> experiment[<name>]: ..." as AZURE2 prints it. */
  static QString checkExperimentLines(const QStringList &lines);
  /// The <name> of an experiment[<name>] line, or "" if the line is not one.
  static QString experimentName(const QString &line);

  /*! What the engine will refuse at startup that can be seen without the
      data: every weight file must be readable, two columns, strictly
      increasing E and w > 0 (ThmWeightTable::Read, the engine's own reader).
      A relative path is taken from projectDir, as the engine takes it from
      the directory of the .azr.  Returns an empty string when all is well. */
  QString validate(const QString &projectDir) const;

  /*! The block lines for these settings, given the lines of the block that
      was read (possibly none).  Comment and blank lines are kept in place; a
      key line whose value did not change is kept verbatim (inline comment
      included), a changed one is rewritten keeping its inline comment, a key
      that went back to its default is dropped, and new keys are appended.
      Experiment lines: an old one still in experimentLines stays in place;
      the new lines of an experiment take the place of its first old line;
      lines of an experiment that is gone are dropped; the lines of new
      experiments are appended last. */
  QStringList compose(const QStringList &oldLines) const;
};

/*!
 * One THM experiment as the workspace edits it: the keys it shows as fields,
 * as written in the file, and every other key=value token kept verbatim
 * (keys a later engine knows and this page does not show).
 */
struct ThmExperimentRecord {
  QString name;
  QString originName;        ///< name when read ("" for a new record)
  QList<int> segments;       ///< <segmentsData> line numbers, ascending
  QString segmentsText;      ///< as written; kept while the set is the same
  QString background = "none";
  bool hasBackgroundKey = false;
  QString beam, target, spectator, beamEnergy;  ///< as written; "" = key absent
  /// lineshape=on (the Coulomb line shape of the spectator); written only
  /// when on, since off is the engine's default.
  bool lineshape = false;
  /// ps= (spectator-momentum window) and psNodes= as written; "" = key absent.
  /// The page composes them from its controls (psText) only when edited.
  QString ps, psNodes;
  QStringList extraTokens;   ///< other key=value tokens, verbatim, in order

  bool hasKinematics() const {
    return !beam.isEmpty() || !target.isEmpty() || !spectator.isEmpty() || !beamEnergy.isEmpty();
  }
  /// A ps window is set (ps= other than delta).
  bool hasWindow() const { return !ps.isEmpty() && ps != "delta"; }
  /// Same content (originName aside; segment sets compared, not their text).
  bool sameAs(const ThmExperimentRecord &o) const;
  /// The record as one experiment[...] line.
  QString line() const;

  /// 1,2,3,5 -> "1-3,5".
  static QString segmentsListText(const QList<int> &segments);
  /// "1,2,5-7" -> {1,2,5,6,7} with the engine's parser; false if refused.
  static bool expandSegments(const QString &text, QList<int> &out);

  /*! The records of a set of experiment lines (already checked by
      ThmSettings::parse), merged by name in order of first appearance. */
  static QList<ThmExperimentRecord> read(const QStringList &lines);
  /*! The experiment lines for `records`, given the lines they were read from
      and the records as read: a record that did not change keeps its lines
      verbatim and in place, a changed one is written as one line at the place
      of its first, a removed one is dropped, new ones are appended. */
  static QStringList compose(const QStringList &oldLines, const QList<ThmExperimentRecord> &oldRecords,
                             const QList<ThmExperimentRecord> &records);
};

#endif
