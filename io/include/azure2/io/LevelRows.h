#ifndef AZURE2_IO_LEVELROWS_H
#define AZURE2_IO_LEVELROWS_H

#include <istream>
#include <string>
#include <vector>

namespace azure2 {
namespace io {

/*!
 * One line of the <levels> section: one channel of one level, with the pair it belongs to.  The grammar: 31 numbers,
 *
 *   levelJ levelPi levelE levelFix aa ir 2s 2l levelID isActive channelFix gamma j1 pi1 j2 pi2 e2 m1 m2 z1 z2
 *   entranceSepE sepE j3 pi3 e3 pType chRad g1 g2 ecMultMask
 *
 * then up to three optional ones, told apart by their number: 1 = the light particle's excitation energy e1 (MeV);
 * 2 = the THM binding energy (MeV) and the width input flag (0 = physical width or ANC, 1 = reduced-width amplitude
 * in MeV^1/2); 3 = the binding energy, the flag and e1.  The file writes the channel spin and l doubled; the row holds
 * s = (2s)/2 and l = (2l)/2 (integer division).  Energies in MeV (E_x of the compound nucleus for levelE), masses in
 * u, the channel radius in fm, gamma in the file's unit (eV for a width, fm^-1/2 for an ANC, MeV^1/2 with the flag).
 * pyazr's writer has the field list in Python (pyazr/azrfile.py).
 */
class LevelRow {
 public:
  /// Reads one row from `stream`.  A missing or non-numeric field sets failbit, as do more than three optional fields
  /// or an optional field that is not a number; a 31-field line leaves the stream at eof without failbit.
  explicit LevelRow(std::istream &stream);

  double levelJ() const { return levelJ_; }
  int levelPi() const { return levelPi_; }
  double levelE() const { return levelE_; }
  int levelFix() const { return levelFix_; }
  int aa() const { return aa_; }
  int ir() const { return ir_; }  ///< the pair key
  double s() const { return s_; }
  int l() const { return l_; }
  int levelID() const { return levelID_; }
  int isActive() const { return isActive_; }
  int channelFix() const { return channelFix_; }
  double gamma() const { return gamma_; }
  double j1() const { return j1_; }
  int pi1() const { return pi1_; }
  double j2() const { return j2_; }
  int pi2() const { return pi2_; }
  double e2() const { return e2_; }  ///< the heavy particle's excitation energy (MeV)
  double m1() const { return m1_; }
  double m2() const { return m2_; }
  int z1() const { return z1_; }
  int z2() const { return z2_; }
  double entranceSepE() const { return entranceSepE_; }  ///< not used
  double sepE() const { return sepE_; }
  int j3() const { return j3_; }      ///< not used
  int pi3() const { return pi3_; }    ///< not used
  double e3() const { return e3_; }   ///< not used
  int pType() const { return pType_; }  ///< 0 particle pair, 10 photon, 20 beta
  double chRad() const { return chRad_; }
  double g1() const { return g1_; }
  double g2() const { return g2_; }
  double ecMultMask() const { return ecMultMask_; }  ///< external-capture multipole mask of the pair
  double e1() const { return e1_; }                  ///< the light particle's excitation energy (MeV), 0 when absent
  double bindingE() const { return bindingE_; }      ///< THM binding energy (MeV), 0 when absent
  int gammaIsRWA() const { return gammaIsRWA_; }     ///< 1: gamma is a reduced-width amplitude (MeV^1/2)

 private:
  double levelJ_;
  int levelPi_;
  double levelE_;
  int levelFix_;
  int aa_;
  int ir_;
  double s_;
  int l_;
  int levelID_;
  int isActive_;
  int channelFix_;
  double gamma_;
  double j1_;
  int pi1_;
  double j2_;
  int pi2_;
  double e2_;
  double m1_;
  double m2_;
  int z1_;
  int z2_;
  double entranceSepE_;
  double sepE_;
  int j3_;
  int pi3_;
  double e3_;
  int pType_;
  double chRad_;
  double g1_;
  double g2_;
  unsigned int ecMultMask_;
  double bindingE_;
  int gammaIsRWA_;
  double e1_;
};

/// The rows of a <levels> section, in file order.
struct LevelSection {
  std::vector<LevelRow> rows;  ///< every row up to the first malformed line
  bool malformed = false;      ///< a line did not read as a row (badLine); the rows after it are not read
  std::string badLine;
  bool complete = false;       ///< the section ended with its </levels> line
};

/*!
 * Reads the <levels> section as Config::SectionStream gives it: the first line is the <levels> tag, then one row per
 * non-blank line (blanks and tabs only count as blank) up to </levels>.  A last line without its newline is not read.
 */
LevelSection ReadLevelSection(std::istream &in);

}  // namespace io
}  // namespace azure2

#endif
