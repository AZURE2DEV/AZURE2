#ifndef NUCLINE_H
#define NUCLINE_H

#include <iostream>
#include <istream>

/// A class to read and store a line from a nuclear input file.

/*!
 * The NucLine class reads and stores a line from a formatted
 * nuclear input file.
 */

class NucLine {
 public:
  /*!
   * Constructor fills the NucLine object from an input stream.
   */
  NucLine(std::istream &stream) {
    stream >> levelJ_ >> levelPi_ >> levelE_ >> levelFix_ >> aa_ >> ir_ >> s_ >> l_ >> levelID_ >> isActive_ >> channelFix_ >> gamma_ >> j1_ >> pi1_ >> j2_ >> pi2_ >> e2_ >> m1_ >> m2_ >> z1_ >> z2_ >> entranceSepE_ >> sepE_ >> j3_ >> pi3_ >> e3_ >> pType_ >> chRad_ >> g1_ >> g2_ >> ecMultMask_;

    // Optional trailing columns.  dev's 32nd column (2026-11-05) is the light particle's excitation energy e1
    // (MeV); thm's two (the Trojan Horse binding energy in MeV and the width input flag, 0 = physical width/ANC,
    // 1 = reduced width in MeV^(1/2), particle channels only) came first and are columns 32 and 33 of every THM
    // file.  The number of trailing fields tells them apart: 1 = e1 (dev), 2 = binding energy and flag (thm),
    // 3 = binding energy, flag and e1.  A 31-field line leaves the stream's state as it was before either
    // existed (eof, not fail).  A trailing token that is not a number sets failbit (as dev), more than three
    // trailing fields too.
    e1_ = 0.;
    bindingE_ = 0.0;
    gammaIsRWA_ = 0;
    if (!stream.fail() && !stream.eof()) {
      double v[3];
      int n = 0;
      stream >> std::ws;
      while (!stream.eof() && n < 4) {
        double x;
        if (!(stream >> x)) break;  // failbit set: reported by the caller
        if (n < 3) v[n] = x;
        n++;
        if (!stream.eof()) stream >> std::ws;  // not at eof: a newer libstdc++'s ws sets failbit there
      }
      if (!stream.fail()) {
        if (n == 1) {
          e1_ = v[0];
        } else if (n == 2 || n == 3) {
          bindingE_ = v[0];
          gammaIsRWA_ = (int)v[1];
          if (n == 3) e1_ = v[2];
        } else if (n > 3) {
          stream.setstate(std::ios::failbit);
        }
      }
    }
    s_ /= 2.;
    l_ /= 2;
  };
  /*!
   * Returns the spin of the corresponding level for the channel line.
   */
  double levelJ() const { return levelJ_; };
  /*!
   * Returns the parity of the corresponding level for the channel line.
   */
  int levelPi() const { return levelPi_; };
  /*!
   * Returns the excitation energy of the corresponding level for
   * the channel line.
   */
  double levelE() const { return levelE_; };
  /*!
   * Returns non-zero if the corresponding level for the channel
   * line is not to be varied in the fit.
   */
  int levelFix() const { return levelFix_; };
  /*!
   * Returns the entrance key of the corresponding particle pair
   * for the channel line.
   */
  int aa() const { return aa_; };
  /*!
   * Returns the exit key of the corresponding particle pair for
   * the channel line.
   */
  int ir() const { return ir_; };
  /*!
   * Returns the channel spin of the channel.
   */
  double s() const { return s_; };
  /*!
   * Returns the orbital angular momentum of the channel.
   */
  int l() const { return l_; };
  /*!
   * Returns an indexing variable used by the graphical setup program.
   */
  int levelID() const { return levelID_; };
  /*!
   * Returns non-zero if the corresponding level for the channel line
   * is to be used in the calculation.
   */
  int isActive() const { return isActive_; };
  /*!
   * Returns non-zero if the reduced width amplitude for the channel is not
   * to be varied.
   */
  int channelFix() const { return channelFix_; };
  /*!
   * Returns the initial reduced width or physical amplitude for the channel.
   */
  double gamma() const { return gamma_; };
  /*!
   * Returns the spin of the light particle in the corresponding pair.
   */
  double j1() const { return j1_; };
  /*!
   * Returns the parity of the light particle in the corresponding pair.
   */
  int pi1() const { return pi1_; };
  /*!
   * Returns the spin of the heavy particle in the corresponding pair.
   */
  double j2() const { return j2_; };
  /*!
   * Returns the parity of the heavy particle in the corresponding pair.
   */
  int pi2() const { return pi2_; };
  /*!
   * Returns the excitation energy of the heavy particle in the
   * corresponding pair.
   */
  double e2() const { return e2_; };
  /*!
   * Returns the mass of the light particle in the corresponding pair.
   */
  double m1() const { return m1_; };
  /*!
   * Returns the mass of the heavy particle in the corresponding pair.
   */
  double m2() const { return m2_; };
  /*!
   * Returns the charge of the light particle in the corresponding pair.
   */
  int z1() const { return z1_; };
  /*!
   * Returns the charge of the heavy particle in the corresponding pair.
   */
  int z2() const { return z2_; };
  /*!
   * This function is depriciated and not used.
   */
  double entranceSepE() const { return entranceSepE_; };
  /*!
   * Returns the separation energy for the corresponding particle pair.
   */
  double sepE() const { return sepE_; };
  /*!
   * This function is depriciated and not used.
   */
  int j3() const { return j3_; };
  /*!
   * This function is depriciated and not used.
   */
  int pi3() const { return pi3_; };
  /*!
   * This function is depriciated and not used.
   */
  double e3() const { return e3_; };
  /*!
   * Returns 0 for particle-particle and 10 for particle-gamma
   * types in the corresponding pair.
   */
  int pType() const { return pType_; };
  /*!
   * Returns the channel radius for the corresponding particle pair.
   */
  double chRad() const { return chRad_; };
  /*!
   * Returns the g-factor for the light particle in the corresponding pair.
   */
  double g1() const { return g1_; };
  /*!
   * Returns the g-factor for the heavy particle in the corresponding pair.
   */
  double g2() const { return g2_; };
  /*!
   * Returns the external capture multiplicity mask for the corresponding pair.
   */
  double ecMultMask() const { return ecMultMask_; };
  /*!
   * Returns the excitation energy (MeV) of the light particle in the
   * corresponding pair: the last optional column (32nd, or 34th after the THM columns), 0 when absent.
   */
  double e1() const { return e1_; };
  /*!
   * Returns the THM binding energy (MeV) of the transferred particle in the
   * Trojan-Horse nucleus for the corresponding pair (0 if not specified).
   */
  double bindingE() const { return bindingE_; };
  /*!
   * Returns 1 if the gamma column of this channel line is already a reduced
   * width amplitude (MeV^(1/2)) rather than a physical partial width/ANC.
   */
  int gammaIsRWA() const { return gammaIsRWA_; };

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

#endif
