#ifndef NUCLINE_H
#define NUCLINE_H

#include <iostream>

///A class to read and store a line from a nuclear input file.

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
    stream >> levelJ_ >> levelPi_ >> levelE_ >> levelFix_ >> aa_ >> ir_ 
	   >> s_ >> l_ >> levelID_ >> isActive_ >> channelFix_ >> gamma_ >> j1_ >> pi1_
	   >> j2_ >> pi2_ >> e2_ >> m1_ >> m2_ >> z1_ >> z2_
	   >> entranceSepE_ >> sepE_ >> isupos_ >> secondaryDecayL_ >> Ic_ 
	   >> pType_ >> chRad_ >> g1_ >> g2_ >> ecMultMask_;
    s_/=2.;
    l_/=2;
  };
  /*!
   * Returns the spin of the corresponding level for the channel line.
   */
  double levelJ() const {return levelJ_;};
  /*!
   * Returns the parity of the corresponding level for the channel line.
   */
  int levelPi() const {return levelPi_;};
  /*!
   * Returns the excitation energy of the corresponding level for 
   * the channel line.
   */
  double levelE() const {return levelE_;};
  /*!
   * Returns non-zero if the corresponding level for the channel 
   * line is not to be varied in the fit.
   */
  int levelFix() const {return levelFix_;};
  /*!
   * Returns the entrance key of the corresponding particle pair 
   * for the channel line.
   */
  int aa() const {return aa_;};
  /*!
   * Returns the exit key of the corresponding particle pair for 
   * the channel line.
   */
  int ir() const {return ir_;};
  /*! 
   * Returns the channel spin of the channel.
   */
  double s() const {return s_;};
  /*!
   * Returns the orbital angular momentum of the channel.
   */
  int l() const {return l_;};
  /*!
   * Returns an indexing variable used by the graphical setup program.
   */
  int levelID() const {return levelID_;};
  /*!
   * Returns non-zero if the corresponding level for the channel line
   * is to be used in the calculation.
   */
  int isActive() const {return isActive_;};
  /*!
   * Returns non-zero if the reduced width amplitude for the channel is not
   * to be varied. 
   */
  int channelFix() const {return channelFix_;};
  /*!
   * Returns the initial reduced width or physical amplitude for the channel.
   */
  double gamma() const {return gamma_;};
  /*!
   * Returns the spin of the light particle in the corresponding pair.
   */
  double j1() const {return j1_;};
  /*!
   * Returns the parity of the light particle in the corresponding pair.
   */
  int pi1() const {return pi1_;};
  /*!
   * Returns the spin of the heavy particle in the corresponding pair.
   */
  double j2() const {return j2_;};
  /*!
   * Returns the parity of the heavy particle in the corresponding pair.
   */
  int pi2() const {return pi2_;};
  /*!
   * Returns the excitation energy of the heavy particle in the 
   * corresponding pair.
   */
  double e2() const {return e2_;};
  /*!
   * Returns the mass of the light particle in the corresponding pair.
   */
  double m1() const {return m1_;};
  /*!
   * Returns the mass of the heavy particle in the corresponding pair.
   */
  double m2() const {return m2_;};
  /*!
   * Returns the charge of the light particle in the corresponding pair.
   */
  int z1() const {return z1_;};
  /*!
   * Returns the charge of the heavy particle in the corresponding pair.
   */
  int z2() const {return z2_;};
  /*!
   * This function is depriciated and not used.
   */
  double entranceSepE() const {return entranceSepE_;};
  /*!
   * Returns the separation energy for the corresponding particle pair.
   */
  double sepE() const {return sepE_;};
  /*!
   * Returns 0 for particle-particle and 10 for particle-gamma
   * and 20 for beta delayed particle emmission 
   * types in the corresponding pair.
   */
  int pType() const {return pType_;};
  /*!
   * Returns the channel radius for the corresponding particle pair. 
   */
  double chRad() const {return chRad_;};
  /*!
   * Returns the g-factor for the light particle in the corresponding pair. 
   */
  double g1() const {return g1_;};
  /*!
   * Returns the g-factor for the heavy particle in the corresponding pair. 
   */
  double g2() const {return g2_;};
  /*!
   * Returns the external capture multiplicity mask for the corresponding pair. 
   */
  double ecMultMask() const {return ecMultMask_;};
  /*!
   * Returns 1 if angular distrabutions need to be calculated for unobserved
   * primary, observed secondary decays.
   */
  int isUPOS() const {return 1;}; //isupos_
  /*!
   * Returns the angular momentum of the decay for unobserved particle
   * observed secondary decay
   */
//  int secondaryDecayL() const {return secondaryDecayL_;}; //secondaryDecayL_
  /*!
   * Returns the spin of the final state of the decay for unobserved particle
   * observed secondary decay
   */
//  double Ic() const {return Ic_;}; //Ic_

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
  int pType_;
  double chRad_;
  double g1_;
  double g2_;
  unsigned int ecMultMask_;
  int isupos_;
  int secondaryDecayL_;
  double Ic_;
};

#endif
