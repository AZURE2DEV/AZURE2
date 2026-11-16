// What the compound nucleus writes: its check files and parameters.out / relativistic.out (LevelReport.h).  Reads the
// nucleus through const access only.
#include "LevelReport.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "AChannel.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Constants.h"
#include "JGroup.h"
#include "PPair.h"
#include "PathwayTable.h"
#include "azure2/physics/LevelParameters.h"

namespace {

/*!
 * A check goes to checks/<file> when the file mask has its bit, else to the screen when the screen mask has it.
 * `body` writes the table; if the file cannot be opened (and the screen mask is clear), one line says so.
 */
template <class Body>
void WriteCheck(const Config &configure, unsigned int bit, const char *file, const char *what, Body body) {
  std::streambuf *sbuffer = nullptr;
  std::filebuf fbuffer;
  if (configure.fileCheckMask & bit) {
    std::string outfile = configure.checkdir + file;
    fbuffer.open(configure.OpenPath(outfile).c_str(), std::ios::out);
    sbuffer = &fbuffer;
  } else if (configure.screenCheckMask & bit)
    sbuffer = configure.outStream.rdbuf();
  std::ostream out(sbuffer);
  if (((configure.fileCheckMask & bit) && fbuffer.is_open()) || (configure.screenCheckMask & bit))
    body(out);
  else
    configure.outStream << "Could not write " << what << " check file." << std::endl;
  out.flush();
  if (fbuffer.is_open()) fbuffer.close();
}

}  // namespace

void WriteCompoundNucleusCheck(const CompoundSystem &compound, const LevelInputs &inputs, const Config &configure) {
  WriteCheck(configure, Config::CHECK_COMPOUND_NUCLEUS, "compoundnucleus.chk", "compound nucleus", [&](std::ostream &out) {
    out << std::endl
        << "************************************" << std::endl
        << "*          Particle Pairs          *" << std::endl
        << "************************************" << std::endl;
    for (int i = 1; i <= compound.NumPairs(); i++) {
      const PPair *pair = compound.GetPair(i);
      out << "Pair Number: " << i << "  Pair Key: " << pair->GetPairKey();
      if (compound.NumNuclei() > 1) out << "  Compound Nucleus: " << compound.NucleusLabel(pair->GetNucleus());
      out << std::endl;
      out << std::setw(30) << "Light Particle J: " << pair->GetJ(1) << std::endl
          << std::setw(30) << "Light Particle Parity: " << pair->GetPi(1) << std::endl
          << std::setw(30) << "Light Particle Z: " << pair->GetZ(1) << std::endl
          << std::setw(30) << "Light Particle M: " << pair->GetM(1) << std::endl
          << std::setw(30) << "Light Particle g: " << pair->GetG(1) << std::endl
          << std::setw(30) << "Heavy Particle J: " << pair->GetJ(2) << std::endl
          << std::setw(30) << "Heavy Particle Parity: " << pair->GetPi(2) << std::endl
          << std::setw(30) << "Heavy Particle Z: " << pair->GetZ(2) << std::endl
          << std::setw(30) << "Heavy Particle M: " << pair->GetM(2) << std::endl
          << std::setw(30) << "Heavy Particle g: " << pair->GetG(2) << std::endl
          << std::setw(30) << "Seperation Energy [MeV]: " << pair->GetSepE() << std::endl
          << std::setw(30) << "Excitation Energy [MeV]: " << pair->GetHeavyExE() << std::endl;
      // the light particle's excitation (32nd column) only when it is set, so older outputs are unchanged
      if (pair->GetLightExE() != 0.) out << std::setw(30) << "Light Excitation [MeV]: " << pair->GetLightExE() << std::endl;
      out << std::setw(30) << "Channel Radius: " << pair->GetPairRadius() << std::endl;
      if (pair->GetPType() == 0)
        out << std::setw(30) << "Pair Type: " << "Particle,Particle" << std::endl;
      else if (pair->GetPType() == 10)
        out << std::setw(30) << "Pair Type: " << "Particle,Gamma" << std::endl;
      else if (pair->GetPType() == 20)
        out << std::setw(30) << "Pair Type: " << "Beta Decay" << std::endl;
      else
        out << std::setw(30) << "Pair Type: Unknown" << std::endl;
    }
    out << std::endl
        << "************************************" << std::endl
        << "*              Levels              *" << std::endl
        << "************************************" << std::endl
        << std::setw(11) << "J Group #"
        << std::setw(5) << "J"
        << std::setw(4) << "Pi"
        << std::setw(9) << "Level #"
        << std::setw(14) << "Energy [MeV]"
        << std::setw(11) << "Channel #"
        << std::setw(3) << "l"
        << std::setw(5) << "s"
        << std::setw(8) << "Pair #"
        << std::setw(11) << "Width"
        << std::setw(11) << "Rad. Type" << std::endl;
    for (int i = 1; i <= compound.NumJGroups(); i++) {
      const JGroup *group = compound.GetJGroup(i);
      for (int ii = 1; ii <= group->NumLevels(); ii++) {
        for (int iii = 1; iii <= group->NumChannels(); iii++) {
          out << std::setw(11) << i
              << std::setw(5) << group->GetJ()
              << std::setw(4) << group->GetPi()
              << std::setw(9) << ii
              << std::setw(14) << inputs.E(i, ii)
              << std::setw(11) << iii
              << std::setw(3) << group->GetChannel(iii)->GetL()
              << std::setw(5) << group->GetChannel(iii)->GetS()
              << std::setw(8) << group->GetChannel(iii)->GetPairNum()
              << std::setw(11) << inputs.Gamma(i, ii, iii)
              << std::setw(11) << group->GetChannel(iii)->GetRadType() << std::endl;
        }
      }
      out << std::endl;
    }
  });
}

void WritePathwaysCheck(const CompoundSystem &compound, const Config &configure) {
  WriteCheck(configure, Config::CHECK_PATHWAYS, "pathways.chk", "pathways", [&](std::ostream &out) {
    out << std::endl
        << "************************************" << std::endl
        << "*    Internal Reaction Pathways    *" << std::endl
        << "************************************" << std::endl
        << std::setw(17) << "Entrance Pair #"
        << std::setw(9) << "Decay #"
        << std::setw(14) << "Decay Pair #"
        << std::setw(12) << "K Group #"
        << std::setw(16) << "Entrance Ch. s"
        << std::setw(13) << "Decay Ch. s"
        << std::setw(11) << "M Group #"
        << std::setw(11) << "J Group #"
        << std::setw(16) << "Entrance Ch. #"
        << std::setw(12) << "Exit Ch. #" << std::endl;
    const azure2::PathwayTable &table = compound.Pathways();
    for (int i = 1; i <= compound.NumPairs(); i++) {
      const azure2::EntrancePathways &P = table.Entrance(i);
      for (int ii = 1; ii <= (int)P.exits.size(); ii++) {
        const azure2::ExitRow &decay = P.exits[ii - 1];
        for (int iii = 1; iii <= decay.num_groups; iii++) {
          const azure2::SpinGroupRow &kGroup = P.Group(decay, iii);
          for (int iiii = 1; iiii <= kGroup.num_internal; iiii++) {
            const azure2::InternalRow &mGroup = P.Internal(kGroup, iiii);
            out << std::setw(17) << i
                << std::setw(9) << ii
                << std::setw(14) << decay.pair
                << std::setw(12) << iii
                << std::setw(16) << kGroup.s
                << std::setw(13) << kGroup.sp
                << std::setw(11) << iiii
                << std::setw(11) << mGroup.group
                << std::setw(16) << mGroup.ch
                << std::setw(12) << mGroup.chp << std::endl;
          }
        }
        out << std::endl;
      }
    }
    out << std::endl
        << "************************************" << std::endl
        << "*    External Reaction Pathways    *" << std::endl
        << "************************************" << std::endl
        << std::setw(17) << "Entrance Pair #"
        << std::setw(9) << "Decay #"
        << std::setw(14) << "Decay Pair #"
        << std::setw(12) << "K Group #"
        << std::setw(16) << "Entrance Ch. s"
        << std::setw(13) << "Decay Ch. s"
        << std::setw(11) << "M Group #"
        << std::setw(11) << "Mult."
        << std::setw(11) << "J_i Value"
        << std::setw(11) << "J_f Value"
        << std::setw(11) << "l_i Value"
        << std::setw(11) << "l_f Value"
        << std::setw(13) << "Type"
        << std::setw(13) << "Ch. Decay #"
        << std::setw(11) << "Ch. K #"
        << std::setw(11) << "Ch. M #"
        << std::setw(11) << "Int. Ch #"
        << std::endl;
    for (int i = 1; i <= compound.NumPairs(); i++) {
      const azure2::EntrancePathways &P = table.Entrance(i);
      for (int ii = 1; ii <= (int)P.exits.size(); ii++) {
        const azure2::ExitRow &decay = P.exits[ii - 1];
        for (int iii = 1; iii <= decay.num_groups; iii++) {
          const azure2::SpinGroupRow &kGroup = P.Group(decay, iii);
          for (int iiii = 1; iiii <= kGroup.num_external; iiii++) {
            const azure2::ExternalRow &theECMGroup = P.External(kGroup, iiii);
            const JGroup *theECJGroup = compound.GetJGroup(theECMGroup.final_group);
            out << std::setw(17) << i
                << std::setw(9) << ii
                << std::setw(14) << decay.pair
                << std::setw(12) << iii
                << std::setw(16) << kGroup.s
                << std::setw(13) << theECJGroup->GetChannel(theECMGroup.final_channel)->GetS()
                << std::setw(11) << iiii
                << std::setw(10) << theECMGroup.rad_type << theECMGroup.mult
                << std::setw(11) << theECMGroup.J
                << std::setw(11) << theECJGroup->GetJ()
                << std::setw(11) << theECMGroup.l
                << std::setw(11) << theECJGroup->GetChannel(theECMGroup.final_channel)->GetL();
            if (theECMGroup.channel_capture)
              out << std::setw(13) << "Channel"
                  << std::setw(13) << theECMGroup.cc_exit
                  << std::setw(11) << theECMGroup.cc_group
                  << std::setw(11) << theECMGroup.cc_pathway
                  << std::setw(11) << theECMGroup.photon_channel << std::endl;
            else
              out << std::setw(15) << "Hard Sphere" << std::endl;
          }
        }
        out << std::endl;
      }
    }
  });
}

void WriteBoundaryConditionsCheck(const CompoundSystem &compound, const Config &configure) {
  WriteCheck(configure, Config::CHECK_BOUNDARY_CONDITIONS, "boundaryconditions.chk", "boundary conditions",
             [&](std::ostream &out) {
               out << std::endl
                   << "************************************" << std::endl
                   << "*        Boundary Conditions       *" << std::endl
                   << "************************************" << std::endl;
               out << std::setw(10) << "J Group #"
                   << std::setw(10) << "Channel #"
                   << std::setw(20) << "Boundary Condition"
                   << std::endl;
               for (int j = 1; j <= compound.NumJGroups(); j++) {
                 const JGroup *theJGroup = compound.GetJGroup(j);
                 if (!theJGroup->IsInRMatrix()) continue;
                 for (int ch = 1; ch <= theJGroup->NumChannels(); ch++)
                   out << std::setw(10) << j
                       << std::setw(10) << ch
                       << std::setw(20) << theJGroup->GetChannel(ch)->GetBoundaryCondition()
                       << std::endl;
               }
             });
}

void WriteAngularDistributionsCheck(const CompoundSystem &compound, const Config &configure) {
  WriteCheck(configure, Config::CHECK_ANGULAR_DISTS, "angulardistributions.chk", "angular distributions",
             [&](std::ostream &out) {
               out << std::endl
                   << "************************************" << std::endl
                   << "*       Angular Distributions       *" << std::endl
                   << "************************************" << std::endl;
               out << std::setw(10) << "ir"
                   << std::setw(10) << "k"
                   << std::setw(10) << "L"
                   << std::setw(10) << "m1"
                   << std::setw(10) << "m2"
                   << std::setw(10) << "z1z2"
                   << std::setw(10) << "type"
                   << std::endl;
               // The interferences of each exit, a blank line after each (k, L) group.
               const azure2::PathwayTable &table = compound.Pathways();
               for (int aa = 1; aa <= compound.NumPairs(); aa++) {
                 const azure2::EntrancePathways &P = table.Entrance(aa);
                 for (const azure2::ExitRow &theDecay : P.exits) {
                   for (int i = 0; i < theDecay.num_interferences; i++) {
                     const azure2::FlatInterference &f = P.interferences[theDecay.first_interference + i];
                     out << std::setw(10) << theDecay.pair
                         << std::setw(10) << f.k
                         << std::setw(10) << f.l
                         << std::setw(10) << f.m1
                         << std::setw(10) << f.m2
                         << std::setw(10) << f.z
                         << std::setw(10) << std::string{f.r1, f.r2}
                         << std::endl;
                     const bool lastOfGroup = i + 1 == theDecay.num_interferences ||
                                              P.interferences[theDecay.first_interference + i + 1].k != f.k ||
                                              P.interferences[theDecay.first_interference + i + 1].l != f.l;
                     if (lastOfGroup) out << std::endl;
                   }
                 }
               }
             });
}

void WriteParametersOut(const CompoundSystem &compound, const LevelState &state, const ObservedLevelTable &table,
                        const Config &configure) {
  namespace phys = azure2::physics;
  char filename[256];
  snprintf(filename, sizeof(filename), "%sparameters.out", configure.outputdir.c_str());
  std::ofstream out;
  out.open(configure.OpenPath(filename));
  if (out) {
    out << "PHYSICAL LEVEL PARAMETERS (BOUNDARY CONDITION SET TO SHIFT AT E_LEVEL)" << std::endl;
    out << std::endl;
    for (int j = 1; j <= compound.NumJGroups(); j++) {
      const JGroup *theJGroup = compound.GetJGroup(j);
      for (int la = 1; la <= theJGroup->NumLevels(); la++) {
        out.precision(1);
        out << std::fixed;
        const ObservedLevel &theLevel = table.Level(j, la);
        out << "J = " << std::setw(3) << theJGroup->GetJ();
        if (theJGroup->GetPi() == -1)
          out << '-';
        else
          out << '+';
        out.precision(4);
        out << "  E_level = " << std::setw(8) << theLevel.E << " MeV"
            << "  ITERATIONS = " << std::setw(5) << theLevel.iterations;
        if (compound.NumNuclei() > 1) out << "  NUCLEUS " << compound.NucleusLabel(theJGroup->GetNucleus());
        if (compound.GetIsospin().Active()) {
          const std::string iso = compound.GetIsospin().MemberLabel(j, la);
          if (!iso.empty()) out << "  " << iso;
        }
        out << std::endl;
        for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
          const AChannel *theChannel = theJGroup->GetChannel(ch);
          const PPair *exitPair = compound.GetPair(theChannel->GetPairNum());
          out << "  R = " << std::setw(2) << exitPair->GetPairKey();
          if (theChannel->GetRadType() == 'P')
            out << "  l = " << std::setw(3) << theChannel->GetL();
          else if (theChannel->GetRadType() == 'F')
            out << "  Fermi Beta Decay ";
          else if (theChannel->GetRadType() == 'G')
            out << "   G-T Beta Decay  ";
          else
            out << "  L = " << std::setw(2) << theChannel->GetRadType() << theChannel->GetL();
          out.precision(1);
          if (theChannel->GetRadType() != 'G' && theChannel->GetRadType() != 'F')
            out << "  s = " << std::setw(4) << theChannel->GetS();
          out.precision(6);
          // The observed quantities of TransformOut (physics/ observed_levels): W and the reported value by kind.
          const double W = theLevel.W[ch - 1];
          switch (static_cast<phys::ObservedKind>(theLevel.kind[ch - 1])) {
            case phys::ObservedKind::ANC:
              out << "  C  = " << std::setw(12) << sqrt(fabs(W)) << " fm^(-1/2)";
              break;
            case phys::ObservedKind::MagneticMoment:
              out << "  mu = " << std::setw(12) << theLevel.value[ch - 1] << " nm       ";
              break;
            case phys::ObservedKind::QuadrupoleMoment:
              out << "  Q  = " << std::setw(12) << theLevel.value[ch - 1] << " b        ";
              break;
            case phys::ObservedKind::Beta:
              out << "  B  = " << std::setw(12) << W << "          ";
              break;
            case phys::ObservedKind::Width:
              if (fabs(W) >= 1e-3)
                out << "  G  = " << std::setw(12) << fabs(W) * 1e3 << " keV      ";
              else if (fabs(W) >= 1e-6)
                out << "  G  = " << std::setw(12) << fabs(W) * 1e6 << " eV       ";
              else
                out << "  G  = " << std::setw(12) << fabs(W) * 1e9 << " meV      ";
              break;
          }
          out << "  g_int = " << std::setw(12) << theLevel.gamma[ch - 1];
          if (theChannel->GetRadType() != 'G' && theChannel->GetRadType() != 'F')
            out << " MeV^(1/2) ";
          else
            out << "           ";
          out << "  g_ext = " << std::setw(20) << theLevel.external[ch - 1];
          if (theChannel->GetRadType() != 'G' && theChannel->GetRadType() != 'F') out << " MeV^(1/2) ";
          out << std::endl;
        }
        out << std::endl;
      }
    }
    out << compound.GetGammaDecays().Report(compound, state, table);
  } else
    configure.outStream << "Could not save parameters.out file." << std::endl;
  if (configure.relativistic) WriteRelativisticOut(compound, table, configure);
  // isospin.out: the isospin levels at the same parameters (IsospinTies.h).
  compound.GetIsospin().Write(configure, compound, state);
}

/*!
 * relativistic.out (runs with relativistic=1): the physical level parameters of parameters.out in the level-energy
 * convention of the LANL code EDA (ENDF KRL=1, LA-UR-19-31084).  AZURE2 keeps the R-matrix energy linear in the
 * excitation energy E_x = sqrt(s) - M; EDA uses E(s) = (s - m0^2) / (2 m0), m0 the mass of a reference pair.  With
 * T0 = E_x - S0 the reference pair's kinetic energy, E(s) = T0 + T0^2 / (2 m0), written here on the excitation scale
 * (+ S0).  Widths follow the energy derivative dE(s)/dE_x = 1 + T0/m0 at the level (the penetrabilities are the same
 * in both): Gamma' = Gamma (1 + T0/m0), gamma' = gamma sqrt(1 + T0/m0).  This is exact at each level energy; away
 * from it the two energy dependences differ at second order in T0/m0.
 */
void WriteRelativisticOut(const CompoundSystem &compound, const ObservedLevelTable &table, const Config &configure) {
  char filename[256];
  snprintf(filename, sizeof(filename), "%srelativistic.out", configure.outputdir.c_str());
  std::ofstream out(configure.OpenPath(filename));
  if (!out) {
    configure.outStream << "Could not save relativistic.out file." << std::endl;
    return;
  }
  // the reference pair of each compound nucleus: the configured key if it belongs to it, else its first particle pair
  auto reference = [&](int nucleus) -> const PPair * {
    const PPair *first = nullptr;
    for (int p = 1; p <= compound.NumPairs(); p++) {
      const PPair *pp = compound.GetPair(p);
      if (pp->GetPType() != 0 || pp->GetNucleus() != nucleus) continue;
      if (pp->GetPairKey() == configure.relativisticReferencePair) return pp;
      if (!first) first = pp;
    }
    return first;
  };
  out << "# Relativistic kinematics (relativistic=1).  Level parameters of parameters.out in the E(s) convention" << std::endl
      << "# of the LANL code EDA (ENDF KRL=1): E(s) = (s - m0^2)/(2 m0) for the reference pair (mass m0, threshold S0)," << std::endl
      << "# here E_s = S0 + T0 + T0^2/(2 m0) with T0 = E_x - S0; Gamma_s = Gamma (1 + T0/m0), g_s = g sqrt(1 + T0/m0)." << std::endl
      << "# relativisticReferencePair=<key> in <config> chooses the reference pair (default: the first particle pair)." << std::endl;
  out << std::fixed;
  for (int j = 1; j <= compound.NumJGroups(); j++) {
    const JGroup *theJGroup = compound.GetJGroup(j);
    const PPair *ref = reference(theJGroup->GetNucleus());
    if (!ref) continue;
    const double m0 = (ref->GetM(1) + ref->GetM(2)) * uconv;
    const double s0 = ref->GetSepE() + ref->GetExE();
    for (int la = 1; la <= theJGroup->NumLevels(); la++) {
      const ObservedLevel &theLevel = table.Level(j, la);
      const double ex = theLevel.E;
      const double t0 = ex - s0;
      const double d = 1.0 + t0 / m0;
      out.precision(1);
      out << "J = " << std::setw(4) << theJGroup->GetJ() << (theJGroup->GetPi() == -1 ? '-' : '+');
      out.precision(6);
      out << "  E_x = " << std::setw(11) << ex << " MeV  E_s = " << std::setw(11) << s0 + t0 + t0 * t0 / (2.0 * m0)
          << " MeV  (E_s - E_x = " << std::setw(9) << 1e3 * t0 * t0 / (2.0 * m0) << " keV, reference R = "
          << ref->GetPairKey() << ", dE_s/dE_x = " << std::setprecision(8) << d << ")" << std::endl;
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        const AChannel *theChannel = theJGroup->GetChannel(ch);
        if (theChannel->GetRadType() != 'P') continue;
        const PPair *exitPair = compound.GetPair(theChannel->GetPairNum());
        out.precision(1);
        out << "    R = " << std::setw(2) << exitPair->GetPairKey() << "  l = " << std::setw(2) << theChannel->GetL()
            << "  s = " << std::setw(4) << theChannel->GetS();
        out.precision(6);
        const double localEnergy = ex - exitPair->GetSepE() - exitPair->GetExE();
        if (localEnergy < 0.0) {
          out << "  (bound in this channel)";
        } else {
          out << "  G = " << std::setw(13) << fabs(theLevel.W[ch - 1]) * 1e3 << " keV  G_s = " << std::setw(13)
              << fabs(theLevel.W[ch - 1]) * d * 1e3 << " keV";
        }
        out << "  g = " << std::setw(13) << theLevel.gamma[ch - 1] << "  g_s = " << std::setw(13)
            << theLevel.gamma[ch - 1] * sqrt(d) << " MeV^(1/2)" << std::endl;
      }
      out << std::endl;
    }
  }
}
