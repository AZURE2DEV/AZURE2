#include "PairsTab.h"
#include "LevelsTab.h"

namespace {

std::vector<QString> setInfoStrings(QString infoString1 = "",
                                    QString infoString2 = "",
                                    QString infoString3 = "",
                                    QString infoString4 = "",
                                    QString infoString5 = "") {
  std::vector<QString> stringVector;
  if (infoString1 != "") stringVector.push_back(infoString1);
  if (infoString2 != "") stringVector.push_back(infoString2);
  if (infoString3 != "") stringVector.push_back(infoString3);
  if (infoString4 != "") stringVector.push_back(infoString4);
  if (infoString5 != "") stringVector.push_back(infoString5);
  return stringVector;
};

}  // namespace

const std::vector<QString> PairsTab::infoText =
    setInfoStrings("<ul><li><p>This tab is used to define the particle pairs of the reaction. These are the two particles that fuse together to form or result from the decay of the compound nucleus. The theory is limited to two particle interactions. &gamma;-ray and &beta;-decays are currently only supported as decay pairs and must be distinguished by the different drop down selection when the pair is created. Since these decay pairs are specialized, some of the information fields are automatically filled in and may not be edited. &gamma;-ray decays are limited to bound states.</p></ul>"
                   "<ul><li><p>Use the &#43; or &#45; buttons on the lower left corner to add or delete a decay pair respectively.</p></ul>"
                   "<ul><li><p>When a particle pair is created it is assigned a number that is displayed on the left hand side of the particle pair row. These numbers are referenced when creating segments on the <i>Segments</i> tab.</p></ul>"
                   "<ul><li><p>The tabs are meant to be filled out starting with this tab and moving to the right.  The information provided in this tab is the basis for several automatic calculations that are performed in the subsequent tabs, especially the <i>Levels and Channels</i> tab. If changes are later made to values in this tab after subsequent information has been filled out in other tabs, those changes will be applied automatically, possibly resulting in significant changes in the other tabs.</p></ul>");

const std::vector<QString> LevelsTab::infoText =
    setInfoStrings("<ul><li><p>This tab is used to create the different levels of the compound nucleus considered in the calculation. Levels can be created by clicking the &#43; button on the lower left corner. A level can be deleted by highlighting the level and then clicking the &#45; button. Levels may be above or below the particle separation energies. If a level is below a particle separation energy, the program automatically recognizes this, changing the label from a partial width to an ANC for the appropriate channels.</p></ul>"
                   "<ul><li><p>The energetically and momentum allowed reaction channels are determined automatically by the code for each level based on the information provided in the <i>Particle Pairs</i> tab and the spin-party and energy of the entered levels. If a change is made to values in the <i>Particle Pairs</i> tab, these changes are automatically applied to the levels calculated in this tab.</p></ul>"
                   "<ul><li><p>Levels can be included (excluded) from the calculation by checking (unchecking) the box on the far left for each level. A level can have its energy fixed by checking the box marked under the <b>Fix?</b> column. Partial widths or ANCs can be fixed by checking the <b>Fix?</b> box under the Selected Channels window.</p></ul>"
                   "<ul><li><p>&gamma;-ray decays are limited to bound states.</p></ul>"
                   "<ul><li><p>For &gamma;-ray calculations, the user should define a level corresponding to the each &gamma;-ray particle pair defined in the <i>Particle Pairs</i> tab. This allows the user to set the value of the ANCs and the &gamma;-ray decay widths for each of these bound states. The energies of these levels need to be the same as the excitation energies of the &gamma;-ray particle pairs.  </p></ul>"
                   "<ul><li><p>Relative interferences for channels can be specified by changing the sign of the corresponding partial width or ANC.</p></ul>"
                   "<ul><li><p>As one of the general features of <i>R</i>-matrix theory, the number of levels must be truncated both in number and in spin parity to some finite amount. One result of this is an ambiguity in the number of angular momentum terms that are summed in the hard sphere phase terms that are included in scattering and in the external capture calculations. To define these angular momenta, the code only uses those that are present in the channels of the levels defined in this tab. The user may therefore have to include dummy levels in the calculation. These dummy levels simply need to be created at an arbitrary energy with the spin-parities that are not present in the real levels under consideration. The energy can be fixed and the partial widths may be set to zero. The affects should be investigated for each allowed spin parity combination.</p></ul>"
                   "<ul><li><p>The user needs to be careful that the maximum orbital momentum is high enough to have at least one channels for each particle pair for each level if it allowed or the code will crash. </p></ul>");

