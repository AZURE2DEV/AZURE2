#include "PairsTab.h"
#include "LevelsTab.h"
#include "PotentialPanel.h"

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


// The Potential panel of the Particle Pairs tab (docs/upgrade/qt_gui.md section 19; the physics:
// docs/upgrade/potential_background.md).
const QString PotentialPanel::helpText = QObject::tr(
    "<p>A pair with a <b>potential</b> record takes its channel functions from the radial equation with that nuclear "
    "potential (a Woods-Saxon or Gaussian central term, the Coulomb potential of a charged sphere): the potential acts "
    "between the channel radius a<sub>c</sub> and the matching radius.  This panel shows the record of the pair "
    "selected above, partial wave by partial wave.</p>"
    "<ul><li><p><b>Potential on/off.</b>  Off keeps the record in the file as comment lines (<code>#off ...</code>): the "
    "engine ignores it, and on again restores it exactly as it was.</p></li>"
    "<li><p><b>Mode of a wave.</b>  <i>Central</i>: no <code>wave</code> line; the central potential acts beyond "
    "a<sub>c</sub>.  <i>Exterior only</i>: a <code>wave</code> line with the wave's own depth (and Gaussian width), "
    "acting beyond a<sub>c</sub>.  <i>Exterior + interior</i>: also the R-function of the potential inside "
    "a<sub>c</sub> as the background of the wave's R-matrix element (<code>background interior boundary B</code>); "
    "the fitted levels sit on top of it.  <i>None</i>: point-Coulomb channel functions in this wave; its depth only "
    "gives the bound-state shape of interior capture (<code>only bound_state</code>, allowed only when a wave has "
    "<code>capture interior</code>).  B<sub>bg</sub> is the boundary condition of the background; <i>Subtract n</i> "
    "removes the pole of the box state with n nodes (a state a level of the project represents); <i>Capture "
    "interior</i> adds the interior direct capture of the background.  What is refused is refused by the engine, "
    "with its words, under the table.</p></li>"
    "<li><p><b>In the calculation.</b>  The engine makes a J<sup>&pi;</sup> group only for the J<sup>&pi;</sup> of an "
    "active level, or, with <i>hard-sphere scattering in every J<sup>&pi;</sup></i> (<code>includeAllJPi</code>) on, "
    "for every J<sup>&pi;</sup> up to <code>maxOrbitalMomentum</code>.  A wave whose J<sup>&pi;</sup> has neither is "
    "left out, and its potential does nothing: the panel says so in red.</p></li>"
    "<li><p><b>Poles.</b>  The box states of the selected wave: the eigenvalues E<sub>p</sub> of the potential inside "
    "a<sub>c</sub> with a<sub>c</sub>u'/u = B<sub>bg</sub> (B = 0 for a wave without background), their residues "
    "&gamma;<sub>p</sub><sup>2</sup> and node numbers.  The interior R-function is "
    "&Sigma;<sub>p</sub> &gamma;<sub>p</sub><sup>2</sup>/(E<sub>p</sub> &minus; E): these are the poles a background "
    "brings, and they move with a<sub>c</sub> and B.  Deep ones are the potential's (Pauli-forbidden) bound states.  "
    "Beside them, the project's levels of the same J<sup>&pi;</sup>: a level near an unsubtracted pole counts that "
    "state twice.</p></li>"
    "<li><p><b>Tuning.</b>  The scattering length a and effective range r of the potential alone (its own phase "
    "shifts, the Coulomb-modified effective-range expansion for a charged pair and l = 0), with the number of its bound "
    "states.  <i>Depth for a</i> finds the wave's depth that gives a chosen a with the same number of bound states (a "
    "falls from +&infin; to &minus;&infin; between the depths at which states are bound, so there is one such depth). "
    "<i>Depth and range for a, r</i> finds the depth and the Gaussian width r<sub>0</sub> that give both, on the same "
    "branch, and says so when the branch has no such pair (it shows the range of r it found).  The scattering length "
    "is the more robust input: it is the intercept of the expansion; the effective range, its slope, depends more on "
    "the potential's shape and on the energy window, so fix a first and use r as a second, softer condition.  "
    "<i>Depth to bind</i> finds the depth within &plusmn;30 % that binds the wave at a level's energy, with the "
    "potential's own ANC beside the level's.  <i>Use this depth</i> writes it into the wave; Apply writes the record.</p>"
    "</li></ul>");
