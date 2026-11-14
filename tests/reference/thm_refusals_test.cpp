// The refusals of the engine's THM set-up say what is wrong, what is allowed and how to go on: each case below
// builds one wrong <thm> experiment (or option) and checks that the refusal names the way on, as the text the set-up
// prints (CheckThmExperiments, CheckThmDwVertexOptions, CheckThmCoulombConsistency, ParseThmExperimentLine,
// ThmOptionEdits, ThmResolveReaction).  Before these parts were written, every case but the last two printed the
// rule alone.
#include <iostream>
#include <string>
#include <vector>

#include "ThmBlockEdits.h"
#include "ThmExperiment.h"

namespace {

int failures = 0;

void check(const std::string &name, const std::string &why, const std::string &part) {
  const bool ok = why.find(part) != std::string::npos;
  std::cout << (ok ? "  ok    " : "  FAIL  ") << name;
  if (!ok) std::cout << "  -- '" << part << "' not in: " << why;
  std::cout << std::endl;
  if (!ok) failures++;
}

// The refusal of CheckThmExperiments for the experiment lines given; each experiment gets data set `segment`.
std::string Checked(const std::vector<std::string> &lines, int segment = 1) {
  std::vector<ThmExperiment> xs;
  for (const std::string &l : lines) {
    const std::string why = ParseThmExperimentLine(l, xs);
    if (!why.empty()) return "parse: " + why;
  }
  for (ThmExperiment &x : xs)
    if (segment > 0) x.segments.push_back(segment);
  return CheckThmExperiments(xs);
}

}  // namespace

int main() {
  const std::string kin = " beam=15N target=d spectator=n Ebeam=60";
  check("partial kinematics: the keys missing", Checked({"experiment[A] beam=15N"}),
        "give target=, spectator=, Ebeam= or remove the others");
  check("lineshape=on without kinematics", Checked({"experiment[A] lineshape=on"}), "or set lineshape=off");
  check("a ps window without kinematics", Checked({"experiment[A] ps=hulthen:0-40"}), "or set ps=delta");
  check("psNodes= without a window", Checked({"experiment[A] psNodes=8"}), "give the window or remove psNodes=");
  check("distortion=coulomb without kinematics", Checked({"experiment[A] distortion=coulomb"}),
        "or set distortion=none");
  check("opticalAA= without distortion=optical",
        Checked({"experiment[A]" + kin + " distortion=coulomb opticalAA=plane"}), "set it, or remove them");
  check("spectatorAngle= without a distortion", Checked({"experiment[A]" + kin + " spectatorAngle=10"}),
        "set one, or remove spectatorAngle=");
  check("vertexModel=dw without a distortion", Checked({"experiment[A]" + kin + " vertexModel=dw"}),
        "or vertexModel=pw");
  check("distortionRef= with vertexModel=dw",
        Checked({"experiment[A]" + kin + " distortion=coulomb distortionRef=0.7 vertexModel=dw"}),
        "remove distortionRef= or set vertexModel=pw");
  check("spectatorAngleNodes= without a window", Checked({"experiment[A]" + kin + " spectatorAngleNodes=8"}),
        "give the window, or remove spectatorAngleNodes=");
  check("an experiment no data set names", Checked({"experiment[A] background=linear"}, 0),
        "give a data set the field 'thm A'");
  check("a data set in two experiments", Checked({"experiment[A] background=linear", "experiment[B] background=const"}),
        "a data set belongs to one experiment");
  check("lineshape=on with cbackground=", Checked({"experiment[A]" + kin + " lineshape=on cbackground=1/2+:2"}),
        "remove cbackground= or set lineshape=off");
  std::vector<ThmExperiment> xs;
  check("an unknown key: where the keys are listed", ParseThmExperimentLine("experiment[A] foo=1", xs),
        "pyazr.kinds()[\"thm\"]");
  check("a malformed value: where the keys are listed", ParseThmExperimentLine("experiment[A] Ebeam=-1", xs),
        "pyazr.kinds()[\"thm\"]");
  check("entranceL=coherent with vertexModel=dw", CheckThmDwVertexOptions(true, 0.0, 1),
        "remove entranceL=coherent or set vertexModel=pw");
  check("spectatorEnergy with vertexModel=dw", CheckThmDwVertexOptions(false, 0.3, 2),
        "remove spectatorEnergy[2]= or set it to 0");
  xs.clear();
  ParseThmExperimentLine("experiment[A]" + kin + " vertexModel=dw distortion=coulomb", xs);
  check("vertexModel=dw with coulombIntegral=1", CheckThmCoulombConsistency(xs, true, nullptr),
        "set coulombIntegral=0, or vertexModel=pw");
  std::vector<ThmLineEdit> edits;
  std::string why;
  ThmOptionEdits({}, "weight[1]", "w.dat", ThmChange::Set, edits, why);
  check("a weight= option in format 3", why, "set that field of the record");
  ThmOptionEdits({}, "vertex", "x", ThmChange::Remove, edits, why);
  check("removing an option the block has not", why, "Project.thm() lists them");
  std::cout << (failures ? "FAIL: " + std::to_string(failures) + " case(s)" : std::string("all ok")) << std::endl;
  return failures ? 1 : 0;
}
