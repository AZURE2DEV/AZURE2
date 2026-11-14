// A run request turned into the run's Config (include/RunPreparation.h).
#include "RunPreparation.h"

#include "Config.h"
#include "ProjectSetup.h"
#include "RunDiagnostics.h"
#include "RunMenu.h"
#include "azure2/io/LevelRows.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

namespace {

/*!
 * Whether the run has external capture: a level with an external-capture multipolarity feeds the exit of one of
 * `exits` (the active data sets', or the rate's), or one of them is the total capture (-1).  Sets
 * Config::USE_EXTERNAL_CAPTURE (cleared again, with a warning, under Reich-Moore); false when the <levels> cannot be
 * read (said on the output).
 */
bool checkExternalCapture(Config &configure, const std::vector<std::pair<int, int>> &segPairs) {
  std::istringstream in = SectionText(configure, "levels");
  if (!in) {
    configure.outStream << "Cannot read nuclear information. Check configuration file." << std::endl;
    return false;
  }
  std::string line = "";
  while (line != "<levels>" && !in.eof()) getline(in, line);
  if (line != "<levels>") {
    configure.outStream << "Cannot read nuclear information. Check configuration file." << std::endl;
    return false;
  }
  line = "";
  while (line != "</levels>" && !in.eof() &&
         !(configure.paramMask & Config::USE_EXTERNAL_CAPTURE)) {
    getline(in, line);
    bool empty = true;
    for (unsigned int i = 0; i < line.size(); ++i)
      if (line[i] != ' ' && line[i] != '\t') {
        empty = false;
        break;
      }
    if (empty == true) continue;
    if (line != "</levels>" && !in.eof()) {
      std::istringstream stm;
      stm.str(line);
      azure2::io::LevelRow tempNucLine(stm);
      if ((stm.rdstate() & (std::stringstream::failbit | std::stringstream::badbit))) {
        configure.outStream << "Problem reading nuclear information. Check configuration file.  The line:\n  "
                            << line << std::endl;
        return false;
      }
      if (tempNucLine.ecMultMask() != 0) {
        for (std::size_t i = 0; i < segPairs.size(); i++) {
          if (tempNucLine.ir() == segPairs[i].second || segPairs[i].second == -1) {
            configure.paramMask |= Config::USE_EXTERNAL_CAPTURE;
            break;
          }
        }
      }
    }
  }
  in.clear();
  if ((configure.paramMask & Config::USE_EXTERNAL_CAPTURE) && (configure.paramMask & Config::USE_RMC_FORMALISM)) {
    configure.outStream << "WARNING: External capture is not compatible with Reich-Moore.  Ignoring external capture." << std::endl;
    configure.paramMask &= ~Config::USE_EXTERNAL_CAPTURE;
  }
  return true;
}

/// The calculation of `command`: the bits of the Config that select it.
void SelectCalculation(RunCommand command, Config &configure) {
  switch (command) {
    case RunCommand::kCalculate:
      break;
    case RunCommand::kFit:
      configure.paramMask |= Config::PERFORM_FIT;
      break;
    case RunCommand::kExtrapolate:
      configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
      break;
    case RunCommand::kMinos:
      configure.paramMask |= Config::PERFORM_FIT;
      configure.paramMask |= Config::PERFORM_ERROR_ANALYSIS;
      break;
    case RunCommand::kRate:
      configure.paramMask &= ~Config::CALCULATE_WITH_DATA;
      configure.paramMask |= Config::CALCULATE_REACTION_RATE;
      break;
    case RunCommand::kMcmc:
      configure.paramMask |= Config::PERFORM_MCMC;
      break;
  }
}

/// Whether the run has external capture (Config::USE_EXTERNAL_CAPTURE set from the project and the rate's pairs);
/// stops the run when the <levels> cannot be read.
bool RunHasCapture(Config &configure, const RunRequest &request) {
  std::vector<std::pair<int, int>> exits;
  if (request.command == RunCommand::kRate)
    exits.push_back({request.rateEntrance, request.rateExit});
  else  // the active data sets (or extrapolations, for a run without data) from their specs (ProjectSetup.h)
    exits = ProjectReactionPairs(configure, (configure.paramMask & Config::CALCULATE_WITH_DATA) != 0);
  if (!checkExternalCapture(configure, exits)) ExitRun(configure, 1);
  // the rate reads no amplitude file
  return (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) && request.command != RunCommand::kRate;
}

/// The answers into the Config (the calculation's bits are set already).
void ApplyAnswers(const RunRequest &request, Config &configure) {
  if (request.covarianceBand) configure.paramMask |= Config::CALCULATE_COVARIANCE_BAND;
  if (request.scaleCovariance) configure.paramMask |= Config::SCALE_COVARIANCE_BY_CHI2;
  if (request.command == RunCommand::kMinos) configure.chiVariance = request.deltaChi2;
  if (!request.parameterFile.empty()) {
    configure.paramMask |= Config::USE_PREVIOUS_PARAMETERS;
    configure.paramfile = request.parameterFile;
  }
  if (!request.captureFile.empty()) {
    configure.paramMask |= Config::USE_PREVIOUS_INTEGRALS;
    configure.integralsfile = request.captureFile;
  }
  if (request.command == RunCommand::kRate) {
    RateParams &rate = configure.rateParams;
    rate.entrancePair = request.rateEntrance;
    rate.exitPair = request.rateExit;
    rate.useFile = request.temperaturesFromFile;
    rate.minTemp = rate.maxTemp = rate.tempStep = -1.;  // the grid's "not given"
    if (rate.useFile) {
      rate.temperatureFile = request.temperatureFile;
    } else {
      rate.minTemp = request.temperatureMin;
      rate.maxTemp = request.temperatureMax;
      rate.tempStep = request.temperatureStep;
    }
  }
}

}  // namespace

std::istringstream SectionText(const Config &configure, const std::string &name) {
  std::istringstream in;
  if (std::unique_ptr<std::istream> section = configure.SectionStream(name)) {
    std::ostringstream text;
    text << section->rdbuf();
    in.str(text.str());
  } else {
    in.setstate(std::ios::failbit);
  }
  return in;
}

std::string RequestFileRefusal(const Config &configure, const RunRequest &request) {
  for (const RunQuestion &q : RunQuestions()) {
    if (q.kind != AnswerKind::kFile || (request.*q.text).empty()) continue;
    const std::string &name = request.*q.text;
    if (std::ifstream(configure.OpenPath(name).c_str())) continue;
    return CommandRefusal(std::string(q.option) + "=" + name + ": cannot read " + configure.OpenPath(name) +
                              " (a file an option names is relative to the directory of the project file).",
                          request.command,
                          std::string("Name a file that exists there, or leave ") + q.option + " out.");
  }
  return "";
}

void PrepareRun(Config &configure, RunRequest &request, RunMenu *menu) {
  SelectCalculation(request.command, configure);
  const std::vector<RunQuestion> &questions = RunQuestions();
  const unsigned command = CommandBit(request.command);
  for (std::size_t i = 0; i < questions.size(); i++) {
    const RunQuestion &q = questions[i];
    // The project decides here, after the rate's pairs and before the amplitude file, whether the run has
    // external capture: for every command, as the menu always did.
    if (q.when == AskedWhen::kProjectNeedsCapture && !RunHasCapture(configure, request)) {
      if (!menu && !(request.*q.text).empty()) {
        configure.outStream << "Notice: " << q.option << "=" << request.*q.text
                            << " is not used: this run has no external capture." << std::endl;
        (request.*q.text).clear();
      }
      continue;
    }
    if (!menu || !(q.menuAsks & command) || !QuestionApplies(q, request)) continue;
    menu->Ask(q, request);
    for (const RunRule &rule : RunRules()) {
      if (std::string(rule.last) != q.option || rule.holds(request)) continue;
      menu->Refuse(rule);
      for (std::size_t first = 0; first < questions.size(); first++)
        if (std::string(rule.first) == questions[first].option) i = first - 1;  // asked again from `first`
    }
  }
  ApplyAnswers(request, configure);
}
