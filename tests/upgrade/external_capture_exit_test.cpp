// The exit of the external-capture pathways is found by the photon pair's number (ExternalCaptureExit.h,
// docs/upgrade/cnuc_split.md section 11.9, defect 1), on tests/ec_exit_position/models/neutron_pair_first.azr: pairs
// n + 12C (1, another compound nucleus), p + 12C (2), photon (3).
//   - the photon pair is the second exit of p + 12C: position 2, not 3;
//   - an exit that p + 12C has no pathway to (n + 12C) is refused with a message naming the data segment and the
//     two pair keys (1 and 3), never read at another position.
//
//   external_capture_exit_test <dir>/neutron_pair_first.azr
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

#include "CNuc.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "ExternalCaptureExit.h"
#include "PPair.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {
int failures = 0;
void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}
}  // namespace

int main(int argc, char **argv) {
  namespace fs = std::filesystem;
  if (argc != 2) {
    std::fprintf(stderr, "usage: external_capture_exit_test <dir>/neutron_pair_first.azr\n");
    return 2;
  }
  const fs::path file = fs::absolute(argv[1]);
  fs::current_path(file.parent_path());
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = file.filename().string();
  cfg.paramMask |= Config::USE_EXTERNAL_CAPTURE;
  g_config = &cfg;
  CNuc compound;
  if (cfg.ReadConfigFile() < 0) {
    check("reads the configuration", false, sink.str());
    return 1;
  }
  cfg.fileCheckMask = 0;
  cfg.screenCheckMask = 0;
  cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache();
  if (compound.Fill(cfg) == -1) {
    check("reads the levels", false, sink.str().substr(0, 400));
    return 1;
  }
  const int neutron = compound.GetPairNumFromKey(3), proton = compound.GetPairNumFromKey(1),
            photon = compound.GetPairNumFromKey(2);
  check("pairs numbered n + 12C 1, p + 12C 2, photon 3", neutron == 1 && proton == 2 && photon == 3,
        std::to_string(neutron) + " " + std::to_string(proton) + " " + std::to_string(photon));
  compound.GetPair(proton)->SetEntrance();
  compound.Initialize(cfg);

  int position = -1;
  try {
    position = ExternalCaptureExit(compound, proton, photon, 7);
  } catch (const std::exception &e) {
    check("the photon exit is found", false, e.what());
  }
  check("the photon pair (3) is exit 2 of p + 12C", position == 2, std::to_string(position));

  std::string message;
  try {
    position = ExternalCaptureExit(compound, proton, neutron, 7);
  } catch (const std::runtime_error &e) {
    message = e.what();
  }
  check("an exit without a pathway is refused", !message.empty(), "returned " + std::to_string(position));
  check("the refusal names the data segment and both pair keys",
        message.find("data segment #7") != std::string::npos && message.find("from pair 1 to pair 3") != std::string::npos,
        message);
  std::printf("external_capture_exit: %s (%d failure%s)\n", failures ? "FAIL" : "all checks passed", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
