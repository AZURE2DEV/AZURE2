// Stub implementation for API functions when USE_API is disabled
#include "../include/Config.h"

int start_api(int port, Config &configure) {
  // Do nothing when API is disabled
  return 0;
}