#include <iostream>

#include "AZUREAPI.h"
#include "AZURESocket.h"
#include "Config.h"

int start_api(int port, Config &configure) {
  AZUREAPI *azureApi_ = new AZUREAPI(configure);
  AZURESocket *azureSocket_ = new AZURESocket(port, azureApi_);

  azureSocket_->start();

  return 0;
}
