#ifndef AZURESOCKET_H
#define AZURESOCKET_H

#include <iostream>
#include <string>
#include <cstring>
// Check if Windows or Linux
#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#include <windows.h>
#pragma comment(lib, "ws2_32.lib")
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif

// Windows compatibility definitions for close
#ifdef _WIN32
#define close closesocket
#define socklen_t int
#else
#include <unistd.h>
#endif

#ifdef _WIN32
using socket_t = SOCKET;
#else
using socket_t = int;
#endif

#include "AZUREAPI.h"

const int BUFFER_SIZE = 10000;

class AZURESocket {

private:
    int port_;
    socket_t serverSocket_;
    socket_t clientSocket_;
    struct sockaddr_in serverAddress_;
    struct sockaddr_in clientAddress_;

    AZUREAPI* api_;

public:
    AZURESocket(int port, AZUREAPI* api): port_(port), api_( api ) {
      serverSocket_ = -1;
      clientSocket_ = -1;
    };

    ~AZURESocket(){
      close(serverSocket_);
      close(clientSocket_);
    };

    bool start( );

    bool sendPacket( vector_r response );
    bool sendPacket( std::string );
    bool sendPacket( std::vector<bool> response );

};


#endif