#ifndef AZURESOCKET_H
#define AZURESOCKET_H

#include <iostream>
#include <string>
#include <cstring>
#include <cstdint>
#include <vector>
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

// Upper bound on the number of doubles accepted in a single frame.  This is a
// sanity guard against a corrupt/garbage length prefix causing a runaway
// allocation; it is NOT a hard limit on legitimate payloads in the way the old
// fixed BUFFER_SIZE was.  ~1e8 doubles == 800 MB, far above any real payload.
const std::uint64_t MAX_FRAME_DOUBLES = 100000000ULL;

class AZURESocket {

private:
    int port_;
    socket_t serverSocket_;
    socket_t clientSocket_;
    struct sockaddr_in serverAddress_;
    struct sockaddr_in clientAddress_;

    AZUREAPI* api_;

    // Low-level, loop-until-complete socket helpers.  recv/send may transfer
    // fewer bytes than requested, so we must loop.  Return false on
    // error/disconnect.
    bool recvAll( void* buffer, size_t length );
    bool sendAll( const void* buffer, size_t length );

    // Length-prefixed framing: [uint64 count][count * double].
    // readMessage returns false on a clean disconnect or fatal error.
    bool readMessage( vector_r& out );
    bool writeMessage( const double* values, std::uint64_t count );

    // Dispatch a single decoded request frame.
    void handle( const vector_r& request );

public:
    AZURESocket(int port, AZUREAPI* api): port_(port), api_( api ) {
      serverSocket_ = -1;
      clientSocket_ = -1;
    };

    ~AZURESocket(){
      if( clientSocket_ != -1 ) close(clientSocket_);
      if( serverSocket_ != -1 ) close(serverSocket_);
    };

    bool start( );

    bool sendPacket( const vector_r& response );
    bool sendPacket( const std::string& response );
    bool sendPacket( const std::vector<bool>& response );

};


#endif
