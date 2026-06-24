#include "AZURESocket.h"

#include "Constants.h"

#include <vector>
#include <cstdint>
#include <cerrno>

// ---------------------------------------------------------------------------
// Low-level helpers
// ---------------------------------------------------------------------------

// Receive exactly `length` bytes, looping over short reads.  Returns false on
// peer disconnect (recv == 0) or a fatal error.
bool AZURESocket::recvAll( void* buffer, size_t length ) {
  char* p = static_cast<char*>( buffer );
  size_t total = 0;
  while( total < length ) {
    ssize_t n = recv( clientSocket_, p + total, length - total, 0 );
    if( n > 0 ) {
      total += static_cast<size_t>( n );
    } else if( n == 0 ) {
      return false;            // peer closed the connection
    } else {
#ifndef _WIN32
      if( errno == EINTR ) continue;
#endif
      return false;            // fatal error
    }
  }
  return true;
}

// Send exactly `length` bytes, looping over short writes.
bool AZURESocket::sendAll( const void* buffer, size_t length ) {
  const char* p = static_cast<const char*>( buffer );
  size_t total = 0;
  while( total < length ) {
    ssize_t n = send( clientSocket_, p + total, length - total, 0 );
    if( n > 0 ) {
      total += static_cast<size_t>( n );
    } else {
#ifndef _WIN32
      if( n < 0 && errno == EINTR ) continue;
#endif
      return false;
    }
  }
  return true;
}

// ---------------------------------------------------------------------------
// Length-prefixed framing: [uint64 count][count * double]
// ---------------------------------------------------------------------------

bool AZURESocket::readMessage( vector_r& out ) {
  std::uint64_t count = 0;
  if( !recvAll( &count, sizeof( count ) ) ) return false;

  if( count > MAX_FRAME_DOUBLES ) {
    std::cerr << "AZURESocket: rejecting oversized frame (" << count
              << " doubles)." << std::endl;
    return false;
  }

  out.resize( count );
  if( count > 0 ) {
    if( !recvAll( out.data(), count * sizeof( double ) ) ) return false;
  }
  return true;
}

bool AZURESocket::writeMessage( const double* values, std::uint64_t count ) {
  if( !sendAll( &count, sizeof( count ) ) ) return false;
  if( count > 0 ) {
    if( !sendAll( values, count * sizeof( double ) ) ) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// Public send helpers
// ---------------------------------------------------------------------------

bool AZURESocket::sendPacket( const vector_r& response ) {
  return writeMessage( response.data(), response.size() );
}

bool AZURESocket::sendPacket( const std::vector<bool>& response ) {
  vector_r tmp( response.size() );
  for( size_t i = 0; i < response.size(); ++i ) tmp[i] = response[i] ? 1.0 : 0.0;
  return writeMessage( tmp.data(), tmp.size() );
}

// Strings are transmitted as one double per character so the client decodes
// them with the same framing as every other response.
bool AZURESocket::sendPacket( const std::string& response ) {
  vector_r tmp( response.size() );
  for( size_t i = 0; i < response.size(); ++i )
    tmp[i] = static_cast<double>( static_cast<unsigned char>( response[i] ) );
  return writeMessage( tmp.data(), tmp.size() );
}

// ---------------------------------------------------------------------------
// Server loop
// ---------------------------------------------------------------------------

bool AZURESocket::start() {
  // Create socket
  serverSocket_ = socket(AF_INET, SOCK_STREAM, 0);
  if (serverSocket_ == -1) {
    std::cerr << "Error creating socket." << std::endl;
    return false;
  }

  // Prepare the server address structure
  serverAddress_.sin_family = AF_INET;
  serverAddress_.sin_addr.s_addr = INADDR_ANY;
  serverAddress_.sin_port = htons(port_);

  int iSetOption = 1;
  setsockopt(serverSocket_, SOL_SOCKET, SO_REUSEADDR, (char*)&iSetOption, sizeof(iSetOption));

  // Bind the socket to a specific address and port
  if (bind(serverSocket_, (struct sockaddr *)&serverAddress_, sizeof(serverAddress_)) == -1) {
    std::cerr << "Error binding socket on port " << port_ << "." << std::endl;
    return false;
  }

  // Listen for incoming connections
  if (listen(serverSocket_, 1) == -1) {
    std::cerr << "Error listening on socket." << std::endl;
    return false;
  }

  socklen_t clientAddressSize = sizeof(clientAddress_);

  // Outer loop: accept clients.  Surviving a client disconnect and waiting for
  // a fresh connection makes the server robust to client restarts.
  while (true) {

    clientSocket_ = accept(serverSocket_, (struct sockaddr *)&clientAddress_, &clientAddressSize);
    if (clientSocket_ == -1) {
      std::cerr << "Error accepting connection." << std::endl;
      break;
    }

    // Inner loop: process requests until the client disconnects.
    vector_r request;
    while ( readMessage( request ) ) {
      if( request.empty() ) {
        // Malformed (no command word); ignore and keep listening.
        continue;
      }
      handle( request );
    }

    close(clientSocket_);
    clientSocket_ = -1;
  }

  close(serverSocket_);
  serverSocket_ = -1;

  return true;
}

// ---------------------------------------------------------------------------
// Request dispatch
//
// Frame layout (request):  [cmd, arg0, arg1, ...]
//   request[0]            -> command id
//   request[1 .. end]     -> arguments (index for getters, parameter vector
//                            for calculators)
// ---------------------------------------------------------------------------

void AZURESocket::handle( const vector_r& request ) {

  const int cmd = static_cast<int>( request[0] );

  // Arguments start at index 1.  Helpers below keep the handlers terse.
  const size_t nargs = request.size() - 1;

  // First argument as an integer index (0 if none supplied).
  const int idx = ( nargs >= 1 ) ? static_cast<int>( request[1] ) : 0;

  // Arguments as a parameter vector.
  vector_r params;
  if( nargs > 0 ) params.assign( request.begin() + 1, request.end() );

  switch( cmd ) {

    // Initialize
    case 0: {
      api_->Initialize( );
      sendPacket( std::vector<bool>{ true } );
      break;
    }

    // Calculate segments from params
    case 1: {
      double nSegments = (double)api_->UpdateSegments( params );
      sendPacket( vector_r{ nSegments } );
      break;
    }

    // Send the parameters
    case 2:
      sendPacket( api_->params_values( ) );
      break;

    // Send the parameter name
    case 3:
      sendPacket( api_->params_names( idx ) );
      break;

    // Send all the parameters
    case 4:
      sendPacket( api_->params_all( ) );
      break;

    // Calculate external capture
    case 5: {
      api_->CalculateExternalCapture( );
      sendPacket( std::vector<bool>{ true } );
      break;
    }

    // Get data energies
    case 6:
      sendPacket( api_->data_energies( idx ) );
      break;

    // Get data segments
    case 7:
      sendPacket( api_->data_segments( idx ) );
      break;

    // Get data segments errors
    case 8:
      sendPacket( api_->data_segments_errors( idx ) );
      break;

    // Update data
    case 9:
      sendPacket( vector_r{ (double)api_->UpdateData( ) } );
      break;

    // Set to data mode
    case 10: {
      api_->SetData( );
      sendPacket( std::vector<bool>{ true } );
      break;
    }

    // Set to extrapolation mode
    case 11: {
      api_->SetExtrap( );
      sendPacket( std::vector<bool>{ true } );
      break;
    }

    // Get calculated segments
    case 12:
      sendPacket( api_->calculated_segments( idx ) );
      break;

    // Get calculated energies
    case 13:
      sendPacket( api_->calculated_energies( idx ) );
      break;

    // Change radius
    case 14: {
      // params = [idx, radius]
      double radius = ( nargs >= 2 ) ? request[2] : 0.0;
      api_->SetRadius( idx, radius );
      sendPacket( std::vector<bool>{ true } );
      break;
    }

    // Send the norms
    case 15:
      api_->UpdateNorms( );
      sendPacket( api_->norms( ) );
      break;

    // Send the norms errors
    case 16:
      api_->UpdateNorms( );
      sendPacket( api_->norms_errors( ) );
      break;

    // Get calculated E1 segments
    case 17:
      sendPacket( api_->calculated_segments_e1( idx ) );
      break;

    // Get calculated E2 segments
    case 18:
      sendPacket( api_->calculated_segments_e2( idx ) );
      break;

    // Get data conversion
    case 19:
      sendPacket( api_->data_conv( idx ) );
      break;

    // Get calculated conversion
    case 20:
      sendPacket( api_->calculated_conv( idx ) );
      break;

    // Get data angles
    case 21:
      sendPacket( api_->data_angles( idx ) );
      break;

    // Get calculated angles
    case 22:
      sendPacket( api_->calculated_angles( idx ) );
      break;

    // Get RWA params
    case 23:
      sendPacket( api_->params_values_rwa( ) );
      break;

    // Calculate from RWA params
    case 24:
      sendPacket( vector_r{ (double)api_->UpdateSegmentsRWA( params ) } );
      break;

    // Transform RWA parameters
    case 25:
      sendPacket( api_->TransformRWAParameters( params ) );
      break;

    // Transform all RWA parameters
    case 26:
      sendPacket( api_->TransformAllRWAParameters( params ) );
      break;

    // Calculate chi2 from RWA parameters
    case 27:
      sendPacket( vector_r{ api_->CalculateChi2RWA( params ) } );
      break;

    // Calculate chi2 from physical parameters
    case 28:
      sendPacket( vector_r{ api_->CalculateChi2Physical( params ) } );
      break;

    // Calculate segments from all-RWA params
    case 29:
      sendPacket( vector_r{ (double)api_->UpdateSegmentsAllRWA( params ) } );
      break;

    // Get indices of normalization parameters
    case 30:
      sendPacket( api_->GetNormalizationIndices( ) );
      break;

    // Get indices of energy shift parameters
    case 31:
      sendPacket( api_->GetEnergyShiftIndices( ) );
      break;

    // Get indices of fixed parameters
    case 32:
      sendPacket( api_->params_fixed( ) );
      break;

    // Get data excitation energies
    case 33:
      sendPacket( api_->data_excitation_energies( idx ) );
      break;

    // Get calculated excitation energies
    case 34:
      sendPacket( api_->calculated_excitation_energies( idx ) );
      break;

    // Calculate log-likelihood from RWA params with per-segment error inflation
    case 35:
      sendPacket( vector_r{ api_->CalculateLnLRWA( params ) } );
      break;

    // Calculate log-likelihood from RWA params with a correlated per-segment
    // covariance matrix
    case 36:
      sendPacket( vector_r{ api_->CalculateLnLCovRWA( params ) } );
      break;

    // Get structured per-parameter metadata (level, J, L, S, ...)
    case 37:
      sendPacket( api_->GetParameterInfo( ) );
      break;

    default:
      std::cerr << "AZURESocket: unknown command " << cmd << "." << std::endl;
      // Reply with an empty frame so the client does not block forever.
      sendPacket( vector_r{} );
      break;
  }
}
