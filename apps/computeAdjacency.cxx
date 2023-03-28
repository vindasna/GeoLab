///////////////////////////////////////////////////////////////////////////////
//---------------------------- Libraries ------------------------------------//
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#include <ncurses.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include "computeAdjacency.h"
#include "ioWrapper.h"



///////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ///////////////
///////////////////////////////////////////////////////////////////////////////

int getFlagPosition( int argc, char* argv[], const std::string& flag )
{

  for ( int i = 0 ; i < argc ; i++ )
  {

    std::string arg = argv[ i ] ;
    if ( arg == flag )
    {

      return i ;

    }

  }
  return 0 ;

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  int index_bundle1, index_bundle2, index_thr, index_nbPoints, index_verbose,
                                                                    index_help ;

  std::vector<std::string> possibleFlags{ "-b1", "-b2", "-thr", "-nbPoints",
                                                                  "-v", "-h" } ;

  std::vector<bool> possibleFlagsNeedArgument{ true, true, true, true, true,
                                                                       false } ;

  if ( possibleFlagsNeedArgument.size() != possibleFlags.size() )
  {

    std::cout << "ProjectAtlas CODING ERROR : possibleFlags and "
              << "possibleFlagsNeedArgument must have the same size"
              << std::endl ;
    exit( 1 ) ;

  }



  for ( int k = 1 ; k < argc ; k++ )
  {

    std::string arg = argv[ k ] ;
    auto it = std::find( possibleFlags.begin(), possibleFlags.end(), arg ) ;
    if( it == possibleFlags.end() )
    {

      if ( k == 1 )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for ProjectAtlas command" << std::endl ;
        exit( 1 ) ;

      }

      std::string arg2 = argv[ k - 1 ] ;
      it = std::find( possibleFlags.begin(), possibleFlags.end(), arg2 ) ;
      if( it == possibleFlags.end() )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for ProjectAtlas command" << std::endl ;
        exit( 1 ) ;

      }
      int indexFound = it - possibleFlags.begin() ;
      if ( !possibleFlagsNeedArgument[ indexFound ] )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for ProjectAtlas command" << std::endl ;
        exit( 1 ) ;

      }

    }
    else
    {

      int indexFound = it - possibleFlags.begin() ;
      if ( possibleFlagsNeedArgument[ indexFound ] )
      {

        if ( k == argc - 1 )
        {

          std::cout << "ERROR in arguments : " << arg << " needs an input but "
                    << "none was given ( see -h )" << std::endl ;
          exit( 1 ) ;

        }

        std::string arg2 = argv[ k + 1 ] ;
        it = std::find( possibleFlags.begin(), possibleFlags.end(), arg2 ) ;
        if( it != possibleFlags.end() )
        {

          std::cout << "ERROR in arguments : " << arg << " needs an input but "
                    << "none was given ( see -h )" << std::endl ;
          exit( 1 ) ;

        }

      }

    }

  }

  index_bundle1 = getFlagPosition( argc, argv, "-b1" ) ;
  index_bundle2 = getFlagPosition( argc, argv, "-b2" ) ;
  index_thr = getFlagPosition( argc, argv, "-thr" ) ;
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to compute adjacency bundle 1 to bundle 2 : \n"
              << "-b1 : Input bundle from which to extract bundles \n"
              << "-b2 : Output directory to save the extracted bundles \n"
              << "-thr : threshold for adjacency in mm \n"
              << "-nbPoints : number of points per track \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_bundle1 )
  {

    std::cout << "-b1 argument required..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_bundle2 )
  {

    std::cout << "-b2 argument required..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_thr )
  {

    std::cout << "-thr argument required..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_nbPoints )
  {

    std::cout << "-nbPoints argument required..." << std::endl ;
    exit( 1 ) ;

  }


  /// Getting bundle 1
  char lastChar ;
  std::string inputBundle1Filename( argv[ index_bundle1 + 1 ] ) ;
  lastChar = inputBundle1Filename[ inputBundle1Filename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputBundle1Filename = inputBundle1Filename.substr( 0,
                                             inputBundle1Filename.size() - 1 ) ;

  }

  /// Getting bundle 2
  std::string inputBundle2Filename( argv[ index_bundle2 + 1 ] ) ;
  lastChar = inputBundle2Filename[ inputBundle2Filename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputBundle2Filename = inputBundle2Filename.substr( 0,
                                             inputBundle2Filename.size() - 1 ) ;

  }

  /// Getting threshold
  float thrDistance = std::stof( argv[ index_thr + 1 ] ) ;

  /// Getting number of points per track
  int nbPoints = std::stoi( argv[ index_nbPoints + 1 ] ) ;


  // Reading bundles
  if ( is_file( inputBundle1Filename ) && is_file( inputBundle2Filename ) )
  {

    BundlesData bundle1( inputBundle1Filename.c_str() ) ;
    BundlesData bundle2( inputBundle2Filename.c_str() ) ;

    AtlasBundles tmpAtlas ;
    float bundleAdjacency = tmpAtlas.bundlesAdjacency( bundle1.matrixTracks,
                                                      bundle2.matrixTracks,
                                                      bundle1.curves_count,
                                                      bundle2.curves_count,
                                                      nbPoints,
                                                      thrDistance,
                                                      0 ) ;
    std::cout << bundleAdjacency << std::endl ;

  }


  return 0 ;

}
