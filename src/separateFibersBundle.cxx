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

#include "separateFibersBundle.h"


///////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ///////////////
///////////////////////////////////////////////////////////////////////////////

int getFlagPosition( int argc, char* argv[], const std::string& flag )
{

  for ( int i = 0 ; i < argc ; i++ )
  {

    std::string arg = argv[ i ] ;
    // if ( 0 == arg.find( flag ) )
    if ( arg == flag )
    {

      return i ;

    }

  }
  return 0 ;

}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  int index_input_bundle, index_output_bundle, index_verbose, index_help ;

  index_input_bundle = getFlagPosition( argc, argv, "-i") ;
  index_output_bundle = getFlagPosition( argc, argv, "-o") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-i : Input bundle \n"
              << "-o : Output directory \n"
              << "[-v] : Set verbosity level \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input_bundle )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output_bundle )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_verbose )
  {
    if ( argv[ index_verbose + 1 ] )
    {

      verbose = std::stoi( argv[ index_verbose + 1 ] ) ;

    }
    else
    {

      verbose = 1 ;

    }

  }

  /////////////////////////////////////////////////////////////////////////////

  std::string inputFilename( argv[ index_input_bundle + 1 ] ) ;
  char lastChar = inputFilename[ inputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputFilename = inputFilename.substr( 0, inputFilename.size() - 1 ) ;

  }

  std::string outputDirectory( argv[ index_output_bundle + 1 ] ) ;
  lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectory = outputDirectory + '/' ;

  }


  //   xxxxxxxxxxxxxxxxxxxxxxx Reading input bundle xxxxxxxxxxxxxxxxxxxxxxx   //

  std::string inputBundlesFilename ;
  std::string inputBundlesDataFilename ;
  bool isBundlesFormat = false ;

  if ( inputFilename.find( ".bundlesdata" ) != std::string::npos )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = inputFilename ;
    size_t index = inputFilename.find( ".bundlesdata" ) ;
    inputBundlesFilename.replace( index, 12, ".bundles") ;

    isBundlesFormat = true ;

  }
  else if ( inputFilename.find( ".bundles" ) != std::string::npos )
  {

    inputBundlesFilename = inputFilename ;

    inputBundlesDataFilename = inputFilename ;
    size_t index = inputFilename.find( ".bundles" ) ;
    inputBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

    isBundlesFormat = true ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesDataFormat inputBundlesData ;
  BundlesFormat inputBundleInfo ;



  if ( isBundlesFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << inputBundlesFilename << std::endl ;

    }

    inputBundleInfo.bundlesReading( inputBundlesFilename.c_str(),
                                        verbose ) ;


    if ( verbose )
    {

      std::cout << "Reading : " << inputBundlesDataFilename << std::endl ;

    }
    inputBundlesData.bundlesdataReading( inputBundlesDataFilename.c_str(),
                                        inputBundlesFilename.c_str(),
                                        verbose ) ;

  }

  // Sanity checks //
  int nbPoints = inputBundlesData.pointsPerTrack[ 0 ] ;

  int nbCurves = inputBundlesData.curves_count ;

  for( int fiber = 1 ; fiber < nbCurves ; fiber++ )
  {

    if ( nbPoints != inputBundlesData.pointsPerTrack[ fiber ] )
    {

      std::cout << "ERROR : All the fibers in the bundle must have the same "
                << "number of points..." << std::endl ;
      exit( 1 ) ;

    }

  }


  // Extract fibers
  std::vector<float>& tractogramFibers = inputBundlesData.matrixTracks ;
  std::vector<BundlesDataFormat> separatedFibers ;
  std::vector<BundlesFormat> separatedFibersInfo ;
  for ( int fiberIndex = 0 ; fiberIndex < nbCurves ; fiberIndex++ )
  {

    std::vector<float> fiber( 3 * nbPoints, 0 ) ;
    inputBundlesData.getFiberFromTractogram( tractogramFibers, fiberIndex,
                                                             nbPoints, fiber ) ;

    std::vector<int32_t> fiberPointsPerTrack( 1, nbPoints ) ;
    BundlesDataFormat fiberBundlesDataFormat( fiber, fiberPointsPerTrack, 1 ) ;
    separatedFibers.push_back( fiberBundlesDataFormat ) ;

    BundlesFormat fiberBundlesFormat( inputBundleInfo ) ;
    fiberBundlesFormat.curves_count = 1 ;
    separatedFibersInfo.push_back( fiberBundlesFormat ) ;

  }


  // Saving results
  for ( int fiberIndex = 0 ; fiberIndex < nbCurves ; fiberIndex++ )
  {

    std::string outputBundlesFilename = outputDirectory + "fiber" +
                                     std::to_string( fiberIndex ) + ".bundles" ;
    std::string outputBundlesDataFilename = outputDirectory + "fiber" +
                                 std::to_string( fiberIndex ) + ".bundlesdata" ;

    if ( verbose > 1 )
    {

      std::cout << "Saving file : " << outputBundlesFilename << std::endl ;

    }
    separatedFibersInfo[ fiberIndex ].bundlesWriting(
                                      outputBundlesFilename.c_str(), verbose ) ;

    if ( verbose > 1 )
    {

      std::cout << "Saving file : " << outputBundlesDataFilename << std::endl ;

    }
    separatedFibers[ fiberIndex ].bundlesdataWriting(
                                  outputBundlesDataFilename.c_str(), verbose ) ;


  }


  return 0 ;

}
