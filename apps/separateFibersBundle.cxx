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
#include "ioWrapper.h"


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

  if ( !( is_dir( outputDirectory ) ) )
  {

    mkdir( outputDirectory ) ;

  }


  //   xxxxxxxxxxxxxxxxxxxxxxx Reading input bundle xxxxxxxxxxxxxxxxxxxxxxx   //

  std::string inputBundlesFilename ;
  std::string inputBundlesDataFilename ;
  bool isBundlesFormat = false ;
  bool isTrk = false ;
  bool isTck = false ;

  if ( endswith( inputFilename, ".bundlesdata" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".bundles" ) ;

    isBundlesFormat = true ;
    format = ".bundles" ;

  }
  else if ( endswith( inputFilename, ".bundles" ) )
  {

    inputBundlesFilename = inputFilename ;

    inputBundlesDataFilename = replaceExtension( inputFilename,
                                                              ".bundlesdata" ) ;

    isBundlesFormat = true ;
    format = ".bundles" ;

  }
  else if ( endswith( inputFilename, ".trk" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".minf" ) ;

    isTrk = true ;
    format = ".trk" ;

  }
  else if ( endswith( inputFilename, ".tck" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".minf" ) ;

    isTck = true ;
    format = ".tck" ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesData inputBundlesData( inputBundlesDataFilename.c_str() ) ;
  BundlesMinf inputBundleInfo( inputBundlesFilename.c_str() ) ;


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
  std::vector<BundlesData> separatedFibers ;
  std::vector<BundlesMinf> separatedFibersInfo ;
  for ( int fiberIndex = 0 ; fiberIndex < nbCurves ; fiberIndex++ )
  {

    std::vector<float> fiber( 3 * nbPoints, 0 ) ;
    inputBundlesData.getFiberFromTractogram( tractogramFibers, fiberIndex,
                                                             nbPoints, fiber ) ;

    std::vector<int32_t> fiberPointsPerTrack( 1, nbPoints ) ;
    BundlesData fiberBundlesDataFormat = inputBundlesData ;
    fiberBundlesDataFormat.matrixTracks = fiber ;
    fiberBundlesDataFormat.pointsPerTrack = fiberPointsPerTrack ;
    fiberBundlesDataFormat.curves_count = 1 ;

    if ( format == ".trk" )
    {

      fiberBundlesDataFormat.tracksScalars.resize( 1 ) ;
      fiberBundlesDataFormat.tracksProperties.resize( 1 ) ;


    }

    separatedFibers.push_back( fiberBundlesDataFormat ) ;

    BundlesMinf fiberBundlesFormat( inputBundleInfo ) ;
    fiberBundlesFormat.curves_count = 1 ;
    separatedFibersInfo.push_back( fiberBundlesFormat ) ;

  }


  // Saving results
  for ( int fiberIndex = 0 ; fiberIndex < nbCurves ; fiberIndex++ )
  {

    std::stringstream outputBundlesFilenameOss ;
    outputBundlesFilenameOss << outputDirectory << "fiber" << fiberIndex
                                                                     << format ;
    std::string outputBundlesFilename = outputBundlesFilenameOss.str() ;


    if ( verbose > 1 )
    {

      std::cout << "Saving file : " << outputBundlesFilename << std::endl ;

    }
    separatedFibers[ fiberIndex ].write( outputBundlesFilename.c_str(),
                                           separatedFibersInfo[ fiberIndex ] ) ;


  }


  return 0 ;

}
