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

#include "alignFibersInBundle.h"
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
              << "-o : Output aligned bundle \n"
              << "[-v] : Set verbosity level at 1 \n"
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

  std::string outputFilename( argv[ index_output_bundle + 1 ] ) ;
  lastChar = outputFilename[ outputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

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



  //---------------------------- Aligning fibers -----------------------------//
  BundlesData outputBundlesData( inputBundlesData ) ;
  outputBundlesData.matrixTracks.resize( 3 * nbPoints * nbCurves +
                                                         3 * 2 * nbCurves, 0 ) ;
  outputBundlesData.pointsPerTrack.resize( 2 * nbCurves ) ;
  outputBundlesData.curves_count = 2 * nbCurves ;

  BundlesMinf outputBundleInfo( inputBundleInfo ) ;
  outputBundleInfo.curves_count = 2 * nbCurves ;

  std::vector<float> referenceFiber( 3 * nbPoints, 0 ) ;
  inputBundlesData.getFiberFromTractogram( inputBundlesData.matrixTracks,
                                                 0, nbPoints, referenceFiber ) ;

  std::vector<float> normalVectorReference( 3, 0 ) ;
  inputBundlesData.computeNormalVectorFiberTractogram( referenceFiber,
                                                       normalVectorReference ) ;

  std::vector<float> medialPointReferenceFiber( 3 , 0 ) ;
  int pointReference = inputBundlesData.computeMedialPointFiberWithDistance(
                                                   referenceFiber,
                                                   medialPointReferenceFiber ) ;

  std::vector<float> directionVectorReference( 3, 0 ) ;
  inputBundlesData.computeDirectionVectorFiberTractogram( referenceFiber,
                             normalVectorReference, directionVectorReference ) ;

  std::copy( referenceFiber.begin(), referenceFiber.end(),
             outputBundlesData.matrixTracks.begin() + 3 * nbPoints * 0 +
                                                                   3 * 2 * 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    directionVectorReference[ i ] += medialPointReferenceFiber[ i ] ;

  }
  std::copy( directionVectorReference.begin(),
             directionVectorReference.end(),
             outputBundlesData.matrixTracks.begin() + 3 * nbPoints *
                                                                   ( 0 + 1 ) ) ;

  std::copy( medialPointReferenceFiber.begin(),
             medialPointReferenceFiber.end(),
             outputBundlesData.matrixTracks.begin() + 3 * nbPoints *
                                                         ( 0 + 1 ) + 3 ) ;

  outputBundlesData.pointsPerTrack[ 2 * 0 ] = nbPoints ;
  outputBundlesData.pointsPerTrack[ 2 * 0 + 1 ] = 2 ;


  for ( int fiber = 1 ; fiber < nbCurves ; fiber++ )
  {

    std::vector<float> movingFiber( 3 * nbPoints, 0 ) ;
    inputBundlesData.getFiberFromTractogram( inputBundlesData.matrixTracks,
                                             fiber, nbPoints, movingFiber ) ;
    std::vector<float> registeredFiber( 3 * nbPoints, 0 ) ;
    std::vector<float> normalVectorRegistered( 3, 0 ) ;
    inputBundlesData.registerFiber( referenceFiber,
                                    movingFiber,
                                    nbPoints,
                                    registeredFiber,
                                    normalVectorRegistered ) ;

    std::vector<float> directionVectorRegistered( 3, 0 ) ;
    inputBundlesData.computeDirectionVectorFiberTractogram( registeredFiber,
                           normalVectorRegistered, directionVectorRegistered ) ;
    // inputBundlesData.computeDirectionVectorFiberTractogram(
    //                                                registeredFiber,
    //                                                directionVectorRegistered ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      directionVectorRegistered[ i ] += medialPointReferenceFiber[ i ] ;

    }

    float angleBetweenPlanes = inputBundlesData.computeAngleBetweenPlanes(
                               normalVectorReference, normalVectorRegistered ) ;
    if ( angleBetweenPlanes > 5 )
    {

      std::cout << "Problem aligning planes, got angle between reference and "
                << "moved of "  << angleBetweenPlanes << std::endl ;
      exit( 1 ) ;

    }


    float angleBetweenDirections =
                                 inputBundlesData.computeAngleBetweenDirections(
                         directionVectorReference, directionVectorRegistered ) ;
    if ( angleBetweenDirections > 5 )
    {

      std::cout << "Problem aligning directions, got angle between reference "
                << "and moved of "  << angleBetweenDirections << std::endl ;
      exit( 1 ) ;

    }


    std::copy( registeredFiber.begin(), registeredFiber.end(),
               outputBundlesData.matrixTracks.begin() + 3 * nbPoints * fiber +
                                                               3 * 2 * fiber ) ;

    std::copy( directionVectorRegistered.begin(),
               directionVectorRegistered.end(),
               outputBundlesData.matrixTracks.begin() + 3 * nbPoints *
                                               ( fiber + 1 ) + 3 * 2 * fiber ) ;

    std::copy( medialPointReferenceFiber.begin(),
               medialPointReferenceFiber.end(),
               outputBundlesData.matrixTracks.begin() + 3 * nbPoints *
                                           ( fiber + 1 ) + 3 * 2 * fiber + 3 ) ;

    outputBundlesData.pointsPerTrack[ 2 * fiber ] = nbPoints ;
    outputBundlesData.pointsPerTrack[ 2 * fiber + 1 ] = 2 ;


  }



  // Saving results
  std::string outputBundlesFilename ;
  std::string outputBundlesDataFilename ;

  if ( endswith( outputFilename, ".bundlesdata" ) )
  {

    if ( format != ".bundles" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = replaceExtension( outputFilename, ".bundles" ) ;

  }
  else if ( endswith( outputFilename, ".bundles" ) )
  {

    if ( format != ".bundles" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesFilename = outputFilename ;

    outputBundlesDataFilename = replaceExtension( outputFilename,
                                                              ".bundlesdata" ) ;

  }
  else if ( endswith( outputFilename, ".trk" ) )
  {

    if ( format != ".trk" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = replaceExtension( outputFilename, ".minf" ) ;

  }
  else if ( endswith( outputFilename, ".tck" ) )
  {

    if ( format != ".tck" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = replaceExtension( outputFilename, ".minf" ) ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }

  outputBundlesData.write( outputBundlesDataFilename.c_str(),
                                                            outputBundleInfo ) ;

  return 0 ;

}
