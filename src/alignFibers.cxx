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

#include "alignFibers.h"
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
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              std::vector<float> vector1,
                              std::vector<float> vector2,
                              std::vector<float> vector3,
                              int nbPoints,
                              std::string outputDirectory,
                              std::string fiberName )
{

  int fiberIndex = 0 ;

  std::vector<float> medialPointFiber( 3, 0 ) ;
  inputFiber.computeMedialPointFiberWithDistance(
                                                fiberIndex, medialPointFiber ) ;

  std::vector<float> matrixTracksfiber( 3 * nbPoints + 3 * 2 + 3 * 3, 0 ) ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      matrixTracksfiber[ 3 * point + i ] = inputFiber[ 3 * point + i ] ;

    }

  }

  matrixTracksfiber[ 3 * nbPoints + 3 * 0 + 0 ] = medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 0 + 1 ] = medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 0 + 2 ] = medialPointFiber[ 2 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 1 + 0 ] = vector1[ 0 ]
                                                       + medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 1 + 1 ] = vector1[ 1 ]
                                                       + medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 1 + 2 ] = vector1[ 2 ]
                                                       + medialPointFiber[ 2 ] ;

  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 0 ] = medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 1 ] = medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 2 ] = medialPointFiber[ 2 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 0 ] = vector2[ 0 ]
                                                       + medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 1 ] = vector2[ 1 ]
                                                       + medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 2 ] = vector2[ 2 ]
                                                       + medialPointFiber[ 2 ] ;

  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 0 ] = medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 1 ] = medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 2 ] = medialPointFiber[ 2 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 0 ] = vector3[ 0 ]
                                                       + medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 1 ] = vector3[ 1 ]
                                                       + medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 2 ] = vector3[ 2 ]
                                                       + medialPointFiber[ 2 ] ;

  std::vector<int> pointsPerTrackFiber = { nbPoints, 2, 2, 2 } ;
  int curves_countFiber = 4 ;

  BundlesData referenceFiber = inputFiber ;
  referenceFiber.matrixTracks = matrixTracksfiber ;
  referenceFiber.pointsPerTrack = pointsPerTrackFiber ;
  referenceFiber.curves_count = curves_countFiber ;

  inputFiberInfo.curves_count = curves_countFiber ;
  BundlesMinf referenceFiberInfo( inputFiberInfo ) ;
  std::string referenceFiberFilename = outputDirectory + fiberName + format ;
  referenceFiber.write( referenceFiberFilename.c_str(), referenceFiberInfo ) ;


}

// -------------------------------------------------------------------------- //

void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              std::vector<float>& vector1,
                              std::vector<float>& vector2,
                              int nbPoints,
                              std::string outputDirectory,
                              std::string fiberName )
{

  int fiberIndex = 0 ;

  std::vector<float> medialPointFiber( 3, 0 ) ;
  inputFiber.computeMedialPointFiberWithDistance(
                                               fiberIndex, medialPointFiber ) ;
  // inputFiber.computeGravityCenterFiber( fiberIndex, medialPointFiber ) ;

  std::vector<float> matrixTracksfiber( 3 * nbPoints + 3 * 2 + 3 * 2, 0 ) ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      matrixTracksfiber[ 3 * point + i ] = inputFiber[ 3 * point + i ] ;

    }

  }

  matrixTracksfiber[ 3 * nbPoints + 3 * 0 + 0 ] = medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 0 + 1 ] = medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 0 + 2 ] = medialPointFiber[ 2 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 1 + 0 ] = vector1[ 0 ] +
                                                         medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 1 + 1 ] = vector1[ 1 ] +
                                                         medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 1 + 2 ] = vector1[ 2 ] +
                                                         medialPointFiber[ 2 ] ;

  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 0 ] = medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 1 ] = medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 2 + 2 ] = medialPointFiber[ 2 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 0 ] = vector2[ 0 ] +
                                                         medialPointFiber[ 0 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 1 ] = vector2[ 1 ] +
                                                         medialPointFiber[ 1 ] ;
  matrixTracksfiber[ 3 * nbPoints + 3 * 3 + 2 ] = vector2[ 2 ] +
                                                         medialPointFiber[ 2 ] ;

  std::vector<int> pointsPerTrackFiber = { nbPoints, 2, 2 } ;
  int curves_countFiber = 3 ;

  BundlesData referenceFiber = inputFiber ;
  referenceFiber.matrixTracks = matrixTracksfiber ;
  referenceFiber.pointsPerTrack = pointsPerTrackFiber ;
  referenceFiber.curves_count = curves_countFiber ;


  inputFiberInfo.curves_count = curves_countFiber ;
  BundlesMinf referenceFiberInfo( inputFiberInfo ) ;
  std::string referenceFiberFilename = outputDirectory + fiberName + format ;

  referenceFiber.write( referenceFiberFilename.c_str(), referenceFiberInfo ) ;


}

// -------------------------------------------------------------------------- //

void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              int nbPoints,
                              std::string outputDirectory,
                              std::string fiberName )
{

  int fiberIndex = 0 ;

  std::vector<float> medialPointFiber( 3, 0 ) ;
  inputFiber.computeMedialPointFiberWithDistance(
                                               fiberIndex, medialPointFiber ) ;
  std::vector<float> normalVector( 3, 0 ) ;
  inputFiber.computeNormalVectorFiberTractogram( fiberIndex, normalVector ) ;
  std::vector<float> directionVector( 3, 0 ) ;
  inputFiber.computeDirectionVectorFiberTractogram( fiberIndex, normalVector,
                                                             directionVector ) ;

  computeFiberWithVectors( inputFiber, inputFiberInfo, normalVector,
                                    directionVector ,nbPoints, outputDirectory,
                                                                   fiberName ) ;


}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
float scalarProduct( std::vector<float>& vector1, std::vector<float>& vector2 )
{

  float scalarProductResult = 0 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    scalarProductResult += vector1[ i ] * vector2[ i ] ;

  }

  return scalarProductResult ;

}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////::::::////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void crossProduct( const std::vector<float>& vector1,
                   const std::vector<float>& vector2,
                   std::vector<float>& result )
{

  result[ 0 ] = vector1[ 1 ] * vector2[ 2 ] - vector1[ 2 ] * vector2[ 1 ] ;
  result[ 1 ] = vector1[ 2 ] * vector2[ 0 ] - vector1[ 0 ] * vector2[ 2 ] ;
  result[ 2 ] = vector1[ 0 ] * vector2[ 1 ] - vector1[ 1 ] * vector2[ 0 ] ;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  int index_inpu_fiber_1, index_inpu_fiber_2, index_output, index_verbose,
                                                                    index_help ;

  index_inpu_fiber_1 = getFlagPosition( argc, argv, "-i1") ;
  index_inpu_fiber_2 = getFlagPosition( argc, argv, "-i2") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-i1 : Input fiber 1 \n"
              << "-i2 : Input fiber 2 \n"
              << "-o : Output directory to save the extracted bundles \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_inpu_fiber_1 )
  {

    std::cout << "-i1 argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_inpu_fiber_2 )
  {

    std::cout << "-i2 argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output )
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

  std::string inputFilename1( argv[ index_inpu_fiber_1 + 1 ] ) ;
  char lastChar = inputFilename1[ inputFilename1.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputFilename1 = inputFilename1.substr( 0, inputFilename1.size() - 1 ) ;

  }

  std::string inputFilename2( argv[ index_inpu_fiber_2 + 1 ] ) ;
  lastChar = inputFilename2[ inputFilename2.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputFilename2 = inputFilename2.substr( 0, inputFilename2.size() - 1 ) ;

  }


  std::string outputDirectory( argv[ index_output + 1 ] ) ;
  lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectory = outputDirectory + "/" ;

  }

  if ( !( is_dir( outputDirectory ) ) )
  {

    mkdir( outputDirectory ) ;

  }


  //   xxxxxxxxxxxxxxxxxxxxxx Reading input fibers 1 xxxxxxxxxxxxxxxxxxxxxx   //

  std::string inputBundlesFilename1 ;
  std::string inputBundlesDataFilename1 ;
  bool isBundlesFormat = false ;
  bool isTrk = false ;
  bool isTck = false ;

  if ( endswith( inputFilename1, ".bundlesdata" ) )
  {

    inputBundlesDataFilename1 = inputFilename1 ;

    inputBundlesFilename1 = replaceExtension( inputFilename1, ".bundles" ) ;

    isBundlesFormat = true ;
    isTrk = false ;
    isTck = false ;

    format = ".bundles" ;

  }
  else if ( endswith( inputFilename1, ".bundles" ) )
  {

    inputBundlesFilename1 = inputFilename1 ;

    inputBundlesDataFilename1 = replaceExtension( inputFilename1,
                                                              ".bundlesdata" ) ;

    isBundlesFormat = true ;
    isTrk = false ;
    isTck = false ;

    format = ".bundles" ;

  }
  else if ( endswith( inputFilename1, ".trk" ) )
  {

    inputBundlesDataFilename1 = inputFilename1 ;

    inputBundlesFilename1 = replaceExtension( inputFilename1, ".minf" ) ;

    isBundlesFormat = true ;
    isTrk = true ;
    isTck = false ;

    format = ".trk" ;

  }
  else if ( endswith( inputFilename1, ".tck" ) )
  {

    inputBundlesDataFilename1 = inputFilename1 ;

    inputBundlesFilename1 = replaceExtension( inputFilename1, ".minf" ) ;

    isBundlesFormat = false ;
    isTrk = false ;
    isTck = true ;

    format = ".tck" ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles/.trk/.tck"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesData inputFiber1( inputBundlesDataFilename1.c_str() ) ;
  BundlesMinf inputFiberInfo1( inputBundlesFilename1.c_str() ) ;


  //   xxxxxxxxxxxxxxxxxxxxxx Reading input fibers 2 xxxxxxxxxxxxxxxxxxxxxx   //

  std::string inputBundlesFilename2 ;
  std::string inputBundlesDataFilename2 ;

  if ( endswith( inputFilename2, ".bundlesdata" ) )
  {

    if ( !isBundlesFormat )
    {

      std::cout << "ERROR : both inputs -i1 and -i2 must have the same format"
                << std::endl ;

      exit( 1 ) ;

    }

    inputBundlesDataFilename2 = inputFilename2 ;
    inputBundlesFilename2 = replaceExtension( inputFilename2, ".bundles" ) ;


  }
  else if ( endswith( inputFilename2, ".bundles" ) )
  {

    if ( !isBundlesFormat )
    {

      std::cout << "ERROR : both inputs -i1 and -i2 must have the same format"
                << std::endl ;

      exit( 1 ) ;

    }

    inputBundlesFilename2 = inputFilename2 ;

    inputBundlesDataFilename2 = replaceExtension( inputFilename2,
                                                              ".bundlesdata" ) ;


  }
  else if ( endswith( inputFilename2, ".trk" ) )
  {

    if ( !isTrk )
    {

      std::cout << "ERROR : both inputs -i1 and -i2 must have the same format"
                << std::endl ;

      exit( 1 ) ;

    }

    inputBundlesDataFilename2 = inputFilename2 ;
    inputBundlesFilename2 = replaceExtension( inputFilename2, ".minf" ) ;

  }
  else if ( endswith( inputFilename2, ".tck" ) )
  {

    if ( !isTck )
    {

      std::cout << "ERROR : both inputs -i1 and -i2 must have the same format"
                << std::endl ;

      exit( 1 ) ;

    }

    inputBundlesDataFilename2 = inputFilename2 ;
    inputBundlesFilename2 = replaceExtension( inputFilename2, ".minf" ) ;


  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles/.trk/.tck"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesData inputFiber2( inputBundlesDataFilename2.c_str() ) ;
  BundlesMinf inputFiberInfo2( inputBundlesFilename2.c_str() ) ;



  // Sanity checks //
  int nbPoints1 = inputFiber1.pointsPerTrack[ 0 ] ;
  int nbPoints2 = inputFiber2.pointsPerTrack[ 0 ] ;
  if ( nbPoints1 != nbPoints2 )
  {

    std::cout << "ERROR :  the number of points in each fibers has to be the "
              << "same, got " << nbPoints1 << " for fiber 1 and " << nbPoints2
              << "for fiber 2 " << std::endl ;
    exit( 1 ) ;

  }
  int nbPoints = nbPoints1 ;

  int nbCurves1 = inputFiber1.curves_count ;
  if ( nbCurves1 != 1 )
  {

    std::cout << "ERROR : Only single fibers are supported, got " << nbCurves1
              << " for input 1 " << std::endl ;
    exit( 1 ) ;

  }
  int nbCurves2 = inputFiber2.curves_count ;
  if ( nbCurves2 != 1 )
  {

    std::cout << "ERROR : Only single fibers are supported, got " << nbCurves2
              << " for input 2 " << std::endl ;
    exit( 1 ) ;

  }


  //----------------------------- Aligning fibers ----------------------------//
  // std::vector<float> referenceFiber = inputFiber1.matrixTracks ;
  // std::vector<float> movingFiber = inputFiber2.matrixTracks ;
  //
  //
  // std::vector<float> normalVector1( 3, 0 ) ;
  // inputFiber1.computeNormalVectorFiberTractogram( referenceFiber,
  //                                                              normalVector1 ) ;
  // std::vector<float> directionVector1( 3, 0 ) ;
  // inputFiber1.computeDirectionVectorFiberTractogram( referenceFiber,
  //                                            normalVector1, directionVector1 ) ;
  //

  // std::cout << "Registering fibers ... " << std::endl ;
  // std::vector<float> samePlaneAndDirectionTranslatedFiber2( 3 * nbPoints, 0 ) ;
  // std::vector<float> newNormalVectorFiber2( 3, 0 ) ;
  // inputFiber1.registerFiber( referenceFiber,
  //                            movingFiber,
  //                            nbPoints,
  //                            samePlaneAndDirectionTranslatedFiber2,
  //                            newNormalVectorFiber2,
  //                            verbose ) ;
  // std::cout << "Done" << std::endl ;
  //
  //
  //
  // //----------------------------- Saving results -----------------------------//
  // std::string fiberName1 = "referenceFiber" ;
  // computeFiberWithVectors( inputFiber1,
  //                          inputFiberInfo1,
  //                          nbPoints,
  //                          outputDirectory,
  //                          fiberName1 ) ;
  //
  // BundlesData movedFiber( samePlaneAndDirectionTranslatedFiber2,
  //                               inputFiber2.pointsPerTrack,
  //                               inputFiber2.curves_count ) ;
  // BundlesMinf movedFiberInfo( inputFiberInfo2 ) ;
  // std::string fiberName2 = "movedFiber" ;
  // computeFiberWithVectors( movedFiber,
  //                          movedFiberInfo,
  //                          nbPoints,
  //                          outputDirectory,
  //                          fiberName2 ) ;
  //
  // std::vector<float> normalVectorRegistered( 3, 0 ) ;
  // movedFiber.computeNormalVectorFiberTractogram( 0, normalVectorRegistered ) ;
  // std::vector<float> directionVectorRegistered( 3, 0 ) ;
  // movedFiber.computeDirectionVectorFiberTractogram( 0, normalVectorRegistered,
  //                                                  directionVectorRegistered ) ;
  //
  // float angleBetweenPlanes = movedFiber.computeAngleBetweenPlanes(
  //                                      normalVector1, normalVectorRegistered ) ;
  // float angleBetweenDirections = movedFiber.computeAngleBetweenDirections(
  //                                directionVector1, directionVectorRegistered ) ;
  //
  // std::cout << "Angle between planes : " << angleBetweenPlanes << std::endl ;
  // std::cout << "Angle between directions : " << angleBetweenDirections <<
  //                                                                    std::endl ;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<float> fiber1 = inputFiber1.matrixTracks ;
  std::vector<float> fiber2 = inputFiber2.matrixTracks ;

  std::vector<float> normalVectorFiber1( 3, 0 ) ;
  inputFiber1.computeNormalVectorFiberTractogram( fiber1,
                                                  normalVectorFiber1 ) ;
  std::vector<float> medialPointFiber1( 3, 0 ) ;
  inputFiber1.computeMedialPointFiberWithDistance( fiber1,
                                                   medialPointFiber1 ) ;
  std::vector<float> fiber1XY( 3 * nbPoints, 0 ) ;
  std::vector<float> newNormalVectorFiber1( 3, 0 ) ;
  inputFiber1.registerFiberToPlaneXYAndDirectionX( fiber1,
                                                   normalVectorFiber1,
                                                   medialPointFiber1,
                                                   nbPoints,
                                                   fiber1XY,
                                                   newNormalVectorFiber1 ) ;


  std::vector<float> normalVectorFiber2( 3, 0 ) ;
  inputFiber1.computeNormalVectorFiberTractogram( fiber2,
                                                  normalVectorFiber2 ) ;
  std::vector<float> medialPointFiber2( 3, 0 ) ;
  inputFiber1.computeMedialPointFiberWithDistance( fiber2,
                                                   medialPointFiber2 ) ;
  std::vector<float> fiber2XY( 3 * nbPoints, 0 ) ;
  std::vector<float> newNormalVectorFiber2( 3, 0 ) ;
  inputFiber1.registerFiberToPlaneXYAndDirectionX( fiber2,
                                                   normalVectorFiber2,
                                                   medialPointFiber2,
                                                   nbPoints,
                                                   fiber2XY,
                                                   newNormalVectorFiber2 ) ;
  //----------------------------- Saving results -----------------------------//
  BundlesData movedFiber1 = inputFiber1 ;
  movedFiber1.matrixTracks = fiber1XY ;

  BundlesMinf movedFiber1Info( inputFiberInfo1 ) ;
  std::vector<float> directionVectorFiber1XY( 3, 0 ) ;
  inputFiber1.computeDirectionVectorFiberTractogram( fiber1XY,
                                                     newNormalVectorFiber1,
                                                     directionVectorFiber1XY ) ;
  std::string fiberName1 = "fiber1" ;
  computeFiberWithVectors( movedFiber1,
                           movedFiber1Info,
                           newNormalVectorFiber1,
                           directionVectorFiber1XY,
                           nbPoints,
                           outputDirectory,
                           fiberName1 ) ;


  BundlesData movedFiber2 = inputFiber2 ;
  movedFiber2.matrixTracks = fiber2XY ;
  BundlesMinf movedFiber2Info( inputFiberInfo2 ) ;
  std::vector<float> directionVectorFiber2XY( 3, 0 ) ;
  inputFiber1.computeDirectionVectorFiberTractogram( fiber2XY,
                                                     newNormalVectorFiber2,
                                                     directionVectorFiber2XY ) ;
  std::string fiberName2 = "fiber2" ;
  computeFiberWithVectors( movedFiber2,
                           movedFiber2Info,
                           newNormalVectorFiber2,
                           directionVectorFiber2XY,
                           nbPoints,
                           outputDirectory,
                           fiberName2 ) ;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  return 0 ;

}
