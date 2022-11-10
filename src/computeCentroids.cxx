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

#include "computeCentroids.h"

////////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ////////////////
////////////////////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////////////////////
///////////////// Function to compute medial point fibers //////////////////////
////////////////////////////////////////////////////////////////////////////////
void computeMedialPointsFibersTractogram(
             BundlesDataFormat& inputTractogram,
             std::vector< std::vector< float > >& medialPointsFibersTractogram )
{

  int nbCurves = inputTractogram.curves_count ;

  if ( medialPointsFibersTractogram.size() != nbCurves)
  {

    std::cout << "ERROR in computeMedialPointsFibersTractogram : "
              << " the lenght of medialPointsFibersTractogram has to be the "
              << "same as inputTractogram.curves_count" << std::endl ;
    exit( 1 ) ;

  }

  #pragma omp parallel for
  for ( int fiber = 0 ; fiber < nbCurves ; fiber++ )
  {

    std::vector<float>& medialPoint =  medialPointsFibersTractogram[ fiber ] ;
    if ( medialPoint.size() != 3 )
    {

      std::cout << "ERROR in computeMedialPointsFibersTractogram : the lenght "
                << " of the elemenst of medialPointsFibersTractogram has to "
                << "be 3" << std::endl ;
      exit( 1 ) ;

    }

    inputTractogram.computeMedialPointFiberWithDistance( fiber, medialPoint ) ;

  }


}


////////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ////////////////
////////////////////////////////////////////////////////////////////////////////
double getAverageMinDistanceBetweenFibers(
       BundlesDataFormat& inputTractogram,
       const std::vector< std::vector< float > >& medialPointsFibersTractogram )
{

  int nbCurves = inputTractogram.curves_count ;

  double averageMinDistance = 0 ;

  #pragma omp parallel for
  for ( int fiber1 = 0 ; fiber1 < nbCurves ; fiber1++ )
  {

    std::vector<float> medialPointFiberTractogram1 =
                                        medialPointsFibersTractogram[ fiber1 ] ;

    double minDistance = 10000 ; // In mm

    for ( int fiber2 = 0 ; fiber2 < nbCurves ; fiber2++ )
    {

      if ( fiber2 != fiber1 )
      {

        std::vector<float> medialPointFiberTractogram2 =
                                        medialPointsFibersTractogram[ fiber2 ] ;

        float tmpDistance = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          tmpDistance += pow( medialPointFiberTractogram2[ i ] -
                                         medialPointFiberTractogram1[ i ], 2 ) ;

        }
        tmpDistance = sqrt( tmpDistance ) ;

        if ( tmpDistance < minDistance )
        {

          minDistance = tmpDistance ;

        }

      }

    }

    averageMinDistance += minDistance ;

  }

  averageMinDistance /= nbCurves ;

  return averageMinDistance ;

}


////////////////////////////////////////////////////////////////////////////////
////////// Funtion to compute distance between centers of fibers ///////////////
////////////////////////////////////////////////////////////////////////////////
float computeDistanceCenterFibers(
        const std::vector< std::vector< float > >& medialPointsFibersTractogram,
        int indexFiber1,
        int indexFiber2 )
{

  const std::vector<float>& medialPoint1 = medialPointsFibersTractogram[
                                                                 indexFiber1 ] ;
  const std::vector<float>& medialPoint2 = medialPointsFibersTractogram[
                                                                 indexFiber2 ] ;

  float distanceBetweenCenters = 0 ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    distanceBetweenCenters += pow( medialPoint1[ i ] - medialPoint2[ i ], 2 ) ;

  }
   distanceBetweenCenters = sqrt( distanceBetweenCenters ) ;

   return distanceBetweenCenters ;
}


////////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ////////////////
////////////////////////////////////////////////////////////////////////////////
void computeMeanFiber( const std::vector< std::vector< float > >& fibers,
                       std::vector<float>& meanFiber,
                       int nbCurves,
                       int nbPoints )
{

  if ( meanFiber.size() != 3 * nbPoints )
  {

    std::cout << "ERROR in computeMeanFiber : the size of meanFiber has to be "
              << "3 * nbPoints " << std::endl ;
    exit( 1 ) ;

  }

  for ( int i = 0 ; i < 3 * nbPoints ; i++ )
  {

    meanFiber[ i ] = 0 ;

  }

  if ( fibers.size() != nbCurves )
  {

    std::cout << "ERROR in computeMeanFiber : the size of fibers has to be "
              << "nbCurves " << std::endl ;
    exit( 1 ) ;

  }

  for ( int fiberIndex  = 0 ; fiberIndex < nbCurves ; fiberIndex++ )
  {

    for ( int i = 0 ; i < 3 * nbPoints ; i++)
    {

      meanFiber[ i ] += fibers[ fiberIndex ][ i ] / nbCurves ;

    }

  }


}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Main //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  const auto starting_time = std::chrono::system_clock::now() ;

  const auto start_time = std::chrono::system_clock::now() ;

  int index_input, index_output, index_thr, index_verbose, index_help ;
  index_input = getFlagPosition( argc, argv, "-i") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_thr = getFlagPosition( argc, argv, "-thr") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-i : Input tractogram from which to compute the centroids\n"
              << "-o : Output path for centroids \n"
              << "[-thr] : Threshold distance (Default = average min distance "
              << "between fibers * 1.1 ) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_thr )
  {

    thresholdDistance = std::stof( argv[ index_thr + 1 ] ) ;

    if ( thresholdDistance <= 0 )
    {

      std::cout << "Error argument : thr must be greater than 0"
                << std::endl ;
      exit( 1 ) ;

    }

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

  std::string inputTractogramPath( argv[ index_input + 1 ] ) ;
  char lastChar = inputTractogramPath[ inputTractogramPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputTractogramPath = inputTractogramPath.substr( 0,
                                              inputTractogramPath.size() - 1 ) ;

  }
  std::string inputBundlesDataFilename ;
  std::string inputBundlesFilename ;
  std::string inputTRKFilename ;
  bool isBundlesFormat = false ;
  bool isTRKFormat = false ;

  if ( inputTractogramPath.find( ".bundlesdata" ) != std::string::npos )
  {

    inputBundlesDataFilename = inputTractogramPath ;

    inputBundlesFilename = inputTractogramPath ;
    size_t index = inputTractogramPath.find( ".bundlesdata" ) ;
    inputBundlesFilename.replace( index, 12, ".bundles") ;

    isBundlesFormat = true ;

  }
  else if ( inputTractogramPath.find( ".bundles" ) != std::string::npos )
  {

    inputBundlesFilename = inputTractogramPath ;

    inputBundlesDataFilename = inputTractogramPath ;
    size_t index = inputTractogramPath.find( ".bundles" ) ;
    inputBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

    isBundlesFormat = true ;

  }
  else if ( inputTractogramPath.find( ".trk" ) != std::string::npos )
  {

    inputTRKFilename = inputTractogramPath ;

    isTRKFormat = true ;

  }
  else
  {

    std::cout << "ERROR : Only output format supported is "
              << ".bundles/.bundlesdata " << std::endl ;
    exit( 1 ) ;

  }

  /// ********************************************************************** ///

  std::string outputCentroidsPath( argv[ index_output + 1 ] ) ;
  lastChar = outputCentroidsPath[ outputCentroidsPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    outputCentroidsPath = outputCentroidsPath.substr( 0,
                                              outputCentroidsPath.size() - 1 ) ;

  }
  std::string outputBundlesDataFilename ;
  std::string outputBundlesFilename ;
  std::string outputTRKFilename ;

  if ( outputCentroidsPath.find( ".bundlesdata" ) != std::string::npos )
  {

    outputBundlesDataFilename = outputCentroidsPath ;

    outputBundlesFilename = outputCentroidsPath ;
    size_t index = outputCentroidsPath.find( ".bundlesdata" ) ;
    outputBundlesFilename.replace( index, 12, ".bundles") ;

    isBundlesFormat = true ;

  }
  else if ( outputCentroidsPath.find( ".bundles" ) != std::string::npos )
  {

    outputBundlesFilename = outputCentroidsPath ;

    outputBundlesDataFilename = outputCentroidsPath ;
    size_t index = outputCentroidsPath.find( ".bundles" ) ;
    outputBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

    isBundlesFormat = true ;

  }
  else if ( outputCentroidsPath.find( ".trk" ) != std::string::npos )
  {

    outputTRKFilename = outputCentroidsPath ;

    isTRKFormat = true ;

  }
  else
  {

    std::cout << "ERROR : Only output format supported is "
              << ".bundles/.bundlesdata " << std::endl ;
    exit( 1 ) ;

  }

  //////////////////////////////////////////////////////////////////////////////
  // isBundlesFormat = true ;
  // isTRKFormat = false ;

  //   xxxxxxxxxxxxxxxxxxxxx Reading inputTractogram xxxxxxxxxxxxxxxxxxxxxx   //

  BundlesFormat inputTractogramInfo ;
  BundlesDataFormat inputTractogram ;
  TrkFormat trkData ;


  if ( isBundlesFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << inputBundlesFilename << std::endl ;

    }

    inputTractogramInfo.bundlesReading(
                                       inputBundlesFilename.c_str(), verbose ) ;


    if ( verbose )
    {

      std::cout << "Reading : " << inputBundlesDataFilename << std::endl ;

    }


    inputTractogram.bundlesdataReading( inputBundlesDataFilename.c_str(),
                                       inputBundlesFilename.c_str(),
                                       verbose ) ;

  }
  else if ( isTRKFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << inputTRKFilename << std::endl ;

    }

    trkData.trkReading( inputTRKFilename.c_str(), verbose ) ;
    inputTractogram.matrixTracks = trkData.matrixTracks ;
    inputTractogram.pointsPerTrack = trkData.pointsPerTrack ;
    inputTractogram.curves_count = trkData.curves_count ;

  }
  else
  {

    std::cout << "ERROR : The only supported format are .bundles/.bundlesdata "
              << "and .trk" << std::endl ;
    exit( 1 ) ;

  }


  //////////////////////////////////////////////////////////////////////////////
  // Checking number of points
  int nbCurves = inputTractogram.curves_count ;

  int nbPoints = inputTractogram.pointsPerTrack[ 0 ] ;
  for ( int fiber = 1 ; fiber < nbCurves ; fiber++ )
  {

    if ( nbPoints != inputTractogram.pointsPerTrack[ fiber ] )
    {

      std::cout << "ERROR : the number of points of all fibers in the "
                << "tractogram has to be the same, got " << nbPoints << " for "
                << "fiber 0 and " << inputTractogram.pointsPerTrack[ fiber ]
                << " for fiber" << fiber << std::endl ;
      exit( 1 ) ;

    }

  }

  //////////////////////////////////////////////////////////////////////////////
  // Computing medial points
    std::vector< std::vector< float > > medialPointsFibersTractogram(
                                              nbCurves,
                                              std::vector< float >( 3, 0 ) ) ;

  computeMedialPointsFibersTractogram( inputTractogram,
                                       medialPointsFibersTractogram ) ;



  //////////////////////////////////////////////////////////////////////////////
  // Checking if thresholdDistance is correct

  if ( thresholdDistance <= 0 )
  {

    thresholdDistance = getAverageMinDistanceBetweenFibers(
                                                inputTractogram,
                                                medialPointsFibersTractogram ) ;
    thresholdDistance *= 1.1 ;

  }


  //////////////////////////////////////////////////////////////////////////////
  // Clustering
  std::vector< int > labelsFibers ;
  int nbClusters = 0 ;

  bool stop = false ;
  while ( !stop )
  {

    int tmpNbClusters = 0 ;
    std::vector< int > tmpLabelsFibers( nbCurves, -1 ) ;
    for ( int fiber1 = 0 ; fiber1 < nbCurves ; fiber1++ )
    {

      if ( tmpLabelsFibers[ fiber1 ] < 0 )
      {

        tmpLabelsFibers[ fiber1 ] = tmpNbClusters ;
        tmpNbClusters += 1 ;

      }

      if ( fiber1 + 1 < nbCurves )
      {

        #pragma omp parallel for
        for ( int fiber2 = ( fiber1 + 1 ) ; fiber2 < nbCurves ; fiber2++ )
        {

          int& labelFiber = tmpLabelsFibers[ fiber2 ] ;
          if ( labelFiber < 0 )
          {

            float distanceBetweenCenters = computeDistanceCenterFibers(
                                                   medialPointsFibersTractogram,
                                                   fiber1,
                                                   fiber2 ) ;

            if ( distanceBetweenCenters < thresholdDistance )
            {


              labelFiber = tmpLabelsFibers[ fiber1 ] ;

            }

          }

        }

      }

    }

    if ( tmpNbClusters <= maxNbClusters )
    {

      labelsFibers = tmpLabelsFibers ;
      nbClusters = tmpNbClusters ;
      stop = true ;

    }
    else
    {

      thresholdDistance *= 1.5 ;
    }

  }


  std::vector< std::vector< float > > centroids( nbClusters,
                                     std::vector< float >( 3 * nbPoints, 0 ) ) ;
  #pragma omp parallel for
  for ( int label = 0 ; label < nbClusters ; label++ )
  {

    std::vector< float >& centroidCluster = centroids[ label ] ;

    std::vector< std::vector< float> > labeledFibers ;

    for ( int fiberIndex = 0 ; fiberIndex < nbCurves ; fiberIndex++ )
    {

      if ( labelsFibers[ fiberIndex ] == label )
      {

        std::vector<float> fiber( 3 * nbPoints, 0 ) ;
        inputTractogram.getFiberFromTractogram( inputTractogram.matrixTracks,
                                                fiberIndex,
                                                nbPoints,
                                                fiber ) ;
        labeledFibers.push_back( fiber ) ;


      }

    }

    if ( labeledFibers.size() > 0 )
    {

      int nbFibers = labeledFibers.size() ;
      std::vector<float> meanFiber( 3 * nbPoints, 0 ) ;
      computeMeanFiber( labeledFibers, meanFiber, nbFibers, nbPoints ) ;

      centroidCluster = meanFiber ;

    }

  }

  BundlesDataFormat centroidsBundleData ;
  centroidsBundleData.curves_count = nbClusters ;
  centroidsBundleData.pointsPerTrack.resize( nbClusters, 0 ) ;
  centroidsBundleData.matrixTracks.resize( 3 * nbPoints * nbClusters, 0 ) ;
  for ( int fiberIndex = 0 ; fiberIndex < nbClusters ; fiberIndex++ )
  {

     centroidsBundleData.pointsPerTrack[ fiberIndex ] = nbPoints ;
     std::copy( centroids[ fiberIndex ].begin(),
                centroids[ fiberIndex ].begin() + 3 * nbPoints,
                centroidsBundleData.matrixTracks.begin() + 3 * nbPoints *
                                                                  fiberIndex ) ;

  }

  centroidsBundleData.bundlesdataWriting( outputBundlesDataFilename.c_str(),
                                                                           0 ) ;

  float disimilarity = centroidsBundleData.compareDisimilarityBundles(
                                               centroidsBundleData.matrixTracks,
                                               inputTractogram.matrixTracks,
                                               centroidsBundleData.curves_count,
                                               inputTractogram.curves_count,
                                               nbPoints,
                                               nbPoints ) ;

  float coverage = centroidsBundleData.coverageBundle1toBundle2(
                                                          centroidsBundleData,
                                                          inputTractogram,
                                                          thresholdAdjacency,
                                                          0 ) ;

  float overlap = centroidsBundleData.overlapBundle1toBundle2(
                                                          centroidsBundleData,
                                                          inputTractogram,
                                                          thresholdAdjacency,
                                                          0 ) ;

  //
  float adjacency = centroidsBundleData.bundlesAdjacency(
                                                          centroidsBundleData,
                                                          inputTractogram,
                                                          thresholdAdjacency,
                                                          0 ) ;

  BundlesFormat centroidsBundleInfo = inputTractogramInfo ;
  centroidsBundleInfo.curves_count = centroidsBundleData.curves_count ;
  centroidsBundleInfo.bundlesWriting( outputBundlesFilename.c_str(),
                               disimilarity, coverage, overlap, adjacency, 0 ) ;
  // centroidsBundleInfo.bundlesWriting( outputBundlesFilename.c_str(), 0 ) ;


  const std::chrono::duration< double > duration_computations =
                                 std::chrono::system_clock::now() - start_time ;

  std::cout << "Duration : " << duration_computations.count() << " s \n" ;

  return( 0 ) ;


}
