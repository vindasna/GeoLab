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

// Optimization
#include <Eigen/Core>
#include <Eigen/Dense>
#include <LBFGS.h>

#include "BundlesDataFormat.h"
#include "BundlesFormat.h"

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Operator ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::operator[]( int64_t index )
{

  return ( this->matrixTracks[ index ] ) ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesDataFormat::BundlesDataFormat() {}

BundlesDataFormat::BundlesDataFormat( const char* bundlesdataFilename,
                                      const char* bundlesFilename,
                                      int verbose )
{

  bundlesdataReading( bundlesdataFilename, bundlesFilename, verbose ) ;

}

BundlesDataFormat::BundlesDataFormat( std::vector<float>& matrixTracks,
                                      std::vector<int32_t>& pointsPerTrack,
                                      int curves_count )
{

  this->curves_count = curves_count ;

  this->pointsPerTrack = pointsPerTrack ;

  this->matrixTracks = matrixTracks ;


}

BundlesDataFormat::BundlesDataFormat( const BundlesDataFormat& bundlesData )
{

  this->curves_count = bundlesData.curves_count ;

  this->pointsPerTrack = bundlesData.pointsPerTrack ;

  this->matrixTracks = bundlesData.matrixTracks ;

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesDataFormat::~BundlesDataFormat() {}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::bundlesdataReading( const char* bundlesdataFilename,
                                            const char* bundlesFilename,
                                            int verbose )
{

  BundlesFormat bundlesInfo( bundlesFilename, verbose ) ;

  this->curves_count = bundlesInfo.curves_count ;

  if ( this->curves_count == 0 )
  {

    std::cout << "Problem reading file, the number of traks was not stored "
              << "in the header " << std::endl ;

    std::exit( 1 ) ;

  }
  else
  {

    if ( verbose > 1 )
    {

      std::cout << "Number of curves : " << this->curves_count << std::endl ;

    }


  }

  int64_t sizeMatrixTracks = 0 ;

  this->pointsPerTrack.resize( this->curves_count, 0 ) ;

  std::ifstream file ;
  file.open( bundlesdataFilename, std::ios::binary | std::ios::in ) ;
  if ( file.fail() )
  {

    std::cout << "Problem reading file : " << bundlesdataFilename << std::endl ;
    exit( 1 ) ;

  }

  // Getting number of points per curve
  int64_t offsetBytes = 0 ;
  for( int track = 0 ; track < this->curves_count ; track++ )
  {

    if ( verbose > 1 && ( track % 1000 == 0 ||
                                  ( track + 1 ) == this->curves_count ) )
    {

      printf("\rProcessing tracks [ %10d / %10d ]",
                                        track + 1 , this->curves_count ) ;
      std::cout << "" << std::flush ;

    }

    file.read( reinterpret_cast<char*>( &( this->pointsPerTrack[ track ] ) ),
                                                           sizeof( int32_t ) ) ;

    offsetBytes += sizeof( int32_t ) +
                           this->pointsPerTrack[ track ] * 3 * sizeof( float ) ;

    file.seekg( offsetBytes, std::ios_base::beg ) ;

    sizeMatrixTracks += this->pointsPerTrack[ track ] * 3 ;

  }

  // Getting vector with fiber coordinates
  this->matrixTracks.resize( sizeMatrixTracks, 0 ) ;
  file.seekg( 0, std::ios_base::beg ) ;
  int offsetPoints = 0 ;
  for ( int track = 0 ; track < this->curves_count ; track++ )
  {

    file.seekg( sizeof( int32_t ), std::ios_base::cur ) ;
    for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
    {

      file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 0 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
      file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 1 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
      file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 2 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;

    }

    offsetPoints += 3 * this->pointsPerTrack[ track ] ;

  }

  file.close() ;

  if ( verbose > 1 )
  {

    std::cout << "\n" ;

  }


  // Check for NaN values
  int64_t nbCurves = this->curves_count ;
  int64_t offsetTractogram = 0 ;
  int nbFibersWithNan = 0 ;
  for ( int fiber = 0 ; fiber < nbCurves ; fiber++ )
  {

    if ( verbose > 1 && ( fiber % 1000 == 0 ||
                                  ( fiber + 1 ) == this->curves_count ) )
    {

      printf("\rChecking for NaN values in fibers : [ %10d / %10d ]",
                                        fiber + 1 , this->curves_count ) ;
      std::cout << "" << std::flush ;

    }

    int nbPoints = this->pointsPerTrack[ fiber ] ;
    bool fiberResampled = false ;
    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        if ( isnan( this->matrixTracks[ offsetTractogram + 3 * point + i ] ) )
        {

          resampleFiberWithNan( this->matrixTracks ,
                                fiber,
                                nbPoints,
                                verbose + 1 ) ;

          fiberResampled = true ;
          fibersWithNans.push_back( fiber ) ;
          nbFibersWithNan++ ;

          break ;

        }

      }

      if ( fiberResampled )
      {

        break ;

      }

    }

    offsetTractogram += 3 * nbPoints ;

  }

  if ( verbose > 1 )
  {

    if ( nbFibersWithNan )
    {

      std::cout << "\nThere are " << nbFibersWithNan << " with NaN values \n" ;

    }
    else
    {

      std::cout << "\nThere are no fibers with NaN values" << std::endl ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::bundlesdataWriting( const char* bundlesdataFilename,
                                            int verbose )
{

  if ( verbose > 1 )
  {

    std::cout << "Writing " << bundlesdataFilename << std::endl ;

  }



  std::ofstream file ;
  file.open( bundlesdataFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << bundlesdataFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  int n_curves = this->curves_count ;


  int offsetPoints = 0 ;


  for ( int track = 0 ; track < n_curves ; track++ )
  {

    file.write( reinterpret_cast<char*>( &this->pointsPerTrack[ track ] ),
                                                           sizeof( int32_t ) ) ;

    for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
    {

      file.write( reinterpret_cast<char*>( &this->matrixTracks[ 0 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>( &this->matrixTracks[ 1 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>( &this->matrixTracks[ 2 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;

    }

    offsetPoints += 3 * this->pointsPerTrack[ track ] ;

  }

  file.close() ;

}


////////////////////////////////////////////////////////////////////////////////
// void BundlesFormat::toTRK( const BundlesFormat& bundlesInfo,
//                            TrkFormat& trkData )
// {
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     trkData.dim[ i ] = bundlesInfo.resolution[ i ] ;
//     trkData.voxel_size[ i ] = bundlesInfo.size[ i ] ;
//
//   }
//
//   trkData.curves_count = this->curves_count ;
//   trkData.matrixTracks = this->matrixTracks ;
//   trkData.pointsPerTrack = this->pointsPerTrack ;
//
// }



////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::computeLengthFiber(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int nbPoints )
{

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  float lengthFiber = 0 ;
  for ( int point = 1 ; point < nbPoints ; point++ )
  {

    float tmpDistancePoints = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

       tmpDistancePoints += pow(
               tractogramFibers[ 3 * ( point - 1 ) + i + offsetTractogram ] -
               tractogramFibers[ 3 * point + i + offsetTractogram ], 2 ) ;

    }

    tmpDistancePoints = sqrt( tmpDistancePoints ) ;

    lengthFiber += tmpDistancePoints ;

  }

  return lengthFiber ;

}

//----------------------------------------------------------------------------//
float BundlesDataFormat::computeLengthFiber( const std::vector<float>& fiber )
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;

  float lengthFiber = computeLengthFiber( fiber, fiberIndex, nbPoints ) ;

  return lengthFiber ;

}

//----------------------------------------------------------------------------//
float BundlesDataFormat::computeLengthFiber( int fiberIndex )
{

  int nbPoints = this->pointsPerTrack[ fiberIndex ] ;
  std::vector<float>& tractogramFibers = this->matrixTracks ;

  float lengthFiber = computeLengthFiber( tractogramFibers, fiberIndex,
                                                                    nbPoints ) ;

  return lengthFiber ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::resampleFiberEquidistant(
                                           const std::vector<float>& inputFiber,
                                           std::vector<float>& outputFiber,
                                           int nbPointsToResample )
{

  int nbPointsInputFiber = ( int )( inputFiber.size() / 3 ) ;

  float lengthFiber = computeLengthFiber( inputFiber ) ;

  float stepSize = lengthFiber / ( nbPointsToResample - 1 ) ;

  outputFiber.resize( 3 * nbPointsToResample, 0 ) ;

  std::vector<float> pointToTranslate( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    int point = 0 ;
    outputFiber[ 3 * point + i ] = inputFiber[ 3 * point + i ] ;
    pointToTranslate[ i ] = inputFiber[ 3 * point + i ] ;

  }


  int pointInInputFiber = 1 ;
  for ( int point = 1 ; point < nbPointsToResample ; point++ )
  {

    float translationDistance = 0 ;
    while ( translationDistance < stepSize )
    {

      std::vector<float> translationVector( 3, 0 ) ;
      std::vector<float> point1( 3, 0 ) ;
      std::vector<float> point2( 3, 0 ) ;

      for ( int i = 0 ; i < 3 ; i++)
      {

        point1[ i ] = pointToTranslate[ i ] ;
        point2[ i ] = inputFiber[ 3 * pointInInputFiber + i ] ;

      }
      vectorFromPoints( point2, point1, translationVector ) ;
      float distanceBetweenPoints = normVector( translationVector ) ;
      normalizeVector( translationVector ) ;

      float residualDistance = stepSize - translationDistance ;

      if ( distanceBetweenPoints < residualDistance )
      {

        if ( pointInInputFiber == nbPointsInputFiber - 1 )
        {

          translatePoint( pointToTranslate, translationVector,
                                          residualDistance, pointToTranslate ) ;
          translationDistance += residualDistance ;

        }
        else
        {

          translatePoint( pointToTranslate, translationVector,
                                     distanceBetweenPoints, pointToTranslate ) ;
          translationDistance += distanceBetweenPoints ;
          pointInInputFiber++ ;

        }

      }
      else
      {

        translatePoint( pointToTranslate, translationVector,
                                       residualDistance, pointToTranslate ) ;

        translationDistance += residualDistance ;

      }

    }

    for ( int i = 0 ; i < 3 ; i++ )
    {

      outputFiber[ 3 * point + i ] = pointToTranslate[ i ] ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::computeMedialPointFiberTractogram(
                                const std::vector<float>& tractogramFibers,
                                int fiberIndex,
                                int nbPoints,
                                std::vector<float>& medialPointFiberTractogram )
{

  int64_t offsetTractogram = nbPoints * 3 * fiberIndex ;

  if ( nbPoints % 2 != 0 )
  {

    int64_t tmpIndexTract = 3 * ( nbPoints - 1 ) / 2 + offsetTractogram ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointFiberTractogram[ i ] = tractogramFibers[ tmpIndexTract + i ] ;

    }

  }
  else
  {

    int64_t tmpIndexTract = 3 * nbPoints / 2 + offsetTractogram ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointFiberTractogram[ i ] = ( tractogramFibers[ tmpIndexTract + i ]
                        + tractogramFibers[ tmpIndexTract + i + 3 ] ) / 2 ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeMedialPointFiberTractogram(
                                int fiberIndex,
                                std::vector<float>& medialPointFiberTractogram )
{

  std::vector<float>& tractogramFibers = this->matrixTracks ;

  int nbPoints = this->pointsPerTrack[ fiberIndex ] ;

  computeMedialPointFiberTractogram( tractogramFibers, fiberIndex, nbPoints,
                                                  medialPointFiberTractogram ) ;

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeMedialPointFiberTractogram(
                                const std::vector<float>& fiber,
                                std::vector<float>& medialPointFiberTractogram )
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;
  computeMedialPointFiberTractogram( fiber, fiberIndex, nbPoints,
                                                  medialPointFiberTractogram ) ;


}

////////////////////////////////////////////////////////////////////////////////
int BundlesDataFormat::computeMedialPointFiberWithDistance(
                                const std::vector<float>& tractogramFibers,
                                int fiberIndex,
                                int nbPoints,
                                std::vector<float>& medialPointFiberTractogram )
{

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  float lengthFiber = computeLengthFiber( tractogramFibers, fiberIndex,
                                                                    nbPoints ) ;

  float dist = 0 ;
  int point = 0 ;
  while ( dist < lengthFiber / 2 )
  {

    point++ ;

    float tmpDist = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      tmpDist += pow( tractogramFibers[ 3 * point + i + offsetTractogram ] -
             tractogramFibers[ 3 * ( point - 1 ) + i + offsetTractogram ], 2 ) ;
    }

    tmpDist = sqrt( tmpDist ) ;

    dist += tmpDist ;

  }

  std::vector<float> translationVector( 3 , 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    int index1 = 3 * ( point - 1 ) + i + offsetTractogram ;
    int index2 = 3 * point + i + offsetTractogram ;

    // Vector going from index2 to index1 because we are going to search the
    // center going from dist > length / 2 until dist < length by changing
    // point2 where dist > length ( and not point1 where dist < length / 2 )
    translationVector[ i ] = tractogramFibers[ index1 ] -
                                                    tractogramFibers[ index2 ] ;

  }

  float normTranslationVector = scalarProduct( translationVector,
                                                           translationVector ) ;
  normTranslationVector = sqrt( normTranslationVector ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translationVector[ i ] /= normTranslationVector ; // Unitary vector

  }

  normTranslationVector = scalarProduct( translationVector,
                                                           translationVector ) ;
  normTranslationVector = sqrt( normTranslationVector ) ;

  float step = 0.1 ;
  int numberOfSteps = 0 ;
  while ( dist > lengthFiber / 2 )
  {

    dist -= step ;
    numberOfSteps++ ;

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    int index = 3 * point + i + offsetTractogram ;
    medialPointFiberTractogram[ i ] = tractogramFibers[ index ] + step *
                                        numberOfSteps * translationVector[ i ] ;
                                        // It's + step and not - step because
                                        // it's a vector translation

  }

  return ( point - 1 ) ; // Return point where dist < length / 2

}

//----------------------------------------------------------------------------//
int BundlesDataFormat::computeMedialPointFiberWithDistance(
                                int fiberIndex,
                                std::vector<float>& medialPointFiberTractogram )
{

  std::vector<float>& tractogramFibers = this->matrixTracks ;

  int nbPoints = this->pointsPerTrack[ fiberIndex ] ;

  int point = computeMedialPointFiberWithDistance( tractogramFibers, fiberIndex,
                                        nbPoints, medialPointFiberTractogram ) ;

  return point ; // It is point and no point - 1 because the function use to
                 // compute point returns already point - 1

}

//----------------------------------------------------------------------------//
int BundlesDataFormat::computeMedialPointFiberWithDistance(
                                const std::vector<float>& fiber,
                                std::vector<float>& medialPointFiberTractogram )
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;
  int point = computeMedialPointFiberWithDistance( fiber, fiberIndex, nbPoints,
                                                  medialPointFiberTractogram ) ;

  return point ; // It is point and no point - 1 because the function use to
                 // compute point returns already point - 1


}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::computeGravityCenterFiber(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int nbPoints,
                                     std::vector<float>& gravityCenterFiber )
{

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterFiber[ i ] = 0 ;

  }

  int sampling = nbPoints ;

  float dis = computeLengthFiber( tractogramFibers, fiberIndex, nbPoints ) ;

  float step = dis / sampling ;

  std::vector<float> pointLeft( 3, 0 ) ;
  std::vector<float> pointRight( 3 , 0 ) ;
  std::vector<float> pointStep( 3, 0 ) ;
  std::vector<float> direction( 3, 0 ) ;

  int pointLeftIndex = 0 ;
  int pointRightIndex = 1 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    pointLeft[ i ] = tractogramFibers[ 3 * pointLeftIndex + i +
                                                            offsetTractogram ] ;
    pointRight[ i ] = tractogramFibers[ 3 * pointRightIndex + i +
                                                            offsetTractogram ] ;

  }

  vectorFromPoints( pointRight, pointLeft, direction ) ;
  float segmentLength = normVector( direction ) ;
  normalizeVector( direction ) ;

  translatePoint( pointLeft, direction, step, pointStep ) ;

  if ( segmentLength < step )
  {

    float translationDistance = 0 ;

    while ( segmentLength < step )
    {

      pointLeftIndex++ ;
      pointRightIndex++ ;

      if ( pointRightIndex > ( nbPoints - 1 ) )
      {

        dis = -1 ;
        break ;

      }

      pointLeft = pointRight ;
      for ( int i = 0 ; i < 3 ; i++ )
      {

        pointRight[ i ] = tractogramFibers[ 3 * pointRightIndex + i +
                                                            offsetTractogram ] ;

      }

      vectorFromPoints( pointRight, pointLeft, direction ) ;
      translationDistance = step - segmentLength ;
      segmentLength += normVector( direction ) ;

    }

    translatePoint( pointLeft, direction, translationDistance, pointStep ) ;


  }

  int count = 0 ;


  while ( dis > 0 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      gravityCenterFiber[ i ] += pointStep[ i ] ;

    }

    pointLeft = pointStep ;

    vectorFromPoints( pointRight, pointLeft, direction ) ;
    segmentLength = normVector( direction ) ;
    normalizeVector( direction ) ;

    translatePoint( pointLeft, direction, step, pointStep ) ;

    if ( segmentLength < step )
    {

      float translationDistance = 0 ;

      while ( segmentLength < step )
      {

        pointLeftIndex++ ;
        pointRightIndex++ ;

        if ( pointRightIndex > ( nbPoints - 1 ) )
        {

          dis = -1 ;
          break ;

        }

        pointLeft = pointRight ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          pointRight[ i ] = tractogramFibers[ 3 * pointRightIndex + i +
                                                            offsetTractogram ] ;

        }

        vectorFromPoints( pointRight, pointLeft, direction ) ;
        translationDistance = step - segmentLength ;
        segmentLength += normVector( direction ) ;

      }

      translatePoint( pointLeft, direction, translationDistance, pointStep ) ;

    }

    dis -= step ;

    count++ ;

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterFiber[ i ] /= count ;

  }

  // std::cout << "Sampling = " << sampling << "   |   Count = " << count << "\n" ;

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeGravityCenterFiber(
                                        int fiberIndex,
                                        std::vector<float>& gravityCenterFiber )
{

  int nbPoints = this->pointsPerTrack[ fiberIndex ] ;
  std::vector<float>& tractogramFibers = this->matrixTracks ;

  computeGravityCenterFiber( tractogramFibers,
                             fiberIndex,
                             nbPoints,
                             gravityCenterFiber ) ;

}
// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeGravityCenterFiber(
                                        const std::vector<float>& fiber,
                                        std::vector<float>& gravityCenterFiber )
{

  int fiberIndex = 0 ;
  int nbPoints = ( int )( fiber.size() / 3 ) ;


  computeGravityCenterFiber( fiber,
                             fiberIndex,
                             nbPoints,
                             gravityCenterFiber ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::computeNormalVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& normalVector )
{

  Eigen::MatrixXf AMatrix ;
  Eigen::VectorXf BVector ;

  getMatrixAndVectorForLeastSquares( tractogramFibers,
                                     fiberIndex,
                                     nbPoints,
                                     AMatrix,
                                     BVector ) ;

  // Eigen::VectorXf solutionVector = AMatrix.bdcSvd(
  //                 Eigen::ComputeThinU | Eigen::ComputeThinV ).solve( BVector ) ;
  // Eigen::VectorXf solutionVector = AMatrix.colPivHouseholderQr().solve(
  //                                                                    BVector ) ;
  Eigen::VectorXf solutionVector = ( AMatrix.transpose() * AMatrix
                               ).ldlt().solve( AMatrix.transpose() * BVector ) ;


  normalVector[ 0 ] = solutionVector( 0 ) ;
  normalVector[ 1 ] = solutionVector( 1 ) ;
  normalVector[ 2 ] = -1 ;

  normalizeVector( normalVector ) ;


}
// void BundlesDataFormat::computeNormalVectorFiberTractogram(
//                            const std::vector<float>& tractogramFibers,
//                            const std::vector<float>& medialPointFiberTractogram,
//                            int fiberIndex,
//                            int nbPoints,
//                            std::vector<float>& normalVector )
// {
//
//   int64_t offset = 3 * nbPoints * fiberIndex ;
//
//   std::vector<float> vector1( 3, 0 ) ;
//   std::vector<float> vector2( 3, 0 ) ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     vector1[ i ] = tractogramFibers[ 3 * 0 + i + offset ] -
//                                                medialPointFiberTractogram[ i ] ;
//     vector2[ i ] = tractogramFibers[ 3 * ( nbPoints - 1 ) + i + offset ] -
//                                                medialPointFiberTractogram[ i ] ;
//
//   }
//
//   normalizeVector( vector1 ) ;
//   normalizeVector( vector2 ) ;
//
//   crossProduct( vector1, vector2, normalVector ) ;
//
//   normalizeVector( normalVector ) ;
//
//
// }

// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeNormalVectorFiberTractogram(
                                              int fiberIndex,
                                              std::vector<float>& normalVector )
{

  int nbPoints = this->pointsPerTrack[ fiberIndex ];
  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // this->computeMedialPointFiberTractogram( fiberIndex,
  //                                          medialPointFiberTractogram ) ;
  int medialPointIndex = computeMedialPointFiberWithDistance(
                                                  fiberIndex,
                                                  medialPointFiberTractogram ) ;

  computeNormalVectorFiberTractogram( this->matrixTracks,
                                      medialPointFiberTractogram,
                                      fiberIndex,
                                      nbPoints,
                                      normalVector ) ;

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeNormalVectorFiberTractogram(
                                              const std::vector<float>& fiber,
                                              std::vector<float>& normalVector )
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;

  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber, medialPointFiberTractogram ) ;
  int medialPointIndex = computeMedialPointFiberWithDistance(
                                                  fiber,
                                                  medialPointFiberTractogram ) ;

  int fiberIndex = 0 ;
  computeNormalVectorFiberTractogram( fiber,
                                      medialPointFiberTractogram,
                                      fiberIndex,
                                      nbPoints,
                                      normalVector ) ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::computeDirectionVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           const std::vector<float>& normalVector,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& directionVector )
{


  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  std::vector<float> point1( 3, 0 ) ;
  std::vector<float> point2( 3, 0 ) ;

  std::vector<float> vector1( 3, 0 ) ;
  std::vector<float> vector2( 3, 0 ) ;
  std::vector<float> vector3( 3, 0 ) ;


  for ( int i = 0 ; i < 3 ; i++ )
  {

    point1[ i ] = tractogramFibers[ 3 * 0 + i + offsetTractogram ] ;
    point2[ i ] = tractogramFibers[ 3 * ( nbPoints - 1 ) + i +
                                                          offsetTractogram ] ;
    vector1[ i ] = point1[ i ] - medialPointFiberTractogram[ i ] ;
    vector2[ i ] = point2[ i ] - medialPointFiberTractogram[ i ] ;



    // point1[ i ] = tractogramFibers[ 3 * 0 + i + offsetTractogram ] ;
    // point2[ i ] = tractogramFibers[ 3 * ( nbPoints - 1 ) + i +
    //                                                       offsetTractogram ] ;
    // vector1[ i ] = point2[ i ] - point1[ i ] ;

  }

  normalizeVector( vector1 ) ;
  normalizeVector( vector2 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    vector3[ i ] = vector1[ i ] + vector2[ i ] ;

  }


  // crossProduct( vector1, normalVector, vector3 ) ;

  // directionVector = vector3 ;

  projectVectorToPlane( normalVector, vector3, directionVector ) ;

  normalizeVector( directionVector ) ;

}



// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeDirectionVectorFiberTractogram(
                                         int fiberIndex,
                                         const std::vector<float>& normalVector,
                                         std::vector<float>& directionVector )
{

  int nbPoints = this->pointsPerTrack[ 0 ] ;

  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // this->computeMedialPointFiberTractogram( fiberIndex,
  //                                          medialPointFiberTractogram ) ;
  this->computeMedialPointFiberWithDistance( fiberIndex,
                                             medialPointFiberTractogram ) ;

  computeDirectionVectorFiberTractogram( this->matrixTracks,
                                         medialPointFiberTractogram,
                                         normalVector,
                                         fiberIndex,
                                         nbPoints,
                                         directionVector ) ;

}


// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeDirectionVectorFiberTractogram(
                                         const std::vector<float>& fiber,
                                         const std::vector<float>& normalVector,
                                         std::vector<float>& directionVector )
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;

  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // this->computeMedialPointFiberTractogram( fiber,
  //                                          medialPointFiberTractogram ) ;
  int medialPointIndex = computeMedialPointFiberWithDistance( fiber,
                                             medialPointFiberTractogram ) ;

  computeDirectionVectorFiberTractogram( fiber,
                                         medialPointFiberTractogram,
                                         normalVector,
                                         fiberIndex,
                                         nbPoints,
                                         directionVector ) ;

}


////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::computeMDADBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints )
{

  std::vector<float> fiber1( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_1, fiberIndex_1, nbPoints, fiber1 ) ;
  std::vector<float> resampledFiber1( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber1, resampledFiber1, nbPoints ) ;

  std::vector<float> fiber2( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_2, fiberIndex_2, nbPoints, fiber2 ) ;
  std::vector<float> resampledFiber2( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber2, resampledFiber2, nbPoints ) ;

  float dMDA = 0 ;

  ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
  std::vector<float> translation( 3, 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointTractFiber_1[ i ] -
                                                  medialPointTractFiber_2[ i ] ;

  }

  float tmpDirectDistance = 0 ;
  float tmpFlippedDistance = 0 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    tmpDirectDistance += pow( resampledFiber1[ i + 3 * 0 ] -
                                ( resampledFiber2[ i + 3 * 0 ] +
                                                       translation[ i ] ), 2 ) ;

    tmpFlippedDistance += pow( resampledFiber1[
                                i + 3 * 0 ] -
                                ( resampledFiber2[ i + 3 * ( nbPoints - 1 ) ]
                                                     + translation[ i ] ), 2 ) ;

  }

  bool isDirectSens = true ;
  if ( tmpFlippedDistance < tmpDirectDistance )
  {

    isDirectSens = false ;

  }

  // Computing MDA when distance
  dMDA = 0 ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    int k = point ;
    if ( !isDirectSens )
    {

      k = nbPoints - point - 1 ;

    }

    float tmpMDA = 0 ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      tmpMDA += pow( resampledFiber1[ i + 3 * point ] - ( resampledFiber2[
                                         i + 3 * k ] + translation[ i ] ), 2 ) ;

    }

    tmpMDA = sqrt( tmpMDA ) ;

    if ( tmpMDA > dMDA )
    {

      dMDA = tmpMDA ;

    }

  }

  return dMDA ;

}
// float BundlesDataFormat::computeMDADBetweenTwoFibers(
//                               const std::vector<float>& tractogramFibers_1,
//                               const std::vector<float>& tractogramFibers_2,
//                               const std::vector<float>& medialPointTractFiber_1,
//                               const std::vector<float>& medialPointTractFiber_2,
//                               int fiberIndex_1,
//                               int fiberIndex_2,
//                               int nbPoints )
// {
//
//   float dMDA = 0 ;
//
//   int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;
//   int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;
//
//   ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
//   std::vector<float> translation( 3, 0 ) ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     translation[ i ] = medialPointTractFiber_1[ i ] -
//                                                   medialPointTractFiber_2[ i ] ;
//
//   }
//
//   float tmpDirectDistance = 0 ;
//   float tmpFlippedDistance = 0 ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     tmpDirectDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * 0 +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//     tmpFlippedDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * ( nbPoints - 1 ) +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//   }
//
//   bool isDirectSens = true ;
//   if ( tmpFlippedDistance < tmpDirectDistance )
//   {
//
//     isDirectSens = false ;
//
//   }
//
//   // Computing MDA when distance
//   dMDA = 0 ;
//
//   for ( int point = 0 ; point < nbPoints ; point++ )
//   {
//
//     int k = point ;
//     if ( !isDirectSens )
//     {
//
//       k = nbPoints - point - 1 ;
//
//     }
//
//     float tmpMDA = 0 ;
//
//     for ( int i = 0 ; i < 3 ; i++ )
//     {
//
//       tmpMDA += pow( tractogramFibers_1[ i + 3 * point +
//                      offsetTractogram_1 ] - ( tractogramFibers_2[
//                      i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
//                                                                            2 ) ;
//
//     }
//
//     tmpMDA = sqrt( tmpMDA ) ;
//
//     if ( tmpMDA > dMDA )
//     {
//
//       dMDA = tmpMDA ;
//
//     }
//
//   }
//
//   return dMDA ;
//
// }

//----------------------------------------------------------------------------//
float BundlesDataFormat::computeMDADBetweenTwoFibers(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int fiberIndex_1,
                                   int fiberIndex_2,
                                   int nbPoints )
{

  ///////////// Searching the medial point of tractogram 1st fiber /////////////
  int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_1( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_1,
                                       fiberIndex_1,
                                       nbPoints,
                                       medialPointTractFiber_1 ) ;


  ///////////// Searching the medial point of tractogram 2nd fiber /////////////
  int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_2( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_2,
                                       fiberIndex_2,
                                       nbPoints,
                                       medialPointTractFiber_2 ) ;

  float dMDA = computeMDADBetweenTwoFibers( tractogramFibers_1,
                                            tractogramFibers_2,
                                            medialPointTractFiber_1,
                                            medialPointTractFiber_2,
                                            fiberIndex_1,
                                            fiberIndex_2,
                                            nbPoints ) ;

  return dMDA ;


}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::computeMDFBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints )
{


  std::vector<float> fiber1( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_1, fiberIndex_1, nbPoints, fiber1 ) ;
  std::vector<float> resampledFiber1( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber1, resampledFiber1, nbPoints ) ;

  std::vector<float> fiber2( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_2, fiberIndex_2, nbPoints, fiber2 ) ;
  std::vector<float> resampledFiber2( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber2, resampledFiber2, nbPoints ) ;

  float dMDF = 0 ;

  ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
  std::vector<float> translation( 3, 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointTractFiber_1[ i ] -
                                                  medialPointTractFiber_2[ i ] ;

  }

  float tmpDirectDistance = 0 ;
  float tmpFlippedDistance = 0 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    tmpDirectDistance += pow( resampledFiber1[ i + 3 * 0 ] -
                                                ( resampledFiber2[ i + 3 * 0 ] +
                                                       translation[ i ] ), 2 ) ;

    tmpFlippedDistance += pow( resampledFiber1[ i + 3 * 0 ] -
                                 ( resampledFiber2[ i + 3 * ( nbPoints - 1 ) ] +
                                                       translation[ i ] ), 2 ) ;

  }

  bool isDirectSens = true ;
  if ( tmpFlippedDistance < tmpDirectDistance )
  {

    isDirectSens = false ;

  }

  // Computing MDA when distance
  dMDF = 0 ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    int k = point ;
    if ( !isDirectSens )
    {

      k = nbPoints - point - 1 ;

    }

    float tmpMDF = 0 ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      // tmpMDF += pow( resampledFiber1[ i + 3 * point ] - ( resampledFiber2[
      //                                 i + 3 * k ] + translation[ i ] ), 2 ) ;
      tmpMDF += pow( resampledFiber1[ i + 3 * point ] -  resampledFiber2[
                                                              i + 3 * k ], 2 ) ;

    }

    tmpMDF = sqrt( tmpMDF ) ;

    dMDF += tmpMDF ;

  }

  dMDF /= nbPoints ;

  return dMDF ;


}

// float BundlesDataFormat::computeMDFBetweenTwoFibers(
//                               const std::vector<float>& tractogramFibers_1,
//                               const std::vector<float>& tractogramFibers_2,
//                               const std::vector<float>& medialPointTractFiber_1,
//                               const std::vector<float>& medialPointTractFiber_2,
//                               int fiberIndex_1,
//                               int fiberIndex_2,
//                               int nbPoints )
// {
//
//   float dMDF = 0 ;
//
//   int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;
//   int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;
//
//   ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
//   std::vector<float> translation( 3, 0 ) ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     translation[ i ] = medialPointTractFiber_1[ i ] -
//                                                   medialPointTractFiber_2[ i ] ;
//
//   }
//
//   float tmpDirectDistance = 0 ;
//   float tmpFlippedDistance = 0 ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     tmpDirectDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * 0 +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//     tmpFlippedDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * ( nbPoints - 1 ) +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//   }
//
//   bool isDirectSens = true ;
//   if ( tmpFlippedDistance < tmpDirectDistance )
//   {
//
//     isDirectSens = false ;
//
//   }
//
//   // Computing MDA when distance
//   dMDF = 0 ;
//
//   for ( int point = 0 ; point < nbPoints ; point++ )
//   {
//
//     int k = point ;
//     if ( !isDirectSens )
//     {
//
//       k = nbPoints - point - 1 ;
//
//     }
//
//     float tmpMDF = 0 ;
//
//     for ( int i = 0 ; i < 3 ; i++ )
//     {
//
//       // tmpMDF += pow( tractogramFibers_1[ i + 3 * point +
//       //                offsetTractogram_1 ] - ( tractogramFibers_2[
//       //                i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
//       //                                                                      2 ) ;
//       tmpMDF += pow( tractogramFibers_1[ i + 3 * point +
//                      offsetTractogram_1 ] -  tractogramFibers_2[
//                                          i + 3 * k + offsetTractogram_2 ], 2 ) ;
//
//     }
//
//     tmpMDF = sqrt( tmpMDF ) ;
//
//     dMDF += tmpMDF ;
//
//   }
//
//   dMDF /= nbPoints ;
//
//   return dMDF ;
//
//
// }

//----------------------------------------------------------------------------//
float BundlesDataFormat::computeMDFBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints )
{

  ///////////// Searching the medial point of tractogram 1st fiber /////////////
  int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_1( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_1,
                                       fiberIndex_1,
                                       nbPoints,
                                       medialPointTractFiber_1 ) ;


  ///////////// Searching the medial point of tractogram 2nd fiber /////////////
  int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_2( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_2,
                                       fiberIndex_2,
                                       nbPoints,
                                       medialPointTractFiber_2 ) ;

  float dMDF = computeMDFBetweenTwoFibers( tractogramFibers_1,
                                           tractogramFibers_2,
                                           medialPointTractFiber_1,
                                           medialPointTractFiber_2,
                                           fiberIndex_1,
                                           fiberIndex_2,
                                           nbPoints ) ;

  return dMDF ;


}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::compareDisimilarityBundles(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int nbFibersTract_1,
                                   int nbFibersTract_2,
                                   int nbPoints1,
                                   int nbPoints2 )
{

  if ( nbPoints1 != nbPoints2 )
  {

    std::cout << "ERROR : in compareDisimilarityBundles( const float* "
              << "tractogramFibers_1, const float* tractogramFibers_2, int "
              << "nbFibersTract_1, int nbFibersTract_2, int nbPoints1, int "
              << "nbPoints2 ) the number of points of the fibers must be the "
              << "same " << std::endl ;
    exit( 1 ) ;

  }
  int nbPoints = nbPoints1 ;

  float disimilarity = 0 ;

  #pragma omp parallel for reduction( + : disimilarity )
  for ( int fiberIndex_1 = 0 ; fiberIndex_1 < nbFibersTract_1 ; fiberIndex_1++ )
  {

    float minDMDA = 1000 ;

    for ( int fiberIndex_2 = 0 ; fiberIndex_2 < nbFibersTract_2 ; fiberIndex_2++ )
    {

      float dMDA = computeMDADBetweenTwoFibers( tractogramFibers_1,
                                                tractogramFibers_2,
                                                fiberIndex_1,
                                                fiberIndex_2,
                                                nbPoints ) ;

      if ( minDMDA > dMDA )
      {

        minDMDA = dMDA ;

      }

    }

    disimilarity += minDMDA ;

  }

  disimilarity /= nbFibersTract_1 ;

  return disimilarity ;

}

//----------------------------------------------------------------------------//
float BundlesDataFormat::compareDisimilarityBundles(
                                     const std::vector<float>& tractogramFibers,
                                     int nbFibersTract,
                                     int nbPoints )
{

  int nbFibersTract_1 = this->curves_count ;
  int nbFibersTract_2 = nbFibersTract ;
  int nbPoints1 = this->pointsPerTrack[ 0 ] ;
  for ( int fiber = 1 ; fiber < nbFibersTract_2 ; fiber++ )
  {

    if ( this->pointsPerTrack[ fiber ] != nbPoints1 )
    {

      std::cout << "ERROR : in compareDisimilarityBundles( const float* "
                << "tractogramFibers, int nbFibersTract, int nbPoints ) the "
                << "number of points of all the fibers in tractogramFibers "
                << "has to be the same " << std::endl ;
      exit( 1 ) ;
    }

  }

  if ( nbPoints != nbPoints1 )
  {

    std::cout << "ERROR : in compareDisimilarityBundles( const float* "
              << "tractogramFibers, int nbFibersTract, int nbPoints ) the "
              << "number of points of the fibers must be the same "
              << std::endl ;
    exit( 1 ) ;

  }


  float disimilarity = compareDisimilarityBundles( this->matrixTracks,
                                                   tractogramFibers,
                                                   nbFibersTract_1,
                                                   nbFibersTract_2,
                                                   nbPoints1,
                                                   nbPoints ) ;

  return disimilarity ;

}


////////////////////////////////////////////////////////////////////////////////
double BundlesDataFormat::distanceBetweenBundles(
                                              const std::vector<float>& bundle1,
                                              const std::vector<float>& bundle2,
                                              int nbFibersBundle1,
                                              int nbFibersBundle2,
                                              int nbPointsTract_1,
                                              int nbPointsTract_2 )
{

  double disimilarity_1 = compareDisimilarityBundles( bundle1,
                                                      bundle2,
                                                      nbFibersBundle1,
                                                      nbFibersBundle2,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity_2 = compareDisimilarityBundles( bundle2,
                                                      bundle1,
                                                      nbFibersBundle2,
                                                      nbFibersBundle1,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity = 0 ;

  if ( disimilarity_1 > disimilarity_2 )
  {

    disimilarity = disimilarity_1 ;

  }
  else
  {

    disimilarity = disimilarity_2 ;

  }

  return disimilarity ;

}


//----------------------------------------------------------------------------//


double BundlesDataFormat::distanceBetweenBundles(
                                               const std::vector<float>& bundle,
                                               int nbFibersBundle,
                                               int nbPoints )
{

  const std::vector<float>& bundle1 = bundle ;
  int nbFibersBundle1 = nbFibersBundle ;
  int nbPointsTract_1 = nbPoints ;

  const std::vector<float>& bundle2 = this->matrixTracks ;
  int nbFibersBundle2 = this->curves_count ;
  int nbPointsTract_2 = this->pointsPerTrack[ 0 ] ;



  double disimilarity_1 = compareDisimilarityBundles( bundle1,
                                                      bundle2,
                                                      nbFibersBundle1,
                                                      nbFibersBundle2,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity_2 = compareDisimilarityBundles( bundle2,
                                                      bundle1,
                                                      nbFibersBundle2,
                                                      nbFibersBundle1,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity = 0 ;

  if ( disimilarity_1 > disimilarity_2 )
  {

    disimilarity = disimilarity_1 ;

  }
  else
  {

    disimilarity = disimilarity_2 ;

  }

  return disimilarity ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::computeNumberAdjacentFibersBundle1toBundle2(
                                      BundlesDataFormat& bundle1,
                                      BundlesDataFormat& bundle2,
                                      float threshold,
                                      std::vector<int>& nbAdjacentFibersBundle )
{

  if ( threshold <= 0 )
  {
    std::cout << "Warning : the threshold value in bundlesAdjacency method "
              << "must be greater than 0 but got " << threshold
              << "... Using default value of 10 " << std::endl ;

    threshold = 10 ;

  }

  const std::vector<float>& bundlesData1 = bundle1.matrixTracks ;
  int nbFibersBundle1 = bundle1.curves_count ;
  int nbPoints1 = bundle1.pointsPerTrack[ 0 ] ;

  std::vector<float>& bundleData2 = bundle2.matrixTracks ;
  int nbFibersBundle2 = bundle2.curves_count ;
  int nbPoints2 = bundle2.pointsPerTrack[ 0 ] ;

  if ( nbPoints2 != nbPoints1 )
  {

    std::cout << "ERROR : Problem in compareDisimilarityBundles, not the same "
              << "number of points per track in each bundle " << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPoints1 ;

  #pragma omp parallel for
  for ( int fiberIndex1 = 0 ; fiberIndex1 < nbFibersBundle1 ; fiberIndex1++ )
  {

    int64_t offsetTractogram_1 = fiberIndex1 * nbPoints * 3 ;

    // Searching the medial point of tractogram 1 fiber
    std::vector<float> medialPointTractFiber1( 3, 0 ) ;
    bundle1.computeMedialPointFiberWithDistance( fiberIndex1,
                                                 medialPointTractFiber1 ) ;

    float minDMDA = 1000 ;
    float dMDA = 0 ;
    nbAdjacentFibersBundle[ fiberIndex1 ] = 0 ;

    for ( int fiberIndex2 = 0 ; fiberIndex2 < nbFibersBundle2 ; fiberIndex2++ )
    {


      int64_t offsetTractogram_2 = fiberIndex2 * nbPoints * 3 ;

      // Searching the medial point of tractogram 2 fiber
      std::vector<float> medialPointTractFiber2( 3, 0 ) ;
      bundle2.computeMedialPointFiberWithDistance( fiberIndex2,
                                                   medialPointTractFiber2 ) ;


      // Computing translation from tract 1 fiber to tract 2 fiber
      std::vector<float> translation( 3, 0 ) ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        translation[ i ] = medialPointTractFiber1[ i ] -
                                                   medialPointTractFiber2[ i ] ;

      }

      float tmpDirectDistance = 0 ;
      float tmpFlippedDistance = 0 ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        tmpDirectDistance += pow(
                              bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                              ( bundleData2[ i + 3 * 0 + offsetTractogram_2 ] +
                                                       translation[ i ] ), 2 ) ;
        tmpFlippedDistance += pow(
                                bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                                ( bundleData2[ i + 3 * ( nbPoints - 1 ) +
                                offsetTractogram_2 ] + translation[ i ] ), 2 ) ;

      }

      bool isDirectSens = true ;
      if ( tmpFlippedDistance < tmpDirectDistance )
      {

        isDirectSens = false ;

      }

      // Computing MDA when distance between middle points is below threshold
      dMDA = 0 ;

      for ( int point = 0 ; point < nbPoints ; point++ )
      {

        int k = point ;
        if ( !isDirectSens )
        {

          k = nbPoints - point - 1 ;

        }

        float tmpMDA = 0 ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          tmpMDA += pow( bundlesData1[ i + 3 * point +
                         offsetTractogram_1 ] - ( bundleData2[
                         i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
                                                                           2 ) ;

        }

        tmpMDA = sqrt( tmpMDA ) ;

        if ( tmpMDA > dMDA )
        {

          dMDA = tmpMDA ;

        }

      }

      if ( dMDA < threshold )
      {

        nbAdjacentFibersBundle[ fiberIndex1 ] += 1 ;

      }

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::computeNumberAdjacentFibersBundle1toBundle2(
                                      const std::vector<float>& bundle1,
                                      const std::vector<float>& bundle2,
                                      int nbFibersBundle1,
                                      int nbFibersBundle2,
                                      int nbPoints,
                                      float threshold,
                                      std::vector<int>& nbAdjacentFibersBundle )
{

  if ( threshold <= 0 )
  {
    std::cout << "Warning : the threshold value in bundlesAdjacency method "
              << "must be greater than 0 but got " << threshold
              << "... Using default value of 10 " << std::endl ;

    threshold = 10 ;

  }

  if ( nbAdjacentFibersBundle.size() < nbFibersBundle1 )
  {

    nbAdjacentFibersBundle.resize( nbFibersBundle1, 0 ) ;

  }

  const std::vector<float>& bundlesData1 = bundle1 ;

  const std::vector<float>& bundleData2 = bundle2 ;

  #pragma omp parallel for
  for ( int fiberIndex1 = 0 ; fiberIndex1 < nbFibersBundle1 ; fiberIndex1++ )
  {

    int64_t offsetTractogram_1 = fiberIndex1 * nbPoints * 3 ;

    // Searching the medial point of tractogram 1 fiber
    std::vector<float> medialPointTractFiber1( 3, 0 ) ;
    this->computeMedialPointFiberWithDistance( bundle1, fiberIndex1, nbPoints,
                                                 medialPointTractFiber1 ) ;

    float minDMDA = 1000 ;
    float dMDA = 0 ;
    nbAdjacentFibersBundle[ fiberIndex1 ] = 0 ;

    for ( int fiberIndex2 = 0 ; fiberIndex2 < nbFibersBundle2 ; fiberIndex2++ )
    {


      int64_t offsetTractogram_2 = fiberIndex2 * nbPoints * 3 ;

      // Searching the medial point of tractogram 2 fiber
      std::vector<float> medialPointTractFiber2( 3, 0 ) ;
      this->computeMedialPointFiberWithDistance( bundle2, fiberIndex2, nbPoints,
                                                   medialPointTractFiber2 ) ;


      // Computing translation from tract 1 fiber to tract 2 fiber
      std::vector<float> translation( 3, 0 ) ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        translation[ i ] = medialPointTractFiber1[ i ] -
                                                   medialPointTractFiber2[ i ] ;

      }

      float tmpDirectDistance = 0 ;
      float tmpFlippedDistance = 0 ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        tmpDirectDistance += pow(
                              bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                              ( bundleData2[ i + 3 * 0 + offsetTractogram_2 ] +
                                                       translation[ i ] ), 2 ) ;
        tmpFlippedDistance += pow(
                                bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                                ( bundleData2[ i + 3 * ( nbPoints - 1 ) +
                                offsetTractogram_2 ] + translation[ i ] ), 2 ) ;

      }

      bool isDirectSens = true ;
      if ( tmpFlippedDistance < tmpDirectDistance )
      {

        isDirectSens = false ;

      }

      // Computing MDA when distance between middle points is below threshold
      dMDA = 0 ;

      for ( int point = 0 ; point < nbPoints ; point++ )
      {

        int k = point ;
        if ( !isDirectSens )
        {

          k = nbPoints - point - 1 ;

        }

        float tmpMDA = 0 ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          tmpMDA += pow( bundlesData1[ i + 3 * point +
                         offsetTractogram_1 ] - ( bundleData2[
                         i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
                                                                           2 ) ;

        }

        tmpMDA = sqrt( tmpMDA ) ;

        if ( tmpMDA > dMDA )
        {

          dMDA = tmpMDA ;

        }

      }

      if ( dMDA < threshold )
      {

        nbAdjacentFibersBundle[ fiberIndex1 ] += 1 ;

      }

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::coverageBundle1toBundle2( BundlesDataFormat& bundle1,
                                                   BundlesDataFormat& bundle2,
                                                   float threshold,
                                                   int verbose )
{

  int nbFibersBundle = bundle1.curves_count ;

  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersBundle1toBundle2( bundle1,
                                                     bundle2,
                                                     threshold,
                                                     nbAdjacentFibersBundle ) ;

  float coverage = 0 ;
  // #pragma omp parallel for reduction( + : coverage )
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersBundle[ i ] > 0 )
    {

      coverage += 1.0 ;

    }

  }

  coverage /= nbFibersBundle ;

  return coverage ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::overlapBundle1toBundle2( BundlesDataFormat& bundle1,
                                                  BundlesDataFormat& bundle2,
                                                  float threshold,
                                                  int verbose )
{

  int nbFibersBundle = bundle1.curves_count ;

  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersBundle1toBundle2( bundle1,
                                                     bundle2,
                                                     threshold,
                                                     nbAdjacentFibersBundle ) ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
  // #pragma omp parallel for reduction( + : nbAdjacentFibersRecognizedBundle, nbAdjacentFibersAtlasBundle )
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersBundle[ i ] > 0 )
    {

      nbAdjacentFibersRecognizedBundle += 1 ;
      nbAdjacentFibersAtlasBundle += nbAdjacentFibersBundle[ i ] ;

    }

  }

  if ( nbAdjacentFibersAtlasBundle == 0 )
  {

    if ( verbose > 1 )
    {

      std::cout << "Warning : there are not adjacent fibers of bundle 1 "
                << "to atlas, overlap is not defined \n" ;

    }

    return 0 ;

  }

  // Overlap is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::bundlesAdjacency( BundlesDataFormat& bundle1,
                                           BundlesDataFormat& bundle2,
                                           float threshold,
                                           int verbose )
{

  float coverageBundle1toBundle2 =
                           this->coverageBundle1toBundle2( bundle1,
                                                           bundle2,
                                                           threshold,
                                                           verbose ) ;

  float coverageBundle2toBundle1 =
                           this->coverageBundle1toBundle2( bundle2,
                                                           bundle1,
                                                           threshold,
                                                           verbose ) ;

  float bundlesAdjacencyMeasure = ( coverageBundle1toBundle2 +
                                    coverageBundle2toBundle1 ) / 2 ;

  return bundlesAdjacencyMeasure ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::computeAngleBetweenVectors(
                                             const std::vector<float>& vector1,
                                             const std::vector<float>& vector2 )
{

  float dotProduct = scalarProduct( vector1, vector2 ) ;
  float normVector1 = normVector( vector1 ) ;
  float normVector2 = normVector( vector2 ) ;

  float angle = acos( dotProduct / ( normVector1 * normVector2 ) ) *
                                                                  180 / M_PI ;

  return angle ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::computeAngleBetweenPlanes(
                                             const std::vector<float>& vector1,
                                             const std::vector<float>& vector2 )
{

  float angle = computeAngleBetweenVectors( vector1, vector2 ) ;

  if ( angle > 90 )
  {

    angle = 180 - angle ;

  }

  // If is angle = nan then angle = 0
  if ( angle != angle )
  {

    angle = 0 ;

  }

  return angle ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::computeAngleBetweenDirections(
                                             const std::vector<float>& vector1,
                                             const std::vector<float>& vector2 )
{

  float angle = computeAngleBetweenVectors( vector1, vector2 ) ;

  float vector1DotVector2 = scalarProduct( vector1, vector2 ) ;


  if ( vector1DotVector2 < 0 && angle < 90 )
  {

    angle = 180 - angle ;

  }

  // If is angle = nan then angle = 0
  if ( angle != angle )
  {

    angle = 0 ;

  }

  return angle ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::projectVectorToPlane(
                                        const std::vector<float>& normalToPlane,
                                        const std::vector<float>& inputVector,
                                        std::vector<float>& projectedVector )
{

  // Normalisation of normal vector
  float normalVectorNorm = sqrt( pow( normalToPlane[ 0 ], 2 ) +
                                 pow( normalToPlane[ 1 ], 2 ) +
                                 pow( normalToPlane[ 2 ], 2 ) ) ;

  std::vector<float> unitNormalVector( 3, 0 ) ;
  if ( normalVectorNorm == 1 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitNormalVector[ i ] = normalToPlane[ i ] ;

    }

  }
  else
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitNormalVector[ i ] = normalToPlane[ i ] / normalVectorNorm ;

    }

  }

  // Scalar product of normalToPlane by inputVector
  float unitNormalVectorDotInputVector = scalarProduct(
                                               unitNormalVector, inputVector ) ;

  // Projection
  for ( int i = 0 ; i < 3 ; i++ )
  {

    projectedVector[ i ] = inputVector[ i ] - unitNormalVectorDotInputVector *
                                                         unitNormalVector[ i ] ;

  }

}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::computeRotationMatrixFromVectorAndAngle(
                                            const std::vector<float>& vector,
                                            float angle, // in rad
                                            std::vector<float>& rotationMatrix )
{

  /////////////////////// Normalisation of normal vector ///////////////////////
  float vectorNorm = sqrt( pow( vector[ 0 ], 2 ) + pow( vector[ 1 ], 2 ) +
                                                       pow( vector[ 2 ], 2 ) ) ;

  std::vector<float> unitVector( 3, 0 ) ;
  if ( vectorNorm == 1 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitVector[ i ] = vector[ i ] ;

    }
  }
  else
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitVector[ i ] = vector[ i ] / vectorNorm ;

    }

  }

  float cosAngle = cos( angle ) ;
  float sinAngle = sin( angle ) ;

  float m00 = cosAngle + pow( unitVector[ 0 ], 2 ) * ( 1 - cosAngle ) ;
  float m01 = unitVector[ 0 ] * unitVector[ 1 ] * ( 1 - cosAngle ) -
                                                    unitVector[ 2 ] * sinAngle ;
  float m02 = unitVector[ 0 ] * unitVector[ 2 ] * ( 1 - cosAngle ) +
                                                    unitVector[ 1 ] * sinAngle ;
  float m10 = unitVector[ 0 ] * unitVector[ 1 ] * ( 1 - cosAngle ) +
                                                    unitVector[ 2 ] * sinAngle ;
  float m11 = cosAngle + pow( unitVector[ 1 ], 2 ) * ( 1 - cosAngle ) ;
  float m12 = unitVector[ 1 ] * unitVector[ 2 ] * ( 1 - cosAngle ) -
                                                    unitVector[ 0 ] * sinAngle ;
  float m20 = unitVector[ 0 ] * unitVector[ 2 ] * ( 1 - cosAngle ) -
                                                    unitVector[ 1 ] * sinAngle ;
  float m21 = unitVector[ 1 ] * unitVector[ 2 ] * ( 1 - cosAngle ) +
                                                    unitVector[ 0 ] * sinAngle ;
  float m22 = cosAngle + pow( unitVector[ 2 ], 2 ) * ( 1 - cosAngle ) ;



  //////////////////// Computing elements of rotationMatrix ////////////////////
  // ( 0, 0 )
  rotationMatrix[ 0 ] = m00 ;
  // ( 0, 1 )
  rotationMatrix[ 1 ] = m01 ;
  // ( 0, 2 )
  rotationMatrix[ 2 ] = m02 ;
  // ( 1, 0 )
  rotationMatrix[ 3 ] = m10 ;
  // ( 1, 1 )
  rotationMatrix[ 4 ] = m11 ;
  // ( 1, 2 )
  rotationMatrix[ 5 ] = m12 ;
  // ( 2, 0 )
  rotationMatrix[ 6 ] = m20 ;
  // ( 2, 1 )
  rotationMatrix[ 7 ] = m21 ;
  // ( 2, 2 )
  rotationMatrix[ 8 ] = m22 ;


}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::applyRotationMatrixToVector(
                                       const std::vector<float>& vector,
                                       const std::vector<float>& rotationMatrix,
                                       std::vector<float>& rotatedVector )
{

  for ( int i = 0 ; i < 3 ; i++ )
  {

    rotatedVector[ i ] = rotationMatrix[ 3 * i + 0 ] * vector[ 0 ] +
                         rotationMatrix[ 3 * i + 1 ] * vector[ 1 ] +
                         rotationMatrix[ 3 * i + 2 ] * vector[ 2 ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::applyRotationMatrixToFiber(
                                     const std::vector<float>& fiber,
                                     const std::vector<float>& rotationMatrix,
                                     const std::vector<float>& medialPointFiber,
                                     int nbPoints,
                                     std::vector<float>& rotatedFiber )
{

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = fiber[ 3 * point + i ] - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::applyRotationMatrixToFiber(
                                       BundlesDataFormat& inputBundlesData,
                                       const std::vector<float>& rotationMatrix,
                                       int fiberIndex,
                                       int nbPoints,
                                       std::vector<float>& rotatedFiber )
{

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // inputBundlesData.computeMedialPointFiberTractogram(
  //                                               fiberIndex, medialPointFiber ) ;
  inputBundlesData.computeMedialPointFiberWithDistance(
                                                fiberIndex, medialPointFiber ) ;

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = inputBundlesData[ offsetTractogram + 3 * point + i ]
                                                       - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::applyRotationMatrixToFiber(
                                     const std::vector<float>& tractogramFibers,
                                     const std::vector<float>& rotationMatrix,
                                     int fiberIndex,
                                     int nbPoints,
                                     std::vector<float>& rotatedFiber )
{

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( tractogramFibers,
  //                                    fiberIndex,
  //                                    nbPoints,
  //                                    medialPointFiber ) ;
  computeMedialPointFiberWithDistance( tractogramFibers,
                                       fiberIndex,
                                       nbPoints,
                                       medialPointFiber ) ;

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = tractogramFibers[ offsetTractogram + 3 * point + i ]
                                                       - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::applyRotationMatrixToFiber(
                                       const std::vector<float>& rotationMatrix,
                                       int fiberIndex,
                                       int nbPoints,
                                       std::vector<float>& rotatedFiber )
{

  applyRotationMatrixToFiber( *this, rotationMatrix, fiberIndex, nbPoints,
                                                                rotatedFiber ) ;


}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::applyRotationMatrixToFiber(
                                       const std::vector<float>& fiber,
                                       const std::vector<float>& rotationMatrix,
                                       int nbPoints,
                                       std::vector<float>& rotatedFiber )
{

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber, medialPointFiber ) ;
  computeMedialPointFiberWithDistance( fiber, medialPointFiber ) ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = fiber[ 3 * point + i ] - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::putFibersInSamePlane(
                                    const std::vector<float>& normalVector1,
                                    const std::vector<float>& normalVector2,
                                    const std::vector<float>& tractogramFibers2,
                                    int fiberIndex2,
                                    int nbPoints,
                                    std::vector<float>& fiber2ToPlane1,
                                    std::vector<float>& newNormalVectorFiber2 )
{

  float angleNormaleVector = computeAngleBetweenPlanes( normalVector1,
                                                        normalVector2 ) ;

  if ( angleNormaleVector < 1 )
  {

    int64_t offsetTractogram = 3 * nbPoints * fiberIndex2 ;

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        fiber2ToPlane1[ 3 * point + i ] = tractogramFibers2[ 3 * point + i +
                                                            offsetTractogram ] ;

      }

    }

    for ( int i = 0 ; i < 3 ; i++ )
    {

      newNormalVectorFiber2[ i ] = normalVector2[ i ] ;

    }

    return ;

  }

  angleNormaleVector *= M_PI / 180 ;

  std::vector<float> rotationAxisSamePlane( 3, 0 );
  crossProduct( normalVector1, normalVector2, rotationAxisSamePlane ) ;

  std::vector<float> rotationMatrixSamePlane( 9, 0 ) ;
  computeRotationMatrixFromVectorAndAngle( rotationAxisSamePlane,
                                           angleNormaleVector,
                                           rotationMatrixSamePlane ) ;

  applyRotationMatrixToFiber( tractogramFibers2,
                              rotationMatrixSamePlane,
                              fiberIndex2,
                              nbPoints,
                              fiber2ToPlane1 ) ;

  applyRotationMatrixToVector( normalVector2,
                               rotationMatrixSamePlane,
                               newNormalVectorFiber2 ) ;

  float residualAnglePlanes = computeAngleBetweenPlanes(
                                                       normalVector1,
                                                       newNormalVectorFiber2 ) ;

  if ( residualAnglePlanes > 1 )
  {

    std::vector<float> rotationMatrixSamePlaneTmp( 9, 0 );
    computeRotationMatrixFromVectorAndAngle( rotationAxisSamePlane,
                                             - angleNormaleVector,
                                             rotationMatrixSamePlaneTmp ) ;
    std::vector<float> fiber2ToPlane1Tmp( 3 * nbPoints, 0 ) ;
    applyRotationMatrixToFiber( tractogramFibers2,
                                rotationMatrixSamePlaneTmp,
                                fiberIndex2,
                                nbPoints,
                                fiber2ToPlane1Tmp ) ;

    std::vector<float> newNormalVectorTmp( 3, 0 ) ;

    applyRotationMatrixToVector( normalVector2, rotationMatrixSamePlaneTmp,
                                                          newNormalVectorTmp ) ;

    float residualAnglePlanesTmp = computeAngleBetweenPlanes(
                                                          normalVector1,
                                                          newNormalVectorTmp ) ;

    // The comparison f != f is only true if f = nan
    if ( residualAnglePlanesTmp < residualAnglePlanes ||
                              residualAnglePlanesTmp != residualAnglePlanesTmp )
    {

      fiber2ToPlane1 = fiber2ToPlane1Tmp ;
      newNormalVectorFiber2 = newNormalVectorTmp ;
      residualAnglePlanes = residualAnglePlanesTmp ;

    }

  }


  if ( residualAnglePlanes > 5 )
  {

    std::cout << "\nERROR : in BundlesDataFormat::putFibersInSamePlane "
              << "could not put fibers in same plan, got a final angle of "
              << residualAnglePlanes << std::endl ;

    exit( 1 ) ;

  }

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::putFibersInSamePlane(
                                     const std::vector<float>& normalVector1,
                                     const std::vector<float>& normalVector2,
                                     const std::vector<float>& fiber2,
                                     int nbPoints,
                                     std::vector<float>& fiber2ToPlane1,
                                     std::vector<float>& newNormalVectorFiber2 )
{

  int fiberIndex2 = 0 ;

  putFibersInSamePlane( normalVector1,
                        normalVector2,
                        fiber2,
                        fiberIndex2,
                        nbPoints,
                        fiber2ToPlane1,
                        newNormalVectorFiber2 ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::putFiberInPlaneXY(
                                     const std::vector<float>& normalVector,
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int nbPoints,
                                     std::vector<float>& fiberToPlaneXY,
                                     std::vector<float>& newNormalVectorFiber )
{

  std::vector<float> zAxisVector{ 0, 0, 1 } ;

  putFibersInSamePlane( zAxisVector,
                        normalVector,
                        tractogramFibers,
                        fiberIndex,
                        nbPoints,
                        fiberToPlaneXY,
                        newNormalVectorFiber ) ;

}

//----------------------------------------------------------------------------//
void BundlesDataFormat::putFiberInPlaneXY(
                                     const std::vector<float>& normalVector,
                                     const std::vector<float>& fiber,
                                     int nbPoints,
                                     std::vector<float>& fiberToPlaneXY,
                                     std::vector<float>& newNormalVectorFiber )
{

  std::vector<float> zAxisVector{ 0, 0, 1 } ;

  putFibersInSamePlane( zAxisVector,
                        normalVector,
                        fiber,
                        nbPoints,
                        fiberToPlaneXY,
                        newNormalVectorFiber ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::putFibersInSameDirection(
                                     const std::vector<float>& normalVector1,
                                     const std::vector<float>& normalVector2,
                                     const std::vector<float>& directionVector1,
                                     const std::vector<float>& fiber2,
                                     int nbPoints,
                                     std::vector<float>& fiber2ToDirection1 )
{

  std::vector<float> medialPointFiber2( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber2, medialPointFiber2 ) ;
  computeMedialPointFiberWithDistance( fiber2, medialPointFiber2 ) ;

  float angleBetweenPlanes = computeAngleBetweenPlanes( normalVector1,
                                                        normalVector2 ) ;
  if ( angleBetweenPlanes > 5 )
  {

    std::cout << "WARNING : in putFibersInSameDirection, fiber 1 and fiber 2 "
              << "MUST be in the same plane but got an angle between their "
              << "normal vectors of " << angleBetweenPlanes << ". The "
              << "computations will not be accurate." << std::endl ;

  }

  std::vector<float> projectedDirectionVectorFiber2( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber2,
                                         normalVector2,
                                         projectedDirectionVectorFiber2 ) ;

  float angleDirectionVectors = computeAngleBetweenDirections(
                                              directionVector1,
                                              projectedDirectionVectorFiber2 ) ;

  if ( angleDirectionVectors < 1 )
  {

    fiber2ToDirection1 = fiber2 ;
    return ;

  }

  angleDirectionVectors *= M_PI / 180 ;

  std::vector<float> rotationAxisSameDirection( 3, 0 ) ;
  crossProduct( directionVector1, projectedDirectionVectorFiber2,
                                                   rotationAxisSameDirection ) ;

  std::vector<float> rotationMatrixSameDirection( 9, 0 ) ;
  computeRotationMatrixFromVectorAndAngle( rotationAxisSameDirection,
                                           angleDirectionVectors,
                                           rotationMatrixSameDirection) ;

  applyRotationMatrixToFiber( fiber2,
                              rotationMatrixSameDirection,
                              nbPoints,
                              fiber2ToDirection1 ) ;


  std::vector<float> newDirectionVector( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber2ToDirection1, normalVector1,
                                                          newDirectionVector ) ;

  float residualAngleDirections = computeAngleBetweenDirections(
                                                          directionVector1,
                                                          newDirectionVector ) ;

  if ( residualAngleDirections > 3 )
  {

    std::vector<float> rotationMatrixSameDirectionTmp( 9, 0 ) ;
    computeRotationMatrixFromVectorAndAngle( rotationAxisSameDirection,
                                             - angleDirectionVectors,
                                             rotationMatrixSameDirectionTmp ) ;
    std::vector<float> fiber2ToDirection1Tmp( 3 * nbPoints, 0 ) ;
    applyRotationMatrixToFiber( fiber2,
                                rotationMatrixSameDirectionTmp,
                                nbPoints,
                                fiber2ToDirection1Tmp ) ;

    std::vector<float> newDirectionVectorTmp( 3, 0 ) ;
    computeDirectionVectorFiberTractogram( fiber2ToDirection1Tmp,
                                        normalVector1, newDirectionVectorTmp ) ;

    float residualAngleDirectionsTmp = computeAngleBetweenDirections(
                                                       directionVector1,
                                                       newDirectionVectorTmp ) ;

    // The comparison f != f is only true if f = nan
    if ( residualAngleDirectionsTmp < residualAngleDirections ||
                      residualAngleDirectionsTmp != residualAngleDirectionsTmp )
    {

      fiber2ToDirection1 = fiber2ToDirection1Tmp ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::registerFiber(
                                    const std::vector<float>& tractogramFibers2,
                                    const std::vector<float>& normalVector1,
                                    const std::vector<float>& normalVector2,
                                    const std::vector<float>& directionVector1,
                                    const std::vector<float>& medialPointFiber1,
                                    const std::vector<float>& medialPointFiber2,
                                    int fiberIndexTractogram2,
                                    int nbPoints,
                                    std::vector<float>& fiber2Tofiber1,
                                    std::vector<float>& newNormalVectorFiber2,
                                    int verbose )
{

  // Compute translation
  std::vector<float> translation( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointFiber1[ i ] - medialPointFiber2[ i ] ;

  }

  // Create translated fiber 2
  std::vector<float> translatedFiber2( 3 * nbPoints, 0 ) ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ;  i < 3 ; i++ )
    {

      int indexTmp = fiberIndexTractogram2 * nbPoints * 3 + 3 * point + i ;
      translatedFiber2[ 3 * point + i ] = tractogramFibers2[ indexTmp ] +
                                                              translation[ i ] ;

    }

  }

  // Computing roation matrix to put fiber 2 in same plane as fiber 1
  std::vector<float> samePlaneTranslatedFiber2( 3 * nbPoints, 0 ) ;
  putFibersInSamePlane( normalVector1,
                        normalVector2,
                        translatedFiber2,
                        nbPoints,
                        samePlaneTranslatedFiber2,
                        newNormalVectorFiber2 ) ;

  // Computing rotation matrix to put fiber 2 in same direction as fiber 1
  putFibersInSameDirection( normalVector1,
                            newNormalVectorFiber2,
                            directionVector1,
                            samePlaneTranslatedFiber2,
                            nbPoints,
                            fiber2Tofiber1 ) ;


  if ( verbose > 1 )
  {

    // Computing residual angles
    std::vector<float> directionMovedFiber2( 3, 0 ) ;
    computeDirectionVectorFiberTractogram( fiber2Tofiber1,
                                           newNormalVectorFiber2,
                                           directionMovedFiber2 ) ;

    float angleNormaleVectorMoved = computeAngleBetweenPlanes(
                                                       normalVector1,
                                                       newNormalVectorFiber2 ) ;

    float angleDirectionVectorsMoved = computeAngleBetweenDirections(
                                                        directionVector1,
                                                        directionMovedFiber2 ) ;

    std::cout << "Residual angles : \n   Angle between planes : "
              << angleNormaleVectorMoved << "\n   Angle between direction : "
              << angleDirectionVectorsMoved << std::endl ;

  }

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::registerFiber(
                                    const std::vector<float>& tractogramFibers1,
                                    const std::vector<float>& tractogramFibers2,
                                    int fiberIndexTractogram1,
                                    int fiberIndexTractogram2,
                                    int nbPoints,
                                    std::vector<float>& fiber2Tofiber1,
                                    std::vector<float>& newNormalVectorFiber2,
                                    int verbose )
{

  std::vector<float> medialPointFiber1( 3, 0 ) ;
  // computeMedialPointFiberTractogram( tractogramFibers1, fiberIndexTractogram1,
  //                                                nbPoints, medialPointFiber1 ) ;
  int medialPoint1Index = computeMedialPointFiberWithDistance( tractogramFibers1,
                          fiberIndexTractogram1, nbPoints, medialPointFiber1 ) ;
  std::vector<float> normalVector1( 3, 0 ) ;
  computeNormalVectorFiberTractogram( tractogramFibers1, medialPointFiber1,
                              fiberIndexTractogram1, nbPoints, normalVector1 ) ;
  // computeNormalVectorFiberTractogram( tractogramFibers1, medialPointFiber1,
  //          medialPoint1Index, fiberIndexTractogram1, nbPoints, normalVector1 ) ;
  std::vector<float> directionVector1( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( tractogramFibers1, medialPointFiber1,
                                         normalVector1, fiberIndexTractogram1,
                                                  nbPoints, directionVector1 ) ;

  std::vector<float> medialPointFiber2( 3, 0 ) ;
  // computeMedialPointFiberTractogram( tractogramFibers2, fiberIndexTractogram2,
  //                                                nbPoints, medialPointFiber2 ) ;
  int medialPoint2Index = computeMedialPointFiberWithDistance( tractogramFibers2,
                          fiberIndexTractogram2, nbPoints, medialPointFiber2 ) ;
  std::vector<float> normalVector2( 3, 0 ) ;
  computeNormalVectorFiberTractogram( tractogramFibers2, medialPointFiber2,
                              fiberIndexTractogram2, nbPoints, normalVector2 ) ;
  // computeNormalVectorFiberTractogram( tractogramFibers2, medialPointFiber2,
  //          medialPoint2Index, fiberIndexTractogram2, nbPoints, normalVector2 ) ;


  registerFiber( tractogramFibers2,
                 normalVector1,
                 normalVector2,
                 directionVector1,
                 medialPointFiber1,
                 medialPointFiber2,
                 fiberIndexTractogram2,
                 nbPoints,
                 fiber2Tofiber1,
                 newNormalVectorFiber2,
                 verbose ) ;


}
// -------------------------------------------------------------------------- //
void BundlesDataFormat::registerFiber(
                                    const std::vector<float>& fiber2,
                                    const std::vector<float>& normalVector1,
                                    const std::vector<float>& normalVector2,
                                    const std::vector<float>& directionVector1,
                                    const std::vector<float>& medialPointFiber1,
                                    const std::vector<float>& medialPointFiber2,
                                    int nbPoints,
                                    std::vector<float>& fiber2Tofiber1,
                                    std::vector<float>& newNormalVectorFiber2,
                                    int verbose )
{

  int fiberIndexTractogram2 = 0 ;

  registerFiber( fiber2,
                 normalVector1,
                 normalVector2,
                 directionVector1,
                 medialPointFiber1,
                 medialPointFiber2,
                 fiberIndexTractogram2,
                 nbPoints,
                 fiber2Tofiber1,
                 newNormalVectorFiber2,
                 verbose ) ;

}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::registerFiber(
                                      const std::vector<float>& fiber1,
                                      const std::vector<float>& fiber2,
                                      int nbPoints,
                                      std::vector<float>& fiber2Tofiber1,
                                      std::vector<float>& newNormalVectorFiber2,
                                      int verbose )
{

  std::vector<float> medialPointFiber1( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber1, medialPointFiber1 ) ;
  computeMedialPointFiberWithDistance( fiber1, medialPointFiber1 ) ;
  std::vector<float> normalVector1( 3, 0 ) ;
  computeNormalVectorFiberTractogram( fiber1, normalVector1 ) ;
  std::vector<float> directionVector1( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber1, normalVector1,
                                                            directionVector1 ) ;

  std::vector<float> medialPointFiber2( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber2, medialPointFiber2 ) ;
  computeMedialPointFiberWithDistance( fiber2, medialPointFiber2 ) ;
  std::vector<float> normalVector2( 3, 0 ) ;
  computeNormalVectorFiberTractogram( fiber2, normalVector2 ) ;

  registerFiber( fiber2,
                 normalVector1,
                 normalVector2,
                 directionVector1,
                 medialPointFiber1,
                 medialPointFiber2,
                 nbPoints,
                 fiber2Tofiber1,
                 newNormalVectorFiber2,
                 verbose ) ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::registerFiberToPlaneXYAndDirectionX(
                               const std::vector<float>& tractogramFibers,
                               const std::vector<float>& normalVector,
                               const std::vector<float>& medialPointFiber,
                               int fiberIndexTractogram,
                               int nbPoints,
                               std::vector<float>& fiberToPlaneXYAndDirectionX,
                               std::vector<float>& newNormalVectorFiber,
                               int verbose )
{

  std::vector<float> xAxisVector{ 1, 0, 0 } ;
  std::vector<float> zAxisVector{ 0, 0, 1 } ;
  std::vector<float> origin{ 0, 0, 0 } ;

  registerFiber( tractogramFibers,
                 zAxisVector,
                 normalVector,
                 xAxisVector,
                 origin,
                 medialPointFiber,
                 fiberIndexTractogram,
                 nbPoints,
                 fiberToPlaneXYAndDirectionX,
                 newNormalVectorFiber,
                 verbose ) ;


}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::registerFiberToPlaneXYAndDirectionX(
                               const std::vector<float>& fiber,
                               const std::vector<float>& normalVector,
                               const std::vector<float>& medialPointFiber,
                               int nbPoints,
                               std::vector<float>& fiberToPlaneXYAndDirectionX,
                               std::vector<float>& newNormalVectorFiber,
                               int verbose )
{

  int fiberIndexTractogram = 0 ;

  registerFiberToPlaneXYAndDirectionX( fiber,
                                       normalVector,
                                       medialPointFiber,
                                       fiberIndexTractogram,
                                       nbPoints,
                                       fiberToPlaneXYAndDirectionX,
                                       newNormalVectorFiber,
                                       verbose ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::getFiberFromTractogram(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int nbPoints,
                                     std::vector<float>& fiber )
{

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      fiber[ 3 * point + i ] = tractogramFibers[ 3 * nbPoints *
                                               fiberIndex + 3 * point + i ] ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::getFiberWithVectors(
                                       const std::vector<float>& fiber,
                                       const std::vector<float>& referenceFiber,
                                       int nbPoints,
                                       std::vector<float>& fiberWithVectors )
{

  std::vector<float> medialPointReferenceFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( referenceFiber,
  //                                                  medialPointReferenceFiber ) ;
  computeMedialPointFiberWithDistance( referenceFiber,
                                                   medialPointReferenceFiber ) ;

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber, medialPointFiber ) ;
  computeMedialPointFiberWithDistance( fiber, medialPointFiber ) ;

  std::vector<float> translation( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointReferenceFiber[ i ] - medialPointFiber[ i ] ;

  }

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      fiberWithVectors[ 3 * point + i ] = fiber[  3 * point + i ] +
                                                              translation[ i ] ;

    }

  }


  std::vector<float> normalVector( 3, 0 ) ;
  computeNormalVectorFiberTractogram( fiber, normalVector ) ;

  std::vector<float> directionVector( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber, normalVector,
                                                             directionVector ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {


    fiberWithVectors[ 3 * nbPoints + i ] = normalVector[ i ] +
                                      medialPointFiber[ i ] + translation[ i ] ;
    fiberWithVectors[ 3 * nbPoints + 3 + i ] = medialPointFiber[ i ] +
                                                              translation[ i ] ;

    fiberWithVectors[ 3 * nbPoints + 3 * 2 + i ] = directionVector[ i ] +
                                      medialPointFiber[ i ] + translation[ i ] ;
    fiberWithVectors[ 3 * nbPoints + 3 * 3 + i ] = medialPointFiber[ i ] +
                                                              translation[ i ] ;

  }


}

// -------------------------------------------------------------------------- //
void BundlesDataFormat::getFiberWithVectors( const std::vector<float>& fiber,
                                             int nbPoints,
                                             std::vector<float>& fiberWithVectors )
{

  getFiberWithVectors( fiber, fiber, nbPoints, fiberWithVectors ) ;

}





////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Optimization functions ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::getMatrixAndVectorForLeastSquares(
                             const std::vector<float>& tractogramFibers,
                             int fiberIndex,
                             int nbPoints,
                             Eigen::MatrixXf& AMatrix,
                             Eigen::VectorXf& BVector )
{

  int offsetTractogram = 3 * nbPoints * fiberIndex ;

  int nbRows = nbPoints;
  int nbColumns = 3 ;
  AMatrix.resize( nbRows, nbColumns ) ;
  BVector.resize( nbRows ) ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    AMatrix( point, 0 ) = tractogramFibers[ 3 * point + 0 + offsetTractogram ] ;
    AMatrix( point, 1 ) = tractogramFibers[ 3 * point + 1 + offsetTractogram ] ;
    AMatrix( point, 2 ) = 1 ;
    BVector( point ) = tractogramFibers[ 3 * point + 2 + offsetTractogram ] ;

  }


}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Algebra functions //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::vectorFromPoints( const std::vector<float>& point1,
                                          const std::vector<float>& point2,
                                          std::vector<float>& vector )
{

  for ( int i = 0 ; i < 3 ; i++ )
  {

    vector[ i ] = point1[ i ] - point2[ i ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::scalarProduct( const std::vector<float>& vector1,
                                        const std::vector<float>& vector2 )
{

  float result = 0 ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    result += vector1[ i ] * vector2[ i ] ;

  }

  return result ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesDataFormat::normVector( const std::vector<float>& vector )
{

  float norm = scalarProduct( vector, vector ) ;
  norm = sqrt( norm ) ;

  return norm ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::normalizeVector( std::vector<float>& vector )
{

  float norm = normVector( vector ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    vector[ i ] /= norm ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::translatePoint(
                             const std::vector<float>& point,
                             const std::vector<float>& unitaryTranslationVector,
                             float translationDistance,
                             std::vector<float>& translatedPoint )
{

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translatedPoint[ i ] = point[ i ] + translationDistance *
                                                 unitaryTranslationVector[ i ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::translatePoint(
                                    const std::vector<float>& point,
                                    const std::vector<float>& translationVector,
                                    std::vector<float>& translatedPoint )
{

  translatePoint( point, translationVector, 1.0, translatedPoint ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::crossProduct( const std::vector<float>& vector1,
                                      const std::vector<float>& vector2,
                                      std::vector<float>& result )
{

  result[ 0 ] = vector1[ 1 ] * vector2[ 2 ] - vector1[ 2 ] * vector2[ 1 ] ;
  result[ 1 ] = vector1[ 2 ] * vector2[ 0 ] - vector1[ 0 ] * vector2[ 2 ] ;
  result[ 2 ] = vector1[ 0 ] * vector2[ 1 ] - vector1[ 1 ] * vector2[ 0 ] ;

}



////////////////////////////////////////////////////////////////////////////////
int BundlesDataFormat::checkIfFiberPointCanBeResampled(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int point,
                                     int nbPoints )
{

  // Returns :
  //           * 0 if the fiber cannot be resample
  //           * 1 if using two consecutive points after the point
  //           * 2 if using two consecutive points before the point
  //           * 3 if using one point before and one point after the point


  if ( nbPoints < 4 )
  {

    std::cout << "ERROR in BundlesDataFormat::checkIfFiberPointCanBeResampled :"
              << " imposible to resample fiber with NaN and less than 4 points "
              << std::endl ;
    return 0 ;

  }

  int64_t offset = 3 * nbPoints * fiberIndex ;

  // Check if the fiber has the points to do the interpolation :
  // - If NaN in first point : two consecutive points after the first point
  // - If NaN is last point : two consecutive points before the last point
  // - Else : one point beore and one point after or two consecutive points
  //          before or after (execpt for second and before last points which
  //          can only be resampled with before and after point or two
  //          consecutive points after the point for the second point, and two
  //          consecutive points before the point for the before last point


  // Test if method 1 is possible
  bool isMethod1 = true ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( point < ( nbPoints - 2 ) )
    {

      if ( isnan( tractogramFibers[ offset + 3 * ( point + 1 ) + i ] ) ||
           isnan( tractogramFibers[ offset + 3 * ( point + 2 ) + i ] ) )
      {

        isMethod1 = false ;
        break ;

      }

    }
    else
    {

      isMethod1 = false ;
      break ;

    }

  }

  if ( isMethod1 )
  {

    return 1 ;

  }
  else
  {

    if ( point == 0 )
    {

      return 0 ;

    }

  }

  // Test if method 2 is possible
  bool isMethod2 = true ;
  for ( int i = 0 ; i < 3 ; i++ )
  {
    if ( point > 2 )
    {

      if ( isnan( tractogramFibers[ offset + 3 * ( point - 1 ) + i ] ) ||
           isnan( tractogramFibers[ offset + 3 * ( point - 2 ) + i ] ) )
      {

        isMethod2 = false ;
        break ;

      }

    }
    else
    {

      isMethod2 = false ;
      break ;

    }

  }

  if ( isMethod2 )
  {

    return 2 ;

  }
  else
  {

    if ( point == ( nbPoints - 1 ) )
    {

      return 0 ;

    }

  }

  // Test if method 3 is possible
  bool isMethod3 = true ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( point > 0 && point < ( nbPoints - 1 ) )
    {

      if ( isnan( tractogramFibers[ offset + 3 * ( point - 1 ) + i ] ) ||
           isnan( tractogramFibers[ offset + 3 * ( point + 1 ) + i ] ) )
      {

        isMethod3 = false ;
        break ;

      }

    }
    else
    {

      isMethod3 = false ;
      break ;

    }

  }

  if ( isMethod3 )
  {

    return 3 ;

  }
  else
  {

    return 0 ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesDataFormat::resampleFiberWithNan(
                                           std::vector<float>& tractogramFibers,
                                           int fiberIndex,
                                           int nbPoints,
                                           int verbose )
{

  if ( nbPoints < 4 )
  {

    std::cout << "ERROR in BundlesDataFormat::resampleFiberWithNan : "
              << "imposible to resample fiber with NaN and less than 4 points "
              << std::endl ;
    return ;

  }

  int64_t offset = 3 * nbPoints * fiberIndex ;

  float lengthFiberSegment = computeLengthFiber( tractogramFibers, fiberIndex,
                                                  nbPoints ) / (nbPoints - 1 ) ;

  // Check if the fiber has the points to do the interpolation :
  // - If NaN in first point : two consecutive points after the first point
  // - If NaN is last point : two consecutive points before the last point
  // - Else : one point before and one point after or two consecutive points
  //          before or after (execpt for second and before last points which
  //          can only be resampled with before and after point or two
  //          consecutive points after the point for the second point, and two
  //          consecutive points before the point for the before last point
  std::vector<int> nanPointIndices ;
  std::vector<int> isPointOkForInterpolation( nbPoints, 0 ) ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    // Checking if point is NaN
    for ( int i = 0 ; i < 3 ; i++ )
    {

      if ( isnan( tractogramFibers[ offset + 3 * point + i ] ) )
      {

        nanPointIndices.push_back( point ) ;
        break ;

      }

    }

    // Checking if points can be resampled
    isPointOkForInterpolation[ point ] = checkIfFiberPointCanBeResampled(
                                                               tractogramFibers,
                                                               fiberIndex,
                                                               point,
                                                               nbPoints ) ;

  }

  int nbNanPoints = nanPointIndices.size() ;
  bool isResamplingPossible = false ;
  for ( int point = 0 ; point < nbNanPoints ; point++ )
  {

    if ( isPointOkForInterpolation[ nanPointIndices[ point ] ] )
    {

      isResamplingPossible = true ;
      break ;

    }

  }

  if ( !isResamplingPossible )
  {

    std::cout << "ERROR in BundlesDataFormat::resampleFiberWithNan : fiber "
              << fiberIndex << " with NaN values cannot be resampled "
              << "setting values to 0 \n" ;

    for ( int _tmpPointIndex = 0 ; _tmpPointIndex < nbPoints ;
                                                              _tmpPointIndex++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        tractogramFibers[ offset + 3 * _tmpPointIndex + i ] =
                                              tractogramFibers[ offset + 3 *
                                                  ( _tmpPointIndex + 1 ) + i ] ;

      }

    }

    return ;

  }



  // Resampling
  int nbResampledPoints = 0 ;
  int nbIterations = 0 ;
  std::vector<bool> nanPointResampled( nbNanPoints ,false ) ;
  while ( nbResampledPoints < nbNanPoints && nbIterations <= nbNanPoints )
  {

    nbIterations++ ;

    for ( int nanPoint = 0 ; nanPoint < nbNanPoints ; nanPoint++ )
    {

      int nanPointIndex = nanPointIndices[ nanPoint ] ;
      if ( !nanPointResampled[ nanPoint ] )
      {

        int isInterpolationPossible = isPointOkForInterpolation[
                                                               nanPointIndex ] ;

        if (  isInterpolationPossible == 1 )
        {

          for ( int i = 0 ; i < 3 ; i++ )
          {

            tractogramFibers[ offset + 3 * nanPointIndex + i ] =
                                                  tractogramFibers[ offset + 3 *
                                                   ( nanPointIndex + 1 ) + i ] ;

            tractogramFibers[ offset + 3 * ( nanPointIndex + 1 ) + i ] =
                          ( tractogramFibers[ offset + 3 * nanPointIndex + i ] +
                            tractogramFibers[ offset + 3 * ( nanPointIndex + 2 )
                                                                   + i ] ) / 2 ;


          }

          nanPointResampled[ nanPoint ] = true ;
          nbResampledPoints++ ;


        }

        if (  isInterpolationPossible == 2 )
        {

          for ( int i = 0 ; i < 3 ; i++ )
          {

            tractogramFibers[ offset + 3 * nanPointIndex + i ] =
                                                tractogramFibers[ offset + 3 * (
                                                     nanPointIndex - 1 ) + i ] ;

            tractogramFibers[ offset + 3 * ( nanPointIndex - 1 ) + i ] =
                          ( tractogramFibers[ offset + 3 * nanPointIndex + i ] +
                            tractogramFibers[ offset + 3 * ( nanPointIndex - 2 )
                                                                   + i ] ) / 2 ;


          }

          nanPointResampled[ nanPoint ] = true ;
          nbResampledPoints++ ;

        }

        if (  isInterpolationPossible == 3 )
        {

          for ( int i = 0 ; i < 3 ; i++ )
          {

            tractogramFibers[ offset + 3 * nanPointIndex + i ] =
                  ( tractogramFibers[ offset + 3 * ( nanPointIndex + 0 ) + i ] +
                    tractogramFibers[ offset + 3 * ( nanPointIndex - 1 ) + i ] )
                                                                           / 2 ;


          }

          nanPointResampled[ nanPoint ] = true ;
          nbResampledPoints++ ;

        }

      }

    }

    // Checking for new points that can be interpolated
    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      isPointOkForInterpolation[ point ] = checkIfFiberPointCanBeResampled(
                                                               tractogramFibers,
                                                               fiberIndex,
                                                               point,
                                                               nbPoints ) ;

    }

  }

  if ( nbResampledPoints < nbNanPoints && verbose > 1 )
  {

    std::cout << "WARNING : Not all NaN points in fiber " << fiberIndex
              << " could be resampled " << std::endl ;

  }

}
