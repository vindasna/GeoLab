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

#include <Eigen/Core>
#include <LBFGS.h>

#include "./BMDDistance.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BMDDistance::BMDDistance( BundlesDataFormat bundle1,
                 BundlesDataFormat bundle2 )
{

  this->reference = bundle1 ;
  this->moving = bundle2 ;

}

// double BMDDistance::operator()( const Eigen::VectorXd& affineTransform,
//                                                      Eigen::VectorXd& gradient )
// {
//
//     float eps = 1e-3 ;
//
//     int nbPoints = reference.pointsPerTrack[ 0 ] ;
//
//     #pragma omp parallel for
//     // for ( int i = 0 ; i < 9 ; i++ )
//     for ( int i = 0 ; i < 12 ; i++ )
//     {
//
//       BundlesDataFormat movingTransformedPlus( moving ) ;
//       Eigen::VectorXd affineFiniteDifferencePlus( affineTransform ) ;
//       affineFiniteDifferencePlus[ i ] += eps ;
//       applyAffineToBundle( moving.matrixTracks,
//                            affineFiniteDifferencePlus,
//                            moving.curves_count,
//                            nbPoints,
//                            movingTransformedPlus.matrixTracks ) ;
//       float finitePlus = computeBMD( reference, movingTransformedPlus ) ;
//
//
//       BundlesDataFormat movingTransformedMinus( moving ) ;
//       Eigen::VectorXd affineFiniteDifferenceMinus( affineTransform ) ;
//       affineFiniteDifferenceMinus[ i ] -= eps ;
//       applyAffineToBundle( moving.matrixTracks,
//                            affineFiniteDifferenceMinus,
//                            moving.curves_count,
//                            nbPoints,
//                            movingTransformedMinus.matrixTracks ) ;
//       float finiteMinus = computeBMD( reference, movingTransformedMinus ) ;
//
//       gradient[ i ] = ( finitePlus - finiteMinus ) / ( 2 * eps ) ;
//
//     }
//
//     BundlesDataFormat movingTransformed( moving ) ;
//     applyAffineToBundle( moving.matrixTracks,
//                          affineTransform,
//                          moving.curves_count,
//                          nbPoints,
//                          movingTransformed.matrixTracks ) ;
//     float bmdDistance = computeBMD( reference, movingTransformed ) ;
//
//     // std::cout << "Gradients : " << gradient.transpose() << std::endl ;
//
//     return bmdDistance ;
//
// }
double BMDDistance::operator()( const Eigen::VectorXd& transformCoefficients,
                                                     Eigen::VectorXd& gradient )
{

   int nbTransformCoefficients = transformCoefficients.size() ;
   int nbGradientCoefficients = gradient.size() ;

   if ( nbTransformCoefficients != nbGradientCoefficients )
   {

     std::cout << "ERROR bmdDistance::operator() : the transform and the "
               << "gradient must have the same number of coefficients \n" ;
     exit( 1 ) ;

   }

   double eps = 1e-3 ;

   int nbPoints = reference.pointsPerTrack[ 0 ] ;

   BundlesDataFormat movingTransformed( moving ) ;
   if ( nbTransformCoefficients == 12 )
   {

     applyAffineToBundle( moving.matrixTracks,
                          transformCoefficients,
                          moving.curves_count,
                          nbPoints,
                          movingTransformed.matrixTracks ) ;

   }
   else if ( nbTransformCoefficients == 6 )
   {

     applyRigidToBundle( moving.matrixTracks,
                         transformCoefficients,
                         moving.curves_count,
                         nbPoints,
                         movingTransformed.matrixTracks ) ;

   }
   else
   {

     std::cout << "ERROR bmdDistance::operator() : the only transforms "
               << "supported are affine ( 12 coefficients ) or rigid ( 9 "
               << "coefficietns ) \n" ;
     exit( 1 ) ;

   }

   double bmdDistance = computeBMD( reference, movingTransformed ) ;


   // Compute gradient with foward finite difference
   bool continueLoop = true ;
   float nbConsecutiveZeroGradient = 0 ;
   while ( continueLoop )
   {

     #pragma omp parallel for
     for ( int i = 0 ; i < nbTransformCoefficients ; i++ )
     {

       BundlesDataFormat movingTransformedPlus( moving ) ;
       Eigen::VectorXd finiteDifferencePlus( transformCoefficients ) ;
       finiteDifferencePlus[ i ] += eps ;
       if ( nbTransformCoefficients == 12 )
       {

         applyAffineToBundle( moving.matrixTracks,
                              finiteDifferencePlus,
                              moving.curves_count,
                              nbPoints,
                              movingTransformedPlus.matrixTracks ) ;

       }
       else if ( nbTransformCoefficients == 6 )
       {

         applyRigidToBundle( moving.matrixTracks,
                             finiteDifferencePlus,
                             moving.curves_count,
                             nbPoints,
                             movingTransformedPlus.matrixTracks ) ;

       }
       else
       {

         std::cout << "ERROR bmdDistance::operator() : the only transforms "
                   << "supported are affine ( 12 coefficients ) or rigid ( 6 "
                   << "coefficietns ), got transform with "
                   << nbTransformCoefficients << " coefficients \n" ;
         exit( 1 ) ;

       }

       double finitePlus = computeBMD( reference, movingTransformedPlus ) ;

       gradient[ i ] = ( finitePlus - bmdDistance ) / ( eps ) ;

     }

     float normGradientFiniteDifference = normGradient( gradient ) ;
     if ( normGradientFiniteDifference != 0 || nbConsecutiveZeroGradient > 3 )
     {

       continueLoop = false ;

     }
     else
     {

       eps *= 10 ;
       nbConsecutiveZeroGradient++ ;

     }


   }


   return bmdDistance ;

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
float BMDDistance::normGradient( Eigen::VectorXd& gradient )
{

  int sizeGradientVector = gradient.size() ;

    float norm = 0 ;
    for ( int i = 0 ; i < sizeGradientVector ; i++ )
    {

      norm += pow( gradient[ i ], 2 ) ;

    }
    norm = sqrt( norm ) ;

    // std::cout << "Norm gradient : " << norm << std::endl ;

    return( norm ) ;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void BMDDistance::applyAffineToBundle(
                                      const std::vector<float>& bundle,
                                      const Eigen::VectorXd& affineCoefficients,
                                      int nbCurves,
                                      int nbPoints,
                                      std::vector<float>& outputBundle )
{

  // affineCoefficeints = [ a1, a2, ..., a12 ] is a vector form of the affine
  // matrix --->  |  a0   a1   a2   a3  |
  //              |  a4   a5   a6   a7  |
  //              |  a8   a9  a10  a11 |
  //              |  0    0    0    1   |

  #pragma omp parallel for
  for ( int fiber = 0 ; fiber < nbCurves ; fiber++ )
  {

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      int64_t index = 3 * nbPoints * fiber + 3 * point ;

      outputBundle[ index + 0 ] = bundle[ index + 0 ] * affineCoefficients[ 0 ]
                                + bundle[ index + 1 ] * affineCoefficients[ 1 ]
                                + bundle[ index + 2 ] * affineCoefficients[ 2 ]
                                + affineCoefficients[ 3 ] ;

      outputBundle[ index + 1 ] = bundle[ index + 0 ] * affineCoefficients[ 4 ]
                                + bundle[ index + 1 ] * affineCoefficients[ 5 ]
                                + bundle[ index + 2 ] * affineCoefficients[ 6 ]
                                + affineCoefficients[ 7 ] ;

      outputBundle[ index + 2 ] = bundle[ index + 0 ] * affineCoefficients[ 8 ]
                                + bundle[ index + 1 ] * affineCoefficients[ 9 ]
                                + bundle[ index + 2 ] * affineCoefficients[ 10 ]
                                + affineCoefficients[ 11 ] ;

    }

  }

}
void BMDDistance::applyAffineToBundle(
                                   const std::vector<float>& bundle,
                                   const std::vector<float>& affineCoefficients,
                                   int nbCurves,
                                   int nbPoints,
                                   std::vector<float>& outputBundle )
{

  int sizeAffineCoefficients = affineCoefficients.size() ;

  Eigen::VectorXd affineCoefficientsEigen = Eigen::VectorXd::Zero(
                                                      sizeAffineCoefficients ) ;
  for ( int i = 0 ; i < sizeAffineCoefficients ; i++ )
  {

    affineCoefficientsEigen[ i ] = affineCoefficients[ i ] ;

  }

  applyAffineToBundle( bundle,
                       affineCoefficientsEigen,
                       nbCurves,
                       nbPoints,
                       outputBundle ) ;

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void BMDDistance::applyRigidToBundle(
                                      const std::vector<float>& bundle,
                                      const Eigen::VectorXd& rigidCoefficients,
                                      int nbCurves,
                                      int nbPoints,
                                      std::vector<float>& outputBundle )
{

  // rigidCoefficients -> first 3 are rotation angles, last 3 are translation
  // coefficients

  std::vector<float> transformationMatrix( 16, 0 ) ;

  // Translation
  transformationMatrix[ 3 ] = rigidCoefficients[ 3 ] ;  // x
  transformationMatrix[ 7 ] = rigidCoefficients[ 4 ] ;  // y
  transformationMatrix[ 11 ] = rigidCoefficients[ 5 ] ; // z

  // Last line
  transformationMatrix[ 12 ] = 0 ;  // x
  transformationMatrix[ 13 ] = 0 ;  // y
  transformationMatrix[ 14 ] = 0 ; // z
  transformationMatrix[ 15 ] = 1 ; // z

  // Setting rotation to identity
  transformationMatrix[ 0 ] = 1 ;
  transformationMatrix[ 5 ] = 1 ;
  transformationMatrix[ 10 ] = 1 ;

  // X Rotation
  std::vector<float> transformationMatrix1( 16, 0 ) ;
  std::vector<float> rotationX( 16, 0 ) ;
  float angle = rigidCoefficients[ 0 ] ;
  // float angle = rigidCoefficients[ 0 ] * M_PI / 180 ;
  float cosAngle = cos( angle ) ;
  float sinAngle = sin( angle ) ;
  rotationX[ 0 ] = 1 ;
  rotationX[ 5 ] = cosAngle ;
  rotationX[ 6 ] = sinAngle ;
  rotationX[ 9 ] = -sinAngle ;
  rotationX[ 10 ] = cosAngle ;
  rotationX[ 15 ] = 1 ;
  matrix1DotMatrix2( transformationMatrix,
                     rotationX,
                     transformationMatrix1 ) ;

  // Y Rotation
  std::vector<float> transformationMatrix2( 16, 0 ) ;
  std::vector<float> rotationY( 16, 0 ) ;
  angle = rigidCoefficients[ 1 ] ;
  // angle = rigidCoefficients[ 1 ] * M_PI / 180 ;
  cosAngle = cos( angle ) ;
  sinAngle = sin( angle ) ;
  rotationY[ 0 ] = cosAngle ;
  rotationY[ 2 ] = -sinAngle ;
  rotationY[ 5 ] = 1 ;
  rotationY[ 8 ] = sinAngle ;
  rotationY[ 10 ] = cosAngle ;
  rotationY[ 15 ] = 1 ;
  matrix1DotMatrix2( transformationMatrix1,
                     rotationY,
                     transformationMatrix2 ) ;

  // Z Rotation
  std::vector<float> transformationMatrix3( 16, 0 ) ;
  std::vector<float> rotationZ( 16, 0 ) ;
  angle = rigidCoefficients[ 2 ] ;
  // angle = rigidCoefficients[ 2 ] * M_PI / 180 ;
  cosAngle = cos( angle ) ;
  sinAngle = sin( angle ) ;
  rotationZ[ 0 ] = cosAngle ;
  rotationZ[ 1 ] = -sinAngle ;
  rotationZ[ 4 ] = sinAngle ;
  rotationZ[ 5 ] = cosAngle ;
  rotationZ[ 10 ] = 1 ;
  rotationZ[ 15 ] = 1 ;
  matrix1DotMatrix2( transformationMatrix2,
                     rotationZ,
                     transformationMatrix3 ) ;


  // transformationMatrix3 = [ a0, a1, ..., a12 ] is a vector form of the
  // matrix --->  |  a0   a1   a2   a3  |
  //              |  a4   a5   a6   a7  |
  //              |  a8   a9   a10  a11 |
  //              |  0    0    0    1   |

  #pragma omp parallel for
  for ( int fiber = 0 ; fiber < nbCurves ; fiber++ )
  {

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      int64_t index = 3 * nbPoints * fiber + 3 * point ;

      outputBundle[ index + 0 ] =
                                bundle[ index + 0 ] * transformationMatrix3[ 0 ]
                              + bundle[ index + 1 ] * transformationMatrix3[ 1 ]
                              + bundle[ index + 2 ] * transformationMatrix3[ 2 ]
                              + transformationMatrix3[ 3 ] ;

      outputBundle[ index + 1 ] =
                                bundle[ index + 0 ] * transformationMatrix3[ 4 ]
                              + bundle[ index + 1 ] * transformationMatrix3[ 5 ]
                              + bundle[ index + 2 ] * transformationMatrix3[ 6 ]
                              + transformationMatrix3[ 7 ] ;

      outputBundle[ index + 2 ] =
                               bundle[ index + 0 ] * transformationMatrix3[ 8 ]
                             + bundle[ index + 1 ] * transformationMatrix3[ 9 ]
                             + bundle[ index + 2 ] * transformationMatrix3[ 10 ]
                             + transformationMatrix3[ 11 ] ;

    }

  }

}

void BMDDistance::applyRigidToBundle(
                                   const std::vector<float>& bundle,
                                   const std::vector<float>& rigidCoefficients,
                                   int nbCurves,
                                   int nbPoints,
                                   std::vector<float>& outputBundle )
{

  int sizeAffineCoefficients = rigidCoefficients.size() ;

  Eigen::VectorXd rigidCoefficientsEigen = Eigen::VectorXd::Zero(
                                                      sizeAffineCoefficients ) ;
  for ( int i = 0 ; i < sizeAffineCoefficients ; i++ )
  {

    rigidCoefficientsEigen[ i ] = rigidCoefficients[ i ] ;

  }

  applyRigidToBundle( bundle,
                      rigidCoefficientsEigen,
                      nbCurves,
                      nbPoints,
                      outputBundle ) ;

}
////////////////////////////////////////////////////////////////////////////////
double BMDDistance::computeMDF( const std::vector<float>& tractogramFibers1,
                               const std::vector<float>& tractogramFibers2,
                               const std::vector<float>& medialPointFiber1,
                               const std::vector<float>& medialPointFiber2,
                               int fiberIndex1,
                               int fiberIndex2,
                               int nbPoints )
{

  int offsetTractogram1 = 3 * nbPoints * fiberIndex1 ;
  int offsetTractogram2 = 3 * nbPoints * fiberIndex2 ;

  float distance1 = 0 ;
  float distance2 = 0 ;

  std::vector<float> translation( 3, 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointFiber1[ i ] - medialPointFiber2[ i ] ;

    distance1 += pow( tractogramFibers2[ offsetTractogram2 + 3 * 0 + i ] +
                      translation[ i ] - tractogramFibers1[ offsetTractogram1 +
                                                             3 * 0 + i ] , 2 ) ;
    distance2 += pow( tractogramFibers2[ offsetTractogram2 + 3 * 0 + i ] +
                      translation[ i ] - tractogramFibers1[ offsetTractogram1 +
                                              3 * ( nbPoints - 1 ) + i ] , 2 ) ;
  }

  bool isDirectSens = true ;
  if ( distance1 > distance2 )
  {

    isDirectSens = false ;

  }

  float mdfDistance = 0 ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    int k = point ;
    if ( !isDirectSens )
    {

      k = nbPoints - point - 1 ;

    }

    float tmpDistance = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      tmpDistance += pow( tractogramFibers2[ offsetTractogram2 + 3 * point + i ]
                     - tractogramFibers1[ offsetTractogram1 + 3 * k + i ], 2 ) ;

    }

    tmpDistance = sqrt( tmpDistance ) ;

    mdfDistance += tmpDistance ;

  }

  mdfDistance /= nbPoints ;

  return( mdfDistance ) ;

}
////////////////////////////////////////////////////////////////////////////////
double BMDDistance::computeBMD( BundlesDataFormat& bundle1,
                               BundlesDataFormat& bundle2 )
{

  int nbPoints1 = bundle1.pointsPerTrack[ 0 ] ;
  int nbCurves1 = bundle1.curves_count ;
  for ( int i = 0 ; i < nbCurves1 ; i++ )
  {

    if ( nbPoints1 != bundle1.pointsPerTrack[ i ] )
    {

      std::cout << "ERROR in computeBMD : the number of points for each fiber "
      << "in bundle1 has to be the same " << std::endl;
      exit( 1 ) ;

    }

  }


  int nbPoints2 = bundle2.pointsPerTrack[ 0 ] ;
  int nbCurves2= bundle2.curves_count ;
  for ( int i = 0 ; i < nbCurves2 ; i++ )
  {

    if ( nbPoints2 != bundle2.pointsPerTrack[ i ] )
    {

      std::cout << "ERROR in computeBMD : the number of points for each fiber "
      << "in bundle2 has to be the same " << std::endl;
      exit( 1 ) ;

    }

  }


  if ( nbPoints1 != nbPoints2 )
  {

    std::cout << "ERROR in computeBMD : the number of points per fiber in "
    << "each bundle has to be the same " << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPoints1 ;

  // Computing medial points
  std::vector<float> zeroVector( 3, 0 ) ;
  std::vector< std::vector<float> > medialPointsFiber1( nbCurves1,
                                                                  zeroVector ) ;
  #pragma omp parallel for
  for ( int fiber = 0 ; fiber < nbCurves1 ; fiber++ )
  {

    std::vector<float>& medialPoint = medialPointsFiber1[ fiber ] ;

    bundle1.computeMedialPointFiberWithDistance( bundle1.matrixTracks,
                                                 fiber,
                                                 nbPoints,
                                                 medialPoint ) ;

  }

  std::vector< std::vector<float> > medialPointsFiber2( nbCurves2,
                                                                  zeroVector ) ;
  #pragma omp parallel for
  for ( int fiber = 0 ; fiber < nbCurves2 ; fiber++ )
  {

    std::vector<float>& medialPoint = medialPointsFiber2[ fiber ] ;

    bundle2.computeMedialPointFiberWithDistance( bundle2.matrixTracks,
                                                 fiber,
                                                 nbPoints,
                                                 medialPoint ) ;

  }

  // Distance bundle1 to bundle2
  std::vector<float> minBundle1( nbCurves1, 100000 ) ;
  #pragma omp parallel for
  for ( int fiber1 = 0 ; fiber1 < nbCurves1 ; fiber1++ )
  {

    std::vector<float> medialPointFiber1 = medialPointsFiber1[ fiber1 ] ;

    for ( int fiber2 = 0 ; fiber2 < nbCurves2 ; fiber2++ )
    {

      std::vector<float> medialPointFiber2 = medialPointsFiber2[ fiber2 ] ;

      double distance = computeMDF( bundle1.matrixTracks,
                                   bundle2.matrixTracks,
                                   medialPointFiber1,
                                   medialPointFiber2,
                                   fiber1,
                                   fiber2,
                                   nbPoints ) ;

      if ( distance < minBundle1[ fiber1 ] )
      {

        minBundle1[ fiber1 ] = distance ;

      }


    }


  }

  float minDistanceBundle1 = 0 ;
  for ( int fiber = 0 ; fiber < nbCurves1 ; fiber++ )
  {

    minDistanceBundle1 += minBundle1[ fiber ] ;

  }
  minDistanceBundle1 /= nbCurves1 ;



  // Distance bundle2 to bundle1
  std::vector<float> minBundle2( nbCurves2, 100000 ) ;
  #pragma omp parallel for
  for ( int fiber2 = 0 ; fiber2 < nbCurves2 ; fiber2++ )
  {

    std::vector<float> medialPointFiber2 = medialPointsFiber2[ fiber2 ] ;

    for ( int fiber1 = 0 ; fiber1 < nbCurves1 ; fiber1++ )
    {

      std::vector<float> medialPointFiber1 = medialPointsFiber1[ fiber1 ] ;

      double distance = computeMDF( bundle1.matrixTracks,
                                   bundle2.matrixTracks,
                                   medialPointFiber1,
                                   medialPointFiber2,
                                   fiber1,
                                   fiber2,
                                   nbPoints ) ;

      if ( distance < minBundle2[ fiber2 ] )
      {

        minBundle2[ fiber2 ] = distance ;

      }


    }


  }

  float minDistanceBundle2 = 0 ;
  for ( int fiber = 0 ; fiber < nbCurves2 ; fiber++ )
  {

    minDistanceBundle2 += minBundle2[ fiber ] ;

  }
  minDistanceBundle2 /= nbCurves2 ;

  // Computing BMD
  double bmdDistance = pow( minDistanceBundle1 + minDistanceBundle2, 2 ) / 4 ;


  return bmdDistance ;

}

////////////////////////////////////////////////////////////////////////////////
void BMDDistance::matrix1DotMatrix2( const std::vector<float>& matrix1,
                                     const std::vector<float>& matrix2,
                                     std::vector<float>& outputMatrix )
{

  for ( int i = 0 ; i < 4 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++ )
    {

      outputMatrix[ 4 * i + j ] = matrix1[ 4 * i + 0 ] * matrix2[ 4 * 0 + j ] +
                                  matrix1[ 4 * i + 1 ] * matrix2[ 4 * 1 + j ] +
                                  matrix1[ 4 * i + 2 ] * matrix2[ 4 * 2 + j ] +
                                  matrix1[ 4 * i + 3 ] * matrix2[ 4 * 3 + j ] ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BMDDistance::matrixMDotVector( const std::vector<float>& matrix,
                                    const std::vector<float>& vector,
                                    std::vector<float>& outputVector )
{

  for ( int i = 0 ; i < 4 ; i++ )
  {

    outputVector[ i ] = matrix[ 4 * i + 0 ] * vector[ 0 ] +
                        matrix[ 4 * i + 1 ] * vector[ 1 ] +
                        matrix[ 4 * i + 2 ] * vector[ 2 ] +
                        matrix[ 4 * i + 3 ] * vector[ 3 ] ;

  }

}
