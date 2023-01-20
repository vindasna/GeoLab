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

#include "applyTransformBundle.h"
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




//----------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////////
/////////////////////// Function to read affine transform //////////////////////
////////////////////////////////////////////////////////////////////////////////
void readTransform( const std::string& transformFilename,
                    std::vector<float>& affineTransform )
{

  if ( endswith( transformFilename, ".mat" ) )
  {

    std::fstream file ;
    file.open( transformFilename, std::ios::in ) ;
    if ( file.is_open() )
    {

      std::string line ;
      int i = 0 ;
      while ( std::getline( file, line ) )
      {

        std::stringstream ss( line ) ;

        affineTransform[ i ] = std::stof( ss.str() ) ;
        i++ ;

      }

      file.close() ;

    }

  }
  else
  {

    std::cout << "ERROR in transformation : The only format supported is .mat"
              << std::endl ;
    exit( 1 ) ;

  }


}


////////////////////////////////////////////////////////////////////////////////
/////////////////// Function to apply transform to tractogram //////////////////
////////////////////////////////////////////////////////////////////////////////
void applyAffineToBundle( const std::vector<float>& bundle,
                          const std::vector<float>& affineCoefficients,
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



////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Compute matrix product ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<float>> matrixProduct(
                              const std::vector<std::vector<float>>& matrix1,
                              const std::vector<std::vector<float>>& matrix2 )
{

  int linSize1 = matrix1.size() ;

  if ( linSize1 < 1 )
  {

    std::cout << "ERROR in matrixProduct : matrix1 must be non empty"
                                                                  << std::endl ;
    exit( 1 ) ;

  }

  int colSize1 = matrix1[ 0 ].size() ;
  for ( int i = 1 ; i < linSize1 ; i++ )
  {

    if ( matrix1[ i ].size() != colSize1 )
    {

      std::cout << "ERROR in matrixProduct : all lines in matrix1 must have "
                                    "the same number of elements" << std::endl ;
      exit( 1 ) ;


    }

  }

  int linSize2 = matrix2.size() ;
  if ( linSize2 != colSize1 )
  {

    std::cout << "ERROR in matrixProduct : matrix2 must have the same number of"
              << " lines as the number of columns of matrix1" << std::endl ;
    exit( 1 ) ;

  }

  int colSize2 = matrix2[ 0 ].size() ;
  for ( int i = 1 ; i < colSize2 ; i++ )
  {

    if ( matrix1[ i ].size() != colSize2 )
    {

      std::cout << "ERROR in matrixProduct : all lines in matrix2 must have "
                                    "the same number of elements" << std::endl ;
      exit( 1 ) ;


    }

  }

  std::vector<std::vector<float>> outMatrix( linSize1,
                                           std::vector<float>( colSize2, 0 ) ) ;

  for ( int i = 0 ; i < linSize1 ; i++ )
  {

    for ( int j = 0 ; j < colSize1 ; j++ )
    {

      for ( int k = 0 ; k < colSize2 ; k++ )
      {

        outMatrix[ i ][ k ] += matrix1[ i ][ j ] * matrix2[ j ][ k ] ;

      }

    }

  }

  return( outMatrix ) ;

}

////////////////////////////////////////////////////////////////////////////////
//////////////////////// Compute adjacent of 4x4 matrix ////////////////////////
////////////////////////////////////////////////////////////////////////////////
void computeAdjugateMatrix4x4( const std::vector<std::vector<float>>& matrix,
                               std::vector<std::vector<float>>& adjugate )
{

  adjugate.resize( 4, std::vector<float>( 4 ) ) ;

  // Sanity checks
  if ( matrix.size() != 4 )
  {

    std::cout << "ERRROR in computeAdjugateMatrix4x4 : The input matrix must be "
              << "a 4x4 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 4 )
    {

      std::cout << "ERRROR in computeAdjugateMatrix4x4 : The input matrix must be "
                << "a  matrix" << std::endl ;
      exit( 1 ) ;

    }

  }


  // Computing adjugate matrix
  #pragma omp parallel for
  for ( int i = 0 ; i < 4 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++ )
    {

      std::vector<std::vector<float>> minorMatrix ;
      minorMatrix.resize( 3, std::vector<float>( 3 ) ) ;

      int lineIndex = 0 ;
      int columnIndex = 0 ;
      for ( int n = 0 ; n < 4 ; n++ )
      {

        for ( int p = 0 ; p < 4 ; p++ )
        {

          if ( n != j && p != i )
          {

            int tmpN = n ;
            if ( n > j )
            {

              tmpN  = n - 1 ;

            }

            int tmpP = p ;
            if ( p > i )
            {

              tmpP  = p - 1 ;

            }

            minorMatrix[ tmpN ][ tmpP ] = matrix[ n ][ p ] ;

          }

        }

      }

      float sign = 0 ;
      if ( ( ( i + j ) % 2 == 0 ) )
      {

        sign = 1 ;

      }
      else
      {

        sign = -1 ;

      }

      adjugate[ i ][ j ] = computeDeterminantMatrix( minorMatrix ) ;
      adjugate[ i ][ j ] *= sign ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////// Compute determinant  matrix //////////////////////////
////////////////////////////////////////////////////////////////////////////////
float computeDeterminantMatrix( const  std::vector<std::vector<float>>& matrix )
{

  float determinant = 0 ;

  // Sanity checks
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != matrix.size() )
    {

      std::cout << "ERRROR in computeDeterminantMatrix : The input matrix "
                << "lust be square " << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( matrix.size() == 2 )
  {

    determinant = matrix[ 0 ][ 0 ] * matrix[ 1 ][ 1 ] -
                                           matrix[ 1 ][ 0 ] * matrix[ 0 ][ 1 ] ;

    return determinant ;


  }

  // Compute determinant
  int indexForCofactor = 0 ;
  // #pragma omp parallel for // Race condition
  for ( int i = 0 ; i < matrix.size() ; i++ )
  {

    for ( int j = 0 ; j < matrix.size() ; j++ )
    {

      if ( j != indexForCofactor )
      {

        continue ;

      }

      std::vector<std::vector<float>> minorMatrix ;
      minorMatrix.resize( matrix.size() - 1,
                                  std::vector<float>( matrix.size() - 1, 0 ) ) ;


      for ( int n = 0 ; n < matrix.size() ; n++ )
      {

        for ( int p = 0 ; p < matrix.size() ; p++ )
        {

          if ( n != i && p != j )
          {

            int tmpN = n ;
            if ( n > i )
            {

              tmpN  = n - 1 ;

            }

            int tmpP = p ;
            if ( p > j )
            {

              tmpP  = p - 1 ;

            }

            minorMatrix[ tmpN ][ tmpP ] = matrix[ n ][ p ] ;

          }

        }

      }

      float tmpSign = -1 ;
      if ( ( indexForCofactor + i ) % 2 == 0 )
      {

        tmpSign = 1 ;

      }

      // determinant += ( tmpSign ) * matrix[ indexForCofactor ][ j ] *
      //                               computeDeterminantMatrix3x3( minorMatrix ) ;
      determinant += ( tmpSign ) * matrix[ i ][ indexForCofactor ] *
                                       computeDeterminantMatrix( minorMatrix ) ;

    }

  }

  return determinant ;

}



////////////////////////////////////////////////////////////////////////////////
///////////////////////// Compute inverse 4x4 matrix ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
void computeInverseMatrix4x4( const std::vector<std::vector<float>>& matrix,
                              std::vector<std::vector<float>>& inverseMatrix )
{

  // Sanity checks
  if ( matrix.size() != 4 )
  {

    std::cout << "ERRROR in computeInverseMatrix4x4 : The input matrix must be "
              << "a 4x4 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 4 )
    {

      std::cout << "ERRROR in computeInverseMatrix4x4 : The input matrix must "
                << "be a 4x4 matrix" << std::endl ;
      exit( 1 ) ;

    }

  }

  // float determinant = computeDeterminantMatrix4x4( matrix ) ;
  float determinant = computeDeterminantMatrix( matrix ) ;
  if ( abs( determinant ) < 1e-30 )
  {

    std::cout << "ERROR in computeInverseMatrix4x4 : determinant is 0, singular"
              << " matrix are not invertible" << std::endl ;
    exit( 1 ) ;

  }

  std::vector<std::vector<float>> adjugateMatrix ;
  computeAdjugateMatrix4x4( matrix, adjugateMatrix ) ;

  inverseMatrix.resize( 4, std::vector<float>( 4 ) ) ;

  for ( int i = 0 ; i < 4 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++ )
    {

      inverseMatrix[ i ][ j ] = 1 / determinant * adjugateMatrix[ i ][ j ] ;

    }

  }

}



////////////////////////////////////////////////////////////////////////////////
/////////////////// Transform affineTransform into matrix //////////////////////
////////////////////////////////////////////////////////////////////////////////
void affineTransformToMatrix( const std::vector<float>& inAffineTransform,
                              std::vector<std::vector<float>>& outMatrix )
{

  outMatrix.resize( 4, std::vector<float>( 4, 0 ) ) ;
  outMatrix[ 3 ][ 3 ] = 1 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++ )
    {

      outMatrix[ i ][ j ] = inAffineTransform[ 4 * i + j ] ;


    }

  }


}


////////////////////////////////////////////////////////////////////////////////
/////////////////// Transform matrix into affineTransform //////////////////////
////////////////////////////////////////////////////////////////////////////////
void matrixToAffineTransform( const std::vector<std::vector<float>>& inMatrix,
                              std::vector<float>& outAffineTransform )
{

  outAffineTransform.resize( 12 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++ )
    {

      outAffineTransform[ 4 * i + j ] = inMatrix[ i ][ j ] ;


    }

  }


}




////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  int index_input_bundle, index_transform, index_output_bundle, index_inverse,
                                                     index_verbose, index_help ;

  index_input_bundle = getFlagPosition( argc, argv, "-i") ;
  index_transform = getFlagPosition( argc, argv, "-t") ;
  index_output_bundle = getFlagPosition( argc, argv, "-o") ;
  index_inverse = getFlagPosition( argc, argv, "-inv") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to apply affine transform to tractogram : \n"
              << "-i : Moving bundle \n"
              << "-t : Transform \n"
              << "-o : Output bundle \n"
              << "[-inv ] : apply inverse of input transform \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input_bundle )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_transform )
  {

    std::cout << "-t argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output_bundle )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_inverse )
  {

    isInverse = true ;

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

  std::string movingFilename( argv[ index_input_bundle + 1 ] ) ;
  char lastChar = movingFilename[ movingFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    movingFilename = movingFilename.substr( 0, movingFilename.size() - 1 ) ;

  }

  std::string transformFilename( argv[ index_transform + 1 ] ) ;
  lastChar = transformFilename[ transformFilename.size()  - 1 ] ;
  if ( lastChar == '/' )
  {

    transformFilename = transformFilename.substr( 0,
                                                transformFilename.size() - 1 ) ;

  }

  std::string outputFilename( argv[ index_output_bundle + 1 ] ) ;
  lastChar = outputFilename[ outputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

  }


  //   xxxxxxxxxxxxxxxxxxxxxxx Reading moving bundle xxxxxxxxxxxxxxxxxxxxxx   //

  std::string movingBundlesFilename ;
  std::string movingBundlesDataFilename ;
  bool isBundlesFormat = false ;
  bool isTrk = false ;
  bool isTck = false ;

  if ( endswith( movingFilename, ".bundles" ) ||
                                    endswith( movingFilename, ".bundlesdata" ) )
  {

    movingBundlesDataFilename = replaceExtension( movingFilename,
                                                              ".bundlesdata" ) ;

    movingBundlesFilename = replaceExtension( movingFilename, ".bundles" ) ;

    isBundlesFormat = true ;
    format = ".bundles" ;

  }
  else if ( endswith( movingFilename, ".trk" )  )
  {

    movingBundlesDataFilename = movingFilename ;

    movingBundlesFilename = replaceExtension( movingFilename, ".minf" ) ;

    isTrk = true ;
    format = ".trk" ;

  }
  else if ( endswith( movingFilename, ".tck" )  )
  {

    movingBundlesDataFilename = movingFilename ;

    movingBundlesFilename = replaceExtension( movingFilename, ".minf" ) ;

    isTck = true ;
    format = ".tck" ;

  }
  else
  {

    std::cout << "The only tractogram format supported are .bundles/.trk/.tck"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesData movingBundlesData( movingBundlesDataFilename.c_str() ) ;
  BundlesMinf movingBundleInfo( movingBundlesFilename.c_str() ) ;


  //   xxxxxxxxxxxxxxxxxxxxxxxx Reading transform xxxxxxxxxxxxxxxxxxxxxxxx   //
  readTransform( transformFilename, affineTransform ) ;

  if ( isInverse )
  {

    std::vector<std::vector<float>> affineMatrix ;
    affineTransformToMatrix( affineTransform, affineMatrix ) ;

    std::vector<std::vector<float>> inverseAffineMatrix ;
    computeInverseMatrix4x4( affineMatrix,
                             inverseAffineMatrix ) ;

    matrixToAffineTransform( inverseAffineMatrix, affineTransform ) ;

    std::vector<std::vector<float>> productMatrix = matrixProduct(
                                           affineMatrix, inverseAffineMatrix ) ;

    std::cout << "affineMatrix : " << std::endl ;
    for ( int i = 0 ; i < affineMatrix.size() ; i++ )
    {

      for ( int j = 0 ; j < affineMatrix[ i ].size() ; j++ )
      {

        std::cout << affineMatrix[ i ][ j ] << "    " ;

      }
      std::cout << "\n" ;

    }
    std::cout << "\n" ;

    std::cout << "inverseAffineMatrix : " << std::endl ;
    for ( int i = 0 ; i < inverseAffineMatrix.size() ; i++ )
    {

      for ( int j = 0 ; j < inverseAffineMatrix[ i ].size() ; j++ )
      {

        std::cout << inverseAffineMatrix[ i ][ j ] << "    " ;

      }
      std::cout << "\n" ;

    }
    std::cout << "\n" ;

    std::cout << "productMatrix : " << std::endl ;
    for ( int i = 0 ; i < productMatrix.size() ; i++ )
    {

      for ( int j = 0 ; j < productMatrix[ i ].size() ; j++ )
      {

        std::cout << productMatrix[ i ][ j ] << "    " ;

      }
      std::cout << "\n" ;

    }
    std::cout << "\n" ;


  }

  if ( verbose )
  {

    std::cout << "Transform : \n" ;
    for ( int i = 0 ; i < 12 ; i++ )
    {

      if ( i == 3 || i == 7 || i == 11 )
      {

        std::cout << affineTransform[ i ] << std::endl ;

      }
      else
      {

        std::cout << affineTransform[ i ] << "    " ;

      }

    }

  }


  //------------------------------ Sanity checks -----------------------------//
  int nbPointsMoving = movingBundlesData.pointsPerTrack[ 0 ] ;
  int nbCurvesMoving = movingBundlesData.curves_count ;
  for( int fiber = 1 ; fiber < nbCurvesMoving ; fiber++ )
  {

    if ( nbPointsMoving != movingBundlesData.pointsPerTrack[ fiber ] )
    {

      std::cout << "ERROR : All the fibers in the moving bundle must have the "
                << "same number of points, got " << nbPointsMoving
                << " for point 0 and "
                << movingBundlesData.pointsPerTrack[ fiber ]
                << " for point " << fiber << std::endl ;
      exit( 1 ) ;

    }

  }


  int nbPoints = nbPointsMoving ;

  if ( verbose )
  {

    std::cout << "Sanity cheks : OK " << std::endl ;

  }

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  if ( verbose )
  {

    std::cout << "Applying transform... " ;

  }
  BundlesData movedBundlesData( movingBundlesData ) ;
  applyAffineToBundle( movingBundlesData.matrixTracks,
                       affineTransform,
                       movingBundlesData.curves_count,
                       nbPoints,
                       movedBundlesData.matrixTracks ) ;

  BundlesMinf movedBundleInfo( movingBundleInfo ) ;

  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }

  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //

  // Saving results
  if ( verbose )
  {

    std::cout << "Saving result... " ;

  }
  std::string outputBundlesFilename ;
  std::string outputBundlesDataFilename ;

  if ( endswith( outputFilename, ".bundles" ) || endswith( outputFilename,
                                                              ".bundlesdata" ) )
  {

    if ( format != ".bundles" )
    {

      std::cout << "ERROR : input and output must be the same format"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesDataFilename = replaceExtension( outputFilename,
                                                              ".bundlesdata" ) ;

    outputBundlesFilename = replaceExtension( outputFilename, ".bundles" ) ;

  }
  else if ( endswith( outputFilename, ".trk" ) )
  {

    if ( format != ".trk" )
    {

      std::cout << "ERROR : input and output must be the same format"
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

      std::cout << "ERROR : input and output must be the same format"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = replaceExtension( outputFilename, ".minf" ) ;

  }
  else
  {

    std::cout << "The only tractogram format supported are .bundles/.trk/.tck"
              << std::endl ;
    exit( 1 ) ;

  }


  movedBundlesData.write( outputBundlesDataFilename.c_str(), movedBundleInfo ) ;

  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }

  return 0 ;

}
