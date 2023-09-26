///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------- Libraries ----------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <ncurses.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <zlib.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

#include <stdexcept>

#include "ConvertBundleFormat.h"
#include "ioWrapper.h"
#include "niftiImage.h"




///////////////////////////////////////////////////////////////////////////////
//////////////// .bundlesdata to .trk with reference image  ///////////////////
///////////////////////////////////////////////////////////////////////////////
void bundlesData2Trk( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{

  std::string inputFileStr = inputFile ;
  if ( !endswith( inputFileStr, ".bundles" ) &&
                                     !endswith( inputFileStr, ".bundlesdata" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->bundlesData2Trk : input file must "
                  << "be .bundles/.bundlesdata" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string outputFileStr = outputFile ;
  if ( !endswith( outputFileStr, ".trk" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->bundlesData2Trk : output file must "
                  << "be .trk" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string refFileStr = refFile ;

  std::cout << "Reading : " << inputFile << std::endl ;
  BundlesData bundlesdata( inputFile ) ;
  BundlesMinf bundlesInfo( inputFile ) ;

  if ( is_file( refFileStr ) )
  {

    std::cout << "Reading : " << refFile << std::endl ;
    NiftiImage referenceAnatomy( refFile ) ;

    bundlesInfo.vox_to_ras = referenceAnatomy.vox_to_ras ;
    bundlesInfo.resolution = referenceAnatomy.resolution ;
    bundlesInfo.size = referenceAnatomy.size ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->bundlesData2Trk : reference image "
                  << refFileStr << " does not exists\n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  bundlesdata.isBundles = false ;
  bundlesdata.isTrk = true ;

  bundlesInfo.isBundles = false ;
  bundlesInfo.isTrk = true ;

  // Because we do not want to produce a .minf when converting .bundles->.trk
  bundlesInfo.haveMinf = false ;

  std::cout << "Writing : " << outputFile << std::endl ;
  bundlesdata.write( outputFile, bundlesInfo ) ;


}


///////////////////////////////////////////////////////////////////////////////
//////////////////////// Function trk -> bundlesdata //////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Trk2BundlesData( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{


  std::string inputFileStr = inputFile ;
  if ( !endswith( inputFileStr, ".trk" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Trk2BundlesData : input file must "
                  << "be .trk" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string outputFileStr = outputFile ;
  if ( !endswith( outputFileStr, ".bundles" ) &&
                                    !endswith( outputFileStr, ".bundlesdata" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Trk2BundlesData : output file must "
                  << "be .bundles/.bundlesdata" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }


  std::string refFileStr = refFile ;

  std::cout << "Reading : " << inputFile << std::endl ;
  BundlesData bundlesdata( inputFile ) ;
  BundlesMinf bundlesInfo( inputFile ) ;

  if ( is_file( refFileStr ) )
  {

    std::cout << "Reading : " << refFile << std::endl ;
    NiftiImage referenceAnatomy( refFile ) ;

    bundlesInfo.vox_to_ras = referenceAnatomy.vox_to_ras ;
    bundlesInfo.resolution = referenceAnatomy.resolution ;
    bundlesInfo.size = referenceAnatomy.size ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Trk2BundlesData : reference image "
                  << refFileStr << " does not exists\n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  bundlesInfo.fillDefaultBundles() ;

  bundlesdata.isBundles = true ;
  bundlesdata.isTrk = false ;

  bundlesInfo.isBundles = true ;
  bundlesInfo.isTrk = false ;

  // Because we want to produce a .bundles when converting .trk->.bundles
  bundlesInfo.haveMinf = true ;

  std::cout << "Writing : " << outputFile << std::endl ;
  bundlesdata.write( outputFile, bundlesInfo ) ;

}



///////////////////////////////////////////////////////////////////////////////
//////////////////////// Function trk -> bundlesdata //////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Tck2BundlesData( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{

  std::string inputFileStr = inputFile ;
  if ( !endswith( inputFileStr, ".tck" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Tck2BundlesData : input file must "
                  << "be .tck" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string outputFileStr = outputFile ;
  if ( !endswith( outputFileStr, ".bundles" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Tck2BundlesData : output file must "
                  << "be .bundles/.bundlesdata" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string refFileStr = refFile ;

  std::cout << "Reading : " << inputFile << std::endl ;
  BundlesData bundlesdata( inputFile ) ;
  BundlesMinf bundlesInfo( inputFile ) ;

  if ( is_file( refFileStr ) )
  {

    std::cout << "Reading : " << refFile << std::endl ;
    NiftiImage referenceAnatomy( refFile ) ;

    bundlesInfo.vox_to_ras = referenceAnatomy.vox_to_ras ;
    bundlesInfo.resolution = referenceAnatomy.resolution ;
    bundlesInfo.size = referenceAnatomy.size ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Tck2BundlesData : reference image "
                  << refFileStr << " does not exists\n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }


  bundlesInfo.fillDefaultBundles() ;


  // ------------------------- Writing .bundlesdata ------------------------- //
  std::vector<std::vector<float>> ras_to_vox ;
  computeInverseVoxToRas( bundlesInfo.vox_to_ras, ras_to_vox ) ;

  std::vector<float> correctedTracks ;
  // Applying vox_to_ras to tck fibers
  applyVoxToRasToTractogram( bundlesdata.matrixTracks, ras_to_vox,
                                                             correctedTracks ) ;

  bundlesdata.matrixTracks = correctedTracks ;
  correctedTracks = std::vector<float>() ; // Free memory

  std::vector<float> flippedTracks ;
  flipTractogram( bundlesdata.matrixTracks,
                  bundlesInfo.resolution,
                  bundlesInfo.size,
                  flippedTracks,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;

  bundlesdata.matrixTracks = flippedTracks ;
  flippedTracks = std::vector<float>() ; // Free memory

  bundlesInfo.isTck = false ;
  bundlesInfo.isBundles = true ;

  bundlesdata.isTck = false ;
  bundlesdata.isBundles = true ;

  // Because we do not want to produce a .minf when converting .tck->.bundles
  bundlesInfo.haveMinf = true ;

  std::cout << "Writing : " << outputFile << std::endl ;
  bundlesdata.write( outputFile, bundlesInfo ) ;


  // bundlesInfo.isTck = false ;
  // bundlesInfo.isBundles = true ;
  //
  // bundlesdata.isTck = false ;
  // bundlesdata.isBundles = true ;
  //
  // std::cout << "Writing : " << outputFile << std::endl ;
  // bundlesdata.write( outputFile, bundlesInfo ) ;

}



///////////////////////////////////////////////////////////////////////////////
//////////////////////// Function bundlesdata -> tck //////////////////////////
///////////////////////////////////////////////////////////////////////////////
void BundlesData2Tck( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{

  std::string inputFileStr = inputFile ;
  if ( !endswith( inputFileStr, ".bundles" ) &&
                                     !endswith( inputFileStr, ".bundlesdata" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->BundlesData2Tck : input file must "
                  << "be .bundles/.bundlesdata" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string outputFileStr = outputFile ;
  if ( !endswith( outputFileStr, ".tck" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->BundlesData2Tck : output file must "
                  << "be .tck" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string refFileStr = refFile ;

  std::cout << "Reading : " << inputFile << std::endl ;
  BundlesData bundlesdata( inputFile ) ;
  BundlesMinf bundlesInfo( inputFile ) ;

  if ( is_file( refFileStr ) )
  {

    bundlesInfo.fillDefaultTck() ;

    std::cout << "Reading : " << refFile << std::endl ;
    NiftiImage referenceAnatomy( refFile ) ;

    bundlesInfo.vox_to_ras = referenceAnatomy.vox_to_ras ;
    bundlesInfo.resolution = referenceAnatomy.resolution ;
    bundlesInfo.size = referenceAnatomy.size ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->BundlesData2Tck : reference image "
                  << refFileStr << " does not exists\n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }


  // ---------------------------- Writing .tck ----------------------------- //
  std::vector<float> flippedMatrixTracks ;
  flipTractogram( bundlesdata.matrixTracks,
                  bundlesInfo.resolution,
                  bundlesInfo.size,
                  flippedMatrixTracks,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;

  bundlesdata.matrixTracks = flippedMatrixTracks ;
  flippedMatrixTracks = std::vector<float>() ; // Free memory


  // Applying vox_to_ras to tck fibers
  std::vector<float> correctedTracks ;
  applyVoxToRasToTractogram( bundlesdata.matrixTracks, bundlesInfo.vox_to_ras,
                                                             correctedTracks ) ;


  bundlesdata.matrixTracks = correctedTracks ;
  correctedTracks = std::vector<float>() ; // Free mememory

  bundlesdata.isBundles = false ;
  bundlesdata.isTck = true ;

  bundlesInfo.isBundles = false ;
  bundlesInfo.isTck = true ;

  // Because we do not want to produce a .minf when converting .bundles->.tck
  bundlesInfo.haveMinf = false ;

  std::cout << "Writing : " << outputFile << std::endl ;
  bundlesdata.write( outputFile, bundlesInfo ) ;

  // bundlesdata.isBundles = false ;
  // bundlesdata.isTck = true ;
  //
  // bundlesInfo.isBundles = false ;
  // bundlesInfo.isTck = true ;
  //
  // std::cout << "Writing : " << outputFile << std::endl ;
  // bundlesdata.write( outputFile, bundlesInfo ) ;

}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// Function .trk -> tck /////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Trk2Tck( const char* inputFile,
              const char* outputFile,
              const char* refFile,
              int flip_x,
              int flip_y,
              int flip_z )
{

  std::string inputFileStr = inputFile ;
  if ( !endswith( inputFileStr, ".trk" ) &&
                                     !endswith( inputFileStr, ".bundlesdata" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Trk2Tck : input file must "
                  << "be .trk" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string outputFileStr = outputFile ;
  if ( !endswith( outputFileStr, ".tck" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->BundlesData2Tck : output file must "
                  << "be .tck" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string refFileStr = refFile ;

  std::cout << "Reading : " << inputFile << std::endl ;
  BundlesData bundlesdata( inputFile ) ;
  BundlesMinf bundlesInfo( inputFile ) ;

  if ( is_file( refFileStr ) )
  {

    bundlesInfo.fillDefaultTck() ;

    std::cout << "Reading : " << refFile << std::endl ;
    NiftiImage referenceAnatomy( refFile ) ;

    bundlesInfo.vox_to_ras = referenceAnatomy.vox_to_ras ;
    bundlesInfo.resolution = referenceAnatomy.resolution ;
    bundlesInfo.size = referenceAnatomy.size ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Trk2Tck : reference image "
                  << refFileStr << " does not exists\n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }


  // ---------------------------- Writing .tck ----------------------------- //
  std::vector<float> flippedMatrixTracks ;
  flipTractogram( bundlesdata.matrixTracks,
                  bundlesInfo.resolution,
                  bundlesInfo.size,
                  flippedMatrixTracks,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;

  bundlesdata.matrixTracks = flippedMatrixTracks ;
  flippedMatrixTracks = std::vector<float>() ; // Free memory


  // Applying vox_to_ras to tck fibers
  std::vector<float> correctedTracks ;
  applyVoxToRasToTractogram( bundlesdata.matrixTracks, bundlesInfo.vox_to_ras,
                                                             correctedTracks ) ;


  bundlesdata.matrixTracks = correctedTracks ;
  correctedTracks = std::vector<float>() ; // Free mememory

  bundlesdata.isTrk = false ;
  bundlesdata.isTck = true ;

  bundlesInfo.isTrk = false ;
  bundlesInfo.isTck = true ;

  // Because we do not want to produce a .minf when converting .bundles->.tck
  bundlesInfo.haveMinf = false ;

  std::cout << "Writing : " << outputFile << std::endl ;
  bundlesdata.write( outputFile, bundlesInfo ) ;

  // bundlesdata.isBundles = false ;
  // bundlesdata.isTck = true ;
  //
  // bundlesInfo.isBundles = false ;
  // bundlesInfo.isTck = true ;
  //
  // std::cout << "Writing : " << outputFile << std::endl ;
  // bundlesdata.write( outputFile, bundlesInfo ) ;

}



///////////////////////////////////////////////////////////////////////////////
//////////////////////////// Function tck -> trk //////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Tck2Trk( const char* inputFile,
              const char* outputFile,
              const char* refFile,
              int flip_x,
              int flip_y,
              int flip_z )
{

  std::string inputFileStr = inputFile ;
  if ( !endswith( inputFileStr, ".tck" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Tck2Trk : input file must "
                  << "be .tck" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string outputFileStr = outputFile ;
  if ( !endswith( outputFileStr, ".trk" ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Tck2Trk : output file must "
                  << "be .trk" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::string refFileStr = refFile ;

  std::cout << "Reading : " << inputFile << std::endl ;
  BundlesData bundlesdata( inputFile ) ;
  BundlesMinf bundlesInfo( inputFile ) ;

  if ( is_file( refFileStr ) )
  {

    std::cout << "Reading : " << refFile << std::endl ;
    NiftiImage referenceAnatomy( refFile ) ;

    bundlesInfo.vox_to_ras = referenceAnatomy.vox_to_ras ;
    bundlesInfo.resolution = referenceAnatomy.resolution ;
    bundlesInfo.size = referenceAnatomy.size ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ConvertBundleFormat->Tck2Trk : reference image "
                  << refFileStr << " does not exists\n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }


  bundlesInfo.fillDefaultTrk() ;


  // ------------------------- Writing .bundlesdata ------------------------- //
  std::vector<std::vector<float>> ras_to_vox ;
  computeInverseVoxToRas( bundlesInfo.vox_to_ras, ras_to_vox ) ;

  std::vector<float> correctedTracks ;
  // Applying vox_to_ras to tck fibers
  applyVoxToRasToTractogram( bundlesdata.matrixTracks, ras_to_vox,
                                                             correctedTracks ) ;

  bundlesdata.matrixTracks = correctedTracks ;
  correctedTracks = std::vector<float>() ; // Free memory

  std::vector<float> flippedTracks ;
  flipTractogram( bundlesdata.matrixTracks,
                  bundlesInfo.resolution,
                  bundlesInfo.size,
                  flippedTracks,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;

  bundlesdata.matrixTracks = flippedTracks ;
  flippedTracks = std::vector<float>() ; // Free memory

  bundlesInfo.isTck = false ;
  bundlesInfo.isTrk = true ;

  bundlesdata.isTck = false ;
  bundlesdata.isTrk = true ;

  // Because we do not want to produce a .minf when converting .tck->.bundles
  bundlesInfo.haveMinf = false ;

  std::cout << "Writing : " << outputFile << std::endl ;
  bundlesdata.write( outputFile, bundlesInfo ) ;


  // bundlesInfo.isTck = false ;
  // bundlesInfo.isBundles = true ;
  //
  // bundlesdata.isTck = false ;
  // bundlesdata.isBundles = true ;
  //
  // std::cout << "Writing : " << outputFile << std::endl ;
  // bundlesdata.write( outputFile, bundlesInfo ) ;

}



///////////////////////////////////////////////////////////////////////////////
///////////// Function to apply vox_to_ras to point matrices //////////////////
///////////////////////////////////////////////////////////////////////////////
void applyVoxToRasToTractogram(
                              const std::vector<float>& tractogram,
                              const std::vector<std::vector<float>>& vox_to_ras,
                              std::vector<float>& newTractogram )
{

  int64_t sizeTractogram = tractogram.size() ;
  newTractogram.resize( sizeTractogram ) ;

  int64_t nbPointsTractogram = ( int64_t )( sizeTractogram / 3 ) ;

  #pragma omp parallel for
  for ( int pointIndex = 0 ; pointIndex < nbPointsTractogram ; pointIndex++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      newTractogram[ 3 * pointIndex + i ] = 0 ;
      for ( int j = 0 ; j < 3 ; j++ )
      {

        newTractogram[ 3 * pointIndex + i ] += tractogram[ 3 * pointIndex + j ]
                                                        * vox_to_ras[ i ][ j ] ;
      }
      newTractogram[ 3 * pointIndex + i ] += vox_to_ras[ i ][ 3 ] ;

    }

  }

}


///////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute adjacent of 3x3 matrix ////////////////////////
///////////////////////////////////////////////////////////////////////////////
void computeAdjugateMatrix( const std::vector<std::vector<float>>& matrix,
                            std::vector<std::vector<float>>& adjugate )
{

  adjugate.resize( 3 ) ;
  adjugate[ 0 ].resize( 3 ) ;
  adjugate[ 1 ].resize( 3 ) ;
  adjugate[ 2 ].resize( 3 ) ;

  // Sanity checks
  if ( matrix.size() != 3 )
  {

    std::cout << "ERRROR in computeAdjugateMatrix : The input matrix must be "
              << "a 3x3 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 3 )
    {

      std::cout << "ERRROR in computeAdjugateMatrix : The input matrix must be "
                << "a 3x3 matrix" << std::endl ;
      exit( 1 ) ;

    }

  }

  // Transposing matrix
  std::vector<std::vector<float>> matrixTranspose ;
  matrixTranspose.resize( 3 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    matrixTranspose[ i ].resize( 3 ) ;
    for ( int j = 0 ; j < 3 ; j++ )
    {
      matrixTranspose[ i ][ j ] = matrix[ j ][ i ] ;

    }

  }

  // Computing adjugate matrix
  // #pragma omp parallel for
  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      std::vector<std::vector<float>> minorMatrix ;
      minorMatrix.resize( 2 ) ;
      minorMatrix[ 0 ].resize( 2 ) ;
      minorMatrix[ 1 ].resize( 2 ) ;

      int lineIndex = 0 ;
      int columnIndex = 0 ;
      bool stop = false ;
      for ( int n = 0 ; n < 3 ; n++ )
      {

        for ( int p = 0 ; p < 3 ; p++ )
        {

          if ( n != i && p != j )
          {

            minorMatrix[ lineIndex ][ columnIndex ] =
                                                     matrixTranspose[ n ][ p ] ;
            if ( columnIndex == 1 )
            {

              if ( lineIndex < 1 )
              {

                lineIndex++ ;
                columnIndex = 0 ;

              }
              else
              {

                stop = true ;
                break ;

              }

            }
            else
            {

              columnIndex++ ;

            }

          }

        }

        if ( stop )
        {

          break ;

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

      adjugate[ i ][ j ] = minorMatrix[ 0 ][ 0 ] * minorMatrix[ 1 ][ 1 ] -
                           minorMatrix[ 1 ][ 0 ] * minorMatrix[ 0 ][ 1 ] ;
      adjugate[ i ][ j ] *= sign ;

    }

  }

}

///////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute determinant 3x3 matrix ////////////////////////
///////////////////////////////////////////////////////////////////////////////
float computeDeterminantMatrix( const std::vector<std::vector<float>>& matrix )
{

  // Sanity checks
  if ( matrix.size() != 3 )
  {

    std::cout << "ERRROR in computeDeterminantMatrix : The input matrix must "
              << "be a 3x3 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 3 )
    {

      std::cout << "ERRROR in computeDeterminantMatrix : The input matrix must "
                << "be a 3x3 matrix" << std::endl ;
      exit( 1 ) ;

    }

  }


  // Computing determinant
  float determinant = 0 ;
  int j = 0 ; // Select first row to compute determinant
  // #pragma omp parallel for
  for ( int i = 0 ; i < 3 ; i++ )
  {

    std::vector<std::vector<float>> minorMatrix ;
    minorMatrix.resize( 2 ) ;
    minorMatrix[ 0 ].resize( 2 ) ;
    minorMatrix[ 1 ].resize( 2 ) ;

    int lineIndex = 0 ;
    int columnIndex = 0 ;
    bool stop = false ;
    for ( int n = 0 ; n < 3 ; n++ )
    {

      for ( int p = 0 ; p < 3 ; p++ )
      {

        if ( n != i && p != j )
        {

          minorMatrix[ lineIndex ][ columnIndex ] = matrix[ n ][ p ] ;
          if ( columnIndex == 1 )
          {

            if ( lineIndex < 1 )
            {

              lineIndex++ ;
              columnIndex = 0 ;

            }
            else
            {

              stop = true ;
              break ;

            }

          }
          else
          {

            columnIndex++ ;

          }

        }

      }

      if ( stop )
      {

        break ;

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

    determinant += sign * matrix[ i ][ j ] * (
                               minorMatrix[ 0 ][ 0 ] * minorMatrix[ 1 ][ 1 ] -
                               minorMatrix[ 1 ][ 0 ] * minorMatrix[ 0 ][ 1 ] ) ;

  }

  return determinant ;

}


///////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute inverse 3x3 matrix ////////////////////////
///////////////////////////////////////////////////////////////////////////////
void computeInverseMatrix( const std::vector<std::vector<float>>& matrix,
                           std::vector<std::vector<float>>& inverseMatrix )
{

  // Sanity checks
  if ( matrix.size() != 3 )
  {

    std::cout << "ERRROR in computeInverseMatrix : The input matrix must be a "
              << "3x3 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 3 )
    {

      std::cout << "ERRROR in computeInverseMatrix : The input matrix must be "
                << "a 3x3 matrix" << std::endl ;
      exit( 1 ) ;

    }

  }

  float determinant = computeDeterminantMatrix( matrix ) ;
  if ( abs( determinant ) < 1e-30 )
  {

    std::cout << "ERROR in computeInverseMatrix : determinant is 0, singular "
              << "matrix are not invertible" << std::endl ;
    exit( 1 ) ;

  }

  std::vector<std::vector<float>> adjugateMatrix ;
  computeAdjugateMatrix( matrix, adjugateMatrix ) ;

  inverseMatrix.resize( 3 ) ;
  inverseMatrix[ 0 ].resize( 3, 0 ) ;
  inverseMatrix[ 1 ].resize( 3, 0 ) ;
  inverseMatrix[ 2 ].resize( 3, 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      inverseMatrix[ i ][ j ] = 1 / determinant * adjugateMatrix[ i ][ j ] ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute inverse vox_to_ras /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void computeInverseVoxToRas( const std::vector<std::vector<float>>& vox_to_ras,
                             std::vector<std::vector<float>>& ras_to_vox )
{

  // Sanity checks
  if ( vox_to_ras.size() != 4 )
  {

    std::cout << "ERRROR in computeInverseVoxToRas : The input matrix must be "
              << "a 4x4 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( int i = 0 ; i < 4 ; i++ )
  {

    if ( vox_to_ras[ i ].size() != 4 )
    {

      std::cout << "ERRROR in computeInverseVoxToRas : The input matrix must "
                << "be a 4x4 matrix" << std::endl ;
      exit( 1 ) ;

    }

    if ( i == 3 )
    {

      if ( vox_to_ras[ i ][ 0 ] != 0.0f || vox_to_ras[ i ][ 1 ] != 0.0f ||
           vox_to_ras[ i ][ 2 ] != 0.0f || vox_to_ras[ i ][ 3 ] != 1.0f )
        {

          std::cout << "ERRROR in computeInverseVoxToRas : The input matrix is "
          << "not a 4x4 vox_to_ras matrix" << std::endl ;
          exit( 1 ) ;

        }

    }

  }

  // Resizeing ras_to_vox
  ras_to_vox.resize( 4 ) ;
  ras_to_vox[ 0 ].resize( 4, 0 ) ;
  ras_to_vox[ 1 ].resize( 4, 0 ) ;
  ras_to_vox[ 2 ].resize( 4, 0 ) ;
  ras_to_vox[ 3 ].resize( 4, 0 ) ;

  // Getting 3x3 subTransform
  std::vector<std::vector<float>> subTransform ;
  subTransform.resize( 3 ) ;
  subTransform[ 0 ].resize( 3 ) ;
  subTransform[ 1 ].resize( 3 ) ;
  subTransform[ 2 ].resize( 3 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      subTransform[ i ][ j ] = vox_to_ras[ i ][ j ] ;

    }

  }

  // Computing inverse 3x3 subTransform
  std::vector<std::vector<float>> inverseSubTransform ;
  computeInverseMatrix( subTransform, inverseSubTransform ) ;

  // Getting translation
  std::vector<float> translation( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = vox_to_ras[ i ][ 3 ] ;

  }

  // Compute inverseSubTransform.dot( - traslation )
  std::vector<float> translation2( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      translation2[ i ] += - inverseSubTransform[ i ][ j ] * translation[ j ] ;

    }

  }


  // Filling ras_to_vox
  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      ras_to_vox[ i ][ j ] = inverseSubTransform[ i ][ j ] ;

    }

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    ras_to_vox[ i ][ 3 ] = translation2[ i ] ;

  }

  ras_to_vox[ 3 ][ 3 ] = 1.0f ;

}


///////////////////////////////////////////////////////////////////////////////
////////////////////////////// Flip tractogrtam ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void flipTractogram( BundlesData& inputTractogram,
                     BundlesMinf& inputTractogramInfo,
                     BundlesData& outputTractogram,
                     int flip_x,
                     int flip_y,
                     int flip_z )
{

  int nbFibers = inputTractogramInfo.curves_count ;
  int64_t nbPointsTractogram = ( int64_t )( inputTractogram.matrixTracks.size()
                                                                         / 3 ) ;

  if ( nbFibers != inputTractogram.pointsPerTrack.size() )
  {

    std::cout << "ERROR : the number of fibers in the .bundles file is "
             << "different from the number of fibers in .bundlesdata file \n" ;
    exit( 1 ) ;

  }

  for ( int fiber = 0 ; fiber < nbFibers ; fiber++ )
  {

    int nbPoints = inputTractogram.pointsPerTrack[ fiber ] ;
    if ( nbPoints < 0 )
    {

      std::cout << "ERROR : the number of points for each fiber has to be "
                << "strictly positive " << std::endl ;
      exit( 1 ) ;

    }

  }


  outputTractogram.pointsPerTrack = inputTractogram.pointsPerTrack ;
  outputTractogram.curves_count = inputTractogram.curves_count ;
  outputTractogram.matrixTracks.resize( inputTractogram.matrixTracks.size() ) ;

  float translation[ 3 ] = { 0 } ;
  int flip[ 3 ] = { 0 } ;
  flip[ 0 ] = flip_x ;
  flip[ 1 ] = flip_y ;
  flip[ 2 ] = flip_z ;


  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( flip[ i ] == -1 )
    {

      int tmpMiddle = ceilf( inputTractogramInfo.size[ i ] / 2 ) + 1 ;

      if ( inputTractogramInfo.size[ i ] % 2 != 0 )
      {

        translation[ i ] = tmpMiddle * inputTractogramInfo.resolution[ i ] ;

      }
      else
      {

        translation[ i ] = ( ( 2 * tmpMiddle - 1 ) / 2 ) *
                                          inputTractogramInfo.resolution[ i ] ;

      }

    }

  }

  std::vector< float >& tractogramFibers = inputTractogram.matrixTracks ;
  for ( int point = 0 ; point < nbPointsTractogram ; point++ )
  {

    outputTractogram.matrixTracks[ 3 * point + 0 ]  =
                                      flip_x * tractogramFibers[ 3 * point + 0 ]
                                                        + 2 * translation[ 0 ] ;

    outputTractogram.matrixTracks[ 3 * point + 1 ]  =
                                      flip_y * tractogramFibers[ 3 * point + 1 ]
                                                        + 2 * translation[ 1 ] ;

    outputTractogram.matrixTracks[ 3 * point + 2 ]  =
                                      flip_z * tractogramFibers[ 3 * point + 2 ]
                                                        + 2 * translation[ 2 ] ;


  }

}



// -------------------------------------------------------------------------- //

void flipTractogram( const std::vector<float>& inputTractogram,
                     const std::vector<float>& resolution,
                     const std::vector<short int>& size,
                     std::vector<float>& outputTractogram,
                     int flip_x,
                     int flip_y,
                     int flip_z )
{

  int64_t nbPointsTractogram = ( int64_t )( inputTractogram.size() / 3 ) ;

  outputTractogram.resize( inputTractogram.size() ) ;

  float translation[ 3 ] = { 0 } ;
  int flip[ 3 ] = { 0 } ;
  flip[ 0 ] = flip_x ;
  flip[ 1 ] = flip_y ;
  flip[ 2 ] = flip_z ;


  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( flip[ i ] == -1 )
    {

      int tmpMiddle = ceilf( size[ i ] / 2 ) + 1 ;

      if ( size[ i ] % 2 != 0 )
      {

        translation[ i ] = tmpMiddle * resolution[ i ] ;

      }
      else
      {

        translation[ i ] = ( ( 2 * tmpMiddle - 1 ) / 2 ) * resolution[ i ] ;

      }

    }

  }


  for ( int point = 0 ; point < nbPointsTractogram ; point++ )
  {

    outputTractogram[ 3 * point + 0 ]  =
                                      flip_x * inputTractogram[ 3 * point + 0 ]
                                                        + 2 * translation[ 0 ] ;

    outputTractogram[ 3 * point + 1 ]  =
                                      flip_y * inputTractogram[ 3 * point + 1 ]
                                                        + 2 * translation[ 1 ] ;

    outputTractogram[ 3 * point + 2 ]  =
                                      flip_z * inputTractogram[ 3 * point + 2 ]
                                                        + 2 * translation[ 2 ] ;


  }

}


///////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ///////////////
///////////////////////////////////////////////////////////////////////////////

int getFlagPosition( int argc, char* argv[], const std::string& flag )
{

  std::string cmd ;
  for ( int i = 0 ; i < argc ; i++ )
  {

    std::string arg = argv[ i ] ;
    if ( 0 == arg.find( flag ) )
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

  int index_input, index_output, index_reference, index_x, index_y, index_z,
                                                    index_verbose, index_help ;
  index_input = getFlagPosition( argc, argv, "-i" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_reference = getFlagPosition( argc, argv, "-r" ) ;
  index_x = getFlagPosition( argc, argv, "-x" ) ;
  index_y = getFlagPosition( argc, argv, "-y" ) ;
  index_z = getFlagPosition( argc, argv, "-z" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : " << std::endl
              << "-i : Input file name " << std::endl
              << "-o : Out file name " << std::endl
              << "-r : Reference image of the bundles" << std::endl
              << "[-x] : Flip around x brain axis " << std::endl
              << "[-y] : Flip around y brain axis " << std::endl
              << "[-z] : Flip around z brain axis " << std::endl
              << "[-v] : set verbosity level at 1 " << std::endl
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

  if ( !index_reference )
  {

    std::cout << "-r argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_x )
  {

    flip_x = -1 ;

  }
  else
  {

    flip_x = 1 ;

  }

  if ( index_y )
  {

    flip_y = -1 ;

  }
  else
  {

    flip_y = 1 ;

  }

  if ( index_z )
  {

    flip_z = -1 ;

  }
  else
  {

    flip_z = 1 ;

  }

  if ( index_verbose )
  {

    verbose = 1 ;

  }

  // Checking format input reference
  char* refFilename ;
  if ( index_reference )
  {

    refFilename = argv[ index_reference + 1 ] ;
    std::string refString = refFilename ;

    if ( !endswith( refString, ".nii" ) )
    {

      std::string outMessage = "ConvertBundleFormat : Only supported format " \
                               "for reference image is .nii \n" ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

    if ( !is_file( refString ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ConvertBundleFormat : reference image "
                    << refString << " does not exists\n" ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }


  // Different scenarios depending on the inputs
  std::string inputFilename ( argv[ index_input + 1 ] ) ;
  std::string outputFilename( argv[ index_output + 1 ] ) ;

  if ( !endswith( inputFilename, ".bundles" ) &&
       !endswith( inputFilename, ".bundlesdata" ) &&
       !endswith( inputFilename, ".trk" ) &&
       !endswith( inputFilename, ".tck" ) )
  {

    std::string outMessage = "ConvertBundleFormat: Only supported formats " \
                             "for input are .bundles/.bundlesdata, .trk, " \
                             ".tck \n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }

  if ( !endswith( outputFilename, ".bundles" ) &&
       !endswith( outputFilename, ".bundlesdata" ) &&
       !endswith( outputFilename, ".trk" ) &&
       !endswith( outputFilename, ".tck" ) )
  {

    std::stringstream outMessageOss ;
    std::string outMessage = "ConvertBundleFormat: Only supported formats " \
                             "for output ares .bundles/.bundlesdata, .trk, " \
                             ".tck \n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }

  if ( ( endswith( inputFilename, ".bundles" ) ||
         endswith( inputFilename, ".bundlesdata" ) ) &&
       ( endswith( outputFilename, ".bundles" ) ||
              endswith( outputFilename, ".bundlesdata" ) ) )
  {

    std::cout << "Please choose a different format for input and output, "
              << "currently they are both .bundles" << std::endl ;
    exit( 1 ) ;

  }
  else if ( endswith( inputFilename, ".trk" ) &&
       endswith( outputFilename, ".trk" ) )
  {

    std::cout << "Please choose a different format for input and output, "
              << "currently they are both .trk" << std::endl ;
    exit( 1 ) ;

  }
  else if ( endswith( (std::string)inputFilename, ".tck" ) &&
       endswith( (std::string)outputFilename, ".tck" ) )
  {

    std::cout << "Please choose a different format for input and output, "
              << "currently they are both .trk" << std::endl ;
    exit( 1 ) ;

  }



  if ( endswith( inputFilename, ".bundles" ) && endswith( outputFilename,
                                                                      ".trk" ) )
  {


    std::cout << "Conversion .bundles/.bundles -> .trk ... \n" ;

    bundlesData2Trk( inputFilename.c_str(),
                     outputFilename.c_str(),
                     refFilename,
                     flip_x,
                     flip_y,
                     flip_z ) ;

    return 0 ;


  }


  if ( endswith( inputFilename, ".trk" ) && endswith( outputFilename,
                                                                  ".bundles" ) )
  {


    std::cout << "Conversion .trk -> .bundles/.bundlesdata ... \n" ;

    Trk2BundlesData( inputFilename.c_str(),
                     outputFilename.c_str(),
                     refFilename,
                     flip_x,
                     flip_y,
                     flip_z ) ;

    return 0 ;


  }


  if ( endswith( inputFilename, ".tck" ) && endswith( outputFilename,
                                                                  ".bundles" ) )
  {


    std::cout << "Conversion .tck -> .bundles/.bundlesdata ... \n" ;

    Tck2BundlesData( inputFilename.c_str(),
                     outputFilename.c_str(),
                     refFilename,
                     flip_x,
                     flip_y,
                     flip_z ) ;

    return 0 ;

  }


  if ( endswith( inputFilename, ".bundles" ) && endswith( outputFilename,
                                                                      ".tck" ) )
  {


    std::cout << "Conversion .bundles/.bundlesdata -> .tck ... \n" ;

    BundlesData2Tck( inputFilename.c_str(),
                     outputFilename.c_str(),
                     refFilename,
                     flip_x,
                     flip_y,
                     flip_z ) ;

    return 0 ;

  }


  if ( endswith( inputFilename, ".trk" ) && endswith( outputFilename, ".tck" ) )
  {


    std::cout << "Conversion .trk -> .tck ... \n" ;

    Trk2Tck( inputFilename.c_str(),
             outputFilename.c_str(),
             refFilename,
             flip_x,
             flip_y,
             flip_z ) ;

    return 0 ;

  }


  if ( endswith( inputFilename, ".tck" ) && endswith( outputFilename, ".trk" ) )
  {


    std::cout << "Conversion .tck -> .trk ... \n" ;

    Tck2Trk( inputFilename.c_str(),
             outputFilename.c_str(),
             refFilename,
             flip_x,
             flip_y,
             flip_z ) ;

    return 0 ;

  }

}
