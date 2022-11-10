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


#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include "processAtlasInformationMDF.h"

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



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------//
//----------------------------- Read .nii header -----------------------------//
//----------------------------------------------------------------------------//
niiFormat readNiiHeader( const char* niiFilename )
{

  niiFormat niiHeader ;

  std::string niiString( niiFilename ) ;

  if ( niiString.find( ".nii" ) != std::string::npos &&
                                niiString.find( ".gz" ) == std::string::npos )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << niiFilename << std::endl ;

    }

    std::FILE *f = std::fopen( niiFilename, "rb" ) ;
    if ( !f )
    {

      std::cout << "Problem opening file : " << niiFilename << std::endl ;
      exit( 1 ) ;

    }

    std::fread( &niiHeader.sizeof_hdr, sizeof( niiHeader.sizeof_hdr ), 1, f ) ;
    std::fread( &niiHeader.data_type, sizeof( niiHeader.data_type ), 1, f ) ;
    std::fread( &niiHeader.db_name, sizeof( niiHeader.db_name ), 1, f ) ;
    std::fread( &niiHeader.extents, sizeof( niiHeader.extents ), 1, f ) ;
    std::fread( &niiHeader.session_error, sizeof( niiHeader.session_error ),
                                                                       1, f ) ;
    std::fread( &niiHeader.regular, sizeof( niiHeader.regular ), 1, f ) ;
    std::fread( &niiHeader.dim_info, sizeof( niiHeader.dim_info ), 1, f ) ;
    std::fread( &niiHeader.dim, sizeof( niiHeader.dim ), 1, f );
    std::fread( &niiHeader.intent_p1, sizeof( niiHeader.intent_p1 ), 1, f ) ;
    std::fread( &niiHeader.intent_p2, sizeof( niiHeader.intent_p2 ), 1, f ) ;
    std::fread( &niiHeader.intent_p3, sizeof( niiHeader.intent_p3 ), 1, f ) ;
    std::fread( &niiHeader.intent_code, sizeof( niiHeader.intent_code ),
                                                                      1, f ) ;
    std::fread( &niiHeader.datatype, sizeof( niiHeader.datatype ), 1, f ) ;
    std::fread( &niiHeader.bitpix, sizeof( niiHeader.bitpix ), 1, f ) ;
    std::fread( &niiHeader.slice_start, sizeof( niiHeader.slice_start ),
                                                                      1, f ) ;
    std::fread( &niiHeader.pixdim, sizeof( niiHeader.pixdim ), 1, f ) ;
    std::fread( &niiHeader.scl_slope, sizeof( niiHeader.scl_slope ), 1, f ) ;
    std::fread( &niiHeader.scl_inter, sizeof( niiHeader.scl_inter ), 1, f ) ;
    std::fread( &niiHeader.slice_end, sizeof( niiHeader.slice_end ), 1, f ) ;
    std::fread( &niiHeader.slice_code, sizeof( niiHeader.slice_code ), 1, f ) ;
    std::fread( &niiHeader.xyzt_units, sizeof( niiHeader.xyzt_units ), 1, f ) ;
    std::fread( &niiHeader.cal_max, sizeof( niiHeader.cal_max ), 1, f ) ;
    std::fread( &niiHeader.cal_min, sizeof( niiHeader.cal_min ), 1, f ) ;
    std::fread( &niiHeader.slice_duration, sizeof( niiHeader.slice_duration ),
                                                                        1, f ) ;
    std::fread( &niiHeader.toffset, sizeof( niiHeader.toffset ), 1, f ) ;
    std::fread( &niiHeader.glmax, sizeof( niiHeader.glmax ), 1, f ) ;
    std::fread( &niiHeader.gmin, sizeof( niiHeader.gmin ), 1, f ) ;
    std::fread( &niiHeader.descrip, sizeof( niiHeader.descrip ), 1, f ) ;
    std::fread( &niiHeader.aux_file, sizeof( niiHeader.aux_file ), 1, f ) ;
    std::fread( &niiHeader.qform_code, sizeof( niiHeader.qform_code ), 1, f ) ;
    std::fread( &niiHeader.sform_code, sizeof( niiHeader.sform_code ), 1, f ) ;
    std::fread( &niiHeader.quatern_b, sizeof( niiHeader.quatern_b ), 1, f ) ;
    std::fread( &niiHeader.quatern_c, sizeof( niiHeader.quatern_c ), 1, f ) ;
    std::fread( &niiHeader.quatern_d, sizeof( niiHeader.quatern_d ), 1, f ) ;
    std::fread( &niiHeader.qoffset_x, sizeof( niiHeader.qoffset_x ), 1, f ) ;
    std::fread( &niiHeader.qoffset_y, sizeof( niiHeader.qoffset_y ), 1, f ) ;
    std::fread( &niiHeader.qoffset_z, sizeof( niiHeader.qoffset_z ), 1, f ) ;
    std::fread( &niiHeader.srow_x, sizeof( niiHeader.srow_x ), 1, f ) ;
    std::fread( &niiHeader.srow_y, sizeof( niiHeader.srow_y ), 1, f ) ;
    std::fread( &niiHeader.srow_z, sizeof( niiHeader.srow_z ), 1, f ) ;
    std::fread( &niiHeader.intent_name, sizeof( niiHeader.intent_name ), 1,
                                                                           f ) ;
    std::fread( &niiHeader.magic, sizeof( niiHeader.magic ), 1, f ) ;

    std::fclose ( f ) ;

    if ( verbose )
    {

      std::cout << "Done" << std::endl ;

    }

  }
  else
  {

    std::cout << "Error file format : the only format supported for the "
              << "reference image is .nii" << std::endl ;
    exit( 1 ) ;

  }

  return niiHeader ;

}


//----------------------------------------------------------------------------//
//-------------------------- Find radius of bundle ---------------------------//
//----------------------------------------------------------------------------//
void computeCenterAtlasBundleFibers(
                             BundlesDataFormat& atlasBundleData,
                             std::vector<float>& medialPointsAtlasBundleFibers )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    // Searching the medial point of atlas fiber
    std::vector<float> medialPointAtlasBundleFiber( 3, 0 ) ;
    atlasBundleData.computeMedialPointFiberWithDistance(
                                                 atlasBundleFiberIndex,
                                                 medialPointAtlasBundleFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointsAtlasBundleFibers[ 3 * atlasBundleFiberIndex + i ] =
                                              medialPointAtlasBundleFiber[ i ] ;

    }

  }

}

//----------------------------------------------------------------------------//
//----------------------- Compute average fiber bundle -----------------------//
//----------------------------------------------------------------------------//
void computeAverageFiberBundle(
                        BundlesDataFormat& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::vector<float>& medialPointAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;


  std::vector<float> referenceFiber( 3 * nbPoints, 0 ) ;
  int fiberIndex = 0 ;
  atlasBundleData.getFiberFromTractogram( atlasBundleData.matrixTracks,
                                          fiberIndex,
                                          nbPoints,
                                          referenceFiber ) ;



  std::vector<float> centerReferenceFiber( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    centerReferenceFiber[ i ] = medialPointsAtlasBundleFibers[ 3 * 0 + i ] ;

  }


  // #pragma omp parallel for reduction( +:averageFiber[:nbElementsAverageFiber ] )
  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    int offsetAtlas = 3 * nbPoints * atlasBundleFiberIndex ;

    std::vector<float> translation( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      translation[ i ] = centerReferenceFiber[ i ] -
               medialPointsAtlasBundleFibers[ 3 * atlasBundleFiberIndex + i ] ;

    }

    float directDistance = 0 ;
    float indirectDistance = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      directDistance += pow( referenceFiber[ 3 * 0 + i ] - ( atlasBundleData[ 3
                             * 0 + offsetAtlas + i ] + translation[ i ] ), 2 ) ;
      indirectDistance += pow( referenceFiber[ 3 * 0 + i ] - ( atlasBundleData[
            3 * ( nbPoints - 1 ) + offsetAtlas + i ] + translation[ i ] ), 2 ) ;

    }

    bool isDirectSens = true ;
    if ( directDistance > indirectDistance )
    {

      isDirectSens = false ;

    }

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      int k = point ;
      if ( !isDirectSens )
      {

        k = nbPoints - point - 1 ;

      }

      for ( int i = 0 ; i < 3 ; i++ )
      {

        averageFiber[ 3 * point + i ] +=
                                    atlasBundleData[ 3 * k + offsetAtlas + i ] ;

      }

    }

  }


  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      averageFiber[ 3 * point + i ] /= nbFibersAtlasBundle ;

    }

  }



  // Searching the medial point of average atlas bundle fiber
  atlasBundleData.computeMedialPointFiberWithDistance(
                                                      averageFiber,
                                                      medialPointAtlasBundle ) ;

}


//----------------------------------------------------------------------------//
//-------------------------- Find radius of bundle ---------------------------//
//----------------------------------------------------------------------------//
void computeDistancesToCenterBundle(
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        const std::vector<float>& medialPointAtlasBundle,
                        int nbFibersAtlasBundle,
                        std::vector<float>& distancesToCenterAtlasBundle )
{

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    float distanceToCenter = 0 ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      distanceToCenter += pow( medialPointAtlasBundle[ i ] -
           medialPointsAtlasBundleFibers[ 3 * atlasBundleFiberIndex + i ], 2 ) ;


    }

    distanceToCenter = sqrt( distanceToCenter ) ;

    distancesToCenterAtlasBundle[ atlasBundleFiberIndex ] = distanceToCenter ;

  }

}


//----------------------------------------------------------------------------//
//------------- Compute normal vector of fibers in atlas bundle --------------//
//----------------------------------------------------------------------------//
void computeNormalVectorFibersAtlasBundle(
                                       BundlesDataFormat& atlasBundleData,
                                       std::vector<float>& normalVectorsBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    std::vector<float> normalVector( 3, 0 ) ;
    atlasBundleData.computeNormalVectorFiberTractogram( atlasBundleFiberIndex,
                                                        normalVector ) ;

    normalVectorsBundle[ 3 * atlasBundleFiberIndex + 0 ] = normalVector[ 0 ] ;
    normalVectorsBundle[ 3 * atlasBundleFiberIndex + 1 ] = normalVector[ 1 ] ;
    normalVectorsBundle[ 3 * atlasBundleFiberIndex + 2 ] = normalVector[ 2 ] ;


  }

}

//----------------------------------------------------------------------------//
//--------------------------- Compute angle bundle ---------------------------//
//----------------------------------------------------------------------------//
void computeAnglesBundle( BundlesDataFormat& atlasBundleData,
                          const std::vector<float>& normalVectorsBundle,
                          std::vector<float>& anglesAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  int indexAnglesAtlasBundle = 0 ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex_1 = 0 ;
                             atlasBundleFiberIndex_1 < nbFibersAtlasBundle ;
                                                    atlasBundleFiberIndex_1++ )
  {

    std::vector<float> normalVectorBundleFiber_1( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      normalVectorBundleFiber_1[ i ] =
                      normalVectorsBundle[ 3 * atlasBundleFiberIndex_1 + i ] ;

    }

    if ( atlasBundleFiberIndex_1 + 1 < nbFibersAtlasBundle )
    {

      for ( int atlasBundleFiberIndex_2 = atlasBundleFiberIndex_1 + 1 ;
                                 atlasBundleFiberIndex_2 < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex_2++ )
      {

        if ( atlasBundleFiberIndex_2 != atlasBundleFiberIndex_1 )
        {

          std::vector<float> normalVectorBundleFiber_2( 3, 0 ) ;
          for ( int i = 0 ; i < 3 ; i++ )
          {

            normalVectorBundleFiber_2[ i ] =
                        normalVectorsBundle[ 3 * atlasBundleFiberIndex_2 + i ] ;

          }

          float angle = atlasBundleData.computeAngleBetweenVectors(
                                                   normalVectorBundleFiber_1,
                                                   normalVectorBundleFiber_2 ) ;

          indexAnglesAtlasBundle = atlasBundleFiberIndex_1 *
                                                 ( 2 * nbFibersAtlasBundle - 1 -
                                                 atlasBundleFiberIndex_1 ) / 2 +
                                                 atlasBundleFiberIndex_2 -
                                                   atlasBundleFiberIndex_1 - 1 ;

          if ( angle > 90 )
          {

            angle = 180 - angle ;

          }

          anglesAtlasBundle[ indexAnglesAtlasBundle ] = angle ;

        }

      }

    }

  }

}


//----------------------------------------------------------------------------//
// --------------------- Compute direction angles bundle -------------------- //
//----------------------------------------------------------------------------//

void computeDirectionAnglesBundle(
                                BundlesDataFormat& atlasBundleData,
                                const std::vector<float>& normalVectorsBundle,
                                int nbPoints,
                                std::vector<float>& directionAnglesAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  int indexAnglesAtlasBundle = 0 ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex_1 = 0 ;
                             atlasBundleFiberIndex_1 < nbFibersAtlasBundle ;
                                                    atlasBundleFiberIndex_1++ )
  {

    std::vector<float> normalVector1( 3, 0 ) ;
    for ( int i =0 ; i < 3 ; i++ )
    {

      normalVector1[ i ] = normalVectorsBundle[
                                             3 * atlasBundleFiberIndex_1 + i ] ;

    }

    std::vector<float> fiber1( 3 * nbPoints, 0 ) ;
    std::vector<float> newNormalVectorFiber1( 3, 0 ) ;
    atlasBundleData.putFiberInPlaneXY( normalVector1,
                                       atlasBundleData.matrixTracks,
                                       atlasBundleFiberIndex_1,
                                       nbPoints,
                                       fiber1,
                                       newNormalVectorFiber1 ) ;

    std::vector<float> directionVectorFiber1( 3, 0 ) ;
    atlasBundleData.computeDirectionVectorFiberTractogram(
                                                       fiber1,
                                                       newNormalVectorFiber1,
                                                       directionVectorFiber1 ) ;


    if ( atlasBundleFiberIndex_1 + 1 < nbFibersAtlasBundle )
    {

      for ( int atlasBundleFiberIndex_2 = atlasBundleFiberIndex_1 + 1 ;
                                 atlasBundleFiberIndex_2 < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex_2++ )
      {

        if ( atlasBundleFiberIndex_2 != atlasBundleFiberIndex_1 )
        {

          std::vector<float> normalVector2( 3, 0 ) ;
          for ( int i =0 ; i < 3 ; i++ )
          {

            normalVector2[ i ] = normalVectorsBundle[
                                             3 * atlasBundleFiberIndex_2 + i ] ;

          }

          std::vector<float> fiber2( 3 * nbPoints, 0 ) ;
          std::vector<float> newNormalVectorFiber2( 3, 0 ) ;
          atlasBundleData.putFiberInPlaneXY( normalVector2,
                                             atlasBundleData.matrixTracks,
                                             atlasBundleFiberIndex_2,
                                             nbPoints,
                                             fiber2,
                                             newNormalVectorFiber2 ) ;

          std::vector<float> directionVectorFiber2( 3, 0 ) ;
          atlasBundleData.computeDirectionVectorFiberTractogram(
                                                       fiber2,
                                                       newNormalVectorFiber2,
                                                       directionVectorFiber2 ) ;


          float angleBetweenPlanes = atlasBundleData.computeAngleBetweenPlanes(
                                                       newNormalVectorFiber1,
                                                       newNormalVectorFiber2 ) ;

          if ( angleBetweenPlanes > 5 )
          {

            std::cout << "\nERROR : could not align fibers, got minimum angle "
                      << "between planes of " << angleBetweenPlanes << "\n" ;
            exit( 1 ) ;

          }



          // Computing angle between directions
          float angle = atlasBundleData.computeAngleBetweenDirections(
                                                       directionVectorFiber1,
                                                       directionVectorFiber2 ) ;

          indexAnglesAtlasBundle = atlasBundleFiberIndex_1 *
                                                 ( 2 * nbFibersAtlasBundle - 1 -
                                                 atlasBundleFiberIndex_1 ) / 2 +
                                                  atlasBundleFiberIndex_2 -
                                                   atlasBundleFiberIndex_1 - 1 ;

          // if ( angle > 90 )
          // {
          //
          //   angle = 180 - angle ;
          //
          // }

          directionAnglesAtlasBundle[ indexAnglesAtlasBundle ] = angle ;

        }

      }

    }

  }

}

//----------------------------------------------------------------------------//
//--------------------- Compute direction angles bundle ----------------------//
//----------------------------------------------------------------------------//
void computeShapeAnglesBundle(
                        BundlesDataFormat& atlasBundleData,
                        int nbPoints,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        std::vector<float>& shapeAnglesAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    int offsetAtlas = 3 * nbPoints * atlasBundleFiberIndex ;

    std::vector<float> point1( 3, 0 ) ;
    std::vector<float> point2( 3, 0 ) ;

    std::vector<float> vector1( 3, 0 ) ;
    std::vector<float> vector2( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      point1[ i ] = atlasBundleData[ offsetAtlas + i ] ;
      point2[ i ] = atlasBundleData[ offsetAtlas + 3 * ( nbPoints - 1 ) + i ] ;

      vector1[ i ] = point1[ i ] - medialPointsAtlasBundleFibers[
                                               3 * atlasBundleFiberIndex + i ] ;

      vector2[ i ] = point2[ i ] - medialPointsAtlasBundleFibers[
                                               3 * atlasBundleFiberIndex + i ] ;

    }

    float angle = atlasBundleData.computeAngleBetweenVectors( vector1,
                                                                     vector2 ) ;

    shapeAnglesAtlasBundle[ atlasBundleFiberIndex ] = angle ;

  }

}

//----------------------------------------------------------------------------//
//------------------- Compute min similarity measure (dMDA) ------------------//
//----------------------------------------------------------------------------//
void computeAverageDisimilarity(
                        BundlesDataFormat& atlasBundleData,
                        int nbPoints,
                        std::vector<float>& disimilaritiesAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  int indexDisimilarityAtlasBundle = 0 ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex_1 = 0 ;
                             atlasBundleFiberIndex_1 < nbFibersAtlasBundle ;
                                                    atlasBundleFiberIndex_1++ )
  {

    if ( atlasBundleFiberIndex_1 + 1 < nbFibersAtlasBundle )
    {

      for ( int atlasBundleFiberIndex_2 = atlasBundleFiberIndex_1 + 1 ;
                    atlasBundleFiberIndex_2 < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex_2++ )
      {

        // Computing MDF distance
        float dMDF = atlasBundleData.computeMDFBetweenTwoFibers(
                                                   atlasBundleData.matrixTracks,
                                                   atlasBundleData.matrixTracks,
                                                   atlasBundleFiberIndex_1,
                                                   atlasBundleFiberIndex_2,
                                                   nbPoints ) ;


        indexDisimilarityAtlasBundle = atlasBundleFiberIndex_1 *
                                               ( 2 * nbFibersAtlasBundle - 1 -
                                                atlasBundleFiberIndex_1 ) / 2 +
                                                  atlasBundleFiberIndex_2 -
                                                  atlasBundleFiberIndex_1 - 1 ;

        disimilaritiesAtlasBundle[ indexDisimilarityAtlasBundle ] = dMDF ;

      }

    }

  }

}
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  int index_atlas, index_reference, index_averageFibers, index_useMDF,
                                                     index_verbose, index_help ;
  index_atlas =   getFlagPosition( argc, argv, "-a") ;
  index_reference =   getFlagPosition( argc, argv, "-r") ;
  index_averageFibers = getFlagPosition( argc, argv, "-af") ;
  index_useMDF = getFlagPosition( argc, argv, "-useMDF" ) ;
  index_verbose =   getFlagPosition( argc, argv, "-v") ;
  index_help =   getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-a : Directory with the atlas (one file per bundle) \n"
              << "[-r] : Reference .nii image where the atlas is \n"
              << "[-af] : Output directory where to save the average fibers of "
              << "the bundles \n"
              << "[-useMDF] : NOT IMPLEMENTED\n"
              << "[-v] : set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  std::string atlasDirectory ;
  if ( !index_atlas )
  {

    std::cout << "-a argument required ..." << std::endl ;
    exit( 1 ) ;

  }
  else
  {

    atlasDirectory = argv[ index_atlas + 1 ] ;
    char lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      atlasDirectory = atlasDirectory + "/" ;

    }

  }

  std::string niiFilename ;
  if ( index_reference )
  {

    niiFilename = argv[ index_reference + 1 ] ;
    char lastChar = niiFilename[ niiFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      niiFilename = niiFilename.substr( 0, niiFilename.size() - 1 ) ;

    }

  }

  if ( index_useMDF )
  {

    useMDFDistance = true ;

  }
  else
  {

    useMDFDistance = false ;

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

  // Cheking format reference image
  if ( index_reference )
  {

    if ( niiFilename.find( ".nii") == std::string::npos ||
                            niiFilename.find( ".nii.gz") != std::string::npos )
    {

      std::cout << "Error format : the only supported format for the reference "
                << "image is .nii " << std::endl ;
      exit( 1 ) ;

    }

  }

  std::string averageFibersDirectory ;
  if ( index_averageFibers )
  {

    averageFibersDirectory = argv[ index_averageFibers + 1 ] ;

    char lastChar = averageFibersDirectory[ averageFibersDirectory.size()
                                                                      - 1 ] ;
    if ( lastChar != '/' )
    {

      averageFibersDirectory = averageFibersDirectory + "/" ;

    }

  }



  //   xxxxxxxxxxxxxxxxxxxxxxxxxx Reading Atlas xxxxxxxxxxxxxxxxxxxxxxxxxx   //

  int sizeBytesFullAtlas = 0 ;
  AtlasBundles atlasData ;

  std::vector< std::string > atlasBundlesFilenames ;
  std::vector< std::string > atlasBundlesDataFilenames ;
  std::string tmpBundlesDataFilename ;
  std::string tmpBundlesFilename ;

  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                    atlasDirectory.c_str() ) )
  {

    tmpBundlesDataFilename = file.path() ;
    tmpBundlesFilename = tmpBundlesDataFilename ;
    std::string key (".bundlesdata") ;

    if ( tmpBundlesDataFilename.find( ".bundlesdata" ) != std::string::npos )
    {

      std::size_t found = tmpBundlesFilename.rfind( key ) ;
      tmpBundlesFilename.replace( found, key.length(), ".bundles" ) ;

      atlasBundlesDataFilenames.push_back( tmpBundlesDataFilename ) ;
      atlasBundlesFilenames.push_back( tmpBundlesFilename ) ;

    }
    else if ( tmpBundlesDataFilename.find( ".trk" ) != std::string::npos )
    {

      std::cout << "Format error : Input tractogram is .bundles but .trk "
                << "files found in the atlas directory. Projection between "
                << "different formats is not supported yet \n" ;
      exit( 1 ) ;

    }

  }

  int nbBundlesDataFiles = atlasBundlesDataFilenames.size() ;
  int nbBundlesFiles = atlasBundlesFilenames.size() ;

  if ( nbBundlesDataFiles != nbBundlesFiles )
  {

    std::cout << "There has to be a .bundles file for each .bundlesdata : \n"
              << "Number of .bundles files : " << nbBundlesFiles
              << "\nNumber of .bundlesdata files : " << nbBundlesDataFiles
              << std::endl ;
    exit( 1 ) ;

  }

  atlasData.bundles.resize( nbBundlesFiles ) ;
  atlasData.bundlesData.resize( nbBundlesFiles ) ;

  if ( verbose )
  {

    std::cout << "Reading atlas : " << atlasDirectory << std::endl ;
  }

  for ( int k = 0 ; k < nbBundlesFiles ; k++ )
  {

    if ( verbose > 1 )
    {

      printf("\rReading bundle [ %d / %d ]   |   %s", k + 1 , nbBundlesFiles,
                                          atlasBundlesFilenames[ k ].c_str() ) ;
      std::cout << "" << std::flush ;

    }
    else if ( verbose )
    {

      printf("\rReading bundle [ %d / %d ]", k + 1 , nbBundlesFiles ) ;
      std::cout << "" << std::flush ;

    }

    std::string bundlesFilename = atlasBundlesFilenames[ k ] ;
    std::string bundlesdataFilename = atlasBundlesDataFilenames[ k ] ;


    BundlesFormat bundleInfo( bundlesFilename.c_str(), 0 ) ;

    atlasData.bundles[ k ] = bundleInfo ;

    BundlesDataFormat bundleData( bundlesdataFilename.c_str(),
                                  bundlesFilename.c_str(),
                                  0 ) ;


    // atlasData.bundlesData.push_back( bundleData ) ;
    atlasData.bundlesData[ k ] = bundleData  ;

    for ( int fiber = 0 ; fiber < atlasData.bundlesData[ k ].curves_count ;
                                                                      fiber++ )
    {

      sizeBytesFullAtlas +=
                       3 * atlasData.bundlesData[ k ].pointsPerTrack[ fiber ] *
                                                               sizeof( float ) ;

    }

  }

  if ( verbose > 0 )
  {

    std::cout << "\n" ;

  }


  //////////////////////////////// Sanity cheks ////////////////////////////////
  int nbPoints = atlasData.bundlesData[ 0 ].pointsPerTrack[ 0 ] ;
  int nbBundlesAtlas = atlasData.bundlesData.size() ;


  for ( int bundleIndex = 0 ; bundleIndex < nbBundlesAtlas ; bundleIndex++ )
  {

    BundlesDataFormat& bundleFibers = atlasData.bundlesData[ bundleIndex ] ;

    int nbFibersBundle = bundleFibers.curves_count ;

    for ( int fiber = 1 ; fiber < nbFibersBundle ; fiber++ )
    {

      if ( bundleFibers.pointsPerTrack[ fiber ] != nbPoints )
      {

        std::cout << "Error atlas : the number of points in each fiber "
                  << "of the atlas has to be the same, got "
                  << bundleFibers.pointsPerTrack[ fiber ]
                  << " and " << nbPoints << " for fiber " << fiber
                  << " in bundle " << bundleIndex << " and for fiber 0 in "
                  << "bundle 0" << std::endl ;
        exit( 1 ) ;

      }

    }

  }



  ////////////////////////////// Analysing Atlas ///////////////////////////////

  // Looping over bundles
  for ( int atlasBundleIndex = 0 ; atlasBundleIndex < nbBundlesAtlas ;
                                                      atlasBundleIndex++ )
  {

    BundlesFormat& atlasBundleInfo = atlasData.bundles[ atlasBundleIndex ] ;
    BundlesDataFormat& atlasBundleData = atlasData.bundlesData[
                                                            atlasBundleIndex ] ;

    if ( verbose )
    {

      printf( "\rProcessing atlas bundles : [ %d  /  %d ]",
                                        atlasBundleIndex + 1, nbBundlesAtlas ) ;
      fflush( stdout ) ;

    }

    int nbFibersAtlasBundle = atlasBundleData.curves_count ;

    // Compute center points of atlas bundle fibers
    std::vector<float> medialPointsAtlasBundleFibers(
                                                  3 * nbFibersAtlasBundle, 0 ) ;

    computeCenterAtlasBundleFibers( atlasBundleData,
                                    medialPointsAtlasBundleFibers ) ;

    // Compute average atlas bundle fiber
    std::vector<float> averageFiber( 3 * nbPoints, 0 ) ;
    std::vector<float> medialPointAtlasBundle( 3, 0 ) ;
    computeAverageFiberBundle( atlasBundleData,
                               medialPointsAtlasBundleFibers,
                               nbPoints,
                               averageFiber,
                               medialPointAtlasBundle ) ;

    atlasBundleInfo.centerBundle[ 0 ] = medialPointAtlasBundle[ 0 ] ;
    atlasBundleInfo.centerBundle[ 1 ] = medialPointAtlasBundle[ 1 ] ;
    atlasBundleInfo.centerBundle[ 2 ] = medialPointAtlasBundle[ 2 ] ;

    if ( index_averageFibers )
    {

      std::string bundleName = getFilenameNoExtension(
                                  atlasBundlesFilenames[ atlasBundleIndex ] ) ;

      std::string averageFiberBundlesFilename = averageFibersDirectory +
                                         "average_" + bundleName + ".bundles" ;

      std::string averageFiberBundlesDataFilename = averageFibersDirectory +
                                     "average_" + bundleName + ".bundlesdata" ;

      BundlesFormat averageFiberBundlesInfo( atlasBundleInfo ) ;
      averageFiberBundlesInfo.curves_count = 1 ;

      BundlesDataFormat averageFiberBundlesData ;
      averageFiberBundlesData.curves_count = 1 ;

      averageFiberBundlesData.pointsPerTrack.push_back( nbPoints ) ;

      averageFiberBundlesData.matrixTracks = averageFiber ;

      averageFiberBundlesInfo.bundlesWriting(
                                            averageFiberBundlesFilename.c_str(),
                                            0 ) ;

      averageFiberBundlesData.bundlesdataWriting(
                                        averageFiberBundlesDataFilename.c_str(),
                                        0 ) ;

    }


    // Compute radius atlas bundle
    std::vector<float> distancesToCenterAtlasBundle( nbFibersAtlasBundle, 0 ) ;
    computeDistancesToCenterBundle( medialPointsAtlasBundleFibers,
                                    medialPointAtlasBundle,
                                    nbFibersAtlasBundle,
                                    distancesToCenterAtlasBundle ) ;

    // Since accumulate is used in vector of floats, it is necesary to put 0.0f
    // instead of 0 or to not put anything at all at the end
    float averageRadius = std::accumulate(
                                    distancesToCenterAtlasBundle.begin(),
                                    distancesToCenterAtlasBundle.end(), 0.0f ) /
                                    nbFibersAtlasBundle ;
    atlasBundleInfo.averageRadius = averageRadius ;
    float minRadius = *(std::min_element( distancesToCenterAtlasBundle.begin(),
                                        distancesToCenterAtlasBundle.end() )) ;
    atlasBundleInfo.minRadius = minRadius ;
    float maxRadius = *( std::max_element( distancesToCenterAtlasBundle.begin(),
                                        distancesToCenterAtlasBundle.end() ) ) ;
    atlasBundleInfo.maxRadius = maxRadius ;



    // Compute angle atlas bundle
    std::vector<float> normalVectorsBundle( 3 * nbFibersAtlasBundle, 0 ) ;
    computeNormalVectorFibersAtlasBundle( atlasBundleData,
                                          normalVectorsBundle ) ;

    int32_t numberAnglesAtlasBundle = ( nbFibersAtlasBundle - 1 ) *
                                              ( nbFibersAtlasBundle ) / 2 ;
    std::vector<float> anglesAtlasBundle( numberAnglesAtlasBundle, 0 ) ;
    computeAnglesBundle( atlasBundleData,
                         normalVectorsBundle,
                         anglesAtlasBundle ) ;

    float averageAngle = std::accumulate( anglesAtlasBundle.begin(),
                                          anglesAtlasBundle.end(), 0.0f ) /
                                          numberAnglesAtlasBundle ;
    atlasBundleInfo.averageAngle = averageAngle ;
    float minAngle = *( std::min_element( anglesAtlasBundle.begin(),
                                       anglesAtlasBundle.end() ) ) ;
    atlasBundleInfo.minAngle = minAngle ;
    float maxAngle = *( std::max_element( anglesAtlasBundle.begin(),
                                       anglesAtlasBundle.end() ) ) ;
    atlasBundleInfo.maxAngle = maxAngle ;



    // Compute direction angles atlas bundles
    std::vector<float> directionAnglesAtlasBundle( numberAnglesAtlasBundle,
                                                                           0 ) ;
    computeDirectionAnglesBundle( atlasBundleData,
                                  normalVectorsBundle,
                                  nbPoints,
                                  directionAnglesAtlasBundle ) ;

    float averageDirectionAngle = std::accumulate(
                                      directionAnglesAtlasBundle.begin(),
                                      directionAnglesAtlasBundle.end(), 0.0f ) /
                                      numberAnglesAtlasBundle ;
    atlasBundleInfo.averageDirectionAngle = averageDirectionAngle ;
    float minDirectionAngle = *( std::min_element(
                                          directionAnglesAtlasBundle.begin(),
                                          directionAnglesAtlasBundle.end() ) ) ;
    atlasBundleInfo.minDirectionAngle = minDirectionAngle ;
    float maxDirectionAngle = *( std::max_element(
                                          directionAnglesAtlasBundle.begin(),
                                          directionAnglesAtlasBundle.end() ) ) ;
    atlasBundleInfo.maxDirectionAngle = maxDirectionAngle ;



    // Compute the angle between the vectors formed by central point - first
    // point and central point - last points
    std::vector<float> shapeAnglesAtlasBundle( nbFibersAtlasBundle, 0 ) ;
    computeShapeAnglesBundle( atlasBundleData,
                              nbPoints,
                              medialPointsAtlasBundleFibers,
                              shapeAnglesAtlasBundle ) ;

    float averageShapeAngle = std::accumulate(
                                          shapeAnglesAtlasBundle.begin(),
                                          shapeAnglesAtlasBundle.end(), 0.0f ) /
                                          nbFibersAtlasBundle ;
    atlasBundleInfo.averageShapeAngle = averageShapeAngle ;
    float minShapeAngle = *( std::min_element( shapeAnglesAtlasBundle.begin(),
                                            shapeAnglesAtlasBundle.end() ) ) ;
    atlasBundleInfo.minShapeAngle = minShapeAngle ;
    float maxShapeAngle = *( std::max_element( shapeAnglesAtlasBundle.begin(),
                                            shapeAnglesAtlasBundle.end() ) ) ;
    atlasBundleInfo.maxShapeAngle = maxShapeAngle ;



    // Compute lenght fibers atlas bundle
    std::vector<float> lengthsAtlasBundle( nbFibersAtlasBundle, 0 ) ;
    atlasData.computeLengthsAtlasBundleFibers( atlasBundleData,
                                               nbPoints,
                                               lengthsAtlasBundle ) ;

    float averageLength = std::accumulate( lengthsAtlasBundle.begin(),
                                       lengthsAtlasBundle.end(), 0.0f ) /
                                       nbFibersAtlasBundle ;
    atlasBundleInfo.averageLength = averageLength ;
    averageLength /= nbFibersAtlasBundle ;
    float minLength = *( std::min_element( lengthsAtlasBundle.begin(),
                                        lengthsAtlasBundle.end() ) ) ;
    atlasBundleInfo.minLength = minLength ;
    float maxLength = *( std::max_element( lengthsAtlasBundle.begin(),
                                        lengthsAtlasBundle.end() ) ) ;
    atlasBundleInfo.maxLength = maxLength ;



    // Compute disimilarity atlas bundle
    int32_t numberDisimilaritiesAtlasBundle = ( nbFibersAtlasBundle - 1) *
                                              ( nbFibersAtlasBundle ) / 2 ;

    std::vector<float> disimilaritiesAtlasBundle(
                                          numberDisimilaritiesAtlasBundle, 0 ) ;
    computeAverageDisimilarity( atlasBundleData,
                                nbPoints,
                                disimilaritiesAtlasBundle ) ;

    float averageDisimilarity = std::accumulate(
                                       disimilaritiesAtlasBundle.begin(),
                                       disimilaritiesAtlasBundle.end(), 0.0f ) /
                                       numberDisimilaritiesAtlasBundle ;
    atlasBundleInfo.averageDisimilarity = averageDisimilarity ;
    float minDisimilarity = *( std::min_element(
                                           disimilaritiesAtlasBundle.begin(),
                                           disimilaritiesAtlasBundle.end() ) ) ;
    atlasBundleInfo.minDisimilarity = minDisimilarity ;
    float maxDisimilarity = *( std::max_element(
                                           disimilaritiesAtlasBundle.begin(),
                                           disimilaritiesAtlasBundle.end() ) ) ;
    atlasBundleInfo.maxDisimilarity = maxDisimilarity ;


    // Computing density
    float density = atlasBundleInfo.curves_count / (
                                          4 / 3 * M_PI * pow( maxRadius, 3 ) ) ;
    atlasBundleInfo.density = density ;



  }


  /////////////////////////////// Saving results ///////////////////////////////
  if ( verbose )
  {

    std::cout << "\n" ;

  }


  if ( verbose )
  {

    std::cout << "Checking other information in .bundles " << std::endl ;

  }

  // Checking other information in .bundles files
  if ( index_reference )
  {

    niiFormat niiHeader = readNiiHeader( niiFilename.c_str() ) ;
    for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        if ( atlasData.bundles[ bundle ].resolution[ i ] == 0 )
        {

          atlasData.bundles[ bundle ].resolution[ i ] =
                                                   niiHeader.pixdim[ i + 1 ] ;

        }

        if ( atlasData.bundles[ bundle ].size[ i ] == 0 )
        {

          atlasData.bundles[ bundle ].size[ i ] = niiHeader.dim[ i + 1 ] ;

        }


      }

      if ( atlasData.bundles[ bundle ].space_dimension == 0 )
      {

        atlasData.bundles[ bundle ].space_dimension =  3 ;

      }

    }

  }

  // Saving .bundles files

  if ( verbose )
  {

    std::cout << "Saving .bundles " << std::endl ;

  }

  for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
  {

    atlasData.bundles[ bundle ].bundlesWriting(
                                        atlasBundlesFilenames[ bundle ].c_str(),
                                        0 ) ;

  }

  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }

}
