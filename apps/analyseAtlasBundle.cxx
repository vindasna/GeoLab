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

#include "analyseAtlasBundle.h"
#include "ioWrapper.h"

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
////////////////////////////////// Functions ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------//
//-------------------------- Find radius of bundle ---------------------------//
//----------------------------------------------------------------------------//
void computeCenterAtlasBundleFibers(
                             BundlesData& atlasBundleData,
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
                        BundlesData& atlasBundleData,
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
//--------------------- Compute gravity center of bundle ---------------------//
//----------------------------------------------------------------------------//
void computeGravityCenterAtlasBundle(
                                  BundlesData& atlasBundleData,
                                  int nbPoints,
                                  std::vector<float>& gravityCenterAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterAtlasBundle[ i ] = 0 ;

  }

  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    int offsetAtlas = 3 * nbPoints * atlasBundleFiberIndex ;

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        gravityCenterAtlasBundle[ i ] +=
                                    atlasBundleData[ 3 * point + offsetAtlas + i ] ;

      }

    }

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterAtlasBundle[ i ] /= nbFibersAtlasBundle * nbPoints ;

  }

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
                                       BundlesData& atlasBundleData,
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
//--------- Compute distances between medial points fibers in bundle ---------//
//----------------------------------------------------------------------------//
void computeDistancesBetweenMedialPointsBundle(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        std::vector<float>& distancesBetweenMedialPointsBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  int indexDistanceBetweenMedialPointsBundle = 0 ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex_1 = 0 ;
                             atlasBundleFiberIndex_1 < nbFibersAtlasBundle ;
                                                    atlasBundleFiberIndex_1++ )
  {

    std::vector<float> medialPointFiber1( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointFiber1[ i ] = medialPointsAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex_1 + i ] ;

    }

    if ( atlasBundleFiberIndex_1 + 1 < nbFibersAtlasBundle )
    {

      for ( int atlasBundleFiberIndex_2 = atlasBundleFiberIndex_1 + 1 ;
                                 atlasBundleFiberIndex_2 < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex_2++ )
      {

        if ( atlasBundleFiberIndex_2 != atlasBundleFiberIndex_1 )
        {

          std::vector<float> medialPointFiber2( 3, 0 ) ;
          for ( int i = 0 ; i < 3 ; i++ )
          {

            medialPointFiber2[ i ] = medialPointsAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex_2 + i ] ;

          }

          float tmpDistance = 0 ;

          for ( int i = 0 ; i < 3 ; i++ )
          {

            tmpDistance += pow( medialPointFiber1[ i ] -
                                                   medialPointFiber2[ i ], 2 ) ;


          }

          tmpDistance = sqrt( tmpDistance ) ;

          indexDistanceBetweenMedialPointsBundle = atlasBundleFiberIndex_1 *
                                                 ( 2 * nbFibersAtlasBundle - 1 -
                                                 atlasBundleFiberIndex_1 ) / 2 +
                                                 atlasBundleFiberIndex_2 -
                                                   atlasBundleFiberIndex_1 - 1 ;


          distancesBetweenMedialPointsBundle[
                        indexDistanceBetweenMedialPointsBundle ] = tmpDistance ;

        }

      }

    }

  }

}

//----------------------------------------------------------------------------//
//--------------------------- Compute angle bundle ---------------------------//
//----------------------------------------------------------------------------//
void computeAnglesBundle( BundlesData& atlasBundleData,
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
                                BundlesData& atlasBundleData,
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
    atlasBundleData.getFiberFromTractogram( atlasBundleData.matrixTracks,
                                            atlasBundleFiberIndex_1,
                                            nbPoints,
                                            fiber1 ) ;

    std::vector<float> directionVectorFiber1( 3, 0 ) ;
    atlasBundleData.computeDirectionVectorFiberTractogram(
                                                       fiber1,
                                                       normalVector1,
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
          atlasBundleData.getFiberFromTractogram( atlasBundleData.matrixTracks,
                                                  atlasBundleFiberIndex_2,
                                                  nbPoints,
                                                  fiber2 ) ;

          std::vector<float> directionVectorFiber2( 3, 0 ) ;
          atlasBundleData.computeDirectionVectorFiberTractogram(
                                                       fiber2,
                                                       normalVector2,
                                                       directionVectorFiber2 ) ;



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
                        BundlesData& atlasBundleData,
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
                        BundlesData& atlasBundleData,
                        const std::vector<float>& normalVectorsBundle,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        bool useMeanForMDAD,
                        std::vector<float>& disimilaritiesAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  int indexDisimilarityAtlasBundle = 0 ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex_1 = 0 ;
                             atlasBundleFiberIndex_1 < nbFibersAtlasBundle ;
                                                    atlasBundleFiberIndex_1++ )
  {

    std::vector<float> normalVector1( 3 , 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      normalVector1[ i ] = normalVectorsBundle[
                                             3 * atlasBundleFiberIndex_1 + i ] ;

    }

    std::vector<float> medialPoint1( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPoint1[ i ] = medialPointsAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex_1 + i ] ;

    }

    std::vector<float> fiber1( 3 * nbPoints, 0 ) ;
    std::vector<float> newNormalVectorFiber1( 3, 0 ) ;
    atlasBundleData.registerFiberToPlaneXYAndDirectionX(
                            atlasBundleData.matrixTracks,
                            normalVector1,
                            medialPoint1,
                            atlasBundleFiberIndex_1,
                            nbPoints,
                            fiber1,
                            newNormalVectorFiber1 ) ;


    if ( atlasBundleFiberIndex_1 + 1 < nbFibersAtlasBundle )
    {

      for ( int atlasBundleFiberIndex_2 = atlasBundleFiberIndex_1 + 1 ;
                    atlasBundleFiberIndex_2 < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex_2++ )
      {

        std::vector<float> normalVector2( 3 , 0 ) ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          normalVector2[ i ] = normalVectorsBundle[
                                           3 * atlasBundleFiberIndex_2 + i ] ;

        }

        std::vector<float> medialPoint2( 3, 0 ) ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          medialPoint2[ i ] = medialPointsAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex_2 + i ] ;

        }

        std::vector<float> fiber2( 3 * nbPoints, 0 ) ;
        std::vector<float> newNormalVectorFiber2( 3, 0 ) ;
        atlasBundleData.registerFiberToPlaneXYAndDirectionX(
                                atlasBundleData.matrixTracks,
                                normalVector2,
                                medialPoint2,
                                atlasBundleFiberIndex_2,
                                nbPoints,
                                fiber2,
                                newNormalVectorFiber2 ) ;



        // Computing MDA distance
        std::vector<float> origin( 3, 0 ) ;
        float dMDA = atlasBundleData.computeMDADBetweenTwoFibersAfterAlignement( 
                                                                  fiber1,
                                                                  fiber2,
                                                                  0,
                                                                  0,
                                                                  useMeanForMDAD,
                                                                  nbPoints ) ;


        indexDisimilarityAtlasBundle = atlasBundleFiberIndex_1 *
                                               ( 2 * nbFibersAtlasBundle - 1 -
                                                atlasBundleFiberIndex_1 ) / 2 +
                                                  atlasBundleFiberIndex_2 -
                                                  atlasBundleFiberIndex_1 - 1 ;

        disimilaritiesAtlasBundle[ indexDisimilarityAtlasBundle ] = dMDA ;

      }

    }

  }

}

//----------------------------------------------------------------------------//
//------------------- Compute min similarity measure (MDF) ------------------//
//----------------------------------------------------------------------------//
void computeAverageDisimilarityMDF(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& normalVectorsBundle,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
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

    std::vector<float> medialPoint1( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPoint1[ i ] = medialPointsAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex_1 + i ] ;

    }

    if ( atlasBundleFiberIndex_1 + 1 < nbFibersAtlasBundle )
    {

      for ( int atlasBundleFiberIndex_2 = atlasBundleFiberIndex_1 + 1 ;
                    atlasBundleFiberIndex_2 < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex_2++ )
      {

        std::vector<float> medialPoint2( 3, 0 ) ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          medialPoint2[ i ] = medialPointsAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex_2 + i ] ;

        }

        // Computing MDA distance
        float dMDF = atlasBundleData.computeMDADBetweenTwoFibers(
                                                   atlasBundleData.matrixTracks,
                                                   atlasBundleData.matrixTracks,
                                                   medialPoint1,
                                                   medialPoint2,
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

  int index_atlas, index_outDir, index_format, index_useMDF, 
                              index_useMeanForMDAD, index_verbose, index_help ;
  index_atlas =   getFlagPosition( argc, argv, "-a") ;
  index_outDir = getFlagPosition( argc, argv, "-o") ;
  index_format = getFlagPosition( argc, argv, "-f") ;
  index_useMDF = getFlagPosition( argc, argv, "-useMDF") ;
  index_useMeanForMDAD = getFlagPosition( argc, argv, "-useMeanMDAD" ) ;
  index_verbose =   getFlagPosition( argc, argv, "-v") ;
  index_help =   getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-a : Directory with the atlas (one file per bundle) \n"
              << "-o : Output directory where to save the analysis \n"
              << "-f : Format of atlas bundles ( optios = [ .bundles, .trk, "
              << ".tck ] ) \n"
              << "[-useMDF] : analyse using MDF \n"
              << "[-useMeanMDAD] : Use mean instead of max for MDAD (default : false ) \n"
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

  std::string outputDirectory ;
  if ( index_outDir )
  {

    outputDirectory = argv[ index_outDir + 1 ] ;

    char lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      outputDirectory = outputDirectory + "/" ;

    }

  }
  else
  {

    std::cout << "-o argument required..." << std::endl ;
    exit( 1 ) ;

  }

  std::string format ;
  if ( index_format )
  {

    format = argv[ index_format + 1 ] ;

    if ( format != ".bundles" && format != ".trk" && format != ".tck" )
    {

      std::cout << "The only supported formats for the atlas bundles are "
                << ".bundles, .trk and .tck" << std::endl ;
      exit( 1 ) ;

    }

  }
  else
  {

    std::cout << "-f argument required..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_useMDF )
  {

    useMDF = true ;

  }

  if ( index_useMeanForMDAD )
  {

    std::string _tmpUseMeanForMDAD( argv[ index_useMeanForMDAD + 1 ] ) ;
    if ( _tmpUseMeanForMDAD == "true" )
    {

      useMeanForMDAD = true ;

    }
    else if ( _tmpUseMeanForMDAD == "false" )
    {

      useMeanForMDAD = false ;

    }
    else
    {

      std::cout << "Argument of -useMeanMDAD must be either \"true\" or \"false\" "
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



  //   xxxxxxxxxxxxxxxxxxxxxxxxxx Reading Atlas xxxxxxxxxxxxxxxxxxxxxxxxxx   //

  std::vector< std::string > atlasBundlesFilenames ;
  std::string tmpBundlesFilename ;

  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                    atlasDirectory.c_str() ) )
  {

    tmpBundlesFilename = file.path() ;

    if ( endswith( tmpBundlesFilename, format) )
    {
;
      atlasBundlesFilenames.push_back( tmpBundlesFilename ) ;

    }

  }

  bool isBundles = false ;
  bool isTrk = false ;
  bool isTck = false ;

  if ( format == ".bundles" )
  {

    isBundles = true ;

  }
  else if ( format == ".trk" )
  {

    isTrk = true ;

  }
  else
  {

    isTck = true ;

  }



  bool isBundlesFormat = true ;
  bool isTRKFormat = false ;
  AtlasBundles atlasData( atlasDirectory.c_str(),
                          isBundles,
                          isTrk,
                          isTck,
                          verbose ) ;



  //   xxxxxxxxxxxxxxxxxxxxxxxxx Analysing Atlas xxxxxxxxxxxxxxxxxxxxxxxxxx   //
  int nbPoints = atlasData.bundlesData[ 0 ].pointsPerTrack[ 0 ] ;


  // Looping over bundles
  int nbBundlesAtlas = atlasData.bundlesData.size() ;
  for ( int atlasBundleIndex = 0 ; atlasBundleIndex < nbBundlesAtlas ;
                                                      atlasBundleIndex++ )
  {

    BundlesData& atlasBundleData = atlasData.bundlesData[ atlasBundleIndex ] ;

    if ( verbose )
    {

      printf( "\rProcessing atlas bundles : [ %d  /  %d ]",
                                        atlasBundleIndex + 1, nbBundlesAtlas ) ;
      fflush( stdout ) ;

    }

    int nbFibersAtlasBundle =
                        atlasData.bundlesData[ atlasBundleIndex ].curves_count ;

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


    std::vector<float> gravityCenterAtlasBundle( 3, 0 ) ;
    computeGravityCenterAtlasBundle( atlasBundleData,
                                     nbPoints,
                                     gravityCenterAtlasBundle ) ;


    // Compute radius atlas bundle
    std::vector<float> distancesToCenterAtlasBundle( nbFibersAtlasBundle, 0 ) ;
    computeDistancesToCenterBundle( medialPointsAtlasBundleFibers,
                                    gravityCenterAtlasBundle,
                                    nbFibersAtlasBundle,
                                    distancesToCenterAtlasBundle ) ;

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

    // Compute direction angles atlas bundles
    std::vector<float> directionAnglesAtlasBundle( numberAnglesAtlasBundle,
                                                                           0 ) ;
    computeDirectionAnglesBundle( atlasBundleData,
                                  normalVectorsBundle,
                                  nbPoints,
                                  directionAnglesAtlasBundle ) ;

    // Compute the angle between the vectors formed by central point - first
    // point and central point - last points
    std::vector<float> shapeAnglesAtlasBundle( nbFibersAtlasBundle, 0 ) ;
    computeShapeAnglesBundle( atlasBundleData,
                              nbPoints,
                              medialPointsAtlasBundleFibers,
                              shapeAnglesAtlasBundle ) ;


    // Compute lenght fibers atlas bundle
    std::vector<float> lengthsAtlasBundle( nbFibersAtlasBundle, 0 ) ;
    atlasData.computeLengthsAtlasBundleFibers( atlasBundleData,
                                                     nbPoints,
                                                     lengthsAtlasBundle ) ;


    // Compute similarity atlas bundle
    int32_t numberDisimilaritiesAtlasBundle = ( nbFibersAtlasBundle - 1) *
                                              ( nbFibersAtlasBundle ) / 2 ;
    std::vector<float> disimilaritiesAtlasBundle(
                                          numberDisimilaritiesAtlasBundle, 0 ) ;
    if ( !useMDF )
    {

      computeAverageDisimilarity( atlasBundleData,
                                  normalVectorsBundle,
                                  medialPointsAtlasBundleFibers,
                                  nbPoints,
                                  useMeanForMDAD,
                                  disimilaritiesAtlasBundle ) ;


    }
    else
    {

      computeAverageDisimilarityMDF( atlasBundleData,
                                     normalVectorsBundle,
                                     medialPointsAtlasBundleFibers,
                                     nbPoints,
                                     disimilaritiesAtlasBundle ) ;

    }

    // Computing distance between medial points
    int32_t numberDistancesBetweenMedialPoints = ( nbFibersAtlasBundle - 1) *
                                                   ( nbFibersAtlasBundle ) / 2 ;
    std::vector<float> distancesBetweenMedialPointsBundle(
                                       numberDistancesBetweenMedialPoints, 0 ) ;
    computeDistancesBetweenMedialPointsBundle(
                                          atlasBundleData,
                                          medialPointsAtlasBundleFibers,
                                          distancesBetweenMedialPointsBundle ) ;


    // Saving analysis
    std::string bundleName = getFilenameNoExtension(
                                atlasBundlesFilenames[ atlasBundleIndex ] ) ;
    std::string analysisFilename = outputDirectory + bundleName + ".ana" ;

    atlasData.analysisWriting( analysisFilename.c_str(),
                               distancesToCenterAtlasBundle,
                               lengthsAtlasBundle,
                               disimilaritiesAtlasBundle,
                               anglesAtlasBundle,
                               directionAnglesAtlasBundle,
                               shapeAnglesAtlasBundle,
                               distancesBetweenMedialPointsBundle,
                               nbFibersAtlasBundle,
                               numberDisimilaritiesAtlasBundle,
                               numberDistancesBetweenMedialPoints,
                               verbose ) ;


  }

  std::cout << "\nDone" << std::endl ;

}
