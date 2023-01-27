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

#include "RecognizedBundles.h"


template< typename T >
void fillArrayPointer( T* array,
                       int64_t nbElements,
                       T x )
{

  for ( int i = 0 ; i < nbElements ; i++ )
  {

    array[ i ] = x ;

  }

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
RecognizedBundles::RecognizedBundles() : AtlasBundles() {}

RecognizedBundles::RecognizedBundles(
                                  int verbose,
                                  float p,
                                  float thrDistance ,
                                  float minimumLenghtFiber,
                                  float maximumLenghtFiber,
                                  float maxAngle,
                                  float maxDirectionAngle,
                                  float minShapeAngle,
                                  float maxShapeAngle,
                                  bool compareRecognizedWithAtlas,
                                  bool useMDFDistance,
                                  bool useMedialPointAverageFiber,
                                  bool useSimpleProjection,
                                  float toleranceP,
                                  float toleranceThr,
                                  float toleranceMaxAngle,
                                  float toleranceMaxDirectionAngle,
                                  float toleranceMinShapeAngle,
                                  float toleranceMaxShapeAngle,
                                  float toleranceLenght,
                                  float toleranceDistanceBetweenMedialPoints,
                                  float thrPercentageSimilarity,
                                  float thrDistanceBetweenMedialPoints,
                                  float minimumNumberFibers,
                                  float thresholdAdjacency,
                                  bool useDefautlP,
                                  bool useDefaultThr,
                                  bool useDefaultMaxAngle,
                                  bool useDefautlMaxDirectionAngle,
                                  bool useDefautlMinShapeAngle,
                                  bool useDefautlMaxShapeAngle,
                                  bool useDefaultThrDistanceBetweenMedialPoints,
                                  bool useDefaultMinLen,
                                  bool useDefaultMaxLen,
                                  bool saveExtractedBundles,
                                  bool saveUnlabeled,
				                          bool useAvgThr,
                                  int nbThreads ) : AtlasBundles()
{

  this->verbose = verbose ;
  this->p = p ;
  this->thrDistance = thrDistance ;
  this->minimumLenghtFiber = minimumLenghtFiber ;
  this->maximumLenghtFiber = maximumLenghtFiber ;
  this->maxAngle = maxAngle ;
  this->maxDirectionAngle = maxDirectionAngle ;
  this->minShapeAngle = minShapeAngle ;
  this->maxShapeAngle = maxShapeAngle ;
  this->compareRecognizedWithAtlas = compareRecognizedWithAtlas ;
  this->useMDFDistance = useMDFDistance ;
  this->useMedialPointAverageFiber = useMedialPointAverageFiber ;
  this->useSimpleProjection = useSimpleProjection ;
  this->toleranceP = toleranceP ;
  this->toleranceThr = toleranceThr ;
  this->toleranceMaxAngle = toleranceMaxAngle ;
  this->toleranceMaxDirectionAngle = toleranceMaxDirectionAngle ;
  this->toleranceMinShapeAngle = toleranceMinShapeAngle ;
  this->toleranceMaxShapeAngle = toleranceMaxShapeAngle ;
  this->toleranceLenght = toleranceLenght ;
  this->toleranceDistanceBetweenMedialPoints =
                                          toleranceDistanceBetweenMedialPoints ;
  this->thrPercentageSimilarity = thrPercentageSimilarity ;
  this->thrDistanceBetweenMedialPoints = thrDistanceBetweenMedialPoints ;
  this->minimumNumberFibers = minimumNumberFibers ;
  this->thresholdAdjacency = thresholdAdjacency ;
  this->useDefautlP = useDefautlP ;
  this->useDefaultThr = useDefaultThr ;
  this->useDefaultMaxAngle = useDefaultMaxAngle ;
  this->useDefautlMaxDirectionAngle = useDefautlMaxDirectionAngle ;
  this->useDefautlMinShapeAngle = useDefautlMinShapeAngle ;
  this->useDefautlMaxShapeAngle = useDefautlMaxShapeAngle ;
  this->useDefaultThrDistanceBetweenMedialPoints =
                                      useDefaultThrDistanceBetweenMedialPoints ;
  this->useDefaultMinLen = useDefaultMinLen ;
  this->useDefaultMaxLen = useDefaultMaxLen ;
  this->saveExtractedBundles = saveExtractedBundles ;
  this->saveUnlabeled = saveUnlabeled ;
  this->useAvgThr = useAvgThr ;

  this->nbThreads = nbThreads ;

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
RecognizedBundles::~RecognizedBundles(){}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void RecognizedBundles::MDADLabeling(
                         BundlesData& atlasBundleData,
                         BundlesData& subjectBundlesData,
                         int atlasBundleIndex,
                         const std::vector<float>& medialPointAtlasBundle,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         const std::vector<float>& normalVectorsAtlasBundle,
                         const std::vector<float>& directionVectorsAtlasBundle,
                         const std::vector<float>& lengthsAtlasBundleFibers,
                         int nbPoints,
                         std::vector<std::vector<int16_t>>& labels )
{

  if ( this->nbThreads )
  {

    omp_set_num_threads( this->nbThreads ) ;

  }

  float p = this->p ;
  float thrDistance = this->thrDistance ;
  float minimumLenghtFiber = this->minimumLenghtFiber ;
  float maximumLenghtFiber = this->maximumLenghtFiber ;
  float maxAngle = this->maxAngle ;
  float maxDirectionAngle = this->maxDirectionAngle ;
  float minShapeAngle = this->minShapeAngle ;
  float maxShapeAngle = this->maxShapeAngle ;
  float thrPercentageSimilarity = this->thrPercentageSimilarity ;
  float thrDistanceBetweenMedialPoints = this->thrDistanceBetweenMedialPoints ;

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;
  int nbFibersTract = subjectBundlesData.curves_count ;

  std::vector<float>& atlasBundleFibers = atlasBundleData.matrixTracks ;
  std::vector<float>& tractogramFibers = subjectBundlesData.matrixTracks ;

  std::vector<int16_t> tmpLabels( nbFibersTract, -1 ) ;

  #pragma omp parallel for shared(tmpLabels) schedule(dynamic)
  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    float angle = 90 ;
    float directionAngle = 180 ;

    float lengthFiber = subjectBundlesData.computeLengthFiber( tractogramFibers,
                                                               fiberIndex,
                                                               nbPoints ) ;
    if ( lengthFiber < minimumLenghtFiber || lengthFiber > maximumLenghtFiber )
    {

      continue ;

    }


    // Searching the medial point of tractogram fiber
    std::vector<float> medialPointTractFiber( 3, 0 ) ;
    subjectBundlesData.computeMedialPointFiberWithDistance(
                                                     fiberIndex,
                                                     medialPointTractFiber ) ;

    // Computing the distance between medial point of bundle and medial
    // point of tractogram fiber to see if the fibers are too far and save
    // in computing time

    float distanceMedialPoints = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      distanceMedialPoints += pow( ( medialPointAtlasBundle[ i ] -
                                             medialPointTractFiber[ i ] ), 2 ) ;

    }
    distanceMedialPoints = sqrt( distanceMedialPoints ) ;

    // Computing MDA when distance between middle points is below threshold
    if ( distanceMedialPoints < p )
    {

      // Loop over fibers in atlas bundle
      int nbAtlasBundleFibersSimilar = 0 ;

      for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex++ )
      {

        std::vector<float> medialPointAtlasFiber{
                                             medialPointAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex + 0 ],
                                             medialPointAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex + 1 ],
                                             medialPointAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex + 2 ] } ;

        //
        float distanceMedialPointsFiberSubjectToAtlas = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          distanceMedialPointsFiberSubjectToAtlas += pow(
                                          ( medialPointAtlasBundleFibers[ i ] -
                                          medialPointTractFiber[ i ] ), 2 ) ;

        }
        distanceMedialPointsFiberSubjectToAtlas = sqrt(
                                     distanceMedialPointsFiberSubjectToAtlas ) ;
        if ( distanceMedialPointsFiberSubjectToAtlas >
                                                thrDistanceBetweenMedialPoints )
        {

          continue ;

        }

        ////////////////////////////////////////////////////////////////////////
        // float lengthAltasFiber = atlasBundleData.computeLengthFiber(
        //                                                   atlasBundleFibers,
        //                                                   atlasBundleFiberIndex,
        //                                                   nbPoints ) ;
        //
        // if ( abs( lengthAltasFiber - lengthFiber ) > thrDistance / 2.0 )
        // {
        //
        //   continue ;
        //
        //
        // }

        ///
        std::vector<float> translation( 3, 0 ) ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          translation[ i ] = medialPointAtlasBundleFibers[ i ] -
                                                    medialPointTractFiber[ i ] ;

        }

        float tmpDirectDistance = 0 ;
        float tmpFlippedDistance = 0 ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          tmpDirectDistance += pow( atlasBundleFibers[
                            3 * nbPoints * atlasBundleFiberIndex + 3 * 0 + i ] -
                            ( tractogramFibers[  3 * nbPoints * fiberIndex +
                                         3 * 0 + i ] + translation[ i ] ), 2 ) ;

          tmpFlippedDistance += pow( atlasBundleFibers[
                            3 * nbPoints * atlasBundleFiberIndex + 3 * 0 + i ] -
                            ( tractogramFibers[  3 * nbPoints * fiberIndex +
                          3 * ( nbPoints - 1 ) + i ] + translation[ i ] ), 2 ) ;

        }

        bool isDirectSens = true ;
        if ( tmpFlippedDistance < tmpDirectDistance )
        {

          isDirectSens = false ;

        }
        ///

        std::vector<float> endPointAtlas1 ;
        std::vector<float> endPointAtlas2 ;
        std::vector<float> endPointTract1 ;
        std::vector<float> endPointTract2 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          endPointAtlas1.push_back( atlasBundleFibers[
                          3 * nbPoints * atlasBundleFiberIndex + 3 * 0 + i ] ) ;
          endPointAtlas2.push_back( atlasBundleFibers[
                          3 * nbPoints * atlasBundleFiberIndex + 3 * (
                                                        nbPoints - 1 ) + i ] ) ;

          endPointTract1.push_back( tractogramFibers[
                                     3 * nbPoints * fiberIndex + 3 * 0 + i ] ) ;
          endPointTract2.push_back( tractogramFibers[
                                     3 * nbPoints * fiberIndex + 3 * (
                                                        nbPoints - 1 ) + i ] ) ;
        }

        float tmpDist1 = 0 ;
        float tmpDist2 = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          if ( isDirectSens )
          {

            tmpDist1 += pow( endPointAtlas1[ i ] - endPointTract1[ i ], 2 ) ;
            tmpDist2 += pow( endPointAtlas2[ i ] - endPointTract2[ i ], 2 ) ;

          }
          else
          {

            tmpDist1 += pow( endPointAtlas1[ i ] - endPointTract2[ i ], 2 ) ;
            tmpDist2 += pow( endPointAtlas2[ i ] - endPointTract1[ i ], 2 ) ;

          }



        }
        tmpDist1 = sqrt( tmpDist1 ) ;
        tmpDist2 = sqrt( tmpDist2 ) ;

        if ( tmpDist1 > thrDistance / 2.0 || tmpDist2 > thrDistance / 2.0  )
        {

          continue ;

        }


        ////////////////////////////////////////////////////////////////////////



        // Computing shape angle
        std::vector<float> point1( 3, 0 ) ;
        std::vector<float> point2( 3, 0 ) ;

        std::vector<float> vector1( 3, 0 ) ;
        std::vector<float> vector2( 3, 0 ) ;
        int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          point1[ i ] = tractogramFibers[ 3 * 0 + i + offsetTractogram ] ;
          point2[ i ] = tractogramFibers[ 3 * ( nbPoints - 1 ) + i
                                                          + offsetTractogram ] ;

          vector1[ i ] = point1[ i ] - medialPointTractFiber[ i ] ;
          vector2[ i ] = point2[ i ] - medialPointTractFiber[ i ] ;

        }

        float shapeAngle = subjectBundlesData.computeAngleBetweenVectors(
                                                        vector1, vector2 ) ;

        if ( shapeAngle > maxShapeAngle || shapeAngle < minShapeAngle )
        {

          continue ;

        }
        //

        float dMDA = 0 ;

        // Computing angle between atlas fiber and tractogram fiber
        std::vector<float> normalVectorTractogramFiber( 3, 0 ) ;
        subjectBundlesData.computeNormalVectorFiberTractogram(
                                               tractogramFibers,
                                               medialPointTractFiber,
                                               fiberIndex,
                                               nbPoints,
                                               normalVectorTractogramFiber ) ;

        std::vector<float> normalVectorsAtlasBundleFiber( 3, 0 ) ;
        normalVectorsAtlasBundleFiber[ 0 ] = normalVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 0 ] ;
        normalVectorsAtlasBundleFiber[ 1 ] = normalVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 1 ] ;
        normalVectorsAtlasBundleFiber[ 2 ] = normalVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 2 ] ;

        angle = subjectBundlesData.computeAngleBetweenVectors(
                                             normalVectorTractogramFiber,
                                             normalVectorsAtlasBundleFiber ) ;


        if ( angle > maxAngle )
        {

          continue ;

        }

        std::vector<float> directionVectorTractogramFiber( 3, 0 ) ;
        subjectBundlesData.computeDirectionVectorFiberTractogram(
                                            tractogramFibers,
                                            medialPointTractFiber,
                                            normalVectorTractogramFiber,
                                            fiberIndex,
                                            nbPoints,
                                            directionVectorTractogramFiber ) ;

        std::vector<float> directionVectorsAtlasBundleFiber( 3, 0 ) ;
        directionVectorsAtlasBundleFiber[ 0 ] = directionVectorsAtlasBundle[
                                          3 * atlasBundleFiberIndex + 0 ] ;
        directionVectorsAtlasBundleFiber[ 1 ] = directionVectorsAtlasBundle[
                                          3 * atlasBundleFiberIndex + 1 ] ;
        directionVectorsAtlasBundleFiber[ 2 ] = directionVectorsAtlasBundle[
                                          3 * atlasBundleFiberIndex + 2 ] ;

        directionAngle = subjectBundlesData.computeAngleBetweenDirections(
                                          directionVectorTractogramFiber,
                                          directionVectorsAtlasBundleFiber ) ;

        if ( directionAngle > maxDirectionAngle )
        {

          continue ;

        }


        std::vector<float> tractogramFiberToAtlasFiber( 3 * nbPoints, 0 ) ;
        std::vector<float> normalVectorTractogramFiberToAtlasFiber( 3, 0 ) ;
        subjectBundlesData.registerFiber(
                                 atlasBundleData.matrixTracks,
                                 subjectBundlesData.matrixTracks,
                                 atlasBundleFiberIndex,
                                 fiberIndex,
                                 nbPoints,
                                 tractogramFiberToAtlasFiber,
                                 normalVectorTractogramFiberToAtlasFiber ) ;


        // Computing MDA distance
        dMDA = subjectBundlesData.computeMDADBetweenTwoFibers(
                                            tractogramFiberToAtlasFiber,
                                            atlasBundleFibers,
                                            medialPointAtlasFiber,
                                            medialPointAtlasFiber,
                                            0,
                                            atlasBundleFiberIndex,
                                            nbPoints ) ;


        if ( dMDA < thrDistance )
        {

          nbAtlasBundleFibersSimilar += 1 ;
          float percentageFibersSimilar = ( float )nbAtlasBundleFibersSimilar
                                                         / nbFibersAtlasBundle ;


	        if ( percentageFibersSimilar > thrPercentageSimilarity || ( dMDA < 1
	        		      && distanceMedialPoints < thrDistanceBetweenMedialPoints
                    && directionAngle < 5
                    && angle < 5 ) )
          {

            tmpLabels[ fiberIndex ] = atlasBundleIndex ;
            break ; // This break is for the for loop of the atlas bundle
                    // fibers, it does not break the for of openMP
                    // hence it is allowed

          }

        }

      }

    }

  }

  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    if ( tmpLabels[ fiberIndex ] > -1 )
    {

      labels[ fiberIndex ].push_back( tmpLabels[ fiberIndex ] ) ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
void RecognizedBundles::MDADLabelingSimple(
                         BundlesData& atlasBundleData,
                         BundlesData& subjectBundlesData,
                         int atlasBundleIndex,
                         const std::vector<float>& medialPointAtlasBundle,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         const std::vector<float>& lengthsAtlasBundleFibers,
                         int nbPoints,
                         std::vector<std::vector<int16_t>>& labels )
{

  if ( this->nbThreads )
  {

    omp_set_num_threads( this->nbThreads ) ;

  }

  float p = this->p ;
  float thrDistance = this->thrDistance ;
  float minimumLenghtFiber = this->minimumLenghtFiber ;
  float maximumLenghtFiber = this->maximumLenghtFiber ;
  float thrPercentageSimilarity = this->thrPercentageSimilarity ;
  float thrDistanceBetweenMedialPoints = this->thrDistanceBetweenMedialPoints ;

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;
  int nbFibersTract = subjectBundlesData.curves_count ;

  std::vector<float>& atlasBundleFibers = atlasBundleData.matrixTracks ;
  std::vector<float>& tractogramFibers = subjectBundlesData.matrixTracks ;

  std::vector<int16_t> tmpLabels( nbFibersTract, -1 ) ;

  #pragma omp parallel for
  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    float lengthFiber = subjectBundlesData.computeLengthFiber( tractogramFibers,
                                                               fiberIndex,
                                                               nbPoints ) ;
    if ( lengthFiber < minimumLenghtFiber || lengthFiber > maximumLenghtFiber )
    {

      continue ;

    }

    // Searching the medial point of tractogram fiber
    std::vector<float> medialPointTractFiber( 3, 0 ) ;
    subjectBundlesData.computeMedialPointFiberWithDistance(
                                                     fiberIndex,
                                                     medialPointTractFiber ) ;

    // Computing the distance between medial point of bundle and medial
    // point of tractogram fiber to see if the fibers are too far and save
    // in computing time

    float distanceMedialPoints = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      distanceMedialPoints += pow( ( medialPointAtlasBundle[ i ] -
                                        medialPointTractFiber[ i ] ), 2 ) ;

    }
    distanceMedialPoints = sqrt( distanceMedialPoints ) ;

    // Computing MDA when distance between middle points is below threshold
    if ( distanceMedialPoints < p )
    {

      // Loop over fibers in atlas bundle
      int nbAtlasBundleFibersSimilar = 0 ;

      for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex++ )
      {

        std::vector<float> medialPointAtlasFiber{
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 0 ],
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 1 ],
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 2 ] } ;
        //
        float distanceMedialPointsFiberSubjectToAtlas = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          distanceMedialPointsFiberSubjectToAtlas += pow(
                                          ( medialPointAtlasBundleFibers[ i ] -
                                          medialPointTractFiber[ i ] ), 2 ) ;

        }
        distanceMedialPointsFiberSubjectToAtlas = sqrt(
                                     distanceMedialPointsFiberSubjectToAtlas ) ;
        if ( distanceMedialPointsFiberSubjectToAtlas >
                                                thrDistanceBetweenMedialPoints )
        {

          continue ;

        }
        //

        float dMDA = 0 ;

        std::vector<float> medialPointFiberAtlasBundle( 3, 0 ) ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          medialPointFiberAtlasBundle[ i ] = medialPointAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex + i ] ;

        }

        // Computing MDA distance
        dMDA = subjectBundlesData.computeMDADBetweenTwoFibers(
                                                  tractogramFibers,
                                                  atlasBundleFibers,
                                                  medialPointTractFiber,
                                                  medialPointFiberAtlasBundle,
                                                  fiberIndex,
                                                  atlasBundleFiberIndex,
                                                  nbPoints ) ;

        if ( dMDA < thrDistance )
        {

          nbAtlasBundleFibersSimilar += 1 ;
          float percentageFibersSimilar = ( float )nbAtlasBundleFibersSimilar /
                                                          nbFibersAtlasBundle ;

          // if ( percentageFibersSimilar > thrPercentageSimilarity ||
          //                                                          dMDA == 0 )
          // if ( percentageFibersSimilar > thrPercentageSimilarity )
	        if ( percentageFibersSimilar > thrPercentageSimilarity || ( dMDA < 1
                                                 && distanceMedialPoints < 1 ) )
          {

            // #pragma omp critical
            // {
            //
            //   labels[ fiberIndex ].push_back( atlasBundleIndex ) ;
            //
            // }
            tmpLabels[ fiberIndex ] = atlasBundleIndex ;
            break ; // This break is for the for loop of the atlas bundle
                    // fibers, it does not break the for of openMP
                    // hence it is allowed

          }

        }

      }

    }

  }

  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    if ( tmpLabels[ fiberIndex ] > -1 )
    {

      labels[ fiberIndex ].push_back( tmpLabels[ fiberIndex ] ) ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
void RecognizedBundles::MDFLabeling(
                         BundlesData& atlasBundleData,
                         BundlesData& subjectBundlesData,
                         int atlasBundleIndex,
                         const std::vector<float>& medialPointAtlasBundle,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         const std::vector<float>& normalVectorsAtlasBundle,
                         const std::vector<float>& directionVectorsAtlasBundle,
                         const std::vector<float>& lengthsAtlasBundleFibers,
                         int nbPoints,
                         std::vector<std::vector<int16_t>>& labels )
{

  if ( this->nbThreads )
  {

    omp_set_num_threads( this->nbThreads ) ;

  }

  float p = this->p ;
  float thrDistance = this->thrDistance ;
  float minimumLenghtFiber = this->minimumLenghtFiber ;
  float maximumLenghtFiber = this->maximumLenghtFiber ;
  float maxAngle = this->maxAngle ;
  float maxDirectionAngle = this->maxDirectionAngle ;
  float minShapeAngle = this->minShapeAngle ;
  float maxShapeAngle = this->maxShapeAngle ;
  float thrPercentageSimilarity = this->thrPercentageSimilarity ;
  float thrDistanceBetweenMedialPoints = this->thrDistanceBetweenMedialPoints ;

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;
  int nbFibersTract = subjectBundlesData.curves_count ;

  std::vector<float>& atlasBundleFibers = atlasBundleData.matrixTracks ;
  std::vector<float>& tractogramFibers = subjectBundlesData.matrixTracks ;

  std::vector<int16_t> tmpLabels( nbFibersTract, -1 ) ;

  #pragma omp parallel for
  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    float lengthFiber = subjectBundlesData.computeLengthFiber( tractogramFibers,
                                                               fiberIndex,
                                                               nbPoints ) ;
    if ( lengthFiber < minimumLenghtFiber || lengthFiber > maximumLenghtFiber )
    {

      continue ;

    }

    // Searching the medial point of tractogram fiber
    std::vector<float> medialPointTractFiber( 3, 0 ) ;
    subjectBundlesData.computeMedialPointFiberWithDistance(
                                                     fiberIndex,
                                                     medialPointTractFiber ) ;

    // Computing the distance between medial point of bundle and medial
    // point of tractogram fiber to see if the fibers are too far and save
    // in computing time

    float distanceMedialPoints = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      distanceMedialPoints += pow( ( medialPointAtlasBundle[ i ] -
                                        medialPointTractFiber[ i ] ), 2 ) ;

    }
    distanceMedialPoints = sqrt( distanceMedialPoints ) ;

    // Computing MDF when distance between middle points is below threshold
    if ( distanceMedialPoints < p )
    {

      // Loop over fibers in atlas bundle
      int nbAtlasBundleFibersSimilar = 0 ;

      for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex++ )
      {

        std::vector<float> medialPointAtlasFiber{
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 0 ],
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 1 ],
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 2 ] } ;
        //
        float distanceMedialPointsFiberSubjectToAtlas = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          distanceMedialPointsFiberSubjectToAtlas += pow(
                                          ( medialPointAtlasBundleFibers[ i ] -
                                          medialPointTractFiber[ i ] ), 2 ) ;

        }
        distanceMedialPointsFiberSubjectToAtlas = sqrt(
                                     distanceMedialPointsFiberSubjectToAtlas ) ;
        if ( distanceMedialPointsFiberSubjectToAtlas >
                                                thrDistanceBetweenMedialPoints )
        {

          continue ;

        }


        //
        std::vector<float> translation( 3, 0 ) ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          translation[ i ] = medialPointAtlasBundleFibers[ i ] -
                                                    medialPointTractFiber[ i ] ;

        }

        float tmpDirectDistance = 0 ;
        float tmpFlippedDistance = 0 ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          tmpDirectDistance += pow( atlasBundleFibers[
                            3 * nbPoints * atlasBundleFiberIndex + 3 * 0 + i ] -
                            ( tractogramFibers[  3 * nbPoints * fiberIndex +
                                         3 * 0 + i ] + translation[ i ] ), 2 ) ;

          tmpFlippedDistance += pow( atlasBundleFibers[
                            3 * nbPoints * atlasBundleFiberIndex + 3 * 0 + i ] -
                            ( tractogramFibers[  3 * nbPoints * fiberIndex +
                          3 * ( nbPoints - 1 ) + i ] + translation[ i ] ), 2 ) ;

        }

        bool isDirectSens = true ;
        if ( tmpFlippedDistance < tmpDirectDistance )
        {

          isDirectSens = false ;

        }
        ///

        std::vector<float> endPointAtlas1 ;
        std::vector<float> endPointAtlas2 ;
        std::vector<float> endPointTract1 ;
        std::vector<float> endPointTract2 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          endPointAtlas1.push_back( atlasBundleFibers[
                          3 * nbPoints * atlasBundleFiberIndex + 3 * 0 + i ] ) ;
          endPointAtlas2.push_back( atlasBundleFibers[
                          3 * nbPoints * atlasBundleFiberIndex + 3 * (
                                                        nbPoints - 1 ) + i ] ) ;

          endPointTract1.push_back( tractogramFibers[
                                     3 * nbPoints * fiberIndex + 3 * 0 + i ] ) ;
          endPointTract2.push_back( tractogramFibers[
                                     3 * nbPoints * fiberIndex + 3 * (
                                                        nbPoints - 1 ) + i ] ) ;
        }

        float tmpDist1 = 0 ;
        float tmpDist2 = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          if ( isDirectSens )
          {

            tmpDist1 += pow( endPointAtlas1[ i ] - endPointTract1[ i ], 2 ) ;
            tmpDist2 += pow( endPointAtlas2[ i ] - endPointTract2[ i ], 2 ) ;

          }
          else
          {

            tmpDist1 += pow( endPointAtlas1[ i ] - endPointTract2[ i ], 2 ) ;
            tmpDist2 += pow( endPointAtlas2[ i ] - endPointTract1[ i ], 2 ) ;

          }



        }
        tmpDist1 = sqrt( tmpDist1 ) ;
        tmpDist2 = sqrt( tmpDist2 ) ;

        if ( tmpDist1 > thrDistance / 2.0 || tmpDist2 > thrDistance / 2.0  )
        {

          continue ;

        }
        //

        float dMDF = 0 ;

        // Computing angle between atlas fiber and tractogram fiber

        std::vector<float> normalVectorTractogramFiber( 3, 0 ) ;
        subjectBundlesData.computeNormalVectorFiberTractogram(
                                               tractogramFibers,
                                               medialPointTractFiber,
                                               fiberIndex,
                                               nbPoints,
                                               normalVectorTractogramFiber ) ;

        std::vector<float> normalVectorsAtlasBundleFiber( 3, 0 ) ;
        normalVectorsAtlasBundleFiber[ 0 ] = normalVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 0 ] ;
        normalVectorsAtlasBundleFiber[ 1 ] = normalVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 1 ] ;
        normalVectorsAtlasBundleFiber[ 2 ] = normalVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 2 ] ;

        float angle = subjectBundlesData.computeAngleBetweenVectors(
                                             normalVectorTractogramFiber,
                                             normalVectorsAtlasBundleFiber ) ;


        // if ( angle < maxAngle || 180 - angle < maxAngle )
        if ( angle < maxAngle )
        {

          std::vector<float> tractogramFiberToAtlasFiberPlane( 3 * nbPoints,
                                                                           0 ) ;
          // std::vector<float> newNormalVectorTractogramFiber( 3, 0 ) ;
          // subjectBundlesData.putFibersInSamePlane(
          //                                   normalVectorsAtlasBundleFiber,
          //                                   normalVectorTractogramFiber,
          //                                   tractogramFibers,
          //                                   fiberIndex,
          //                                   nbPoints,
          //                                   tractogramFiberToAtlasFiberPlane,
          //                                   newNormalVectorTractogramFiber ) ;

          std::vector<float> directionVectorTractogramFiber( 3, 0 ) ;
          // subjectBundlesData.computeDirectionVectorFiberTractogram(
          //                                   tractogramFiberToAtlasFiberPlane,
          //                                   newNormalVectorTractogramFiber,
          //                                   directionVectorTractogramFiber ) ;
          subjectBundlesData.computeDirectionVectorFiberTractogram(
                                              tractogramFibers,
                                              medialPointTractFiber,
                                              normalVectorTractogramFiber,
                                              fiberIndex,
                                              nbPoints,
                                              directionVectorTractogramFiber ) ;

          std::vector<float> directionVectorsAtlasBundleFiber( 3, 0 ) ;
          directionVectorsAtlasBundleFiber[ 0 ] = directionVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 0 ] ;
          directionVectorsAtlasBundleFiber[ 1 ] = directionVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 1 ] ;
          directionVectorsAtlasBundleFiber[ 2 ] = directionVectorsAtlasBundle[
                                            3 * atlasBundleFiberIndex + 2 ] ;



          float directionAngle =
                               subjectBundlesData.computeAngleBetweenDirections(
                                          directionVectorTractogramFiber,
                                          directionVectorsAtlasBundleFiber ) ;

          if ( directionAngle < maxDirectionAngle )
          {

            std::vector<float> tractogramFiberToAtlasFiber( 3 * nbPoints, 0 ) ;
            // subjectBundlesData.putFibersInSameDirection(
            //                                  normalVectorsAtlasBundleFiber,
            //                                  newNormalVectorTractogramFiber,
            //                                  directionVectorsAtlasBundleFiber,
            //                                  tractogramFiberToAtlasFiberPlane,
            //                                  nbPoints,
            //                                  tractogramFiberToAtlasFiber ) ;
            std::vector<float> normalVectorTractogramFiberToAtlasFiber( 3, 0 ) ;
            subjectBundlesData.registerFiber(
                                     atlasBundleData.matrixTracks,
                                     subjectBundlesData.matrixTracks,
                                     atlasBundleFiberIndex,
                                     fiberIndex,
                                     nbPoints,
                                     tractogramFiberToAtlasFiber,
                                     normalVectorTractogramFiberToAtlasFiber ) ;

            std::vector<float> point1( 3, 0 ) ;
            std::vector<float> point2( 3, 0 ) ;

            std::vector<float> vector1( 3, 0 ) ;
            std::vector<float> vector2( 3, 0 ) ;
            for ( int i = 0 ; i < 3 ; i++ )
            {

              point1[ i ] = tractogramFiberToAtlasFiber[ i ] ;
              point2[ i ] = tractogramFiberToAtlasFiber[ 3 * ( nbPoints - 1 )
                                                                       + i ] ;

              vector1[ i ] = point1[ i ] - medialPointTractFiber[ i ] ;

              vector2[ i ] = point2[ i ] - medialPointTractFiber[ i ] ;

            }

            float shapeAngle = subjectBundlesData.computeAngleBetweenVectors(
                                                          vector1, vector2 ) ;

            if ( shapeAngle < maxShapeAngle && shapeAngle > minShapeAngle )
            {

              std::vector<float> medialPointFiberAtlasBundle( 3, 0 ) ;
              for ( int i = 0 ; i < 3 ; i++ )
              {

                medialPointFiberAtlasBundle[ i ] = medialPointAtlasBundleFibers[
                                               3 * atlasBundleFiberIndex + i ] ;

              }

              // Computing MDF distance
              dMDF = subjectBundlesData.computeMDFBetweenTwoFibers(
                                                  tractogramFiberToAtlasFiber,
                                                  atlasBundleFibers,
                                                  medialPointFiberAtlasBundle,
                                                  medialPointFiberAtlasBundle,
                                                  0,
                                                  atlasBundleFiberIndex,
                                                  nbPoints ) ;

            }
            else
            {

              dMDF = 50 ;

            }

          }
          else
          {
             dMDF = 50 ;

          }

        }
        else
        {

          dMDF = 50 ;

        }

        if ( dMDF < thrDistance )
        {

          nbAtlasBundleFibersSimilar += 1 ;
          float percentageFibersSimilar = ( float )nbAtlasBundleFibersSimilar
                                                       / nbFibersAtlasBundle ;

          // if ( percentageFibersSimilar > thrPercentageSimilarity )
          if ( percentageFibersSimilar > thrPercentageSimilarity ||
                                                                   dMDF == 0 )
          {

            // #pragma omp critical
            // {
            //
            //   labels[ fiberIndex ].push_back( atlasBundleIndex ) ;
            //
            // }
            tmpLabels[ fiberIndex ] = atlasBundleIndex ;
            break ; // This break is for the for loop of the atlas bundle
                    // fibers, it does not break the for of openMP
                    // hence it is allowed

          }

        }

      }

    }

  }

  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    if ( tmpLabels[ fiberIndex ] > -1 )
    {

      labels[ fiberIndex ].push_back( tmpLabels[ fiberIndex ] ) ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
void RecognizedBundles::MDFLabelingSimple(
                         BundlesData& atlasBundleData,
                         BundlesData& subjectBundlesData,
                         int atlasBundleIndex,
                         const std::vector<float>& medialPointAtlasBundle,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         const std::vector<float>& lengthsAtlasBundleFibers,
                         int nbPoints,
                         std::vector<std::vector<int16_t>>& labels )
{

  if ( this->nbThreads )
  {

    omp_set_num_threads( this->nbThreads ) ;

  }

  // if ( this->nbThreads <= 0 )
  // {
  //
  //   this->nbThreads = omp_get_num_procs() ;
  //
  // }
  //
  // omp_set_num_threads( this->nbThreads ) ;

  float p = this->p ;
  float thrDistance = this->thrDistance ;
  float minimumLenghtFiber = this->minimumLenghtFiber ;
  float maximumLenghtFiber = this->maximumLenghtFiber ;
  float thrPercentageSimilarity = this->thrPercentageSimilarity ;
  float thrDistanceBetweenMedialPoints = this->thrDistanceBetweenMedialPoints ;

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;
  int nbFibersTract = subjectBundlesData.curves_count ;

  std::vector<float>& atlasBundleFibers = atlasBundleData.matrixTracks ;
  std::vector<float>& tractogramFibers = subjectBundlesData.matrixTracks ;

  std::vector<int16_t> tmpLabels( nbFibersTract, -1 ) ;

  #pragma omp parallel for
  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    float lengthFiber = subjectBundlesData.computeLengthFiber( tractogramFibers,
                                                               fiberIndex,
                                                               nbPoints ) ;
    if ( lengthFiber < minimumLenghtFiber || lengthFiber > maximumLenghtFiber )
    {

      continue ;

    }

    // Searching the medial point of tractogram fiber
    std::vector<float> medialPointTractFiber( 3, 0 ) ;
    subjectBundlesData.computeMedialPointFiberWithDistance(
                                                     fiberIndex,
                                                     medialPointTractFiber ) ;

    // Computing the distance between medial point of bundle and medial
    // point of tractogram fiber to see if the fibers are too far and save
    // in computing time

    float distanceMedialPoints = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      distanceMedialPoints += pow( ( medialPointAtlasBundle[ i ] -
                                        medialPointTractFiber[ i ] ), 2 ) ;

    }
    distanceMedialPoints = sqrt( distanceMedialPoints ) ;

    // Computing MDF when distance between middle points is below threshold
    if ( distanceMedialPoints < p )
    {

      // Loop over fibers in atlas bundle
      int nbAtlasBundleFibersSimilar = 0 ;

      for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                     atlasBundleFiberIndex++ )
      {

        std::vector<float> medialPointAtlasFiber{
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 0 ],
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 1 ],
                                           medialPointAtlasBundleFibers[
                                           3 * atlasBundleFiberIndex + 2 ] } ;
        //
        float distanceMedialPointsFiberSubjectToAtlas = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          distanceMedialPointsFiberSubjectToAtlas += pow(
                                          ( medialPointAtlasBundleFibers[ i ] -
                                          medialPointTractFiber[ i ] ), 2 ) ;

        }
        distanceMedialPointsFiberSubjectToAtlas = sqrt(
                                     distanceMedialPointsFiberSubjectToAtlas ) ;
        if ( distanceMedialPointsFiberSubjectToAtlas >
                                                thrDistanceBetweenMedialPoints )
        {

          continue ;

        }
        //

        float dMDF = 0 ;

        std::vector<float> medialPointFiberAtlasBundle( 3, 0 ) ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          medialPointFiberAtlasBundle[ i ] = medialPointAtlasBundleFibers[
                                             3 * atlasBundleFiberIndex + i ] ;

        }

        // Computing MDF distance
        dMDF = subjectBundlesData.computeMDFBetweenTwoFibers(
                                                  tractogramFibers,
                                                  atlasBundleFibers,
                                                  medialPointTractFiber,
                                                  medialPointFiberAtlasBundle,
                                                  fiberIndex,
                                                  atlasBundleFiberIndex,
                                                  nbPoints ) ;


        if ( dMDF < thrDistance )
        {

          nbAtlasBundleFibersSimilar += 1 ;
          float percentageFibersSimilar = ( float )nbAtlasBundleFibersSimilar
                                                       / nbFibersAtlasBundle ;

          // if ( percentageFibersSimilar > thrPercentageSimilarity )
          if ( percentageFibersSimilar > thrPercentageSimilarity ||
                                                                  dMDF == 0 )
          {

            // #pragma omp critical
            // {
            //
            //   labels[ fiberIndex ].push_back( atlasBundleIndex ) ;
            //
            //
            // }
            tmpLabels[ fiberIndex ] = atlasBundleIndex ;
            break ; // This break is for the for loop of the atlas bundle
                    // fibers, it does not break the for of openMP
                    // hence it is allowed

          }

        }

      }

    }

  }

  for ( int fiberIndex = 0 ; fiberIndex < nbFibersTract ; fiberIndex++ )
  {

    if ( tmpLabels[ fiberIndex ] > -1 )
    {

      labels[ fiberIndex ].push_back( tmpLabels[ fiberIndex ] ) ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
void RecognizedBundles::saveLabels(
                              const char* labelsBinaryFilename,
                              const std::vector<std::vector<int16_t>>& labels  )
{

  std::ofstream file ;
  file.open( labelsBinaryFilename ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << labelsBinaryFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  int n_curves = labels.size() ;

  for ( int track = 0 ; track < n_curves ; track++ )
  {

    int nbLabels = labels[ track ].size() ;
    for ( int _labelIndex = 0 ; _labelIndex < nbLabels ; _labelIndex++ )
    {

      int label = labels[ track ][ _labelIndex ] ;

      file << track << " : " << label << std::endl ;


    }

  }

  file.close() ;

}


////////////////////////////////////////////////////////////////////////////////
void RecognizedBundles::saveLabelsDict(
                                const char* labelsDictFilename,
                                const std::vector< std::string >& bundlesNames )
{

  std::ofstream file( labelsDictFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path " ;
    exit( 1 );

  }

  int nbBundles = bundlesNames.size() ;

  for ( int bundle = 0 ; bundle < nbBundles ; bundle++ )
  {

    file << bundlesNames[ bundle ] << " : " << bundle << std::endl ;

  }

  file.close() ;

}



////////////////////////////////////////////////////////////////////////////////
void RecognizedBundles::projectAtlas( AtlasBundles& atlasData,
                                      BundlesData& subjectBundlesData,
                                      float thresholdAdjacency,
                                      std::string outputDirectory,
                                      std::string labelsName )
{

  verbose = this->verbose ;

  int nbBundlesAtlas = atlasData.bundlesMinf.size() ;

  int nbFibersTractogram = subjectBundlesData.curves_count ;
  int nbPoints = subjectBundlesData.pointsPerTrack[ 0 ] ;


  ///////////////////// Showing parameters for projection //////////////////////
  if ( verbose == 1 )
  {

    std::cout << "Parameters for projection : " << std::endl ;

    if ( !useDefautlP )
    {

      std::cout << "   p : bundle specific ( based on distribution's variance )"
                << std::endl ;

    }
    else
    {

      std::cout << "   p : " << p << std::endl ;

    }

    if ( !useDefaultThr )
    {

      std::cout << "   thrDistance : bundle specific ( based on distribution's "
                << "variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   thrDistance : " << thrDistance << std::endl ;

    }

    if ( !useDefaultMinLen )
    {

      std::cout << "   minimumLenghtFiber : bundle specific ( based on "
                << "distribution's variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   minimumLenghtFiber : " << minimumLenghtFiber
                                              << std::endl ;

    }

    if ( !useDefaultMaxLen )
    {

      std::cout << "   maximumLenghtFiber : bundle specific ( based on "
                << "distribution's variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   maximumLenghtFiber : " << maximumLenghtFiber
                                              << std::endl ;

    }

    if ( !useDefaultMaxAngle )
    {

      std::cout << "   maxAngle : bundle specific ( based on "
                << "distribution's variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   maxAngle : " << maxAngle << std::endl ;

    }

    if ( !useDefautlMaxDirectionAngle )
    {

      std::cout << "   maxDirectionAngle : bundle specific ( based on "
                << "distribution's variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   maxDirectionAngle : " << maxDirectionAngle << std::endl ;

    }

    if ( !useDefautlMinShapeAngle )
    {

      std::cout << "   minShapeAngle : bundle specific  ( based on "
                << "distribution's variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   minShapeAngle : " << minShapeAngle << std::endl ;

    }

    if ( !useDefautlMaxShapeAngle )
    {

      std::cout << "   maxShapeAngle : bundle specific  ( based on "
                << "distribution's variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   maxShapeAngle : " << maxShapeAngle << std::endl ;

    }

    if ( !useDefaultThrDistanceBetweenMedialPoints )
    {

      std::cout << "   thrDistanceBetweenMedialPoints : bundle specific ( based"
                << " on distribution's variance ) " << std::endl ;

    }
    else
    {

      std::cout << "   thrDistanceBetweenMedialPoints : " <<
                                   thrDistanceBetweenMedialPoints << std::endl ;

    }

    if ( useMedialPointAverageFiber )
    {

      std::cout << "   medialPointAtlasBundle : medial point of average Fiber "
                << std::endl ;

    }
    else
    {

      std::cout << "   medialPointAtlasBundle : average of medial points of "
                << "the fibers in the bundle" << std::endl ;

    }

    if ( useSimpleProjection )
    {

      std::cout << "Using simple projection " << std::endl ;

    }

    if ( useMDFDistance )
    {

      std::cout << "Using MDF distance " << std::endl ;

    }
    else
    {

      std::cout << "Using MDA distance " << std::endl ;

    }

    if ( useAvgThr )
    {

      std::cout << "Using average for disimilarity" << std::endl ;

    }
    else
    {

      std::cout << "Using max for disimilarity" << std::endl ;

    }


    std::cout << "Tolerance on parameters : " << std::endl ;
    std::cout << "\ttoleranceP : " << toleranceP << std::endl ;
    std::cout << "\ttoleranceThr : " << toleranceThr << std::endl ;
    std::cout << "\ttoleranceMaxAngle : " << toleranceMaxAngle << std::endl ;
    std::cout << "\ttoleranceMaxDirectionAngle : "
	                            << toleranceMaxDirectionAngle << std::endl ;
    std::cout << "\ttoleranceMinShapeAngle : " << toleranceMinShapeAngle
	                                                          << std::endl ;
    std::cout << "\ttoleranceMaxShapeAngle : " << toleranceMaxShapeAngle
	                                                          << std::endl ;
    std::cout << "\ttoleranceLenght : " << toleranceLenght << std::endl ;
    std::cout << "\ttoleranceDistanceBetweenMedialPointst : "
                          << toleranceDistanceBetweenMedialPoints << std::endl ;

  }

  //////////////////////////////////////////////////////////////////////////////

  if ( verbose )
  {

    std::cout << "Projecting atlas ... " << std::endl ;

  }


  std::vector<std::vector<int16_t>> labels( nbFibersTractogram ) ;


  std::cout << "Number of points : " << nbPoints << "\nNumber of fibers in "
            << "input tractogram : " << nbFibersTractogram << std::endl ;

  ///////////////////////////////// PROJECTION /////////////////////////////////


  const auto before_cpu = std::chrono::system_clock::now() ;
  // Cannot put pragma here because changing atributes of this class for each
  // bundle
  for ( int atlasBundleIndex = 0 ; atlasBundleIndex < nbBundlesAtlas ;
                                                            atlasBundleIndex++ )
  {

    if ( verbose )
    {

      printf( "\rProjecting to atlas bundles : [ %d  /  %d ]",
                                   atlasBundleIndex + 1, nbBundlesAtlas ) ;
      std::cout << "" << std::flush ;

    }

    BundlesData& atlasBundleData = atlasData.bundlesData[ atlasBundleIndex ] ;
    BundlesMinf& atlasBundleInfo = atlasData.bundlesMinf[ atlasBundleIndex ] ;
    int nbFibersAtlasBundle = atlasBundleData.curves_count ;

    // Computing the medials point of atlas fibers and medial point of bundle
    std::vector<float> medialPointAtlasBundle( 3, 0 ) ;
    std::vector<float> medialPointAtlasBundleFibers(
                                                  3 * nbFibersAtlasBundle, 0 ) ;

    atlasData.findCenterBundle( atlasBundleData,
                                nbPoints,
                                medialPointAtlasBundle,
                                medialPointAtlasBundleFibers ) ;

    std::vector<float> normalVectorsAtlasBundle( 3 * nbFibersAtlasBundle, 0 ) ;
    atlasData.computeNormalVectorFibersAtlasBundle(
                                                   atlasBundleData,
                                                   nbPoints,
                                                   medialPointAtlasBundleFibers,
                                                   normalVectorsAtlasBundle ) ;


    std::vector<float> directionVectorsAtlasBundle(
                                                  3 * nbFibersAtlasBundle, 0 ) ;
    atlasData.computeDirectionVectorFibersAtlasBundle(
                                                 atlasBundleData,
                                                 normalVectorsAtlasBundle,
                                                 nbPoints,
                                                 medialPointAtlasBundleFibers,
                                                 directionVectorsAtlasBundle ) ;


    std::vector<float> lengthsAtlasBundleFibers( nbFibersAtlasBundle, 0 ) ;
    atlasData.computeLengthsAtlasBundleFibers( atlasBundleData,
                                               nbPoints,
                                               lengthsAtlasBundleFibers ) ;


    if( useMedialPointAverageFiber )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        medialPointAtlasBundle[ i ] = atlasBundleInfo.centerBundle[ i ] ;

      }

    }

    if ( !this->useDefautlP )
    {

      this->p = atlasBundleInfo.maxRadius * ( 1.0 + toleranceP ) ;
      // this->p = atlasData.bundlesMinf[ atlasBundleIndex ].averageRadius *
      //                                                     ( 1.0 + toleranceP ) ;

    }

    if ( !this->useDefaultThr )
    {

      if ( useAvgThr )
      {

        this->thrDistance = atlasBundleInfo.averageDisimilarity *
                                                        ( 1.0 + toleranceThr ) ;

      }
      else
      {
        this->thrDistance = atlasBundleInfo.maxDisimilarity *
                                                        ( 1.0 + toleranceThr ) ;
      }

    }

    if ( !this->useDefaultMaxAngle )
    {

      this->maxAngle = atlasBundleInfo.maxAngle * ( 1 + toleranceMaxAngle ) ;

    }

    if ( !this->useDefautlMaxDirectionAngle )
    {

      this->maxDirectionAngle = atlasBundleInfo.maxDirectionAngle *
                                          ( 1.0 + toleranceMaxDirectionAngle ) ;

    }

    if ( !this->useDefautlMinShapeAngle )
    {

      this->minShapeAngle = atlasBundleInfo.minShapeAngle *
                                                ( 1 - toleranceMinShapeAngle ) ;

    }

    if ( !this->useDefautlMaxShapeAngle )
    {

      this->maxShapeAngle = atlasBundleInfo.maxShapeAngle *
                                              ( 1.0 + toleranceMaxShapeAngle ) ;

    }


    if ( !this->useDefaultMinLen )
    {

      this->minimumLenghtFiber = atlasBundleInfo.minLength *
                                                    ( 1.0 - toleranceLenght ) ;

      if ( this->minimumLenghtFiber < 5 )
      {

        this->minimumLenghtFiber = 5 ;


      }

    }

    if ( !this->useDefaultMaxLen )
    {

      this->maximumLenghtFiber = atlasBundleInfo.maxLength *
                                                     ( 1.0 + toleranceLenght ) ;

      if ( this->maximumLenghtFiber > 200 )
      {

        this->maximumLenghtFiber = 200 ;


      }

    }

    if ( !useDefaultThrDistanceBetweenMedialPoints )
    {

      this->thrDistanceBetweenMedialPoints =
                           atlasBundleInfo.averageDistanceBetweenMedialPoints *
                           ( 1.0 + toleranceDistanceBetweenMedialPoints ) ;


    }



    if ( verbose > 1 )
    {

      std::cout << "\n\nBundle : "
                << atlasBundleInfo.bundles
                << std::endl ;

      printf( "medialPointAtlasBundle : [ %f, %f, %f ] \n",
              medialPointAtlasBundle[ 0 ],
              medialPointAtlasBundle[ 1 ],
              medialPointAtlasBundle[ 2 ] ) ;

      printf( "thrDistance : %f \n", this->thrDistance ) ;
      printf( "p : %f \n", this->p ) ;
      printf( "minimumLenghtFiber : %f \n", this->minimumLenghtFiber ) ;
      printf( "maximumLenghtFiber : %f \n", this->maximumLenghtFiber ) ;
      printf( "maxAngle : %f \n", this->maxAngle ) ;
      printf( "maxDirectionAngle : %f \n", this->maxDirectionAngle ) ;
      printf( "minShapeAngle : %f \n", this->minShapeAngle ) ;
      printf( "maxShapeAngle : %f \n", this->maxShapeAngle ) ;
      printf( "thrDistanceBetweenMedialPoints : %f \n",
                                        this->thrDistanceBetweenMedialPoints ) ;

    }

    if ( useMDFDistance )
    {

      if ( useSimpleProjection )
      {

        this->MDFLabelingSimple( atlasBundleData,
                                 subjectBundlesData,
                                 atlasBundleIndex,
                                 medialPointAtlasBundle,
                                 medialPointAtlasBundleFibers,
                                 lengthsAtlasBundleFibers,
                                 nbPoints,
                                 labels ) ;

      }
      else
      {

        this->MDFLabeling( atlasBundleData,
                           subjectBundlesData,
                           atlasBundleIndex,
                           medialPointAtlasBundle,
                           medialPointAtlasBundleFibers,
                           normalVectorsAtlasBundle,
                           directionVectorsAtlasBundle,
                           lengthsAtlasBundleFibers,
                           nbPoints,
                           labels ) ;

      }


    }
    else
    {

      if ( useSimpleProjection )
      {

        this->MDADLabelingSimple( atlasBundleData,
                                  subjectBundlesData,
                                  atlasBundleIndex,
                                  medialPointAtlasBundle,
                                  medialPointAtlasBundleFibers,
                                  lengthsAtlasBundleFibers,
                                  nbPoints,
                                  labels ) ;

      }
      else
      {

        this->MDADLabeling( atlasBundleData,
                            subjectBundlesData,
                            atlasBundleIndex,
                            medialPointAtlasBundle,
                            medialPointAtlasBundleFibers,
                            normalVectorsAtlasBundle,
                            directionVectorsAtlasBundle,
                            lengthsAtlasBundleFibers,
                            nbPoints,
                            labels ) ;

      }

    }

  }



  if ( verbose )
  {

    std::cout << "\nLabeling finished" << std::endl ;

  }

  int nbFibersFound_CPU = 0 ;
  // #pragma omp parallel for
  for ( int k = 0 ; k < nbFibersTractogram ; k++ )
  {

    if ( !labels[ k ].empty() )
    {

      nbFibersFound_CPU += 1 ;

    }
    else
    {

      std::vector<int16_t> _tmpVectorNoLabel = { -1 } ;
      labels[ k ] = _tmpVectorNoLabel ;

    }

  }

  std::cout << "\nNumber of fibers labeled CPU : " << nbFibersFound_CPU
                                                                  << std::endl ;


  const std::chrono::duration< double > duration_cpu =
                              std::chrono::system_clock::now() - before_cpu ;

  if ( verbose )
  {

    std::cout << "Duration of projection with CPU : "
              << duration_cpu.count() << " s" << std::endl ;

  }



  if ( verbose )
  {

    std::cout << "\n" ;

  }


 ///////////////////////////////////////////////////////////////////////////////
 //////////////////////////////// Saving results ///////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////
 char lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
 if ( lastChar != '/' )
 {

   outputDirectory = outputDirectory + "/" ;

 }

 std::vector<int> labeledFibersBundleCount( nbBundlesAtlas, 0 ) ;

 for ( int i = 0 ; i < nbFibersTractogram ; i++ )
 {

   // Not necessary to check if vector is empty since we filled the empty
   // vectors with { -1 }
   for ( int k = 0 ; k < labels[ i ].size() ; k++ )
   {

     int labelIndex = labels[ i ][ k ] ;

     if ( labelIndex > -1 )
     {

       labeledFibersBundleCount[ labelIndex ] += 1 ;

     }

   }


 }

 for ( int i = 0 ; i < nbFibersTractogram ; i++ )
 {

   // Not necessary to check if vector is empty since we filled the empty
   // vectors with { -1 }
   int countRemovedElements = 0 ;
   int nbLabelsForFiber = labels[ i ].size() ;
   for ( int k = 0 ; k < nbLabelsForFiber ; k++ )
   {

     int labelIndex = labels[ i ][ k - countRemovedElements ] ;

     if ( labelIndex == -1 )
     {

       continue ;

     }

     if ( labeledFibersBundleCount[ labelIndex ] < minimumNumberFibers )
     {

       labels[ i ].erase( labels[ i ].begin() + ( k - countRemovedElements ) ) ;
       countRemovedElements++ ;

     }

   }

   if ( labels[ i ].empty() )
   {

     std::vector<int16_t> _tmpVectorNoLabel = { -1 } ;
     labels[ i ] = _tmpVectorNoLabel ;

   }

 }


 std::string labelsBinaryFilename = outputDirectory + labelsName + ".txt" ;
 std::cout << "Saving labels in : " << labelsBinaryFilename << std::endl ;
 saveLabels( labelsBinaryFilename.c_str(), labels ) ;

 // std::string labelsDictFilename = outputDirectory + "labels.dict" ;
 std::string labelsDictFilename = outputDirectory + labelsName + ".dict" ;
 std::cout << "Saving labels dict in : " << labelsDictFilename << std::endl ;
 saveLabelsDict( labelsDictFilename.c_str(),
                 atlasData.bundlesNames ) ;



 if ( !saveExtractedBundles )
 {

   return ;

 }


  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }

  std::vector<BundlesData>& labeledFibers = this->bundlesData ;
  labeledFibers.resize( nbBundlesAtlas ) ;
  int labelIndex ;
  int numberFibersUnlabeled = 0 ;


  for ( int i = 0 ; i < nbFibersTractogram ; i++ )
  {

    // Not necessary to check if vector is empty since we filled the empty
    // vectors with { -1 }
    for ( int k = 0 ; k < labels[ i ].size() ; k++ )
    {

      labelIndex = labels[ i ][ k ] ;

      if ( labelIndex > -1 )
      {

        labeledFibers[ labelIndex ].curves_count += 1 ;

      }
      else
      {

        numberFibersUnlabeled++ ;

      }


    }

  }

  if ( verbose > 1 )
  {

    for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
    {

      printf( "Labeled bundle : %s   |   Curve count : %ld \n",
               atlasData.bundlesMinf[ bundle ].bundles.c_str(),
               labeledFibers[ bundle ].curves_count ) ;

    }

  }

  int nbBundlesRecognized = 0 ;
  for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
  {

    if ( labeledFibers[ bundle ].curves_count > minimumNumberFibers )
    {

      nbBundlesRecognized += 1 ;

    }

  }

  if ( verbose )
  {

    printf( "Number of bundles recognized : %d \n", nbBundlesRecognized ) ;


  }

  if ( nbBundlesRecognized )
  {

    std::vector<int> labeledFibers_nans( nbBundlesAtlas, 0 ) ;

    std::vector<int64_t> offsetLabeledFibers( nbBundlesAtlas, 0 ) ;
    std::vector<int64_t> offsetpointsPerTrackLabeledFiber( nbBundlesAtlas, 0 ) ;

    for ( int i = 0 ; i < nbFibersTractogram ; i++ )
    {

      // Not necessary to check if vector is empty since we filled the empty
      // vectors with { -1 }
      for ( int k = 0 ; k < labels[ i ].size() ; k++ )
      {

        labelIndex = labels[ i ][ k ] ;

        if ( labelIndex > -1 )
        {

          int64_t offsetTractogram = 3 * nbPoints * i  ;

          for ( int point = 0 ; point < nbPoints ; point++ )
          {

            for ( int coord = 0 ; coord < 3 ; coord++ )
            {

              if ( isnan( subjectBundlesData.matrixTracks[
                                      offsetTractogram + 3 * point + coord ] ) )
              {

                labeledFibers_nans[ labelIndex ]++ ;

              }

            }

          }

          if ( offsetLabeledFibers[ labelIndex ] == 0 )
          {

            labeledFibers[ labelIndex ].matrixTracks.resize(
                  labeledFibers[ labelIndex ].curves_count * 3 * nbPoints, 0 ) ;
            labeledFibers[ labelIndex ].pointsPerTrack.resize(
                                 labeledFibers[ labelIndex ].curves_count, 0 ) ;

          }


          std::copy(
                    subjectBundlesData.matrixTracks.begin() + offsetTractogram,
                    subjectBundlesData.matrixTracks.begin() + offsetTractogram +
                                                                   3 * nbPoints,
                     labeledFibers[ labelIndex ].matrixTracks.begin() +
                                           offsetLabeledFibers[ labelIndex ] ) ;
          labeledFibers[ labelIndex ].pointsPerTrack[
                              offsetpointsPerTrackLabeledFiber[ labelIndex ] ] =
                                        subjectBundlesData.pointsPerTrack[ i ] ;


          offsetLabeledFibers[ labelIndex ] += nbPoints * 3 ;

          offsetpointsPerTrackLabeledFiber[ labelIndex ] += 1 ;

        }

      }

    }

    for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
    {

      if ( verbose && labeledFibers_nans[ bundle ] )
      {

        std::cout << "WARNING : " << labeledFibers_nans[ bundle ] << " NaN "
                  << "in bundle " << atlasData.bundlesNames[ bundle ] << "\n" ;

      }

    }

  }
  else
  {

    std::cout << "No bundles of the atlas found in the input tractogram \n" ;
    return ;

  }


  if ( verbose )
  {

    std::cout << "Saving extracted bundles in : " << outputDirectory
                                                  << std::endl ;

  }

  int countProcessedRecognizedBundles = 0 ;
  // #pragma omp parallel for reduction( + : countProcessedRecognizedBundles )
  for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
  {

    if ( verbose && ( countProcessedRecognizedBundles % 1 == 0 ||
          ( countProcessedRecognizedBundles + 1 ) == nbBundlesRecognized ) &&
          labeledFibers[ bundle ].curves_count > minimumNumberFibers )
    {

      printf("\rProcessing bundle [ %10d / %10d ]",
                                        countProcessedRecognizedBundles + 1 ,
                                                        nbBundlesRecognized ) ;
      std::cout << "" << std::flush ;

    }

    if ( labeledFibers[ bundle ].curves_count > minimumNumberFibers )
    {


      // std::string labelName = getFilenameNoExtension(
      //                                   atlasBundlesFilenames[ bundle ] ) ;
      std::string labelName = atlasData.bundlesNames[ bundle ] ;

      BundlesMinf outBundlesInfo = atlasData.bundlesMinf[ bundle ] ;

      std::string format = outBundlesInfo.getFormat() ;

      outBundlesInfo.curves_count = labeledFibers[ bundle ].curves_count ;

      if ( format == ".tck" )
      {

        outBundlesInfo.computeHeaderTck(
                                        labeledFibers[ bundle ].curves_count ) ;
        outBundlesInfo.computeTckOffsetBinary() ;

      }
      int sizeOutBundlesFilename =
                               outputDirectory.size() + labelName.size() + 8 ;
      char outBundlesFilename[ sizeOutBundlesFilename ] ;
      strcpy( outBundlesFilename, outputDirectory.c_str() ) ;
      strcat( outBundlesFilename, labelName.c_str() ) ;
      strcat( outBundlesFilename, format.c_str() ) ;

      // In the initialization fo labeledFibers we did not give a format for
      // the labeledFibers[ bundle ], hence using format stored in .minf
      if ( format == ".bundles" )
      {

        labeledFibers[ bundle ].isBundles = true ;
        labeledFibers[ bundle ].isTrk = false ;
        labeledFibers[ bundle ].isTck = false ;

      }
      else if ( format == ".trk" )
      {

        labeledFibers[ bundle ].isBundles = false ;
        labeledFibers[ bundle ].isTrk = true ;
        labeledFibers[ bundle ].isTck = false ;

      }
      else if ( ".tck" )
      {

        labeledFibers[ bundle ].isBundles = false ;
        labeledFibers[ bundle ].isTrk = false ;
        labeledFibers[ bundle ].isTck = true ;

      }
      else
      {

        std::string outMessage = "BundlesMinf::getFormat : attributes " \
                                 "isBundles, isTrk and isTck are all set to " \
                                 "false, not supported format found \n" ;

        throw( std::invalid_argument( outMessage ) ) ;

      }

      labeledFibers[ bundle ].write( outBundlesFilename, outBundlesInfo ) ;

      if ( compareRecognizedWithAtlas )
      {

        BundlesData& tractogramFibers_1 = atlasData.bundlesData[ bundle ] ;
        BundlesData& tractogramFibers_2 = labeledFibers[ bundle ] ;
        int nbFibersTract_1 = atlasData.bundlesData[ bundle ].curves_count ;
        int nbFibersTract_2 = labeledFibers[ bundle ].curves_count ;
        float disimilarity =  atlasData.distanceBetweenBundles(
                                                             tractogramFibers_1,
                                                             tractogramFibers_2,
                                                             nbFibersTract_1,
                                                             nbFibersTract_2,
                                                             nbPoints ) ;

        // float thresholdAdjacency = 5 ; // In mm
        std::vector<int> nbAdjacentFibersRecognizedToAtlasBundles(
                                                        nbFibersTract_2, 0 ) ;
        atlasData.computeNumberAdjacentFibersRecognizedToAtlasBundles(
                                  labeledFibers[ bundle ],
                                  atlasData.bundlesNames[ bundle ],
                                  thresholdAdjacency,
                                  nbAdjacentFibersRecognizedToAtlasBundles ) ;

        std::vector<int> nbAdjacentFibersAtlasToRecognizedBundles(
                                                        nbFibersTract_1, 0 ) ;
        atlasData.computeNumberAdjacentFibersAtlasToRecognizedBundles(
                                  labeledFibers[ bundle ],
                                  atlasData.bundlesNames[ bundle ],
                                  thresholdAdjacency,
                                  nbAdjacentFibersAtlasToRecognizedBundles ) ;

        float coverageRecognizedToAtlas =
                                atlasData.coverageRecognizedToAtlasBundles(
                                  nbAdjacentFibersRecognizedToAtlasBundles ) ;

        float coverageAtlasToRecognized =
                                atlasData.coverageRecognizedToAtlasBundles(
                                  nbAdjacentFibersAtlasToRecognizedBundles ) ;

        float overlap = atlasData.overlapRecognizedToAtlasBundles(
                                  nbAdjacentFibersRecognizedToAtlasBundles ) ;

        float adjacency = atlasData.bundlesAdjacency(
                                                 coverageRecognizedToAtlas,
                                                 coverageAtlasToRecognized ) ;

        outBundlesInfo.write( outBundlesFilename,
                              disimilarity,
                              coverageRecognizedToAtlas,
                              overlap,
                              adjacency,
                              true ) ;

      }


      countProcessedRecognizedBundles += 1 ;

    }

  }


  ////////////////////////// Saving unlabeled fibers ///////////////////////////
  if ( saveUnlabeled )
  {

    if ( numberFibersUnlabeled == 0 )
    {

      std::cout << "WARNING : No unlabeled fibers to save" << std::endl ;
      return ;

    }

    if ( verbose )
    {

      std::cout << "\nSaving unlabeled fibers in : " << outputDirectory
                                                                  << std::endl ;

    }

    BundlesData unlabeledBundlesData ;
    unlabeledBundlesData.curves_count = numberFibersUnlabeled ;
    unlabeledBundlesData.matrixTracks.resize(
                                   numberFibersUnlabeled * 3 * nbPoints, 0 ) ;
    unlabeledBundlesData.pointsPerTrack.resize( numberFibersUnlabeled, 0 ) ;

    int64_t offsetUnlabeledFibers = 0 ;
    int offsetpointsPerTrackUnlabeledFiber = 0 ;


    int nbNans = 0 ;
    for ( int fiber = 0 ; fiber < nbFibersTractogram ; fiber++ )
    {

      // Not necessary to check if vector is empty since we filled the empty
      // vectors with { -1 }
      for ( int k = 0 ; k < labels[ fiber ].size() ; k++ )
      {

        labelIndex = labels[ fiber ][ k ] ;
        if ( labelIndex == -1 )
        {

          int64_t offsetTractogram = 3 * nbPoints * fiber ;
          for ( int point = 0 ; point < nbPoints ; point++ )
          {

            for ( int coord = 0 ; coord < 3 ; coord++ )
            {

              if ( isnan( subjectBundlesData.matrixTracks[
                                      offsetTractogram + 3 * point + coord ] ) )
              {

                nbNans++ ;

              }

            }

          }


          std::copy(
                    subjectBundlesData.matrixTracks.begin() + offsetTractogram,
                    subjectBundlesData.matrixTracks.begin() + offsetTractogram +
                                                                   3 * nbPoints,
                     unlabeledBundlesData.matrixTracks.begin() +
                                                       offsetUnlabeledFibers ) ;

          unlabeledBundlesData.pointsPerTrack[
                                    offsetpointsPerTrackUnlabeledFiber ] =
                                    subjectBundlesData.pointsPerTrack[ fiber ] ;


          offsetUnlabeledFibers += nbPoints * 3 ;

          offsetpointsPerTrackUnlabeledFiber += 1 ;

        }


      }

    }
///////////////////
    if ( verbose && nbNans )
    {

      std::cout << "WARNING : Number of NaN detected in unlabeled fibers : "
                                                        << nbNans << std::endl ;

    }
    ///////////////



    /// Saving
    std::string labelName = "unlabeledFibers" ;

    BundlesMinf unlabeledInfo = atlasData.bundlesMinf[ 0 ] ;

    unlabeledInfo.curves_count = numberFibersUnlabeled ;

    std::string format = unlabeledInfo.getFormat() ;

    if ( format == ".tck" )
    {

      unlabeledInfo.computeHeaderTck( numberFibersUnlabeled ) ;
      unlabeledInfo.computeTckOffsetBinary() ;

    }

    int sizeOutBundlesFilename = outputDirectory.size() + labelName.size() + 8 ;
    char outBundlesFilename[ sizeOutBundlesFilename ] ;
    strcpy( outBundlesFilename, outputDirectory.c_str() ) ;
    strcat( outBundlesFilename, labelName.c_str() ) ;
    strcat( outBundlesFilename, format.c_str() ) ;

    unlabeledBundlesData.write( outBundlesFilename, unlabeledInfo ) ;

  }

  if ( verbose )
  {

    std::cout << "\nDone" << std::endl ;

  }

}
