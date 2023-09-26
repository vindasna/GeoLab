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

#include "AtlasBundles.h"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class RecognizedBundles : public AtlasBundles
{

  public :
  //////////////////////////////// Public Fields ///////////////////////////////
  int verbose = 0 ;
  float p = 5 ; // In mm
  float thrDistance = 10 ; // In mm
  float minimumLenghtFiber = 10 ; // In mm
  float maximumLenghtFiber = 200 ; // In mm
  float maxAngle = 45; // In degrees
  float maxDirectionAngle = 45; // In degrees
  float minShapeAngle = 0; // In degrees
  float maxShapeAngle = 180; // In degrees
  bool compareRecognizedWithAtlas = false ;
  bool useMDFDistance = false ;
  bool useMedialPointAverageFiber = false ;
  bool useSimpleProjection = false ;
  float toleranceP = 0.0 ;
  float toleranceThr = 0.0 ;
  float toleranceMaxAngle = 0.0 ;
  float toleranceMaxDirectionAngle = 0.0 ;
  float toleranceMinShapeAngle = 0.0 ;
  float toleranceMaxShapeAngle = 0.0 ;
  float toleranceLenght = 0.3 ;
  float toleranceDistanceBetweenMedialPoints = 0.3 ;
  float thrPercentageSimilarity = 0.10 ;
  float thrDistanceBetweenMedialPoints = 5 ;
  float minimumNumberFibers = 20 ;
  float thresholdAdjacency = 5 ; // In mm
  bool useDefautlP = true ;
  bool useDefaultThr = true ;
  bool useDefaultMaxAngle = true ;
  bool useDefautlMaxDirectionAngle = true ;
  bool useDefautlMinShapeAngle= true ;
  bool useDefautlMaxShapeAngle = true ;
  bool useDefaultThrDistanceBetweenMedialPoints = false ;
  bool useDefaultMinLen = true ;
  bool useDefaultMaxLen = true ;
  bool saveExtractedBundles = true ;
  bool saveUnlabeled = false ;
  bool useAvgThr = false ;

  int nbThreads = -1 ;

  //////////////////////////////// Constructors ////////////////////////////////
  RecognizedBundles() ;

  RecognizedBundles( int verbose,
                     float p,
                     float thrDistance,
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
                     int nbThreads ) ;

  //////////////////////////////// Destructors /////////////////////////////////
  virtual ~RecognizedBundles() ;


  ////////////////////////////////// Mehtods ///////////////////////////////////
  void labeling( const BundlesData& atlasBundleData,
                     const BundlesData& subjectBundlesData,
                     int atlasBundleIndex,
                     const std::vector<float>& medialPointAtlasBundle,
                     const std::vector<float>& medialPointAtlasBundleFibers,
                     const std::vector<float>& normalVectorsAtlasBundle,
                     const std::vector<float>& directionVectorsAtlasBundle,
                     const std::vector<float>& lengthsAtlasBundleFibers,
                     int nbPoints,
                     bool useMeanForMDAD,
                     std::vector<std::vector<int16_t>>& labels ) ;
  void labeling( const BundlesData& atlasBundleData,
                     const BundlesData& subjectBundlesData,
                     const std::vector<float>& medialPointAtlasBundle,
                     const std::vector<float>& medialPointAtlasBundleFibers,
                     const std::vector<float>& normalVectorsAtlasBundle,
                     const std::vector<float>& directionVectorsAtlasBundle,
                     const std::vector<float>& lengthsAtlasBundleFibers,
                     int nbPoints,
                     bool useMeanForMDAD,
                     std::vector<int>& indexInTractogramRecognized ) ;


  void labelingSimple(
                         const BundlesData& atlasBundleData,
                         const BundlesData& subjectBundlesData,
                         int atlasBundleIndex,
                         const std::vector<float>& medialPointAtlasBundle,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         const std::vector<float>& lengthsAtlasBundleFibers,
                         int nbPoints,
                         bool useMeanForMDAD,
                         std::vector<std::vector<int16_t>>& labels ) ;
  void labelingSimple(
                         const BundlesData& atlasBundleData,
                         const BundlesData& subjectBundlesData,
                         const std::vector<float>& medialPointAtlasBundle,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         const std::vector<float>& lengthsAtlasBundleFibers,
                         int nbPoints,
                         bool useMeanForMDAD,
                         std::vector<int>& indexInTractogramRecognized ) ;


  void saveLabels( const char* labelsBinaryFilename,
                   const std::vector<std::vector<int16_t>>& labels ) ;

  void saveLabelsDict( const char* labelsDictFilename,
                       const std::vector< std::string >& atlasData ) ;


  void projectAtlas( const AtlasBundles& atlasData,
                     const BundlesData& subjectBundlesData,
                     float thresholdAdjacency,
                     std::string outputDirectory,
                     std::string labelsName,
                     std::string comparisonWithAtlasFilename,
                     bool comparisonWithAtlasAppend,
                     bool saveBundlesSeparetly,
                     bool useMeanForMDAD ) ;

  void projectBundle( const BundlesData& atlasBundleData,
                      const BundlesMinf& atlasBundleInfo,
                      const BundlesData& subjectBundlesData,
                      std::vector<int>& indexInTractogramRecognized,
                      float& coverage,
                      float& adjacency,
                      float& overlap,
                      float& disimilarity,
                      bool useMeanForMDAD,
                      bool comparisonWithAtlas ) ;

  void printParams() ;




} ;
