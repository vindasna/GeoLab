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

//#include "BundlesFormat.h"
//#include "BundlesDataFormat.h"
#include "TrkFormat.h"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class AtlasBundles
{

  public :
  //////////////////////////////// Public Fields ///////////////////////////////
  std::vector< BundlesFormat > bundles ;
  std::vector< BundlesDataFormat > bundlesData ;
  std::vector< std::string > bundlesNames ;

  //////////////////////////////// Constructors ////////////////////////////////
  AtlasBundles() ;

  AtlasBundles( const char* atlasDirectory,
                bool isBundlesFormat,
                bool isTRKFormat,
                int verbose ) ;

  AtlasBundles( std::vector< BundlesFormat > bundles,
                std::vector< BundlesDataFormat > bundlesData,
                std::vector< std::string > bundlesNames ) ;


  //////////////////////////////// Destructors /////////////////////////////////
  virtual ~AtlasBundles() ;

  ////////////////////////////////// Mehtods ///////////////////////////////////
  void atlasReading( const char* atlasDirectory,
                     bool isBundlesFormat,
                     bool isTRKFormat,
                     int verbose ) ;

  void atlasWriting( const char* outDirectory,
                     bool isBundlesFormat,
                     bool isTRKFormat,
                     int verbose ) ;

  void analysisWriting( const char* analysisFilename,
                        std::vector<float>& distancesToCenterAtlasBundle,
                        std::vector<float>& lengthsAtlasBundle,
                        std::vector<float>& disimilaritiesAtlasBundle,
                        std::vector<float>& anglesAtlasBundle,
                        std::vector<float>& directionAnglesAtlasBundle,
                        std::vector<float>& shapeAnglesAtlasBundle,
                        std::vector<float>& distancesBetweenMedialPointsBundle,
                        int nbFibersAtlasBundle,
                        int numberDisimilaritiesAtlasBundle,
                        int numberDistancesBetweenMedialPoints,
                        int verbose ) ;

  void buildFullAtlas( BundlesFormat& bundles,
                       BundlesDataFormat& bundlesData ) ;

  void computeLengthsAtlasBundleFibers(
                                BundlesDataFormat& atlasBundleData,
                                int nbPoints,
                                std::vector<float>& lengthsAtlasBundleFibers ) ;


  void findCenterBundle( BundlesDataFormat& atlasBundleData,
                         int nbPoints,
                         std::vector<float>& medialPointAtlasBundle,
                         std::vector<float>& medialPointAtlasBundleFibers ) ;
  void findCenterBundle( int bundleIndex,
                         std::vector<float>& medialPointAtlasBundle,
                         std::vector<float>& medialPointAtlasBundleFibers ) ;


  void computeNormalVectorFibersAtlasBundle(
                         BundlesDataFormat& atlasBundleData,
                         int nbPoints,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         std::vector<float>& normalVectors ) ;


  void computeDirectionVectorFibersAtlasBundle(
                         BundlesDataFormat& atlasBundleData,
                         const std::vector<float>& normalVectorsAtlasBundle,
                         int nbPoints,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         std::vector<float>& directionVectors ) ;

  double compareDisimilarityBundles(
                         BundlesDataFormat& bundle1,
                         BundlesDataFormat& bundle2,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         int nbPoints ) ;
  double compareDisimilarityBundles(
                         BundlesDataFormat& bundle1,
                         BundlesDataFormat& bundle2,
                         int nbPoints ) ;


  double distanceBetweenBundles( BundlesDataFormat& bundle1,
                                 BundlesDataFormat& bundle2,
                                 int nbFibersBundle1,
                                 int nbFibersBundle2,
                                 int nbPoints ) ;
  double distanceBetweenBundles( int bundleIndex,
                                 BundlesDataFormat& bundle,
                                 int nbFibersBundle,
                                 int nbPoints ) ;

  void computeNumberAdjacentFibersBundle1ToBundle2(
                                    const std::vector<float>& bundle1,
                                    const std::vector<float>& bundle2,
                                    int nbFibersBundle1,
                                    int nbFibersBundle2,
                                    int nbPoints, // Same for 2 bundles
                                    float threshold,
                                    std::vector<int>& nbAdjacentFibersBundle ) ;

  void computeNumberAdjacentFibersRecognizedToAtlasBundles(
                                    BundlesDataFormat& bundle,
                                    std::string bundleName,
                                    float threshold,
                                    std::vector<int>& nbAdjacentFibersBundle ) ;

  void computeNumberAdjacentFibersAtlasToRecognizedBundles(
                                    BundlesDataFormat& bundle,
                                    std::string bundleName,
                                    float threshold,
                                    std::vector<int>& nbAdjacentFibersBundle ) ;

  float coverageBundle1ToBundle2( const std::vector<float>& bundle1,
                                  const std::vector<float>& bundle2,
                                  int nbFibersBundle1,
                                  int nbFibersBundle2,
                                  int nbPoints, // Same for 2 bundles
                                  float threshold,
                                  int verbose ) ;
  float coverageBundle1ToBundle2(
                    const std::vector<int>& nbAdjacentFibersBundle1ToBundle2 ) ;

  float coverageRecognizedToAtlasBundles( BundlesDataFormat& bundle,
                                          std::string bundleName,
                                          float threshold,
                                          int verbose ) ;
  float coverageRecognizedToAtlasBundles(
            const std::vector<int>& nbAdjacentFibersRecognizedToAtlasBundles ) ;


  float coverageAtlasToRecognizedBundles( BundlesDataFormat& bundle,
                                          std::string bundleName,
                                          float threshold,
                                          int verbose ) ;
  float coverageAtlasToRecognizedBundles(
            const std::vector<int>& nbAdjacentFibersAtlasToRecognizedBundles ) ;

  float overlapBundle1ToBundle2( const std::vector<float>& bundle1,
                                 const std::vector<float>& bundle2,
                                 int nbFibersBundle1,
                                 int nbFibersBundle2,
                                 int nbPoints, // Same for 2 bundles
                                 float threshold,
                                 int verbose ) ;
  float overlapBundle1ToBundle2(
                    const std::vector<int>& nbAdjacentFibersBundle1ToBundle2 ) ;

  float overlapRecognizedToAtlasBundles( BundlesDataFormat& bundle,
                                         std::string bundleName,
                                         float threshold,
                                         int verbose ) ;
  float overlapRecognizedToAtlasBundles(
            const std::vector<int>& nbAdjacentFibersRecognizedToAtlasBundles ) ;

  float overlapAtlasToRecognizedBundles( BundlesDataFormat& bundle,
                                         std::string bundleName,
                                         float threshold,
                                         int verbose ) ;
  float overlapAtlasToRecognizedBundles(
            const std::vector<int>& nbAdjacentFibersAtlasToRecognizedBundles ) ;

  float bundlesAdjacency(  const std::vector<float>& bundle1,
                           const std::vector<float>& bundle2,
                           int nbFibersBundle1,
                           int nbFibersBundle2,
                           int nbPoints, // Same for 2 bundles
                           float threshold,
                           int verbose ) ;

  float bundlesAdjacency(  BundlesDataFormat& bundle,
                           std::string bundleName,
                           float threshold,
                           int verbose ) ;
  float bundlesAdjacency( float coverageRecognizedToAtlasBundlesMeasure,
                          float coverageAtlasToRecognizedBundlesMeasure ) ;

  int findBundleIndexByName( std::string bundleName ) ;

} ;

std::string getFilenameNoExtension( std::string path ) ;
