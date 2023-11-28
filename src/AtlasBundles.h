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

#include "bundlesData.h"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class AtlasBundles
{

  public :
  //////////////////////////////// Public Fields ///////////////////////////////
  std::vector<BundlesMinf> bundlesMinf ;
  std::vector<BundlesData> bundlesData ;
  std::vector<std::string> bundlesNames ;

  bool isBundles = false ;
  bool isTrk = false ;
  bool isTck = false ;

  //////////////////////////////// Constructors ////////////////////////////////
  AtlasBundles() ;

  AtlasBundles( const char* atlasDirectory,
                bool isBundlesFormat,
                bool isTrkFormat,
                bool isTckFormat,
                int verbose ) ;

  AtlasBundles( std::vector<BundlesMinf> bundlesMinf,
                std::vector<BundlesData> bundlesData,
                std::vector<std::string> bundlesNames,
                bool isBundlesFormat,
                bool isTrkFormat,
                bool isTckFormat ) ;


  //////////////////////////////// Destructors /////////////////////////////////
  virtual ~AtlasBundles() ;

  ////////////////////////////////// Mehtods ///////////////////////////////////
  void read( const char* atlasDirectory,
             bool isBundlesFormat,
             bool isTrkFormat,
             bool isTckFormat,
             int verbose ) ;

  void write( const char* outDirectory,
              bool isBundlesFormat,
              bool isTrkFormat,
              bool isTckFormat,
              int verbose ) const ;

  std::string getFormat() const ;



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
                        int verbose ) const ;

  void buildFullAtlas( BundlesMinf& bundles,
                       BundlesData& bundlesData ) const ;

  void computeLengthsAtlasBundleFibers(
                          const BundlesData& atlasBundleData,
                          int nbPoints,
                          std::vector<float>& lengthsAtlasBundleFibers ) const ;


  void findCenterBundle(
                      const BundlesData& atlasBundleData,
                      int nbPoints,
                      std::vector<float>& medialPointAtlasBundle,
                      std::vector<float>& medialPointAtlasBundleFibers ) const ;
  void findCenterBundle(
                      int bundleIndex,
                      std::vector<float>& medialPointAtlasBundle,
                      std::vector<float>& medialPointAtlasBundleFibers ) const ;


  void computeNormalVectorFibersAtlasBundle(
                         const BundlesData& atlasBundleData,
                         int nbPoints,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         std::vector<float>& normalVectors ) const ;


  void computeDirectionVectorFibersAtlasBundle(
                         const BundlesData& atlasBundleData,
                         const std::vector<float>& normalVectorsAtlasBundle,
                         int nbPoints,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         std::vector<float>& directionVectors ) const ;

  double compareDisimilarityBundles(
                         const BundlesData& bundle1,
                         const BundlesData& bundle2,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         int nbPoints ) const ;
  double compareDisimilarityBundles(
                         const BundlesData& bundle1,
                         const BundlesData& bundle2,
                         int nbPoints ) const ;
  double compareDisimilarityBundles(
                         const BundlesData& bundle1,
                         const BundlesData& bundle2,
                         int nbThreads,
                         int nbPoints ) const ;


  double distanceBetweenBundles( const BundlesData& bundle1,
                                 const BundlesData& bundle2,
                                 int nbPoints ) const ;
  double distanceBetweenBundles( const BundlesData& bundle1,
                                 const BundlesData& bundle2,
                                 int nbThreads,
                                 int nbPoints ) const ;
  double distanceBetweenBundles( int bundleIndex,
                                 const BundlesData& bundle,
                                 int nbFibersBundle,
                                 int nbPoints ) const ;

  void computeNumberAdjacentFibersBundle1ToBundle2(
                              const std::vector<float>& bundle1,
                              const std::vector<float>& bundle2,
                              int nbFibersBundle1,
                              int nbFibersBundle2,
                              int nbPoints, // Same for 2 bundles
                              float threshold,
                              std::vector<int>& nbAdjacentFibersBundle ) const ;
void computeNumberAdjacentFibersBundle1ToBundle2(
                            const std::vector<float>& bundle1,
                            const std::vector<float>& bundle2,
                            int nbFibersBundle1,
                            int nbFibersBundle2,
                            int nbPoints, // Same for 2 bundles
                            float threshold,
                            int nbThreads,
                            std::vector<int>& nbAdjacentFibersBundle ) const ;

  void computeNumberAdjacentFibersRecognizedToAtlasBundles(
                              const BundlesData& bundle,
                              std::string bundleName,
                              float threshold,
                              std::vector<int>& nbAdjacentFibersBundle ) const ;

  void computeNumberAdjacentFibersAtlasToRecognizedBundles(
                              const BundlesData& bundle,
                              std::string bundleName,
                              float threshold,
                              std::vector<int>& nbAdjacentFibersBundle ) const ;

  float coverageBundle1ToBundle2( const std::vector<float>& bundle1,
                                  const std::vector<float>& bundle2,
                                  int nbFibersBundle1,
                                  int nbFibersBundle2,
                                  int nbPoints, // Same for 2 bundles
                                  float threshold,
                                  int verbose ) const ;
  float coverageBundle1ToBundle2(
              const std::vector<int>& nbAdjacentFibersBundle1ToBundle2 ) const ;

  float coverageRecognizedToAtlasBundles( const BundlesData& bundle,
                                          std::string bundleName,
                                          float threshold,
                                          int verbose ) const ;
  float coverageRecognizedToAtlasBundles(
      const std::vector<int>& nbAdjacentFibersRecognizedToAtlasBundles ) const ;


  float coverageAtlasToRecognizedBundles( const BundlesData& bundle,
                                          std::string bundleName,
                                          float threshold,
                                          int verbose ) const ;
  float coverageAtlasToRecognizedBundles(
      const std::vector<int>& nbAdjacentFibersAtlasToRecognizedBundles ) const ;

  float overlapBundle1ToBundle2( const std::vector<float>& bundle1,
                                 const std::vector<float>& bundle2,
                                 int nbFibersBundle1,
                                 int nbFibersBundle2,
                                 int nbPoints, // Same for 2 bundles
                                 float threshold,
                                 int verbose ) const ;
  float overlapBundle1ToBundle2(
              const std::vector<int>& nbAdjacentFibersBundle1ToBundle2 ) const ;

  float overlapRecognizedToAtlasBundles( const BundlesData& bundle,
                                         std::string bundleName,
                                         float threshold,
                                         int verbose ) const ;      
  float overlapRecognizedToAtlasBundles(
               const std::vector<int>& nbAdjacentFibersRecognizedToAtlasBundles,
               int verbose = 0 ) const ;

  float overlapAtlasToRecognizedBundles( const BundlesData& bundle,
                                         std::string bundleName,
                                         float threshold,
                                         int verbose ) const ;
  float overlapAtlasToRecognizedBundles(
      const std::vector<int>& nbAdjacentFibersAtlasToRecognizedBundles ) const ;

  float bundlesAdjacency(  const std::vector<float>& bundle1,
                           const std::vector<float>& bundle2,
                           int nbFibersBundle1,
                           int nbFibersBundle2,
                           int nbPoints, // Same for 2 bundles
                           float threshold,
                           int verbose ) const ;

  float bundlesAdjacency(  const BundlesData& bundle,
                           std::string bundleName,
                           float threshold,
                           int verbose ) const ;
  float bundlesAdjacency( float coverageRecognizedToAtlasBundlesMeasure,
                          float coverageAtlasToRecognizedBundlesMeasure ) const;

  int findBundleIndexByName( std::string bundleName ) const ;

} ;
