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

#include <boost/process.hpp>

#include "BMDDistance.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Variables for parse
int index_input, index_atlas, index_reference, index_output, index_fa,
       index_anc, index_thrCov, index_thrAdj, index_minNbFibers, index_thrSim,
       index_adjCB, index_thrDBMP, index_tolP, index_tolThr, index_tolMaxAngle,
       index_tolMaxDirectionAngle, index_tolMinShapeAngle,
       index_tolMaxShapeAngle, index_tolLenght, index_tolThrCN,
       index_tolDistBetMedPts, index_pa, index_cv, index_cn, index_cc,
       index_nbPoints, index_rb, index_slr, index_cp, index_sp, index_force,
       index_verbose, index_nbThreads, index_help ;

// Variables for projection
int verbose = 0 ;
std::string projectAtlasFile = "ProjectAtlas" ;
std::string convertBundleFormatsFile = "ConvertBundleFormat" ;
std::string computeNeighborhoodFile = "computeNeighborhood" ;
std::string applyTransformBundleFile = "applyTransformBundle" ;
std::string computeCentroidsFilename ;
std::string registerBundlesFile ;
float coverageThreshold = 0.0 ;
float adjacencyThreshold = 0.0 ;
float thrDistanceBetweenMedialPoints = 50.0 ;
int nbPointsPerFiber = 0 ;
bool doSLR = false ;
bool saveBundlesSeparetly = true ;
bool force = false ;
int nbThreads = 0 ;
// int minimumNumberFibers = 20 ;
int minimumNumberFibers = 1 ;
// float thrPercentageSimilarity = 0.05 ;
float thrPercentageSimilarity = 0.05 ;
// float thrPercentageSimilarity = 0.50 ;
// float thrPercentageSimilarity = 0.50 ;

bool isFullAtlas = false ;
bool isAtlasNeighborhoodCentroids = false ;

float toleranceP = 0.0 ;
float toleranceThr = 0.0 ;
float toleranceMaxAngle = 0.0 ;
float toleranceMaxDirectionAngle = 0.0 ;
float toleranceMinShapeAngle = 0.0 ;
float toleranceMaxShapeAngle = 0.0 ;
float toleranceLenght = 0.0 ;
float toleranceDistanceBetweenMedialPoints = 0.0 ;

float toleranceThrComputeNeighborhood = 1.2 ;

float adjacencyForCompareBundles = 5 ; // In mm

bool doClassical = false ;


////////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;

//
inline bool is_file( const std::string& path ) ;
//
inline bool is_dir( const std::string& path ) ;
//
inline bool mkdir( const std::string& path ) ;
//
inline bool rmfile( const std::string& path ) ;
//
inline bool rmdir( const std::string& path ) ;
//
inline bool copy( const std::string& source,
                  const std::string& destination ) ;
//
inline bool rename( const std::string& source,
                    const std::string& destination ) ;
//
inline std::string dirname( const std::string& path ) ;
//
inline std::string replaceExtension( const std::string& path,
                                     const std::string& newExtension ) ;
//
inline bool endswith( const std::string& input,
                      const std::string& substring ) ;
//
int countFilesDirectory( const std::string& path ) ;
//
void checkAtlasDirectory( const std::string& path,
                          std::string& outputDirectory ) ;
//
std::string run_sh_process_timeout( const std::string& command, int timeout ) ;
//
int run_sh_process( const std::string& command ) ;
//
void readingPredictedLabels( const char* predictedLabelsFilename,
                             std::vector<std::vector<int16_t>>& predictedLabels,
                             const std::string& outputDirectory,
                             int nbFibers,
                             int verbose ) ;
//
void saveLabels( const char* labelsBinaryFilename,
                 const std::vector<std::vector<int16_t>>& labels,
                 const std::string& outputDirectory ) ;
//
void readLabelsDict( const char* labelsDictFilename,
                     std::vector<std::string>& bundlesNames,
                     const std::string& outputDirectory ) ;
//
void saveLabelsDict( const char* labelsDictFilename,
                     const std::vector<std::string>& bundlesNames,
                     const std::string& outputDirectory ) ;
//
void readIndexInTractogram( const char* predictedLabelsFilename,
                            std::vector<int64_t>& predictedLabels,
                            const std::string& outputDirectory,
                            int nbFibers,
                            int verbose ) ;
//
std::string getAtlasBunldesPaths(
                                 const std::string& outputDirectory,
                                 const std::string& atlasDirectory,
                                 std::vector<std::string>& atlasBundlesPaths ) ;
//
void getAtlasNeighborhoodCentroids(
                              const std::string& outputDirectory,
                              const std::string& inputDirectory,
                              const std::vector<std::string>& atlasBundlesPaths,
                              std::vector<std::string>& bundlesPaths ) ;
//
void getNeighborhoodFilenames(
                             const std::string& tmpNeighborhoodDir,
                             const std::vector<std::string>& atlasBundlesPaths,
                             std::vector<std::string>& neighborhoodFilenames,
                             const std::string& outputDirectory ) ;
//
void convertBundlesFormat( const std::string& inputBundles,
                           const std::string& outputTrk,
                           const std::string& referenceImage,
                           const std::string& outputDirectory ) ;
//
void saveComparisonMeasuresWithAtlas(
                                    const std::vector<float>& coveragesBundles,
                                    const std::vector<float>& adjacencyBundles,
                                    const std::vector<float>& overlapBundles,
                                    const std::vector<std::string>& labelsDict,
                                    const char* fileName,
                                    const std::string& outputDirectory ) ;
//
float getCoverageWithAtlas( const std::string& bundleFilename ) ;
//
float getAdjacencyWithAtlas( const std::string& bundleFilename ) ;
//
float getOverlapWithAtlas( const std::string& bundleFilename ) ;
//
float getAverageRadiusAtlasBundle( const std::string& bundleFilename ) ;
//
float getAverageDistanceBetweenMedialPoints(
                                           const std::string& bundleFilename ) ;
//
int getNbFibers( const std::string& bundleFilename ) ;
//
void applyRecoBundles( const std::string& movedTractogramNeighborhood,
                       const std::string& atlasBundleFile,
                       const std::string& atlasNeighborhoodFile,
                       const std::string& outputDirectory,
                       const std::string& referenceImage,
                       int nbPointsPerFiber,
                       int verbose ) ;
