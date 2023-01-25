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


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool is_file( const std::string& path ) ;
//
bool is_dir( const std::string& path ) ;
//
bool mkdir( const std::string& path ) ;
//
bool rmfile( const std::string& path ) ;
//
bool rmdir( const std::string& path ) ;
//
bool copy( const std::string& source,
           const std::string& destination ) ;
//
bool copytree( const std::string& source,
                      const std::string& destination,
                      bool forceOverride ) ;
//
bool rename( const std::string& source,
             const std::string& destination ) ;
//
std::string dirname( const std::string& path ) ;
//
std::string basename( const std::string& path ) ;
//
std::string basenameNoExtension( const std::string& path ) ;
//
std::string replaceExtension( const std::string& path,
                              const std::string& newExtension ) ;
//
std::string getFilenameNoExtension( const std::string& path ) ;
//
bool endswith( const std::string& input,
               const std::string& substring ) ;
//
int countFilesDirectory( const std::string& path ) ;
//
std::vector<std::string> getFilesInDirectoryWithExtension(
                                                const std::string& path,
                                                const std::string& extension ) ;
//
void listDir( const std::string& path,
              std::vector<std::string>& dirList ) ;
//
void checkAtlasDirectory( const std::string& path,
                          const std::string& format ) ;
//
std::string run_sh_process_timeout( const std::string& command, int timeout ) ;
//
int run_sh_process( const std::string& command, int verbose = 0 ) ;
//
void readingPredictedLabels( const char* predictedLabelsFilename,
                             std::vector<std::vector<int16_t>>& predictedLabels,
                             int nbFibers ) ;
//
void saveLabels( const char* labelsBinaryFilename,
                 const std::vector<std::vector<int16_t>>& labels ) ;
//
void readLabelsDict( const char* labelsDictFilename,
                     std::vector<std::string>& bundlesNames ) ;
//
void saveLabelsDict( const char* labelsDictFilename,
                     const std::vector<std::string>& bundlesNames ) ;
//
void readLabelsWithDict( const char* labelsDictFilename,
                         const char* labelsBinaryFilename,
                         std::vector<std::vector<std::string>>& labelsByName,
                         int nbFibers ) ;
//
int getLabelFromName( const std::vector<std::string>& bundlesDict,
                      const std::string& bundleName ) ;
//
void readIndexInTractogram( const char* predictedLabelsFilename,
                            std::vector<int64_t>& predictedLabels,
                            int nbFibers ) ;
//
std::string getAtlasBunldesPaths(
                                 const std::string& outputDirectory,
                                 const std::string& atlasDirectory,
                                 const std::string& format,
                                 std::vector<std::string>& atlasBundlesPaths ) ;
//
void getAtlasNeighborhoodCentroids(
                              const std::string& inputDirectory,
                              const std::vector<std::string>& atlasBundlesPaths,
                              const std::string& format,
                              std::vector<std::string>& bundlesPaths ) ;
//
void getNeighborhoodFilenames(
                             const std::string& tmpNeighborhoodDir,
                             const std::vector<std::string>& atlasBundlesPaths,
                             const std::string& format,
                             std::vector<std::string>& neighborhoodFilenames ) ;
//
void convertBundlesFormat( const std::string& inputBundles,
                           const std::string& outputTrk,
                           const std::string& referenceImage,
                           bool force = false,
                           int verbose = 0 ) ;
//
void saveComparisonMeasuresWithAtlas(
                                    const std::vector<float>& coveragesBundles,
                                    const std::vector<float>& adjacencyBundles,
                                    const std::vector<float>& overlapBundles,
                                    const std::vector<std::string>& labelsDict,
                                    const char* fileName ) ;
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
