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

#include "./tools/BMDDistance.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int verbose = 0 ;
std::string projectAtlasFile = "ProjectAtlas" ;
std::string convertBundleFormatsFile = "ConvertBundleFormat" ;
std::string computeNeighborhoodFile = "computeNeighborhood" ;
std::string applyTransformBundleFile = "applyTransformBundle" ;
std::string computeCentroidsClientFilename ;
std::string registerBundlesClientFile ;
std::string openDipyServerClientFile ;
std::string closeDipyServerClientFile ;
std::string analyseAtlasBundleFile ;
float coverageThreshold = 0.0 ;
float adjacencyThreshold = 0.0 ;
int nbPointsPerFiber = 0 ;
bool doSLR = false ;
bool saveBundlesSeparetly = true ;
bool force = false ;
int nbThreads = 0 ;
int k_value = -1 ;
float percentageTest = -1 ;

bool isFullAtlas = false ;
bool isAtlasNeighborhoodCentroids = false ;

std::string atlasDirectory ;
std::string referenceFilename ;
std::string outputDirectory ;
std::string fullAtlasBundlesFilename ;
std::string fullAtlasBundlesDataFilename ;
std::string fullAtlasLabelsFilename ;
std::string fullAtlasLabelsDictFilename ;
std::string fullAtlasMultiLabelsFilename ;
std::string fullAtlasMultiLabelsDictFilename ;
std::string atlasNeighborhoodCentroidsDirectory ;

std::string format ;
////////////////////////////////////////////////////////////////////////////////
void saveScorePerBundle( const char* scoresPerBundlePath,
                         const std::vector<float>& precisionPerBundle,
                         const std::vector<float>& recallPerBundle,
                         const std::vector<float>& accuracyPerBundle,
                         const std::vector<float>& weightsBundles,
                         const std::vector<std::string>& bundleNames,
                         float averagePrecision,
                         float averageRecall,
                         float averageAccuracy ) ;
//
bool generateTrueWithXprobability( float probability ) ;
