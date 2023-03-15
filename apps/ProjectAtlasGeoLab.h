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
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include <boost/process.hpp>

#include "RecognizedBundles.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Variables for parse
int index_input, index_atlas, index_reference, index_output, index_fa, index_an,
       index_anc, index_thrCov, index_thrAdj, index_minNbFibers, index_thrSim,
       index_adjCB, index_thrDBMP, index_tolP, index_tolThr, index_tolMaxAngle,
       index_tolMaxDirectionAngle, index_tolMinShapeAngle,
       index_tolMaxShapeAngle, index_tolLenght, index_tolThrCN,
       index_tolDistBetMedPts, index_pa, index_cv, index_cn, index_cc,
       index_nbPoints, index_rb, index_ods, index_cds, index_slr, index_cp,
       index_sp, index_force, index_verbose, index_nbThreads, index_nbThreadsCN,
                                                    index_keep_tmp, index_help ;

// Variables for projection
int verbose = 0 ;
std::string projectAtlasFile = "ProjectAtlas" ;
std::string convertBundleFormatsFile = "ConvertBundleFormat" ;
std::string computeNeighborhoodFile = "computeNeighborhood" ;
std::string applyTransformBundleFile = "applyTransformBundle" ;
std::string computeCentroidsClientFilename ;
std::string registerBundlesClientFile ;
std::string openDipyServerClientFile ;
std::string closeDipyServerClientFile ;
float coverageThreshold = 0.0 ;
float adjacencyThreshold = 0.0 ;
float thrDistanceBetweenMedialPoints = 50.0 ;
int nbPointsPerFiber = 0 ;
bool doSLR = false ;
bool saveBundlesSeparetly = true ;
bool force = false ;
int nbThreads = 0 ;
int nbThreadsCN = 0 ;
int nbThreadsUsed ;
// int minimumNumberFibers = 20 ;
int minimumNumberFibers = 1 ;
float thrPercentageSimilarity = 0.00001 ;

bool isEsbaAtlas = false ;
bool isFullAtlas = false ;
bool isAtlasNeighborhoodCentroids = false ;
bool isAtlasNeighborhood = false ;

float toleranceP = 0.0 ;
float toleranceThr = 0.0 ;
float toleranceMaxAngle = 1.0 ;
float toleranceMaxDirectionAngle = 1.0 ;
float toleranceMinShapeAngle = 1.0 ;
float toleranceMaxShapeAngle = 1.0 ;
float toleranceLenght = 0.0 ;
float toleranceDistanceBetweenMedialPoints = 1.0 ;

float toleranceThrComputeNeighborhood = 6.0 ;

float adjacencyForCompareBundles = 5.0 ; // In mm

bool doClassical = true ;

bool keepTmpFiles = false ;

bool haveMinf = false ;

std::string default_ESBA_DIR ;


////////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;
//
void applyGeoLab( const std::string& movedTractogramNeighborhood,
                  const std::string& atlasBundleDirectory,
                  const std::string& atlasNeighborhoodFile,
                  const std::string& atlasNeighborhoodCentroidsFile,
                  const std::string& outputDirectory,
                  const std::string& referenceImage,
                  const std::string& format,
                  const int minimumNumberFibers,
                  std::vector<int16_t>& indexInNeighborhoodRecognized,
                  float adjacency_classic,
                  int nbFibersClassic,
                  int nbPointsPerFiber,
                  int portDipyServer,
                  bool& keepClassic,
                  float& coverageGeoLab,
                  float& adjacencyGeoLab,
                  float& overlapGeoLab,
                  float& disimilarityGeoLab,
                  int verbose ) ;
//
int getPortNumberDipyService( std::string& logFilePath ) ;
//
void closeDipyServer( int portDipyServer ) ;
