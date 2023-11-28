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
       index_anc, index_thr, index_thrCov, index_thrAdj, index_minNbFibers, 
       index_thrSim, index_adjCB, index_thrDBMP, index_tolP, index_tolThr, 
       index_tolMaxAngle, index_tolMaxDirectionAngle, index_tolMinShapeAngle,
       index_tolMaxShapeAngle, index_tolLenght, index_tolThrCN,
       index_tolDistBetMedPts, index_pa, index_cv, index_cn, index_cc,
       index_nbPoints, index_rb, index_ods, index_cds, index_slr, index_cp,
       index_mdf, index_simple, index_sp, index_force, index_verbose,  
       index_nbThreads, index_nbThreadsCN, index_keep_tmp, index_time_out, 
       index_useMeanForMDAD, index_help ;

// Variables for projection
int verbose = 0 ;
std::string projectAtlasFile = "ProjectAtlas" ;
std::string convertBundleFormatsFile = "ConvertBundleFormat" ;
std::string computeNeighborhoodFile = "computeNeighborhood" ;
std::string applyTransformBundleFile = "applyTransformBundle" ;
std::string computeCentroidsClientFilename = "clientComputeCentroids.py" ;
std::string registerBundlesClientFile = "clientRegisterBundles.py" ;
std::string openDipyServerClientFile = "dipyServer.py" ;
std::string closeDipyServerClientFile = "clientCloseServer.py" ;
float thrDistance = 10 ; // In mm
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
float toleranceMaxAngle = 0.0 ;
float toleranceMaxDirectionAngle = 0.0 ;
float toleranceMinShapeAngle = 0.0 ;
float toleranceMaxShapeAngle = 0.0 ;
float toleranceLenght = 0.0 ;
float toleranceDistanceBetweenMedialPoints = -0.9 ;

float toleranceThrComputeNeighborhood = 0.0 ;

float adjacencyForCompareBundles = 5.0 ; // In mm

bool doClassical = true ;

bool useMDF = false ;

bool useSimple = false ;

bool keepTmpFiles = false ;

bool haveMinf = false ;

float time_out = 50 ; // In s

bool useMeanForMDAD = true ;

bool useDefaultThr = false ;

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
                  float& time_out,
                  int verbose ) ;
//
int getPortNumberDipyService( std::string& logFilePath ) ;
//
void closeDipyServer( int portDipyServer ) ;
