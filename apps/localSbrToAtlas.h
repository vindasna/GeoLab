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
       index_anc, index_tolMaxShapeAngle, index_tolLenght, index_tolThrCN, 
       index_cn, index_slr, index_cc, index_nbPoints, index_rb, index_ods, 
       index_cds, index_force, index_verbose, index_nbThreads, index_nbThreadsCN, 
       index_keep_tmp, index_time_out, index_help ;

// Variables for projection
int verbose = 0 ;

std::string computeNeighborhoodFile = "computeNeighborhood" ;
std::string computeCentroidsClientFilename ;
std::string registerBundlesClientFile ;
std::string openDipyServerClientFile ;
std::string closeDipyServerClientFile ;
int nbPointsPerFiber = 0 ;
bool doSLR = false ;
bool force = false ;
int nbThreads = 0 ;
int nbThreadsCN = 0 ;
int nbThreadsUsed ;

bool isEsbaAtlas = false ;
bool isFullAtlas = false ;
bool isAtlasNeighborhoodCentroids = false ;
bool isAtlasNeighborhood = false ;

float toleranceThrComputeNeighborhood = 6.0 ;


bool keepTmpFiles = false ;

bool haveMinf = false ;

float time_out = 50 ; // In s


std::string default_ESBA_DIR ;


////////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;
//
void applyLocalSBR( const std::string& movedTractogramNeighborhood,
                    const std::string& atlasBundleDirectory,
                    const std::string& atlasNeighborhoodFile,
                    const std::string& atlasNeighborhoodCentroidsFile,
                    const std::string& outputDirectory,
                    const std::string& referenceImage,
                    const std::string& format,
                    int nbPointsPerFiber,
                    int portDipyServer,
                    float& time_out,
                    int verbose ) ;
//
int getPortNumberDipyService( std::string& logFilePath ) ;
//
void closeDipyServer( int portDipyServer ) ;
