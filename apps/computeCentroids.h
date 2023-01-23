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

#include "AtlasBundles.h

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Variables for parse
int index_input, index_reference, index_output, index_cc, index_ods, index_cds,
                     index_nbPoints, index_format, index_nbThreads, index_force,
                                          index_thr, index_verbose, index_help ;

// Variables for projection
int verbose = 0 ;
std::string inputDirPath ;
std::string outputDirectory ;
std::string referenceFilename ;
std::string computeCentroidsClientFilename ;
std::string openDipyServerClientFile ;
std::string closeDipyServerClientFile ;
std::string format ;
bool force = false ;
int nbPointsPerFiber = 0 ;
int nbThreads = 0 ;
float thrQb = 5 ;
int nbThreadsUsed ;

std::vector<std::string> atlasBundlesDataPaths ;
std::vector<std::string> atlasBundlesInfoPaths ;



////////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;
//
int getPortNumberDipyService( std::string& logFilePath ) ;
//
void closeDipyServer( int portDipyServer ) ;
