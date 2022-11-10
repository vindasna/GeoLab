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

#include "RecognizedBundles.h"


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int verbose = 0 ;
float thresholdDistance = 10 ; // In mm
float minLength = 0 ; // In mm
float maxLength = 200 ; // In mm
float toleranceThr = 1.2 ;
bool useDefaultMaxLength = true ;
bool isBundlesFormat = false ;
bool isTRKFormat = false ;
// int minNbCurvesNeighborhood = 50 ;
int minNbCurvesNeighborhood = 500 ;

////////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;

void saveIndexInTractogram( const char* labelsBinaryFilename,
                                          const std::vector<int64_t>& labels ) ;

void readIndexInTractogram( const char* predictedLabelsFilename,
                             std::vector<int64_t>& predictedLabels,
                             int nbFibers,
                             int verbose ) ;


float computeDistanceFiberToBundle( BundlesDataFormat& atlasBundleData,
                                    const std::vector<float>& fiber ) ;
