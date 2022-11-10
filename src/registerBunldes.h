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
#include <future>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include <Eigen/Core>
#include <LBFGS.h>

#include "BMDDistance.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int verbose = 0 ;
std::string transformType = "rigid" ;
// int maxNbFibers = 50 ;
int maxNbFibers = 100 ;
// int maxNbFibers = 200 ;
// int maxNbFibers = 500 ;
// int maxNbFibers = 1000 ;

bool useDefaultMaxNbFibers = true ;

int nbIterationSearchEpsilon = 10 ;

int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;

float computeDistanceFiberToBundle( BundlesDataFormat& atlasBundleData,
                                    const std::vector<float>& fiber ) ;
