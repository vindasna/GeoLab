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
std::string format ;


////////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;

void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              std::vector<float>& vector1,
                              std::vector<float>& vector2,
                              std::vector<float>& vector3,
                              int nbPoints,
                              std::string outputDirectory,
                              std::string fiberName ) ;

void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              std::vector<float>& vector1,
                              std::vector<float>& vector2,
                              int nbPoints,
                              std::string outputDirectory,
                              std::string fiberName ) ;


void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              int nbPoints,
                              std::string outputDirectory,
                              std::string fiberName ) ;

float scalarProduct( std::vector<float>& vector1,
                     std::vector<float>& vector2 ) ;

void crossProduct( const std::vector<float>& vector1,
                   const std::vector<float>& vector2,
                   std::vector<float>& result ) ;
