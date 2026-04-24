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

int nbThreads = 0 ;

int nbThreadsUsed ;

float tolerance = 0.5 ;

///////////////////////////////////////////////////////////////////////////////
void filterBundle( const BundlesData& atlasBundleData,
                   const BundlesMinf& atlasBundleInfo,
                   float tolerance,
                   std::vector<bool>& fibersToKeep ) ;

void saveFeaturesFibersBundle(
                     const std::vector<std::vector<float>>& fibersFeturesValues,
                     int nbPoints,
                     const char* savePath ) ;