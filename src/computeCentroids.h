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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int verbose = 0 ;
int minNbFibers = 100 ;
int maxNbClusters = 300 ;
float thresholdDistance = 0 ;
float thresholdAdjacency = 5 ; // In mm

void computeMedialPointsFibersTractogram(
           BundlesDataFormat& inputTractogram,
           std::vector< std::vector< float > >& medialPointsFibersTractogram ) ;

double getAverageMinDistanceBetweenFibers(
     BundlesDataFormat& inputTractogram,
     const std::vector< std::vector< float > >& medialPointsFibersTractogram ) ;

float computeDistanceCenterFibers(
        const std::vector< std::vector< float > >& medialPointsFibersTractogram,
        int indexFiber1,
        int indexFiber2 ) ;

void computeMeanFiber( const std::vector< std::vector< float > >& fibers,
                       std::vector<float>& meanFiber,
                       int nbCurves,
                       int nbPoints ) ;
