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
// int minimumNumberFibers = 20 ;
int minimumNumberFibers = 0 ;
bool saveUnlabeled = false ;
bool isSupWMA = false ;
std::vector<std::vector<int16_t>> predictedLabels ;


////////////////////////////////////////////////////////////////////////////////
void fillArrayPointer( float* array,
                       int64_t nbElements,
                       float x ) ;

void readingPredictedLabels( const char* predictedLabelsFilename,
                             std::vector<std::vector<int16_t>>& predictedLabels,
                             int nbFibers,
                             int verbose ) ;

void readingPredictedLabelsSupWMA(
                             const char* predictedLabelsFilename,
                             std::vector<std::vector<int16_t>>& predictedLabels,
                             int nbFibers,
                             int verbose ) ;

void readingDictionaryLabels( const char* dictionaryLabelsFilename,
                              std::vector<std::string>& dictionaryLabels,
                              int verbose ) ;

int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;
