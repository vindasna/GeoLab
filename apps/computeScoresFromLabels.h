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

#include "AtlasBundles.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int verbose = 0 ;
bool force = false ;
int nbThreads = 0 ;

std::string predictedLabelsPath ;
std::string predictedDictPath ;
std::string trueLabelsPath ;
std::string trueDictPath ;
std::string subjectTractogramPath ;
std::string outputDirectory ;
std::string confusionMatrixPath ;
std::string scoresPath ;


////////////////////////////////////////////////////////////////////////////////
void saveConfusionMatrix(
                    const char* confusionMatrixPath,
                    const std::vector<std::vector<int32_t>>& confusionMatrix ) ;
//
void saveScorePerBundle( const char* scoresPerBundlePath,
                         const std::vector<float>& precisionPerBundle,
                         const std::vector<float>& sensitivityPerBundle,
                         const std::vector<float>& accuracyPerBundle,
                         const std::vector<float>& jaccardPerBundle,
                         const std::vector<std::string>& bundleNames,
                         float averagePrecision,
                         float averageSensitivity,
                         float averageAccuracy,
                         float averageJaccard ) ;
