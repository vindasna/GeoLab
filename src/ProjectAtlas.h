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
float p = 5 ; // In mm
float thrDistance = 10 ; // In mm
float minimumLenghtFiber = 10 ; // In mm
float maximumLenghtFiber = 200 ; // In mm
float maxAngle = 45 ; // In degrees
float maxDirectionAngle = 45 ; // In degrees
float minShapeAngle = 0 ; // In degrees
float maxShapeAngle = 180 ; // In degrees
float thrDistanceBetweenMedialPoints = 5 ; // In mm
bool compareRecognizedWithAtlas = false ;
bool useMDFDistance = false ;
bool useMedialPointAverageFiber = false ;
bool useSimpleProjection = false ;
bool useAvgThr = false ;

// float toleranceP = 0.0 ;
// float toleranceThr = 0.0 ;
// float toleranceMaxAngle = -0.5 ;
// float toleranceMaxDirectionAngle = -0.5 ;
// float toleranceMinShapeAngle = 0.0 ;
// float toleranceMaxShapeAngle = 0.0 ;
// float toleranceLenght = 0.3 ;

// float toleranceP = 0.0 ;
// float toleranceThr = 0.7 ;
// float toleranceMaxAngle = 0.7 ;
// float toleranceMaxDirectionAngle = 0.7 ;
// float toleranceMinShapeAngle = 0.9 ;
// float toleranceMaxShapeAngle = 0.9 ;
// float toleranceLenght = 0.3 ;

float toleranceP = 0.0 ;
float toleranceThr = 0.0 ;
float toleranceMaxAngle = 0.0 ;
float toleranceMaxDirectionAngle = 0.0 ;
float toleranceMinShapeAngle = 0.0 ;
float toleranceMaxShapeAngle = 0.0 ;
float toleranceLenght = 0.0 ;
float toleranceDistanceBetweenMedialPoints = 0.0 ;

/*
// float toleranceP = 0.0 ;
// float toleranceP = -0.1 ;
float toleranceP = -0.1 ;
// float toleranceThr = 0.2 ;
float toleranceThr = 0.1 ;
float toleranceMaxAngle = 0.0 ;
float toleranceMaxDirectionAngle = 0.1 ;
float toleranceMinShapeAngle = 0.0 ;
float toleranceMaxShapeAngle = 0.0 ;
float toleranceLenght = 0.0 ;
*/

float thrPercentageSimilarity = 0.05 ;
float minimumNumberFibers = 20 ;
float thresholdAdjacency = 10 ; // In mm
// The name useDefautl* can be confusing because it does not enable using the
// default parameter but it enables to use the parameter if passed as input
bool useDefautlP = false ;
bool useDefaultThr = false ;
bool useDefaultMaxAngle = false ;
bool useDefautlMaxDirectionAngle = false ;
bool useDefautlMinShapeAngle = false ;
bool useDefautlMaxShapeAngle = false ;
bool useDefaultMinLen = false ;
bool useDefaultMaxLen = false ;
bool useDefaultThrDistanceBetweenMedialPoints = false ;
bool saveExtractedBundles = true ;
bool saveUnlabeled = false ;


////////////////////////////////////////////////////////////////////////////////
void fillArrayPointer( float* array,
                       int64_t nbElements,
                       float x ) ;

int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;
