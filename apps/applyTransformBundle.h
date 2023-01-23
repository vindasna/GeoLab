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

#include <Eigen/Core>
#include <LBFGS.h>

#include "RecognizedBundles.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool isInverse = false ;

int verbose = 0 ;

std::string format ;

// std::vector<float> affineTransform{ 0.985597, -0.00997088, -0.0112213, -0.000119365,
//                            -0.0021475, 0.996759, -0.00195425, -2.0062e-05,
//                            -0.000457125, -0.00363345, 0.997179, -3.82056e-06 } ;
std::vector<float> affineTransform{ 1.0, 0.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0, 0.0,
                                    0.0, 0.0, 1.0, 0.0 } ;

int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;

void readTransform( const std::string& transformFilename,
                    std::vector<float>& affineTransform ) ;


// affineCoefficeints = [ a1, a2, ..., a11 ] is a vector form of the affine
// matrix --->  |  a0   a1   a2   a3  |
//              |  a4   a5   a6   a7  |
//              |  a8   a9  a10  a11 |
//              |  0    0    0    1   |
void applyAffineToBundle( const std::vector<float>& bundle,
                          const std::vector<float>& affineCoefficients,
                          int nbCurves,
                          int nbPoints,
                          std::vector<float>& outputBundle ) ;
//
std::vector<std::vector<float>> matrixProduct(
                              const std::vector<std::vector<float>>& matrix1,
                              const std::vector<std::vector<float>>& matrix2 ) ;
//
void computeAdjugateMatrix4x4( const std::vector<std::vector<float>>& matrix,
                            std::vector<std::vector<float>>& adjugate ) ;
//
float computeDeterminantMatrix(
                              const  std::vector<std::vector<float>>& matrix ) ;
//
void computeInverseMatrix4x4( const std::vector<std::vector<float>>& matrix,
                           std::vector<std::vector<float>>& inverseMatrix ) ;
//
void affineTransformToMatrix( const std::vector<float>& inAffineTransform,
                              std::vector<std::vector<float>>& outMatrix ) ;
//
void matrixToAffineTransform( const std::vector<std::vector<float>>& inMatrix,
                              std::vector<float>& outAffineTransform ) ;
