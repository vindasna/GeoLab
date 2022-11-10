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

#include "BMDDistance.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int verbose = 0 ;

// std::vector<float> affineTransform{ 0.985597, -0.00997088, -0.0112213, -0.000119365,
//                            -0.0021475, 0.996759, -0.00195425, -2.0062e-05,
//                            -0.000457125, -0.00363345, 0.997179, -3.82056e-06 } ;
std::vector<float> affineTransform{ 1.0, 0.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0, 0.0,
                                    0.0, 0.0, 1.0, 0.0 } ;

int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;
