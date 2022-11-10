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

#include "./RecognizedBundles.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class BMDDistance
{
private:
    BundlesDataFormat reference ;
    BundlesDataFormat moving ;
public:
    BMDDistance( BundlesDataFormat bundle1,
                 BundlesDataFormat bundle2 ) ;

    double operator()( const Eigen::VectorXd& affineTransform,
                       Eigen::VectorXd& gradient ) ;

    float normGradient( Eigen::VectorXd& gradient ) ;

    void applyAffineToBundle( const std::vector<float>& bundle,
                              const Eigen::VectorXd& affineCoefficients,
                              int nbCurves,
                              int nbPoints,
                              std::vector<float>& outputBundle ) ;
    void applyAffineToBundle( const std::vector<float>& bundle,
                              const std::vector<float>& affineCoefficients,
                              int nbCurves,
                              int nbPoints,
                              std::vector<float>& outputBundle ) ;

    void applyRigidToBundle( const std::vector<float>& bundle,
                             const Eigen::VectorXd& rigidCoefficients,
                             int nbCurves,
                             int nbPoints,
                             std::vector<float>& outputBundle ) ;
    void applyRigidToBundle( const std::vector<float>& bundle,
                             const std::vector<float>& rigidCoefficients,
                             int nbCurves,
                             int nbPoints,
                             std::vector<float>& outputBundle ) ;

    double computeMDF( const std::vector<float>& tractogramFibers1,
                      const std::vector<float>& tractogramFibers2,
                      const std::vector<float>& medialPointFiber1,
                      const std::vector<float>& medialPointFiber2,
                      int fiberIndex1,
                      int fiberIndex2,
                      int nbPoints ) ;

    double computeBMD( BundlesDataFormat& bundle1,
                      BundlesDataFormat& bundle2 ) ;

    void matrix1DotMatrix2( const std::vector<float>& matrix1,
                            const std::vector<float>& matrix2,
                            std::vector<float>& outputMatrix ) ;

    void matrixMDotVector( const std::vector<float>& matrix,
                           const std::vector<float>& vector,
                           std::vector<float>& outputVector ) ;


} ;
