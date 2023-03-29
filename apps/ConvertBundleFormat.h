// Libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <omp.h>
#include <stdint.h>
#include <unistd.h>
#include <vector>
#include <algorithm>


#include "bundlesData.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int flip_x = 1 ;
int flip_y = 1 ;
int flip_z = 1 ;
int verbose = 0 ;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void bundlesData2Trk( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;



void Trk2BundlesData( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;

void Tck2BundlesData( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;


void BundlesData2Tck( const char* inputFile,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;

void Trk2Tck( const char* inputFile,
              const char* outputFile,
              const char* refFile,
              int flip_x,
              int flip_y,
              int flip_z ) ;

void Tck2Trk( const char* inputFile,
              const char* outputFile,
              const char* refFile,
              int flip_x,
              int flip_y,
              int flip_z ) ;


void flipTractogram( BundlesData& inputTractogram,
                     BundlesMinf& inputTractogramInfo,
                     BundlesData& outputTractogram,
                     int flip_x,
                     int flip_y,
                     int flip_z ) ;


void flipTractogram( const std::vector<float>& inputTractogram,
                     const std::vector<float>& resolution,
                     const std::vector<short int>& size,
                     std::vector<float>& outputTractogram,
                     int flip_x,
                     int flip_y,
                     int flip_z ) ;


void applyVoxToRasToTractogram(
                              const std::vector<float>& tractogram,
                              const std::vector<std::vector<float>>& vox_to_ras,
                              std::vector<float>& newTractogram ) ;

void computeAdjugateMatrix( const std::vector<std::vector<float>>& matrix,
                            std::vector<std::vector<float>>& adjugate ) ;

float computeDeterminantMatrix(
                               const std::vector<std::vector<float>>& matrix ) ;

void computeInverseMatrix( const std::vector<std::vector<float>>& matrix,
                           std::vector<std::vector<float>>& inverseMatrix ) ;

void computeInverseVoxToRas( const std::vector<std::vector<float>>& vox_to_ras,
                             std::vector<std::vector<float>>& ras_to_vox ) ;

int getFlagPosition( int argc,
                     char* argv[],
                     const std::string& flag ) ;
