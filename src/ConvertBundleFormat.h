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


#include "BMDDistance.h"


struct trkFormat
{
        char id_string[ 6 ] ;
        short int dim[ 3 ] ;
        float voxel_size[ 3 ] ;
        float origin[ 3 ] ;
        short int n_scalars ;
        char scalar_name[ 10 ][ 20 ] ;
        short int n_properties ;
        char property_name[ 10 ][ 20 ] ;
        float vox_to_ras[ 4 ][ 4 ] ;
        char reserved[ 444 ] ;
        char voxel_order[ 4 ] ;
        char pad2[ 4 ] ;
        float image_orientation_patient[ 6 ] ;
        char pad1[ 2 ] ;
        unsigned char invert_x ;
        unsigned char invert_y ;
        unsigned char invert_z ;
        unsigned char swap_xy ;
        unsigned char swap_yz ;
        unsigned char swap_zx ;
        int n_count ;
        int version ;
        int hdr_size ;
        std::vector< float > matrixTracks ;
        std::vector< int > pointsPerTrack ;
        std::vector< float > tracksScalars ;
        std::vector< float >  tracksProperties ;

        bool isOk ;

 };

 struct tckFormat
 {

   std::vector< float > matrixTracks ;
   std::vector< int32_t > pointsPerTrack ;
   int64_t curves_count = 0 ;
   std::vector< std::string > headerInfo ;
   int sizeDataType = sizeof( float ) ;
   int offsetBinary = 0 ;

 };

struct bundlesMapMinf
{

  std::string byte_order ;
  int curve3d_counts ;
  std::string io_mode ;
  int item_count ;
  std::string label_type ;
  std::string labels ;
  std::string object_type ;
  float resolution[ 3 ] ;
  short int size[ 3 ] ;

};


 struct bundlesMapFormat
 {

   int16_t labels ;
   int32_t curves_count ;
   std::vector< int32_t > pointsPerTrack ;
   std::vector< float > matrixTracks ;

 };

 int flip_x = 1 ;
 int flip_y = 1 ;
 int flip_z = 1 ;
 int verbose = 0 ;

//---------------------------IO files Functions---------------------------------
inline bool is_file( const std::string& path ) ;
trkFormat trkReading( const char* trkFilename ) ;
void trkWriting( const char* trkFilename,
                 trkFormat& trkData ) ;
void printTrkHeader( trkFormat& trkData ) ;

tckFormat tckReading( const char* tckFilename ) ;
void tckWriting( const char* tckFilename,
                 tckFormat& tckData ) ;

BundlesFormat bundlesReading( const char* bundlesFilename ) ;
void bundlesWriting( const char* bundlesFilename,
                     BundlesFormat& bundlesInfo ) ;

BundlesDataFormat bundlesdataReading( const char* bundlesdataFilename,
                                      const char* bundlesFilename ) ;
void bundlesdataWriting( const char* bundlesdataFilename,
                         BundlesDataFormat& bundlesdata ) ;

bundlesMapMinf bundlesMapMinfReading( const char* bundlesMapMinfFilename ) ;
void bundlesMapMinfWriting( const char* bundlesMapMinfFilename,
                            bundlesMapMinf& bundlesMapMinfData ) ;

bundlesMapFormat bundlesMapReading( const char* bundlesMapFilename ) ;
void bundlesMapWriting( const char* bundlesMapFilename,
                        bundlesMapFormat& bundlesMapData ) ;

void bundlesData2Trk( const char* inputBundlesdataFormat,
                      const char* inputBundlesFormat,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;


void Trk2BundlesData( const char* inputFile,
                      const char* refFile,
                      const char* bundlesFilename,
                      const char* bundlesdataFilename,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;

void Tck2BundlesData( const char* inputFile,
                      const char* refFile,
                      const char* bundlesMinfFilename,
                      const char* bundlesFilename,
                      const char* bundlesdataFilename,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;

void BundlesData2Tck( const char* inputBundlesdataFormat,
                      const char* inputBundlesFormat,
                      const char* inputBundlesMinf,
                      const char* refFile,
                      const char* outputTckFilename,
                      int flip_x,
                      int flip_y,
                      int flip_z ) ;


void BundlesMap2BundlesData( const char* inputBundlesMapFilename,
                             const char* inputBundlesMapMinfFilename,
                             const char* bundlesFilename,
                             const char* bundlesdataFilename ) ;

void flipTractogram( BundlesDataFormat& inputTractogram,
                     BundlesFormat& inputTractogramInfo,
                     BundlesDataFormat& outputTractogram,
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
