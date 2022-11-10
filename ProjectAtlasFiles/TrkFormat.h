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

#include "BundlesFormat.h"
#include "BundlesDataFormat.h"


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class TrkFormat
{

  public :
  /////////////////////////////// Public Fields ////////////////////////////////
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
  int curves_count = 0 ;
  int version ;
  int hdr_size ;
  std::vector<float> matrixTracks ;
  std::vector<int32_t> pointsPerTrack ;
  std::vector< std::vector < float > > tracksScalars ;
  std::vector< std::vector < float > > tracksProperties ;
  bool isOk ;

  //////////////////////////////// Constructors ////////////////////////////////
  TrkFormat() ;

  TrkFormat( const char* trkFilename, int verbose ) ;

  TrkFormat( char* id_string,
             short int* dim,
             float* voxel_size,
             float* origin,
             short int n_scalars,
             char scalar_name[ 10 ][ 20 ],
             short int n_properties,
             char property_name[ 10 ][ 20 ],
             float vox_to_ras[ 4 ][ 4 ],
             char* reserved,
             char* voxel_order,
             char* pad2,
             float* image_orientation_patient,
             char* pad1,
             unsigned char invert_x,
             unsigned char invert_y,
             unsigned char invert_z,
             unsigned char swap_xy,
             unsigned char swap_yz,
             unsigned char swap_zx,
             int curves_count,
             int version,
             int hdr_size,
             std::vector<float>& matrixTracks,
             std::vector<int32_t>& pointsPerTrack,
             std::vector< std::vector < float > > tracksScalars,
             std::vector< std::vector < float > > tracksProperties,
             bool isOk ) ;


  //////////////////////////////// Destructor ////////////////////////////////
  virtual ~TrkFormat() ;


  ///////////////////////////////// Methods //////////////////////////////////
  void trkReading( const char* trkFilename, int verbose ) ;

  void trkWriting( const char* trkFilename, int verbose ) ;

  void toBundles( BundlesFormat& bundlesInfo, BundlesDataFormat& bundlesData ) ;

  void printTrkHeader( ) ;



};
