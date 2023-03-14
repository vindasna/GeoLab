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



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class BundlesMinf
{

  private :
  std::vector<float> mniResolution = { 1.0, 1.0, 1.0 } ;
  std::vector<short int> mniSize = { 193, 229, 193 } ;
  std::vector<std::vector<float>> mni_vox_to_ras = { { 1.0, 0.0, 0.0, -98.0 },
                                                     { 0.0, 1.0, 0.0, -134.0 },
                                                     { 0.0, 0.0, 1.0, -72.0 },
                                                     { 0.0, 0.0, 0.0, 1.0 }
                                                   } ;
  float epsilon = 1e-3 ;

  public :
  /////////////////////////////// Public Fields ////////////////////////////////
  // Format of input tractogram
  bool isBundles = false ;
  bool isTrk = false ;
  bool isTck = false ;
  bool haveMinf = false ;

  // Common to all formats
  int64_t curves_count = 0 ;
  std::vector<float> resolution = std::vector<float>( 3, 0 ) ; // voxel_size in .trk
  std::vector<short int> size = std::vector<short int>( 3, 0 ) ; // dim in .trk
  float averageRadius = -1 ;
  float minRadius = -1 ;
  float maxRadius = -1 ;
  float averageAngle = -1 ;
  float minAngle = -1 ;
  float maxAngle = -1 ;
  float averageDirectionAngle = -1 ;
  float minDirectionAngle = -1 ;
  float maxDirectionAngle = -1 ;
  float averageShapeAngle = -1 ;
  float minShapeAngle = -1 ;
  float maxShapeAngle = -1 ;
  float averageLength = -1 ;
  float minLength = -1 ;
  float maxLength = -1 ;
  float averageDisimilarity = -1 ;
  float minDisimilarity = -1 ;
  float maxDisimilarity = -1 ;
  float averageDistanceBetweenMedialPoints = -1 ;
  float minDistanceBetweenMedialPoints = -1 ;
  float maxDistanceBetweenMedialPoints = -1 ;
  std::vector<float> centerBundle = std::vector<float>( 3, -1 ) ;
  float density = -1 ;
  float disimilarityWithAtlas = -1 ;
  float coverageWithAtlas = -1 ;
  float adjacencyWithAtlas = -1 ;
  float overlapWithAtlas = -1 ;
  std::string bundleName ;

  //////////////////////////// For .bundles format /////////////////////////////
  int binary = 0 ;
  std::string bundles = "None" ;
  std::string byte_order = "None" ;
  std::string data_file_name = "None" ;
  std::string format = "None" ;
  std::string io_mode = "\'binary\'" ;
  int item_count = 1 ;
  std::string label_type = "std_string" ;
  std::string labels = "[ '255' ]" ;
  std::string object_type = "None" ;
  int space_dimension = 0 ;


  ////////////////////////////// For .trk format ///////////////////////////////
  std::vector<char> id_string = std::vector<char>( 6 ) ;
  std::vector<float> origin = std::vector<float>( 3, 0 ) ;
  short int n_scalars = 0 ;
  std::vector<std::vector<char>> scalar_name = std::vector<std::vector<char>>(
                                            10, std::vector<char>( 20 ) ) ;
  short int n_properties = 0 ;
  std::vector<std::vector<char>> property_name = std::vector<std::vector<char>>(
                                            10, std::vector<char>( 20 ) ) ;
  std::vector<std::vector<float>> vox_to_ras = std::vector<std::vector<float>>(
                                               4, std::vector<float>( 4, 0 ) ) ;
  std::vector<char> reserved = std::vector<char>( 444 ) ;
  std::vector<char> voxel_order = std::vector<char>( 4 ) ;
  std::vector<char> pad2 = std::vector<char>( 4 ) ;
  std::vector<float> image_orientation_patient = std::vector<float>( 6, 0 ) ;
  std::vector<char> pad1 = std::vector<char>( 2 ) ;
  unsigned char invert_x = 0 ;
  unsigned char invert_y = 0 ;
  unsigned char invert_z = 0 ;
  unsigned char swap_xy = 0 ;
  unsigned char swap_yz = 0 ;
  unsigned char swap_zx = 0 ;
  int version = 2 ;
  int hdr_size = 1000 ;

  ////////////////////////////// For .tck format ///////////////////////////////
  std::vector<std::string> tckHeaderInfo ;
  int tckSizeDataType = sizeof( float ) ;
  int tckOffsetBinary = 0 ;
  std::string dataTypeBinary = "" ;


  //////////////////////////////// Constructors ////////////////////////////////
  BundlesMinf() ;

  BundlesMinf( const char* bundlesFilename ) ;

  BundlesMinf( const BundlesMinf& bundlesInfo ) ;

  ///////////////////////////////// Destructor /////////////////////////////////
  virtual ~BundlesMinf() ;


  ///////////////////////////////// Operators //////////////////////////////////
  bool operator == ( const BundlesMinf& bundleInfo ) const;
  bool operator != ( const BundlesMinf& bundleInfo ) const;

  ////////////////////////////////// Methods ///////////////////////////////////
  void read( const char* bundlesFilename ) ;

  void write( const char* bundlesFilename,
              float disimilarity,
              float coverage,
              float overlap,
              float adjacency,
              bool haveMinf ) const ;

  void write( const char* bundlesFilename,
              float disimilarity,
              float coverage,
              float overlap,
              float adjacency ) const ;

  void write( const char* bundlesFilename, bool haveMinf ) const ;

  void write( const char* bundlesFilename ) const ;

  void fillDefaultTrk() ;

  void fillDefaultTrkMni() ;

  void fillDefaultTck() ;

  void fillDefaultBundles() ;

  void computeTckOffsetBinary() ;

  void computeHeaderTck( int curvesCount ) ; // Only for datatype float

  std::string getFormat() const ;

  std::vector<std::vector<float>> getMniVoxToRas() const ;

  std::vector<float> getMniResolution() const ;

  std::vector<short int> getMniSize() const ;


} ;

////////////////////////////// External Operators //////////////////////////////
std::ostream& operator << ( std::ostream& outStream,
                            const BundlesMinf& bundleInfo ) ;
