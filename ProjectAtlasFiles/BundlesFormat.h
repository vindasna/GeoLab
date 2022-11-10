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

class BundlesFormat
{

  public :
  /////////////////////////////// Public Fields ////////////////////////////////
  int binary = 0 ;
  std::string bundles = "None" ;
  std::string byte_order = "None" ;
  int curves_count = 0 ;
  float radio = 0 ;
  float length = 0 ;
  int nSubjects = 0 ;
  std::string data_file_name = "None" ;
  std::string format = "None" ;
  std::string io_mode = "\'binary\'" ;
  int item_count = 1 ;
  std::string label_type = "std_string" ;
  std::string labels = "[ '255' ]" ;
  std::string object_type = "None" ;
  std::vector<float> resolution = std::vector<float>( 3, 0 ) ;
  std::vector<short int> size = std::vector<short int>( 3, 0 ) ;
  int space_dimension = 0 ;
  float averageRadius = 0 ;
  float minRadius = 0 ;
  float maxRadius = 0 ;
  float averageAngle = 0 ;
  float minAngle = 0 ;
  float maxAngle = 0 ;
  float averageDirectionAngle = 0 ;
  float minDirectionAngle = 0 ;
  float maxDirectionAngle = 0 ;
  float averageShapeAngle = 0 ;
  float minShapeAngle = 0 ;
  float maxShapeAngle = 0 ;
  float averageLength = 0 ;
  float minLength = 0 ;
  float maxLength = 0 ;
  float averageDisimilarity = 0 ;
  float minDisimilarity = 0 ;
  float maxDisimilarity = 0 ;
  float averageDistanceBetweenMedialPoints = 0 ;
  float minDistanceBetweenMedialPoints = 0 ;
  float maxDistanceBetweenMedialPoints = 0 ;
  std::vector<float> centerBundle = std::vector<float>( 3, 0 ) ;
  float density = 0 ;
  float coverageWithAtlas = -1 ;
  float adjacencyWithAtlas = -1 ;
  float overlapWithAtlas = -1 ;


  //////////////////////////////// Constructors ////////////////////////////////
  BundlesFormat() ;

  BundlesFormat( const char* bundlesFilename, int verbose ) ;

  BundlesFormat( int binary,
                 std::string bundles,
                 std::string byte_order,
                 int curves_count,
                 float radio,
                 float length,
                 int nSubjects,
                 std::string data_file_name,
                 std::string format,
                 std::string io_mode,
                 int item_count,
                 std::string label_type,
                 std::string labels,
                 std::string object_type,
                 float* resolution,
                 short int* size,
                 int space_dimension,
                 float averageRadius,
                 float minRadius,
                 float maxRadius,
                 float averageAngle,
                 float minAngle,
                 float maxAngle,
                 float averageDirectionAngle,
                 float minDirectionAngle,
                 float maxDirectionAngle,
                 float averageShapeAngle,
                 float minShapeAngle,
                 float maxShapeAngle,
                 float averageLength,
                 float minLength,
                 float maxLength,
                 float averageDisimilarity,
                 float minDisimilarity,
                 float maxDisimilarity,
                 float averageDistanceBetweenMedialPoints,
                 float minDistanceBetweenMedialPoints,
                 float maxDistanceBetweenMedialPoints,
                 std::vector<float> centerBundle,
                 float density ) ;

    BundlesFormat( const BundlesFormat& bundlesDataInfo ) ;

    //////////////////////////////// Destructor ////////////////////////////////
    virtual ~BundlesFormat() ;

    ///////////////////////////////// Methods //////////////////////////////////
    void bundlesReading( const char* bundlesFilename, int verbose ) ;

    void bundlesWriting( const char* bundlesFilename,
                         int verbose ) ;
    void bundlesWriting( const char* bundlesFilename,
                         float disimilarity,
                         float coverage,
                         float overlap,
                         float adjacency,
                         int verbose ) ;

    // void toTRK( const BundlesDataFormat& bundlesData, TrkFormat& trkData ) ;

    void printInfo() ;

} ;
