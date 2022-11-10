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
#include <algorithm>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include "BundlesFormat.h"

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesFormat::BundlesFormat() {}

BundlesFormat::BundlesFormat( const char* bundlesFilename, int verbose )
{

  bundlesReading( bundlesFilename, verbose ) ;

}

BundlesFormat::BundlesFormat( int binary,
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
                              float density )
{

  this->binary = binary ;
  this->bundles = bundles ;
  this->byte_order = byte_order ;
  this->curves_count = curves_count ;
  this->radio = radio ;
  this->length = length ;
  this->nSubjects = nSubjects ;
  this->data_file_name = data_file_name ;
  this->format = format ;
  this->io_mode = io_mode ;
  this->item_count = item_count ;
  this->label_type = label_type ;
  this->labels = labels ;
  this->object_type = object_type ;
  this->space_dimension = space_dimension ;
  this->averageRadius = averageRadius ;
  this->minRadius = minRadius ;
  this->maxRadius = maxRadius ;
  this->averageAngle = averageAngle ;
  this->minAngle = minAngle ;
  this->maxAngle = maxAngle ;
  this->averageDirectionAngle = averageDirectionAngle ;
  this->minDirectionAngle = minDirectionAngle ;
  this->maxDirectionAngle = maxDirectionAngle ;
  this->averageShapeAngle = averageShapeAngle ;
  this->minShapeAngle = minShapeAngle ;
  this->maxShapeAngle = maxShapeAngle ;
  this->averageLength = averageLength ;
  this->minLength = minLength ;
  this->maxLength = maxLength ;
  this->averageDisimilarity = averageDisimilarity ;
  this->minDisimilarity = minDisimilarity ;
  this->maxDisimilarity = maxDisimilarity ;
  this->averageDistanceBetweenMedialPoints =
                                            averageDistanceBetweenMedialPoints ;
  this->minDistanceBetweenMedialPoints = minDistanceBetweenMedialPoints ;
  this->maxDistanceBetweenMedialPoints = maxDistanceBetweenMedialPoints ;
  this->density = density ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    this->resolution[ i ] = resolution[ i ] ;
    this->size[ i ] = size[ i ] ;
    this->centerBundle[ i ] = centerBundle[ i ] ;

  }

}

BundlesFormat::BundlesFormat( const BundlesFormat& bundlesDataInfo )
{

  this->binary = bundlesDataInfo.binary ;
  this->bundles = bundlesDataInfo.bundles ;
  this->byte_order = bundlesDataInfo.byte_order ;
  this->curves_count = bundlesDataInfo.curves_count ;
  this->radio = bundlesDataInfo.radio ;
  this->length = bundlesDataInfo.length ;
  this->nSubjects = bundlesDataInfo.nSubjects ;
  this->data_file_name = bundlesDataInfo.data_file_name ;
  this->format = bundlesDataInfo.format ;
  this->io_mode = bundlesDataInfo.io_mode ;
  this->item_count = bundlesDataInfo.item_count ;
  this->label_type = bundlesDataInfo.label_type ;
  this->labels = bundlesDataInfo.labels ;
  this->object_type = bundlesDataInfo.object_type ;
  this->space_dimension = bundlesDataInfo.space_dimension ;
  this->averageRadius = bundlesDataInfo.averageRadius ;
  this->minRadius = bundlesDataInfo.minRadius ;
  this->maxRadius = bundlesDataInfo.maxRadius ;
  this->averageAngle = bundlesDataInfo.averageAngle ;
  this->minAngle = bundlesDataInfo.minAngle ;
  this->maxAngle = bundlesDataInfo.maxAngle ;
  this->averageDirectionAngle = bundlesDataInfo.averageDirectionAngle ;
  this->minDirectionAngle = bundlesDataInfo.minDirectionAngle ;
  this->maxDirectionAngle = bundlesDataInfo.maxDirectionAngle ;
  this->averageShapeAngle = bundlesDataInfo.averageShapeAngle ;
  this->minShapeAngle = bundlesDataInfo.minShapeAngle ;
  this->maxShapeAngle = bundlesDataInfo.maxShapeAngle ;
  this->averageLength = bundlesDataInfo.averageLength ;
  this->minLength = bundlesDataInfo.minLength ;
  this->maxLength = bundlesDataInfo.maxLength ;
  this->averageDisimilarity = bundlesDataInfo.averageDisimilarity ;
  this->minDisimilarity = bundlesDataInfo.minDisimilarity ;
  this->maxDisimilarity = bundlesDataInfo.maxDisimilarity ;
  this->averageDistanceBetweenMedialPoints =
                            bundlesDataInfo.averageDistanceBetweenMedialPoints ;
  this->minDistanceBetweenMedialPoints =
                                bundlesDataInfo.minDistanceBetweenMedialPoints ;
  this->maxDistanceBetweenMedialPoints =
                                bundlesDataInfo.maxDistanceBetweenMedialPoints ;
  this->density = bundlesDataInfo.density ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    this->resolution[ i ] = bundlesDataInfo.resolution[ i ] ;
    this->size[ i ] = bundlesDataInfo.size[ i ] ;
    this->centerBundle[ i ] = bundlesDataInfo.centerBundle[ i ] ;

  }

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesFormat::~BundlesFormat() {}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void BundlesFormat::bundlesReading( const char* bundlesFilename, int verbose )
{

  const char delim = ':' ;
  std::string line ;
  std::ifstream minfFile ;
  minfFile.open( bundlesFilename ) ;
  if ( minfFile.fail() )
  {

    std::cout << "Problem reading file : " << bundlesFilename << std::endl ;
    exit( 1 ) ;

  }
  while ( std::getline( minfFile, line ) )
  {

    std::vector< std::string > out ;
    std::stringstream ss( line ) ;
    std::string s ;
    while ( std::getline( ss, s, delim ) )
    {

      std::string sCopy = s ;
      std::replace( sCopy.begin(), sCopy.end(), '[', ' ' ) ;
      if ( sCopy  == s )
      {

        s.erase(std::remove( s.begin(), s.end(), ' ' ), s.end() ) ;
        s.erase(std::remove( s.begin(), s.end(), ',' ), s.end() ) ;
        out.push_back( s ) ;

      }
      else
      {

        out.push_back( s.substr( 0, s.size() - 1 ) ) ;

      }

    }

    if ( out[ 0 ] == "\'binary\'" )
    {

      this->binary = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'bundles\'" )
    {

      this->bundles = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'byte_order\'" )
    {

      this->byte_order = out[ 1 ] ;


    }

    if ( out[ 0 ] == "\'curves_count\'" )
    {

      this->curves_count = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'radio\'" )
    {

      this->radio = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'length\'" )
    {

      this->length = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'nSubjects\'" )
    {

      this->nSubjects = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'data_file_name\'" )
    {

      this->data_file_name = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'format\'" )
    {

      this->format = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'io_mode\'" )
    {

      this->io_mode = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'item_count\'" )
    {

      this->item_count = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'label_type\'" )
    {

      this->label_type = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'labels\'" )
    {

      this->labels = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'object_type\'" )
    {

      this->object_type = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'resolutionX\'" )
    {

      this->resolution[ 0 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'resolutionY\'" )
    {

      this->resolution[ 1 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'resolutionZ\'" )
    {

      this->resolution[ 2 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeX\'" )
    {

      this->size[ 0 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeY\'" )
    {

      this->size[ 1 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeZ\'" )
    {

      this->size[ 2 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'space_dimension\'" )
    {

      this->space_dimension = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageRadius\'" )
    {

      this->averageRadius = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minRadius\'" )
    {

      this->minRadius = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxRadius\'" )
    {

      this->maxRadius = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageAngle\'" )
    {

      this->averageAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minAngle\'" )
    {

      this->minAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxAngle\'" )
    {

      this->maxAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageDirectionAngle\'" )
    {

      this->averageDirectionAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minDirectionAngle\'" )
    {

      this->minDirectionAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxDirectionAngle\'" )
    {

      this->maxDirectionAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageShapeAngle\'" )
    {

      this->averageShapeAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minShapeAngle\'" )
    {

      this->minShapeAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxShapeAngle\'" )
    {

      this->maxShapeAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageLength\'" )
    {

      this->averageLength = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minLength\'" )
    {

      this->minLength = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxLength\'" )
    {

      this->maxLength = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageDisimilarity\'" )
    {

      this->averageDisimilarity = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minDisimilarity\'" )
    {

      this->minDisimilarity = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxDisimilarity\'" )
    {

      this->maxDisimilarity = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageDistanceBetweenMedialPoints\'" )
    {

      this->averageDistanceBetweenMedialPoints = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minDistanceBetweenMedialPoints\'" )
    {

      this->minDistanceBetweenMedialPoints = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxDistanceBetweenMedialPoints\'" )
    {

      this->maxDistanceBetweenMedialPoints = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'centerBundleX\'" )
    {

      this->centerBundle[ 0 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'centerBundleY\'" )
    {

      this->centerBundle[ 1 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'centerBundleZ\'" )
    {

      this->centerBundle[ 2 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'density\'" )
    {

      this->density = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'CoverageWithAtlas\'" )
    {

      this->coverageWithAtlas = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'AdjacencyWithAtlas\'" )
    {

      this->adjacencyWithAtlas = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'OverlapWithAtlas\'" )
    {

      this->overlapWithAtlas = std::stof( out[ 1 ] ) ;

    }

  }

  minfFile.close() ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesFormat::bundlesWriting( const char* bundlesFilename,
                                    int verbose )
{

  if ( verbose > 1 )
  {

    std::cout << "\nWriting " << bundlesFilename << std::endl ;

  }



  std::ofstream file( bundlesFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path "
              << std::endl ;
    exit( 1 );

  }

  file << "attributes = {\n" ;

  file << "    'binary' : "
       << this->binary
       << ",\n" ;

  file << "    'bundles' : "
       << this->bundles
       << ",\n" ;

  file << "    'byte_order' : "
       << this->byte_order
       << ",\n" ;

  file << "    'curves_count' : "
       << this->curves_count
       << ",\n" ;

  file << "    'radio' : "
       << this->radio
       << ",\n" ;

  file << "    'length' : "
       << this->length
       << ",\n" ;

  file << "    'nSubjects' : "
       << this->nSubjects
       << ",\n" ;

  file << "    'data_file_name' : "
       << this->data_file_name
       << ",\n" ;

  file << "    'format' : "
       << this->format
       << ",\n" ;

  file << "    'io_mode' : "
       << this->io_mode
       << ",\n" ;

  file << "    'item_count' : "
       << this->item_count
       << ",\n" ;

  file << "    'label_type' : "
       << this->label_type
       << ",\n" ;

  file << "    'labels' : "
       << this->labels
       << ",\n" ;

  file << "    'object_type' : "
       << this->object_type
       << ",\n" ;

  file << "    'averageRadius' : "
       << this->averageRadius
       << ",\n" ;

  file << "    'minRadius' : "
       << this->minRadius
       << ",\n" ;

  file << "    'maxRadius' : "
       << this->maxRadius
       << ",\n" ;

  file << "    'averageAngle' : "
       << this->averageAngle
       << ",\n" ;

  file << "    'minAngle' : "
       << this->minAngle
       << ",\n" ;

  file << "    'maxAngle' : "
       << this->maxAngle
       << ",\n" ;

  file << "    'averageDirectionAngle' : "
       << this->averageDirectionAngle
       << ",\n" ;

  file << "    'minDirectionAngle' : "
       << this->minDirectionAngle
       << ",\n" ;

  file << "    'maxDirectionAngle' : "
       << this->maxDirectionAngle
       << ",\n" ;

  file << "    'averageShapeAngle' : "
       << this->averageShapeAngle
       << ",\n" ;

  file << "    'minShapeAngle' : "
       << this->minShapeAngle
       << ",\n" ;

  file << "    'maxShapeAngle' : "
       << this->maxShapeAngle
       << ",\n" ;

  file << "    'averageLength' : "
       << this->averageLength
       << ",\n" ;

  file << "    'minLength' : "
       << this->minLength
       << ",\n" ;

  file << "    'maxLength' : "
       << this->maxLength
       << ",\n" ;

  file << "    'averageDisimilarity' : "
       << this->averageDisimilarity
       << ",\n" ;

  file << "    'minDisimilarity' : "
       << this->minDisimilarity
       << ",\n" ;

  file << "    'maxDisimilarity' : "
       << this->maxDisimilarity
       << ",\n" ;

  file << "    'averageDistanceBetweenMedialPoints' : "
       << this->averageDistanceBetweenMedialPoints
       << ",\n" ;

  file << "    'minDistanceBetweenMedialPoints' : "
       << this->minDistanceBetweenMedialPoints
       << ",\n" ;

  file << "    'maxDistanceBetweenMedialPoints' : "
       << this->maxDistanceBetweenMedialPoints
       << ",\n" ;

  file << "    'density' : "
       << this->density
       << ",\n" ;

  file << "    'centerBundleX' : "
       << this->centerBundle[ 0 ]
       << ",\n" ;

  file << "    'centerBundleY' : "
       << this->centerBundle[ 1 ]
       << ",\n" ;

  file << "    'centerBundleZ' : "
       << this->centerBundle[ 2 ]
       << ",\n" ;

  file << "    'resolutionX' : "
       << this->resolution[ 0 ]
       << ",\n" ;

  file << "    'resolutionY' : "
       << this->resolution[ 1 ]
       << ",\n" ;

  file << "    'resolutionZ' : "
       << this->resolution[ 2 ]
       << ",\n" ;

  file << "    'sizeX' : "
       << this->size[ 0 ]
       << ",\n" ;

  file << "    'sizeY' : "
       << this->size[ 1 ]
       << ",\n" ;

  file << "    'sizeZ' : "
       << this->size[ 2 ]
       << ",\n" ;

  file << "    'space_dimension' : "
       << this->space_dimension
       << "\n    }" ;

  file.close() ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesFormat::bundlesWriting( const char* bundlesFilename,
                                    float disimilarity,
                                    float coverage,
                                    float overlap,
                                    float adjacency,
                                    int verbose )
{

  if ( verbose > 1 )
  {

    std::cout << "\nWriting " << bundlesFilename << std::endl ;

  }



  std::ofstream file( bundlesFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path "
              << std::endl ;
    exit( 1 );

  }

  file << "attributes = {\n" ;

  file << "    'binary' : "
       << this->binary
       << ",\n" ;

  file << "    'bundles' : "
       << this->bundles
       << ",\n" ;

  file << "    'byte_order' : "
       << this->byte_order
       << ",\n" ;

  file << "    'curves_count' : "
       << this->curves_count
       << ",\n" ;

  file << "    'radio' : "
       << this->radio
       << ",\n" ;

  file << "    'length' : "
       << this->length
       << ",\n" ;

  file << "    'nSubjects' : "
       << this->nSubjects
       << ",\n" ;

  file << "    'data_file_name' : "
       << this->data_file_name
       << ",\n" ;

  file << "    'format' : "
       << this->format
       << ",\n" ;

  file << "    'io_mode' : "
       << this->io_mode
       << ",\n" ;

  file << "    'item_count' : "
       << this->item_count
       << ",\n" ;

  file << "    'label_type' : "
       << this->label_type
       << ",\n" ;

  file << "    'labels' : "
       << this->labels
       << ",\n" ;

  file << "    'object_type' : "
       << this->object_type
       << ",\n" ;

  file << "    'averageRadius' : "
       << this->averageRadius
       << ",\n" ;

  file << "    'minRadius' : "
       << this->minRadius
       << ",\n" ;

  file << "    'maxRadius' : "
       << this->maxRadius
       << ",\n" ;

  file << "    'averageAngle' : "
       << this->averageAngle
       << ",\n" ;

  file << "    'minAngle' : "
       << this->minAngle
       << ",\n" ;

  file << "    'maxAngle' : "
       << this->maxAngle
       << ",\n" ;

  file << "    'averageDirectionAngle' : "
       << this->averageDirectionAngle
       << ",\n" ;

  file << "    'minDirectionAngle' : "
       << this->minDirectionAngle
       << ",\n" ;

  file << "    'maxDirectionAngle' : "
       << this->maxDirectionAngle
       << ",\n" ;

  file << "    'averageShapeAngle' : "
       << this->averageShapeAngle
       << ",\n" ;

  file << "    'minShapeAngle' : "
       << this->minShapeAngle
       << ",\n" ;

  file << "    'maxShapeAngle' : "
       << this->maxShapeAngle
       << ",\n" ;

  file << "    'averageLength' : "
       << this->averageLength
       << ",\n" ;

  file << "    'minLength' : "
       << this->minLength
       << ",\n" ;

  file << "    'maxLength' : "
       << this->maxLength
       << ",\n" ;

  file << "    'averageDisimilarity' : "
       << this->averageDisimilarity
       << ",\n" ;

  file << "    'minDisimilarity' : "
       << this->minDisimilarity
       << ",\n" ;

  file << "    'maxDisimilarity' : "
       << this->maxDisimilarity
       << ",\n" ;

  file << "    'averageDistanceBetweenMedialPoints' : "
       << this->averageDistanceBetweenMedialPoints
       << ",\n" ;

  file << "    'minDistanceBetweenMedialPoints' : "
       << this->minDistanceBetweenMedialPoints
       << ",\n" ;

  file << "    'maxDistanceBetweenMedialPoints' : "
       << this->maxDistanceBetweenMedialPoints
       << ",\n" ;

  file << "    'density' : "
       << this->density
       << ",\n" ;

  file << "    'disimilarityWithAtlas' : "
       << disimilarity
       << ",\n" ;

  file << "    'CoverageWithAtlas' : "
       << coverage
       << ",\n" ;

  file << "    'OverlapWithAtlas' : "
       << overlap
       << ",\n" ;

  file << "    'AdjacencyWithAtlas' : "
       << adjacency
       << ",\n" ;

  file << "    'centerBundleX' : "
       << this->centerBundle[ 0 ]
       << ",\n" ;

  file << "    'centerBundleY' : "
       << this->centerBundle[ 1 ]
       << ",\n" ;

  file << "    'centerBundleZ' : "
       << this->centerBundle[ 2 ]
       << ",\n" ;

  file << "    'resolutionX' : "
       << this->resolution[ 0 ]
       << ",\n" ;

  file << "    'resolutionY' : "
       << this->resolution[ 1 ]
       << ",\n" ;

  file << "    'resolutionZ' : "
       << this->resolution[ 2 ]
       << ",\n" ;

  file << "    'sizeX' : "
       << this->size[ 0 ]
       << ",\n" ;

  file << "    'sizeY' : "
       << this->size[ 1 ]
       << ",\n" ;

  file << "    'sizeZ' : "
       << this->size[ 2 ]
       << ",\n" ;

  file << "    'space_dimension' : "
       << this->space_dimension
       << "\n    }" ;

  file.close() ;

}


////////////////////////////////////////////////////////////////////////////////
// void BundlesFormat::toTRK( const BundlesDataFormat& bundlesData,
//                            TrkFormat& trkData )
// {
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     trkData.dim[ i ] = this->resolution[ i ] ;
//     trkData.voxel_size[ i ] = this->size[ i ] ;
//
//   }
//
//   trkData.curves_count = bundlesData.curves_count ;
//   trkData.matrixTracks = bundlesData.matrixTracks ;
//   trkData.pointsPerTrack = bundlesData.pointsPerTrack ;
//
// }


////////////////////////////////////////////////////////////////////////////////
void BundlesFormat::printInfo()
{

  std::cout << "binary : " << this->binary << std::endl ;
  std::cout << "bundles : " << this->bundles << std::endl ;
  std::cout << "byte_order : " << this->byte_order << std::endl ;
  std::cout << "curves_count : " << this->curves_count << std::endl ;
  std::cout << "radio : " << this->radio << std::endl ;
  std::cout << "length : " << this->length << std::endl ;
  std::cout << "nSubjects : " << this->nSubjects << std::endl ;
  std::cout << "data_file_name : " << this->data_file_name << std::endl ;
  std::cout << "format : " << this->format << std::endl ;
  std::cout << "io_mode : " << this->io_mode << std::endl ;
  std::cout << "item_count : " << this->item_count << std::endl ;
  std::cout << "label_type : " << this->label_type << std::endl ;
  std::cout << "labels : " << this->labels << std::endl ;
  std::cout << "object_type : " << this->object_type << std::endl ;
  std::cout << "space_dimension : " << this->space_dimension << std::endl ;
  std::cout << "averageRadius : " << this->averageRadius << std::endl ;
  std::cout << "averageAngle : " << this->averageAngle << std::endl ;
  std::cout << "averageDirectionAngle : " << this->averageDirectionAngle
                                                                  << std::endl ;
  std::cout << "averageShapeAngle : " << this->averageShapeAngle << std::endl ;
  std::cout << "averageLength : " << this->averageLength << std::endl ;
  std::cout << "averageDisimilarity : " << this->averageDisimilarity
                                                                  << std::endl ;
  std::cout << "density : " << this->density << std::endl ;


  std::cout << "resolution : [ " << this->resolution[ 0 ] << ", "
                                 << this->resolution[ 1 ] << ", "
                                 << this->resolution[ 2 ] << " ]" << std::endl ;
  std::cout << "size : [ " << this->size[ 0 ] << ", "
                           << this->size[ 1 ] << ", "
                           << this->size[ 2 ] << " ]" << std::endl ;
  std::cout << "centerBundle : [ " << this->centerBundle[ 0 ] << ", "
                                   << this->centerBundle[ 1 ] << ", "
                                   << this->centerBundle[ 2 ] << " ]"
                                   << std::endl ;

}
