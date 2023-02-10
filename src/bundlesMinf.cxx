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

#include <stdexcept>

#include "bundlesMinf.h"
#include "ioWrapper.h"

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesMinf::BundlesMinf() {}

BundlesMinf::BundlesMinf( const char* bundlesFilename )
{

  read( bundlesFilename ) ;

}


BundlesMinf::BundlesMinf( const BundlesMinf& bundlesInfo )
{

  // Format of input tractogram
  this->isBundles = bundlesInfo.isBundles ;
  this->isTck = bundlesInfo.isTck ;
  this->isTrk = bundlesInfo.isTrk ;
  this->haveMinf = bundlesInfo.haveMinf ;

  // Information common to all formats
  this->curves_count = bundlesInfo.curves_count ;
  this->averageRadius = bundlesInfo.averageRadius ;
  this->minRadius = bundlesInfo.minRadius ;
  this->maxRadius = bundlesInfo.maxRadius ;
  this->averageAngle = bundlesInfo.averageAngle ;
  this->minAngle = bundlesInfo.minAngle ;
  this->maxAngle = bundlesInfo.maxAngle ;
  this->averageDirectionAngle = bundlesInfo.averageDirectionAngle ;
  this->minDirectionAngle = bundlesInfo.minDirectionAngle ;
  this->maxDirectionAngle = bundlesInfo.maxDirectionAngle ;
  this->averageShapeAngle = bundlesInfo.averageShapeAngle ;
  this->minShapeAngle = bundlesInfo.minShapeAngle ;
  this->maxShapeAngle = bundlesInfo.maxShapeAngle ;
  this->averageLength = bundlesInfo.averageLength ;
  this->minLength = bundlesInfo.minLength ;
  this->maxLength = bundlesInfo.maxLength ;
  this->averageDisimilarity = bundlesInfo.averageDisimilarity ;
  this->minDisimilarity = bundlesInfo.minDisimilarity ;
  this->maxDisimilarity = bundlesInfo.maxDisimilarity ;
  this->averageDistanceBetweenMedialPoints =
                            bundlesInfo.averageDistanceBetweenMedialPoints ;
  this->minDistanceBetweenMedialPoints =
                                bundlesInfo.minDistanceBetweenMedialPoints ;
  this->maxDistanceBetweenMedialPoints =
                                bundlesInfo.maxDistanceBetweenMedialPoints ;
  this->density = bundlesInfo.density ;
  this->disimilarityWithAtlas = bundlesInfo.disimilarityWithAtlas ;
  this->coverageWithAtlas = bundlesInfo.coverageWithAtlas ;
  this->adjacencyWithAtlas = bundlesInfo.adjacencyWithAtlas ;
  this->overlapWithAtlas = bundlesInfo.overlapWithAtlas ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    this->resolution[ i ] = bundlesInfo.resolution[ i ] ;
    this->size[ i ] = bundlesInfo.size[ i ] ;
    this->centerBundle[ i ] = bundlesInfo.centerBundle[ i ] ;

  }

  // Information for .bundles format
  this->binary = bundlesInfo.binary ;
  this->bundles = bundlesInfo.bundles ;
  this->byte_order = bundlesInfo.byte_order ;
  this->data_file_name = bundlesInfo.data_file_name ;
  this->format = bundlesInfo.format ;
  this->io_mode = bundlesInfo.io_mode ;
  this->item_count = bundlesInfo.item_count ;
  this->label_type = bundlesInfo.label_type ;
  this->labels = bundlesInfo.labels ;
  this->object_type = bundlesInfo.object_type ;
  this->space_dimension = bundlesInfo.space_dimension ;


  // Information for .trk format
  this->id_string = bundlesInfo.id_string ;
  this->origin = bundlesInfo.origin ;
  this->n_scalars = bundlesInfo.n_scalars ;
  this->scalar_name = bundlesInfo.scalar_name ;
  this->n_properties = bundlesInfo.n_properties ;
  this->property_name = bundlesInfo.property_name ;
  this->vox_to_ras = bundlesInfo.vox_to_ras ;
  this->reserved = bundlesInfo.reserved ;
  this->voxel_order = bundlesInfo.voxel_order ;
  this->pad2 = bundlesInfo.pad2 ;
  this->image_orientation_patient = bundlesInfo.image_orientation_patient ;
  this->pad1 = bundlesInfo.pad1 ;
  this->invert_x = bundlesInfo.invert_x ;
  this->invert_y = bundlesInfo.invert_y ;
  this->invert_z = bundlesInfo.invert_z ;
  this->swap_xy = bundlesInfo.swap_xy ;
  this->swap_yz = bundlesInfo.swap_yz ;
  this->swap_zx = bundlesInfo.swap_zx ;
  this->version = bundlesInfo.version ;
  this->hdr_size = bundlesInfo.hdr_size ;

  // Information for .tck format
  this->tckHeaderInfo = bundlesInfo.tckHeaderInfo ;
  this->tckSizeDataType = bundlesInfo.tckSizeDataType ;
  this->tckOffsetBinary = bundlesInfo.tckOffsetBinary ;
  this->dataTypeBinary = bundlesInfo.dataTypeBinary ;

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesMinf::~BundlesMinf() {}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// External Operators ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
std::ostream& operator << ( std::ostream& outStream,
                            const BundlesMinf& bundleInfo )
{

  if ( bundleInfo.isBundles )
  {

    outStream << "Format : .bundlesdata " << std::endl
              << "binary : " << bundleInfo.binary << std::endl
              << "bundles : " << bundleInfo.bundles << std::endl
              << "byte_order : " << bundleInfo.byte_order << std::endl
              << "data_file_name : " << bundleInfo.data_file_name << std::endl
              << "format : " << bundleInfo.format << std::endl
              << "io_mode : " << bundleInfo.io_mode << std::endl
              << "item_count : " << bundleInfo.item_count << std::endl
              << "label_type : " << bundleInfo.label_type << std::endl
              << "labels : " << bundleInfo.labels << std::endl
              << "object_type : " << bundleInfo.object_type << std::endl
              << "space_dimension : " << bundleInfo.space_dimension
              << std::endl ;

  }
  if ( bundleInfo.isTrk )
  {

    outStream << "Format : .trk " << std::endl ;


    outStream << "id_string : " ;
    for ( char c : bundleInfo.id_string )
    {

      outStream << c ;

    }
    outStream << "\n" ;


    outStream << "origin : " ;
    for ( float f : bundleInfo.origin )
    {

      outStream << f << "   " ;

    }
    outStream << "\n" ;


    outStream << "n_scalars : " << bundleInfo.n_scalars << std::endl ;


    outStream << "scalar_name : " ;
    for ( std::vector<char> v : bundleInfo.scalar_name )
    {

      for ( char c : v )
      {

        outStream << c << "   " ;

      }
      outStream << "\n              " ;

    }
    outStream << "\r" ;


    outStream << "n_properties : " << bundleInfo.n_properties << std::endl ;


    outStream << "property_name : " ;
    for ( std::vector<char> v : bundleInfo.property_name )
    {

      for ( char c : v )
      {

        outStream << c << "   " ;

      }
      outStream << "\n              " ;

    }
    outStream << "\r" ;


    outStream << "vox_to_ras : " ;
    for ( std::vector<float> v : bundleInfo.vox_to_ras )
    {

      for ( float f : v )
      {

        outStream << f << "   " ;

      }
      outStream << "\n              " ;

    }
    outStream << "\r" ;


    outStream << "reserved : " ;
    for ( char c : bundleInfo.reserved )
    {

      outStream << c ;

    }
    outStream << "\n" ;


    outStream << "voxel_order : " ;
    for ( char c : bundleInfo.voxel_order )
    {

      outStream << c ;

    }
    outStream << "\n" ;


    outStream << "pad2 : " ;
    for ( char c : bundleInfo.pad2 )
    {

      outStream << c ;

    }
    outStream << "\n" ;


    outStream << "image_orientation_patient : " ;
    for ( float f : bundleInfo.image_orientation_patient )
    {

      outStream << f << "   " ;

    }
    outStream << "\n" ;


    outStream << "pad1 : " ;
    for ( char c : bundleInfo.pad1 )
    {

      outStream << c ;

    }
    outStream << "\n" ;


    outStream  << "invert_x : " << bundleInfo.invert_x << std::endl
               << "invert_y : " << bundleInfo.invert_y << std::endl
               << "invert_z : " << bundleInfo.invert_z << std::endl
               << "swap_xy : " << bundleInfo.swap_xy << std::endl
               << "swap_yz : " << bundleInfo.swap_yz << std::endl
               << "swap_zx : " << bundleInfo.swap_zx << std::endl
               << "version : " << bundleInfo.version << std::endl
               << "hdr_size : " << bundleInfo.hdr_size << std::endl ;


  }



  if ( bundleInfo.isTck )
  {

    outStream << "Format : .tck " << std::endl ;

    outStream << "tckHeaderInfo : " ;
    for ( std::string s : bundleInfo.tckHeaderInfo )
    {

      outStream << s << "\n                " ;

    }
    outStream << "\r" ;


    outStream << "tckSizeDataType : " << bundleInfo.tckSizeDataType << std::endl
              << "tckOffsetBinary : " << bundleInfo.tckOffsetBinary << std::endl
              << "dataTypeBinary : " << bundleInfo.dataTypeBinary << std::endl ;


  }

  //////////////////////////// Common to all formats ///////////////////////////
  outStream << "curves_count : " << bundleInfo.curves_count << std::endl ;

  outStream << "resolution : " ;
  for ( float f : bundleInfo.resolution )
  {

    outStream << f << "   " ;

  }
  outStream << "\n" ;

  outStream << "size : " ;
  for ( short int i : bundleInfo.size )
  {

    outStream << i << "   " ;

  }
  outStream << "\n" ;

  outStream << "averageRadius : " << bundleInfo.averageRadius << std::endl
            << "minRadius : " << bundleInfo.minRadius << std::endl
            << "maxRadius : " << bundleInfo.maxRadius << std::endl
            << "averageAngle : " << bundleInfo.averageAngle << std::endl
            << "minAngle : " << bundleInfo.minAngle << std::endl
            << "maxAngle : " << bundleInfo.maxAngle << std::endl
            << "averageDirectionAngle : " << bundleInfo.averageDirectionAngle
                                                                    << std::endl
            << "minDirectionAngle : " << bundleInfo.minDirectionAngle
                                                                    << std::endl
            << "maxDirectionAngle : " << bundleInfo.maxDirectionAngle
                                                                    << std::endl
            << "averageShapeAngle : " << bundleInfo.averageShapeAngle
                                                                    << std::endl
            << "minShapeAngle : " << bundleInfo.minShapeAngle << std::endl
            << "maxShapeAngle : " << bundleInfo.maxShapeAngle << std::endl
            << "averageLength : " << bundleInfo.averageLength << std::endl
            << "minLength : " << bundleInfo.minLength << std::endl
            << "maxLength : " << bundleInfo.maxLength << std::endl
            << "averageDisimilarity : " << bundleInfo.averageDisimilarity
                                                                    << std::endl
            << "minDisimilarity : " << bundleInfo.minDisimilarity << std::endl
            << "maxDisimilarity : " << bundleInfo.maxDisimilarity << std::endl
            << "averageDistanceBetweenMedialPoints : "
                   << bundleInfo.averageDistanceBetweenMedialPoints << std::endl
            << "minDistanceBetweenMedialPoints : "
                       << bundleInfo.minDistanceBetweenMedialPoints << std::endl
            << "maxDistanceBetweenMedialPoints : "
                     << bundleInfo.maxDistanceBetweenMedialPoints << std::endl ;

  outStream << "centerBundle : " ;
  for ( float f : bundleInfo.centerBundle )
  {

    outStream << f << "   " ;

  }
  outStream << "\n" ;

  outStream << "density : " << bundleInfo.density << std::endl
            << "disimilarityWithAtlas : " << bundleInfo.disimilarityWithAtlas
                                                                    << std::endl
            << "coverageWithAtlas : " << bundleInfo.coverageWithAtlas
                                                                    << std::endl
            << "adjacencyWithAtlas : " << bundleInfo.adjacencyWithAtlas
                                                                    << std::endl
            << "overlapWithAtlas : " << bundleInfo.overlapWithAtlas
                                                                  << std::endl ;


  return( outStream ) ;


}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Operators ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool BundlesMinf::operator == ( const BundlesMinf& bundleInfo ) const
{

  bool isOk = ( ( isBundles == bundleInfo.isBundles )  &&
                ( isTrk == bundleInfo.isTrk ) &&
                ( isTck == bundleInfo.isTck ) &&
                ( curves_count == bundleInfo.curves_count ) &&
                ( std::abs( averageRadius - bundleInfo.averageRadius )
                                                                  < epsilon ) &&
                ( std::abs( minRadius - bundleInfo.minRadius ) < epsilon ) &&
                ( std::abs( maxRadius - bundleInfo.maxRadius ) < epsilon ) &&
                ( std::abs( averageAngle - bundleInfo.averageAngle )
                                                                  < epsilon ) &&
                ( std::abs( minAngle - bundleInfo.minAngle ) < epsilon ) &&
                ( std::abs( maxAngle - bundleInfo.maxAngle ) < epsilon ) &&
                ( std::abs( averageDirectionAngle -
                               bundleInfo.averageDirectionAngle ) < epsilon ) &&
                ( std::abs( minDirectionAngle - bundleInfo.minDirectionAngle )
                                                                  < epsilon ) &&
                ( std::abs( maxDirectionAngle - bundleInfo.maxDirectionAngle )
                                                                  < epsilon ) &&
                ( std::abs( averageShapeAngle - bundleInfo.averageShapeAngle )
                                                                  < epsilon ) &&
                ( std::abs( minShapeAngle - bundleInfo.minShapeAngle )
                                                                  < epsilon ) &&
                ( std::abs( maxShapeAngle - bundleInfo.maxShapeAngle )
                                                                  < epsilon ) &&
                ( std::abs( averageLength - bundleInfo.averageLength )
                                                                  < epsilon ) &&
                ( std::abs( minLength - bundleInfo.minLength ) < epsilon ) &&
                ( std::abs( maxLength - bundleInfo.maxLength ) < epsilon ) &&
                ( std::abs( averageDisimilarity
                               - bundleInfo.averageDisimilarity ) < epsilon ) &&
                ( std::abs( minDisimilarity - bundleInfo.minDisimilarity )
                                                                  < epsilon ) &&
                ( std::abs( maxDisimilarity - bundleInfo.maxDisimilarity )
                                                                  < epsilon ) &&
                ( std::abs( averageDistanceBetweenMedialPoints -
                                 bundleInfo.averageDistanceBetweenMedialPoints )
                                                                  < epsilon ) &&
                ( std::abs( minDistanceBetweenMedialPoints -
                                     bundleInfo.minDistanceBetweenMedialPoints )
                                                                  < epsilon ) &&
                ( std::abs( maxDistanceBetweenMedialPoints -
                                     bundleInfo.maxDistanceBetweenMedialPoints )
                                                                  < epsilon ) &&
                ( std::abs( density - bundleInfo.density ) < epsilon ) &&
                ( std::abs( disimilarityWithAtlas -
                               bundleInfo.disimilarityWithAtlas ) < epsilon ) &&
                ( std::abs( coverageWithAtlas -
                                   bundleInfo.coverageWithAtlas ) < epsilon ) &&
                ( std::abs( adjacencyWithAtlas -
                                  bundleInfo.adjacencyWithAtlas ) < epsilon ) &&
                ( std::abs( overlapWithAtlas -
                                    bundleInfo.overlapWithAtlas ) < epsilon ) &&
                ( binary == bundleInfo.binary ) &&
                ( bundles == bundleInfo.bundles ) &&
                ( byte_order == bundleInfo.byte_order ) &&
                ( data_file_name == bundleInfo.data_file_name ) &&
                ( format == bundleInfo.format ) &&
                ( io_mode == bundleInfo.io_mode ) &&
                ( item_count == bundleInfo.item_count ) &&
                ( label_type == bundleInfo.label_type ) &&
                ( labels == bundleInfo.labels ) &&
                ( object_type == bundleInfo.object_type ) &&
                ( space_dimension == bundleInfo.space_dimension ) &&
                ( n_scalars == bundleInfo.n_scalars ) &&
                ( n_properties == bundleInfo.n_properties ) &&
                ( invert_x == bundleInfo.invert_x ) &&
                ( invert_y == bundleInfo.invert_y ) &&
                ( invert_z == bundleInfo.invert_z ) &&
                ( swap_xy == bundleInfo.swap_xy ) &&
                ( swap_yz == bundleInfo.swap_yz ) &&
                ( swap_zx == bundleInfo.swap_zx ) &&
                ( version == bundleInfo.version ) &&
                ( hdr_size == bundleInfo.hdr_size ) &&
                ( tckSizeDataType == bundleInfo.tckSizeDataType ) &&
                ( tckOffsetBinary == bundleInfo.tckOffsetBinary ) &&
                ( dataTypeBinary == bundleInfo.dataTypeBinary ) ) ;



  if ( isOk )
  {

    for ( int i = 0 ; i < resolution.size() ; i++ )
    {

      if ( std::abs( resolution[ i ] - bundleInfo.resolution[ i ] ) > epsilon )
      {

        return( false ) ;

      }

      if ( size[ i ] != bundleInfo.size[ i ] )
      {

        return( false ) ;

      }

      if ( std::abs( origin[ i ] - bundleInfo.origin[ i ] ) > epsilon )
      {

        return( false ) ;

      }

    }

    for ( int i = 0 ; i < id_string.size() ; i++ )
    {

      if ( id_string[ i ] != bundleInfo.id_string[ i ] )
      {

        return( false ) ;

      }

      if ( std::abs( image_orientation_patient[ i ] -
                         bundleInfo.image_orientation_patient[ i ] ) > epsilon )
      {

        return( false ) ;

      }

    }


    for ( int i = 0 ; i < scalar_name.size() ; i++ )
    {

      for ( int j = 0 ; j < scalar_name[ i ].size() ; j++ )
      {

        if ( scalar_name[ i ][ j ] != bundleInfo.scalar_name[ i ][ j ] )
        {

          return( false ) ;

        }

        if ( property_name[ i ][ j ] != bundleInfo.property_name[ i ][ j ] )
        {

          return( false ) ;

        }

      }

    }

    for ( int i = 0 ; i < vox_to_ras.size() ; i++ )
    {

      for ( int j = 0 ; j < vox_to_ras[ i ].size() ; j++ )
      {

        if ( std::abs( vox_to_ras[ i ][ j ] - bundleInfo.vox_to_ras[ i ][ j ] )
                                                                     > epsilon )
        {

          return( false ) ;

        }

      }

      if ( voxel_order[ i ] != bundleInfo.voxel_order[ i ] )
      {

        return( false ) ;

      }

      if ( pad2[ i ] != bundleInfo.pad2[ i ] )
      {

        return( false ) ;

      }

    }

    for ( int i = 0 ; i < reserved.size() ; i++ )
    {

      if ( reserved[ i ] != bundleInfo.reserved[ i ] )
      {

        return( false ) ;

      }

    }

    if ( pad1[ 0 ] != bundleInfo.pad1[ 0 ] ||
                                             pad1[ 1 ] != bundleInfo.pad1[ 1 ] )
    {

      return( false ) ;

    }

    if ( tckHeaderInfo.size() == bundleInfo.tckHeaderInfo.size() )
    {

      for ( int i = 0 ; i < tckHeaderInfo.size() ; i++ )
      {

        if ( tckHeaderInfo[ i ] != bundleInfo.tckHeaderInfo[ i ] )
        {

          return( false ) ;

        }

      }

    }
    else
    {

      return( false ) ;

    }

  }


  return( isOk ) ;

}


////////////////////////////////////////////////////////////////////////////////
bool BundlesMinf::operator != ( const BundlesMinf& bundleInfo ) const
{

  return( !( *this == bundleInfo ) ) ;

}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// For details on .tck http://trackvis.org/docs/?subsect=fileformat
void BundlesMinf::read( const char* bundlesFilename )
{

  std::string bundlesFilenameStr = bundlesFilename ;
  std::string bundlesDataFilenameStr = bundlesFilename ;

  if ( endswith( bundlesFilenameStr, ".bundles" ) )
  {

    isBundles = true ;
    haveMinf = true ;
    bundlesDataFilenameStr = replaceExtension( bundlesDataFilenameStr,
                                                              ".bundlesdata" ) ;

  }
  else if ( endswith( bundlesFilenameStr, ".bundlesdata" ) )
  {

    isBundles = true ;
    haveMinf = true ;
    bundlesFilenameStr = replaceExtension( bundlesFilenameStr, ".bundles" ) ;

  }
  else if ( endswith( bundlesFilenameStr, ".trk" ) )
  {

    isTrk = true ;
    bundlesFilenameStr = replaceExtension( bundlesFilenameStr, ".minf" ) ;
    if ( !is_file( bundlesFilenameStr ) )
    {

      haveMinf = false ;

    }
    else
    {

      haveMinf = true ;

    }

  }
  else if ( endswith( bundlesFilenameStr, ".tck" ) )
  {

    isTck = true ;
    bundlesFilenameStr = replaceExtension( bundlesFilenameStr, ".minf" ) ;
    if ( !is_file( bundlesFilenameStr ) )
    {

      haveMinf = false ;

    }
    else
    {

      haveMinf = true ;

    }

  }
  else if ( endswith( bundlesFilenameStr, ".minf" ) )
  {

    if ( !is_file( bundlesFilenameStr ) )
    {

      haveMinf = false ;

    }
    else
    {

      haveMinf = true ;

    }

    if ( is_file( replaceExtension( bundlesFilenameStr, ".trk" ) ) )
    {

      isTrk = true ;
      bundlesDataFilenameStr = replaceExtension( bundlesDataFilenameStr,
                                                                       ".trk" ) ;


    }
    else if ( is_file( replaceExtension( bundlesFilenameStr, ".tck" ) ) )
    {

      if ( isTrk )
      {

        std::ostringstream outMessageOss ;
        outMessageOss << "BundlesMinf::read : input given is a .minf but found "
                      << "a .trk and .tck associated with the .minf in the "
                      << "directory : "
                      << bundlesFilenameStr << std::endl ;

        std::string outMessage = outMessageOss.str() ;

        throw( std::invalid_argument( outMessage ) ) ;

      }
      else
      {

        isTck = true ;
        bundlesDataFilenameStr = replaceExtension( bundlesDataFilenameStr,
                                                                      ".tck" ) ;

      }

    }
    else
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesMinf::read : when input is .minf, the "
                    << "corresponding .trk or .tck file must be in the same "
                    << "directory" << bundlesFilenameStr << std::endl ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    std::ostringstream outMessageOss ;
    outMessageOss << "BundlesMinf::read : only formats supported are .bundles/"
                  << ".bundlesdata, .trk, .tck but got "
                  << bundlesFilenameStr << std::endl ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }



  if ( haveMinf )
  {

    const char delim = ':' ;
    std::string line ;
    std::ifstream minfFile ;
    minfFile.open( bundlesFilenameStr.c_str() ) ;
    if ( minfFile.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesMinf::read : Problem reading file : "
                    << bundlesFilenameStr << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }
    while ( std::getline( minfFile, line ) )
    {

      std::vector< std::string > out ;
      std::stringstream ss( line ) ;
      std::string s ;
      while ( std::getline( ss, s, delim ) )
      {

        std::string sCopy = s ;
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

      if ( out[ 0 ] == "\'disimilarityWithAtlas\'" )
      {

        this->disimilarityWithAtlas = std::stof( out[ 1 ] ) ;

      }

      if ( out[ 0 ] == "\'coverageWithAtlas\'" )
      {

        this->coverageWithAtlas = std::stof( out[ 1 ] ) ;

      }

      if ( out[ 0 ] == "\'adjacencyWithAtlas\'" )
      {

        this->adjacencyWithAtlas = std::stof( out[ 1 ] ) ;

      }

      if ( out[ 0 ] == "\'overlapWithAtlas\'" )
      {

        this->overlapWithAtlas = std::stof( out[ 1 ] ) ;

      }

      ///////////////////////////// For .bundles format //////////////////////////
      if ( isBundles )
      {

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


        if ( out[ 0 ] == "\'space_dimension\'" )
        {

          this->space_dimension = std::stoi( out[ 1 ] ) ;

        }

        if ( out[ 0 ] == "\'curves_count\'" )
        {

          this->curves_count = std::stoi( out[ 1 ] ) ;

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

      }

    }

    minfFile.close() ;

    if ( isBundles )
    {

      this->fillDefaultTrk() ;

    }

  }


  ////////////////////////////// For .trk format ///////////////////////////////
  if ( isTrk )
  {

    std::streampos size ;

    std::ifstream file ;
    file.open( bundlesDataFilenameStr.c_str(), std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesMinf::read : Problem reading file : "
                    << bundlesDataFilenameStr << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

    for ( int i = 0 ; i < this->id_string.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->id_string[ i ] ) ),
                                              sizeof( this->id_string[ i ] ) ) ;

    }

    for ( int i = 0 ; i < this->size.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->size[ i ] ) ),
                                                   sizeof( this->size[ i ] ) ) ;

    }

    for ( int i = 0 ; i < this->resolution.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->resolution[ i ] ) ),
                                             sizeof( this->resolution[ i ] ) ) ;

    }

    for ( int i = 0 ; i < this->origin.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->origin[ i ] ) ),
                                                 sizeof( this->origin[ i ] ) ) ;

    }

    file.read( reinterpret_cast<char*>( &( this->n_scalars ) ),
                                                  sizeof( this->n_scalars ) ) ;

    for ( int i = 0 ; i < this->scalar_name.size() ; i++ )
    {

      for ( int j = 0 ; j < this->scalar_name[ i ].size() ; j++ )
      {

        file.read( reinterpret_cast<char*>( &( this->scalar_name[ i ][ j ] ) ),
                                       sizeof( this->scalar_name[ i ][ j ] ) ) ;

      }

    }

    file.read( reinterpret_cast<char*>( &( this->n_properties ) ),
                                               sizeof( this->n_properties ) ) ;

    for ( int i = 0 ; i < this->property_name.size() ; i++ )
    {

      for ( int j = 0 ; j < this->property_name[ i ].size() ; j++ )
      {

        file.read( reinterpret_cast<char*>(
                                     &( this->property_name[ i ][ j ] ) ),
                                     sizeof( this->property_name[ i ][ j ] ) ) ;


      }

    }

    for ( int i = 0 ; i < this->vox_to_ras.size() ; i++ )
    {

      for ( int j = 0 ; j < this->vox_to_ras[ i ].size() ; j++ )
      {

        file.read( reinterpret_cast<char*>( &( this->vox_to_ras[ i ][ j ] ) ),
                                        sizeof( this->vox_to_ras[ i ][ j ] ) ) ;

      }

    }

    for ( int i = 0 ; i < this->reserved.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->reserved[ i ] ) ),
                                               sizeof( this->reserved[ i ] ) ) ;

    }


    for ( int i = 0 ; i < this->voxel_order.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->voxel_order[ i ] ) ),
                                            sizeof( this->voxel_order[ i ] ) ) ;

    }


    for ( int i = 0 ; i < this->pad2.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->pad2[ i ] ) ),
                                                   sizeof( this->pad2[ i ] ) ) ;


    }

    for ( int i = 0 ; i < this->image_orientation_patient.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>(
                              &( this->image_orientation_patient[ i ] ) ),
                              sizeof( this->image_orientation_patient[ i ] ) ) ;

    }

    for ( int i = 0 ; i < this->pad1.size() ; i++ )
    {

      file.read( reinterpret_cast<char*>( &( this->pad1[ i ] ) ),
                                                   sizeof( this->pad1[ i ] ) ) ;


    }

    file.read( reinterpret_cast<char*>( &( this->invert_x ) ),
                                                   sizeof( this->invert_x ) ) ;
    file.read( reinterpret_cast<char*>( &( this->invert_y ) ),
                                                   sizeof( this->invert_y ) ) ;
    file.read( reinterpret_cast<char*>( &( this->invert_z ) ),
                                                   sizeof( this->invert_z ) ) ;
    file.read( reinterpret_cast<char*>( &( this->swap_xy ) ),
                                                    sizeof( this->swap_xy ) ) ;
    file.read( reinterpret_cast<char*>( &( this->swap_yz ) ),
                                                    sizeof( this->swap_yz ) ) ;
    file.read( reinterpret_cast<char*>( &( this->swap_zx ) ),
                                                    sizeof( this->swap_zx ) ) ;

    int tmpCurveCount = 0 ; // Need this because in .trk curves count is
                            // encoded in int but here I use int64_t
    file.read( reinterpret_cast<char*>( &( tmpCurveCount) ),
                                                     sizeof( tmpCurveCount ) ) ;
    this->curves_count = tmpCurveCount ;
    file.read( reinterpret_cast<char*>( &( this->version ) ),
                                                    sizeof( this->version ) ) ;
    file.read( reinterpret_cast<char*>( &( this->hdr_size ) ),
                                                   sizeof( this->hdr_size ) ) ;



    if ( this->curves_count == 0 )
    {

      const char* outMessage = "BundlesMinf::read : Unable to read the "
                               "number of tracks in the header\n" ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    if ( this->hdr_size != 1000 )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesMinf::read : invalid header size, got "
                    << this->hdr_size << " but it should be 1000 in file "
                    << bundlesDataFilenameStr << std::endl ;
      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    file.close() ;

  }

  ////////////////////////////// For .tck format ///////////////////////////////
  if ( isTck )
  {

    std::fstream infile ;
    infile.open( bundlesDataFilenameStr.c_str(), std::ios::in ) ;

    if ( infile.is_open() )
    {

      std::string line ;
      bool notFinished = true ;
      std::string spaceDelimiter = " " ;
      while ( std::getline( infile, line ) && notFinished )
      {

        this->tckHeaderInfo.push_back( line ) ;

        std::vector<std::string> words ;
        std::size_t pos = 0 ;
        while ( ( pos = line.find( spaceDelimiter ) ) != std::string::npos )
        {

          words.push_back( line.substr( 0, pos ) ) ;
          line.erase( 0, pos + spaceDelimiter.length() ) ;

        }

        words.push_back( line ) ;

        if ( words[ 0 ] == "file:" )
        {

          this->tckOffsetBinary = std::stoi( words[ words.size() - 1 ] ) ;

        }

        if ( words[ 0 ] == "datatype:" )
        {

          this->dataTypeBinary = words[ words.size() - 1 ] ;

        }

        if ( words[ 0 ] == "count:" )
        {

          this->curves_count = std::stoi( words[ words.size() - 1 ] ) ;

        }

        if ( ( words[ 0 ] == "END" ) )
        {

          notFinished = false ;

        }

      }

      infile.close() ;

    }
    else
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesMinf::read : Problem reading file : "
                    << bundlesDataFilenameStr << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


    }


    if ( this->dataTypeBinary == "Float32LE" ||
                                           this->dataTypeBinary == "Float32BE" )
    {

      this->tckSizeDataType = sizeof( float ) ;

    }
    else
    {

      this->tckSizeDataType = sizeof( double ) ;

    }

    this->fillDefaultTrk() ;

  }

}


////////////////////////////////////////////////////////////////////////////////
// Details on .tck https://readthedocs.org/projects/mrtrix/downloads/pdf/latest/
// Details on .trk http://trackvis.org/docs/?subsect=fileformat
void BundlesMinf::write( const char* bundlesFilename,
                         float disimilarity,
                         float coverage,
                         float overlap,
                         float adjacency,
                         bool haveMinf ) const
{

  if ( !haveMinf )
  {

    return ;

  }

  bool outIsBundles = false ;
  bool outIsTrk = false ;
  bool outIsTck = false ;

  std::string bundlesFilenameStr = bundlesFilename ;

  if ( endswith( bundlesFilenameStr, ".bundles" ) )
  {

    outIsBundles = true ;

  }
  else if ( endswith( bundlesFilenameStr, ".trk" ) )
  {

    outIsTrk = true ;
    bundlesFilenameStr = replaceExtension( bundlesFilenameStr, ".minf" ) ;

  }
  else if ( endswith( bundlesFilenameStr, ".tck" ) )
  {

    outIsTck = true ;
    bundlesFilenameStr = replaceExtension( bundlesFilenameStr, ".minf" ) ;

  }
  else if ( endswith( bundlesFilenameStr, ".minf" ) )
  {

    outIsTck = true ; // This is enough because the .minf is written the
                      // same way for .tck and .trk
  }
  else
  {

    const char* outMessage = "BundlesMinf::write : only formats supported are "
                                                   ".bundles, .trk and .tck\n" ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  // For .bundle format
  if ( outIsBundles )
  {

    std::ofstream file( bundlesFilenameStr.c_str() ) ;
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
         << ",\n" ;

    file << "    'disimilarityWithAtlas' : "
         << disimilarity
         << ",\n" ;

    file << "    'coverageWithAtlas' : "
         << coverage
         << ",\n" ;

    file << "    'overlapWithAtlas' : "
         << overlap
         << ",\n" ;

    file << "    'adjacencyWithAtlas' : "
         << adjacency
         << "\n    }" ;

    file.close() ;

  }

  // For .trk format
  if ( outIsTrk || outIsTck )
  {

    std::ofstream file( bundlesFilenameStr.c_str() ) ;
    if ( !file )
    {

      std::cout << "Cannot save file, there's a problem with the saving path "
                << std::endl ;
      exit( 1 );

    }

    file << "attributes = {\n" ;

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

    file << "    'disimilarityWithAtlas' : "
         << disimilarity
         << ",\n" ;

    file << "    'coverageWithAtlas' : "
         << coverage
         << ",\n" ;

    file << "    'overlapWithAtlas' : "
         << overlap
         << ",\n" ;

    file << "    'adjacencyWithAtlas' : "
         << adjacency
         << "\n    }" ;


    file.close() ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::write( const char* bundlesFilename,
                         float disimilarity,
                         float coverage,
                         float overlap,
                         float adjacency ) const
{


  write( bundlesFilename,
         disimilarity,
         coverage,
         overlap,
         adjacency,
         this->haveMinf ) ;
}

////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::write( const char* bundlesFilename, bool haveMinf ) const
{

  this->write( bundlesFilename,
               this->disimilarityWithAtlas,
               this->coverageWithAtlas,
               this->overlapWithAtlas,
               this->adjacencyWithAtlas,
               haveMinf ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::write( const char* bundlesFilename ) const
{

  this->write( bundlesFilename,
               this->disimilarityWithAtlas,
               this->coverageWithAtlas,
               this->overlapWithAtlas,
               this->adjacencyWithAtlas,
               this->haveMinf ) ;

}



////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::fillDefaultTrk()
{

  char tmpIdString[ 6 ] ;
  memcpy( tmpIdString, "TRACK", 6 * sizeof( char ) ) ;
  for ( int i = 0 ; i < 6 ; i++ )
  {

    this->id_string[ i ] = tmpIdString[ i ] ;
  }

  char tmpVoxelOrder[ 4 ] ;
  memcpy( tmpVoxelOrder, "LPI", 4 * sizeof( char ) ) ;
  for ( int i = 0 ; i < 4 ; i++ )
  {

    this->voxel_order[ i ] = tmpVoxelOrder[ i ] ;

  }

  for ( int i = 0 ; i < 4 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++ )
    {

      if ( i == j )
      {

        this->vox_to_ras[ i ][ j ] = 1 ;

      }
      else
      {

        this->vox_to_ras[ i ][ j ] = 0 ;

      }

    }

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    this->origin[ i ] = 0 ;

  }

  for ( int i = 0 ; i < 6 ; i++ )
  {

    this->image_orientation_patient[ i ] = 0 ;

  }

  this->n_scalars = 0 ;

  this->n_properties = 0 ;


}

////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::fillDefaultTrkMni()
{

  char tmpIdString[ 6 ] ;
  memcpy( tmpIdString, "TRACK", 6 * sizeof( char ) ) ;
  for ( int i = 0 ; i < 6 ; i++ )
  {

    this->id_string[ i ] = tmpIdString[ i ] ;
  }

  char tmpVoxelOrder[ 4 ] ;
  memcpy( tmpVoxelOrder, "LPI", 4 * sizeof( char ) ) ;
  for ( int i = 0 ; i < 4 ; i++ )
  {

    this->voxel_order[ i ] = tmpVoxelOrder[ i ] ;

  }


  for ( int i = 0 ; i < 3 ; i++ )
  {

    this->origin[ i ] = 0 ;

  }

  for ( int i = 0 ; i < 6 ; i++ )
  {

    this->image_orientation_patient[ i ] = 0 ;

  }


  this->resolution = this->mniResolution ;

  this->size = this->mniSize ;

  this->mni_vox_to_ras = this->mni_vox_to_ras ;

  this->n_scalars = 0 ;

  this->n_properties = 0 ;


}

////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::fillDefaultTck()
{

  std::string matrixTracksString = "mrtrix tracks    " ;

  std::string timeStampString = "timestamp: 1661165504.3768622875" ;

  std::string datatypeString = "datatype: Float32LE" ;

  std::ostringstream countStringOss ;
  countStringOss << "count: " << this->curves_count ;
  std::string countString = countStringOss.str() ;

  std::ostringstream totalCountStringOss ;
  totalCountStringOss << "total_count: " << this->curves_count ;
  std::string totalCountString = totalCountStringOss.str() ;

  std::ostringstream fileStringOss ;
  fileStringOss << "file: . " ;
  std::string fileString = fileStringOss.str() ;

  std::string endString = "END" ;

  int tmpSizeBytes = matrixTracksString.size() + timeStampString.size() +
                            datatypeString.size() + countString.size() +
                               totalCountString.size() + fileString.size() +
                                                      endString.size() + 7 ;

  std::ostringstream tmpSizeBytesOss ;
  tmpSizeBytesOss << tmpSizeBytes ;
  std::string tmpSizeBytesString = tmpSizeBytesOss.str() ;

  tmpSizeBytes += tmpSizeBytesString.size() ;

  fileStringOss << tmpSizeBytes ;
  fileString = fileStringOss.str() ;

  this->tckHeaderInfo.clear() ;

  this->tckHeaderInfo.push_back( matrixTracksString ) ;
  this->tckHeaderInfo.push_back( timeStampString ) ;
  this->tckHeaderInfo.push_back( datatypeString ) ;
  this->tckHeaderInfo.push_back( fileString ) ;
  this->tckHeaderInfo.push_back( countString ) ;
  this->tckHeaderInfo.push_back( totalCountString ) ;
  this->tckHeaderInfo.push_back( endString ) ;

  this->tckOffsetBinary = tmpSizeBytes + 1 ;

  this->dataTypeBinary = "Float32LE" ;

  this-> tckSizeDataType = sizeof( float ) ;


}



////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::fillDefaultBundles()
{

  this->binary = 1 ;
  this->bundles = "[ '255', 0 ]" ;
  this->byte_order = "'DCBA'" ;
  this->data_file_name = "'*.bundlesdata'" ;
  this->format = "'bundles_1.0'" ;
  this->io_mode = "'binary'" ;
  this->item_count = ( int )1 ;
  this->label_type = "'std_string'" ;
  this->labels = "[ '255' ]" ;
  this->object_type = "'BundleMap'" ;
  this->space_dimension = ( int )3 ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::computeTckOffsetBinary()
{

  int tmpSizeBytes = 0 ;
  for ( std::string s : this->tckHeaderInfo )
  {

    tmpSizeBytes += s.size() ;

  }

  tmpSizeBytes += 7 ;

  this->tckOffsetBinary = tmpSizeBytes ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesMinf::computeHeaderTck( int curvesCount )
{

  std::string matrixTracksString = "mrtrix tracks    " ;

  std::string timeStampString = "timestamp: 1661165504.3768622875" ;

  std::string datatypeString = "datatype: Float32LE" ;

  std::ostringstream countStringOss ;
  countStringOss << "count: " << curves_count ;
  std::string countString = countStringOss.str() ;

  std::ostringstream totalCountStringOss ;
  totalCountStringOss << "total_count: " << curves_count ;
  std::string totalCountString = totalCountStringOss.str() ;

  std::ostringstream fileStringOss ;
  fileStringOss << "file: . " ;
  std::string fileString = fileStringOss.str() ;

  std::string endString = "END" ;

  int tmpSizeBytes = matrixTracksString.size() + timeStampString.size() +
                            datatypeString.size() + countString.size() +
                               totalCountString.size() + fileString.size() +
                                                      endString.size() + 7 ;

  std::ostringstream tmpSizeBytesOss ;
  tmpSizeBytesOss << tmpSizeBytes ;
  std::string tmpSizeBytesString = tmpSizeBytesOss.str() ;

  tmpSizeBytes += tmpSizeBytesString.size() ;

  fileStringOss << tmpSizeBytes ;
  fileString = fileStringOss.str() ;

  this->tckHeaderInfo.clear() ;

  this->tckHeaderInfo.push_back( matrixTracksString ) ;
  this->tckHeaderInfo.push_back( timeStampString ) ;
  this->tckHeaderInfo.push_back( datatypeString ) ;
  this->tckHeaderInfo.push_back( fileString ) ;
  this->tckHeaderInfo.push_back( countString ) ;
  this->tckHeaderInfo.push_back( totalCountString ) ;
  this->tckHeaderInfo.push_back( endString ) ;

  this->tckOffsetBinary = tmpSizeBytes + 1 ;

}

////////////////////////////////////////////////////////////////////////////////
std::string BundlesMinf::getFormat() const
{

  if ( this->isBundles )
  {

    return( ".bundles" ) ;

  }
  else if ( this->isTrk )
  {

    return( ".trk" ) ;

  }
  else if ( this->isTck )
  {

    return( ".tck" ) ;

  }
  else
  {

    std::string outMessage = "BundlesMinf::getFormat : attributes isBundles, " \
                             "isTrk and isTck are all set to false, not " \
                             "supported format found \n" ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

}


////////////////////////////////////////////////////////////////////////////////
std::vector<float> BundlesMinf::getMniResolution() const
{

  return( this->mniResolution ) ;

}



////////////////////////////////////////////////////////////////////////////////
std::vector<short int> BundlesMinf::getMniSize() const
{

  return( this->mniSize ) ;

}


////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<float>> BundlesMinf::getMniVoxToRas() const
{

  return( this->mni_vox_to_ras ) ;

}
