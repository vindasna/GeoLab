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

#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <LBFGS.h>

#include "niftiImage.h"
#include "bundlesData.h"
#include "ioWrapper.h"


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesData::BundlesData() {}

BundlesData::BundlesData( const char* bundlesdataFilename )
{

  read( bundlesdataFilename ) ;

}

BundlesData::BundlesData( std::vector<float>& matrixTracks,
                          std::vector<int32_t>& pointsPerTrack,
                          std::vector<int64_t>& fibersWithNans,
                          std::vector<std::vector<float>>& tracksScalars,
                          std::vector<std::vector<float>>& tracksProperties,
                          int64_t curves_count,
                          bool isBundles,
                          bool isTrk,
                          bool isTck )
{

  this->matrixTracks = matrixTracks ;
  this->pointsPerTrack = pointsPerTrack ;
  this->fibersWithNans = fibersWithNans ;
  this->tracksScalars = tracksScalars ;
  this->tracksProperties = tracksProperties ;
  this->curves_count = curves_count ;
  this->isBundles = isBundles ;
  this->isTrk = isTrk ;
  this->isTck = isTck ;

}


BundlesData::BundlesData( const BundlesData& bundlesData )
{

  this->matrixTracks = bundlesData.matrixTracks ;

  this->pointsPerTrack = bundlesData.pointsPerTrack ;

  this->fibersWithNans = bundlesData.fibersWithNans ;

  this->tracksScalars = bundlesData.tracksScalars ;

  this->tracksProperties = bundlesData.tracksProperties ;

  this->curves_count = bundlesData.curves_count ;

  this->isBundles = bundlesData.isBundles ;

  this->isTrk = bundlesData.isTrk ;

  this->isTck = bundlesData.isTck ;


}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BundlesData::~BundlesData() {}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Operators ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool BundlesData::operator == ( const BundlesData& bundleData ) const
{

  if ( curves_count != bundleData.curves_count )
  {


    return( false ) ;

  }

  if ( pointsPerTrack.size() == bundleData.pointsPerTrack.size() )
  {

    for ( int i = 0 ; i < pointsPerTrack.size() ; i++ )
    {

      if ( pointsPerTrack[ i ] != bundleData.pointsPerTrack[ i ] )
      {

        return( false ) ;

      }

    }

  }
  else
  {

    std::cout << pointsPerTrack.size() << std::endl ;
    std::cout << bundleData.pointsPerTrack.size() << std::endl ;

    return( false ) ;

  }

  if ( matrixTracks.size() == bundleData.matrixTracks.size() )
  {

    for ( int i = 0 ; i < matrixTracks.size() ; i++ )
    {

      if ( matrixTracks[ i ] != bundleData.matrixTracks[ i ] )
      {

        return( false ) ;

      }

    }

  }
  else
  {


    return( false ) ;

  }


  if ( tracksScalars.size() == bundleData.tracksScalars.size() )
  {

    for ( int i = 0 ; i < tracksScalars.size() ; i++ )
    {

      if ( tracksScalars[ i ] == bundleData.tracksScalars[ i ] )
      {

        for ( int j = 0 ; j < tracksScalars[ i ].size() ; j++ )
        {

          if ( tracksScalars[ i ][ j ] != bundleData.tracksScalars[ i ][ j ] )
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

  }
  else
  {

    return( false ) ;

  }

  if ( tracksProperties.size() == bundleData.tracksProperties.size() )
  {

    for ( int i = 0 ; i < tracksProperties.size() ; i++ )
    {

      if ( tracksProperties[ i ] == bundleData.tracksProperties[ i ] )
      {

        for ( int j = 0 ; j < tracksProperties[ i ].size() ; j++ )
        {

          if ( tracksProperties[ i ][ j ] !=
                                        bundleData.tracksProperties[ i ][ j ] )
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

  }
  else
  {

    return( false ) ;

  }

  return( true ) ;


}


////////////////////////////////////////////////////////////////////////////////

float BundlesData::operator[]( int64_t index ) const
{

  return ( this->matrixTracks[ index ] ) ;

}
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// For details on .tck http://trackvis.org/docs/?subsect=fileformat
// For details on .trk http://trackvis.org/docs/?subsect=fileformat
void BundlesData::read( const char* bundlesFilename )
{

  std::string bundlesDataFilenameStr = bundlesFilename ;
  std::string bundlesInfoFilenameStr = bundlesFilename ;

  if ( endswith( bundlesDataFilenameStr, ".bundles" ) )
  {

    isBundles = true ;
    bundlesDataFilenameStr = replaceExtension( bundlesDataFilenameStr,
                                                              ".bundlesdata" ) ;

  }
  else if ( endswith( bundlesDataFilenameStr, ".bundlesdata" ) )
  {

    isBundles = true ;
    bundlesInfoFilenameStr = replaceExtension( bundlesInfoFilenameStr,
                                                                  ".bundles" ) ;

  }
  else if ( endswith( bundlesDataFilenameStr, ".trk" ) )
  {

    isTrk = true ;

  }
  else if ( endswith( bundlesDataFilenameStr, ".tck" ) )
  {

    isTck = true ;

  }
  else
  {

    std::ostringstream outMessageOss ;
    outMessageOss << "BundlesData::read : Problem reading file "
                  << bundlesDataFilenameStr << ", the Only formats supported "
                  << " are .bundles, .trk and .tck\n"
                  << "stored in the header \n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  BundlesMinf bundlesInfo( bundlesInfoFilenameStr.c_str() ) ;

  this->curves_count = bundlesInfo.curves_count ;

  if ( this->curves_count <= 0 )
  {

    std::ostringstream outMessageOss ;
    outMessageOss << "BundlesData::read : Problem reading the number of tracks "
                  << "stored in header, got " << this->curves_count << "\n" ;

    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  this->tracksScalars.resize( this->curves_count ) ;
  this->tracksProperties.resize( this->curves_count ) ;

  int64_t sizeMatrixTracks = 0 ;

  // If .bundles format
  if ( isBundles )
  {

    this->pointsPerTrack.resize( this->curves_count, 0 ) ;

    std::ifstream file ;
    file.open( bundlesDataFilenameStr.c_str(), std::ios::binary |
                                                                std::ios::in ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesData::read : Problem reading file "
                    << bundlesDataFilenameStr << std::endl ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    // Getting number of points per curve
    int64_t offsetBytes = 0 ;
    for( int track = 0 ; track < this->curves_count ; track++ )
    {

      file.read( reinterpret_cast<char*>( &( this->pointsPerTrack[ track ] ) ),
                                                             sizeof( int32_t ) ) ;

      offsetBytes += sizeof( int32_t ) +
                             this->pointsPerTrack[ track ] * 3 * sizeof( float ) ;

      file.seekg( offsetBytes, std::ios_base::beg ) ;

      sizeMatrixTracks += this->pointsPerTrack[ track ] * 3 ;

    }

    // Getting vector with fiber coordinates
    this->matrixTracks.resize( sizeMatrixTracks, 0 ) ;
    file.seekg( 0, std::ios_base::beg ) ;
    int offsetPoints = 0 ;
    for ( int track = 0 ; track < this->curves_count ; track++ )
    {

      file.seekg( sizeof( int32_t ), std::ios_base::cur ) ;
      for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
      {

        file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 0 + 3 * point +
                                           offsetPoints ] ) ), sizeof( float ) ) ;
        file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 1 + 3 * point +
                                           offsetPoints ] ) ), sizeof( float ) ) ;
        file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 2 + 3 * point +
                                           offsetPoints ] ) ), sizeof( float ) ) ;

      }

      offsetPoints += 3 * this->pointsPerTrack[ track ] ;

    }

    file.close() ;

    // Check for NaN values
    int64_t nbCurves = this->curves_count ;
    int64_t offsetTractogram = 0 ;
    int nbFibersWithNan = 0 ;
    for ( int fiber = 0 ; fiber < nbCurves ; fiber++ )
    {

      int nbPoints = this->pointsPerTrack[ fiber ] ;
      bool fiberResampled = false ;
      for ( int point = 0 ; point < nbPoints ; point++ )
      {

        for ( int i = 0 ; i < 3 ; i++ )
        {

          if ( isnan( this->matrixTracks[ offsetTractogram + 3 * point + i ] ) )
          {

            resampleFiberWithNan( this->matrixTracks ,
                                  fiber,
                                  nbPoints ) ;

            fiberResampled = true ;
            this->fibersWithNans.push_back( fiber ) ;
            nbFibersWithNan++ ;

            break ;

          }

        }

        if ( fiberResampled )
        {

          break ;

        }

      }

      offsetTractogram += 3 * nbPoints ;

    }

    // if ( nbFibersWithNan )
    // {
    //
    //   std::cout << "\nThere are " << nbFibersWithNan << " with NaN values \n" ;
    //
    // }
    // else
    // {
    //
    //   std::cout << "\nThere are no fibers with NaN values" << std::endl ;
    //
    // }

  }


  // If .trk format
  if ( isTrk )
  {

    this->pointsPerTrack.resize( this->curves_count, 0 ) ;

    std::ifstream file ;
    file.open( bundlesDataFilenameStr, std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesData::read : Problem reading file "
                    << bundlesDataFilenameStr << std::endl ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    int64_t sizeMatrixTracks = 0 ;

    // Getting number of points per curve
    int64_t offsetBytes = 1000 ;
    file.seekg( offsetBytes, std::ios_base::beg ) ;
    for( int track = 0 ; track < this->curves_count ; track++ )
    {

      file.read( reinterpret_cast<char*>( &( this->pointsPerTrack[ track ] ) ),
                                                           sizeof( int32_t ) ) ;

      offsetBytes += sizeof( int32_t ) + this->pointsPerTrack[ track ] *
                                 ( 3 + bundlesInfo.n_scalars ) * sizeof( float )
                                  + bundlesInfo.n_properties * sizeof( float ) ;

      file.seekg( offsetBytes, std::ios_base::beg ) ;

      sizeMatrixTracks += this->pointsPerTrack[ track ] * 3 ;

    }

    // Getting vector with fiber coordinates
    this->tracksScalars.resize( this->curves_count ) ;
    this->tracksProperties.resize( this->curves_count ) ;
    this->matrixTracks.resize( sizeMatrixTracks ) ;
    file.seekg( 1000,  std::ios_base::beg ) ;
    int64_t offsetPoints = 0 ;
    for ( int track = 0 ; track < this->curves_count ; track++ )
    {

      std::vector<float>& trackScalars =  this->tracksScalars[ track ] ;
      std::vector<float>& trackProperties = this->tracksProperties[ track ] ;

      trackScalars.resize( bundlesInfo.n_scalars *
                                               this->pointsPerTrack[ track ] ) ;
      trackProperties.resize( bundlesInfo.n_properties ) ;

      file.seekg( sizeof( int32_t ), std::ios_base::cur ) ;

      for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
      {

        file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 0 +
                             3 * point + offsetPoints ] ) ), sizeof( float ) ) ;
        file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 1 +
                             3 * point + offsetPoints ] ) ), sizeof( float ) ) ;
        file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 2 +
                             3 * point + offsetPoints ] ) ), sizeof( float ) ) ;


        for ( int scalar = 0 ; scalar < bundlesInfo.n_scalars ; scalar++ )
        {

          file.read( reinterpret_cast<char*>( &( trackScalars[ scalar +
                        bundlesInfo.n_scalars * point ] ) ), sizeof( float ) ) ;

        }

      }

      for ( int property = 0 ; property < bundlesInfo.n_properties ;
                                                                    property++ )
      {

        file.read( reinterpret_cast<char*>( &( trackProperties[ property ] ) ),
                                                             sizeof( float ) ) ;

      }


      offsetPoints += 3 * this->pointsPerTrack[ track ] ;


    }

    file.close() ;

  }


  // If .tck format
  if ( isTck )
  {

    std::ifstream file ;
    file.open( bundlesDataFilenameStr, std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesData::read : Problem reading file "
                    << bundlesDataFilenameStr << std::endl ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    file.seekg( bundlesInfo.tckOffsetBinary, std::ios_base::beg ) ;

    bool notEndFile = true ;
    uint64_t currentFiber = 0 ;
    int nbPoints = 0 ;
    while ( notEndFile )
    {

      float track[ 3 ] ;
      file.read( reinterpret_cast<char*>( &( track ) ), 3 *
                                                 bundlesInfo.tckSizeDataType ) ;

      if ( std::isinf( track[ 0 ] ) && std::isinf( track[ 1 ] ) &&
                                                      std::isinf( track[ 2 ] ) )
      {

        notEndFile = false ;

      }
      else if ( std::isnan( track[ 0 ] ) && std::isnan( track[ 1 ] ) &&
                                                      std::isnan( track[ 2 ] ) )
      {

        this->pointsPerTrack.push_back( nbPoints ) ;
        nbPoints = 0 ;

        currentFiber++ ;

      }
      else
      {

        this->matrixTracks.push_back( track[ 0 ] ) ;
        this->matrixTracks.push_back( track[ 1 ] ) ;
        this->matrixTracks.push_back( track[ 2 ] ) ;

        nbPoints++ ;

      }

    }

    file.close() ;

    if ( this->curves_count != this->pointsPerTrack.size() )
    {

      std::stringstream outMessageOss ;

      outMessageOss << "BundlesData::read : Problem reading "
                    << bundlesDataFilenameStr
                    << "   ->   the number of fibers stored in the header is "
                    << "different from the number of fibers in the binary data,"
                    << " got " << this->curves_count << " from header and "
                    << this->pointsPerTrack.size() << " from binary "
                    << "( size pointsPerTrack vector ) \n" ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;


    }

  }

}

////////////////////////////////////////////////////////////////////////////////
// For details on .tck http://trackvis.org/docs/?subsect=fileformat
// For details on .trk http://trackvis.org/docs/?subsect=fileformat
void BundlesData::write( const char* bundlesFilename,
                         const BundlesMinf& bundleInfo ) const
{

  bool outIsBundles = false ;
  bool outIsTrk = false ;
  bool outIsTck = false ;

  std::string bundlesFilenameStr = bundlesFilename ;
  std::string bundlesFilenameInfoStr = bundlesFilename ;

  if ( endswith( bundlesFilenameStr, ".bundles" ) )
  {

    outIsBundles = true ;
    bundlesFilenameStr = replaceExtension( bundlesFilenameStr,
                                                              ".bundlesdata" ) ;

  }
  if ( endswith( bundlesFilenameStr, ".bundlesdata" ) )
  {

    outIsBundles = true ;
    bundlesFilenameInfoStr = replaceExtension( bundlesFilenameInfoStr,
                                                                  ".bundles" ) ;

  }
  else if ( endswith( bundlesFilenameStr, ".trk" ) )
  {

    outIsTrk = true ;
    bundlesFilenameInfoStr = replaceExtension( bundlesFilenameInfoStr,
                                                                     ".minf" ) ;

  }
  else if ( endswith( bundlesFilenameStr, ".tck" ) )
  {

    outIsTck = true ;
    bundlesFilenameInfoStr = replaceExtension( bundlesFilenameInfoStr,
                                                                     ".minf" ) ;

  }
  else
  {

    const char* outMessage = "BundlesData::write : only formats supported are "
                                                   ".bundles, .trk and .tck\n" ;

    throw( std::invalid_argument( outMessage ) ) ;

  }


  // Checking if BundlesData and BundlesMinf in same format
  if ( bundleInfo.isBundles != this->isBundles ||
       bundleInfo.isTrk != this->isTrk ||
       bundleInfo.isTck != this->isTck )
  {

    std::cout << "WARNING : in BundlesData::write, target BundlesData and "
              << "input BundlesMinf are not in the same format" << std::endl ;

  }


  // If .bundles format
  if ( outIsBundles )
  {

    std::ofstream file ;
    file.open( bundlesFilenameStr, std::ios::binary | std::ios::out ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesData::write : Problem reading file "
                    << bundlesFilenameStr << std::endl ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    int n_curves = this->curves_count ;


    int offsetPoints = 0 ;


    for ( int track = 0 ; track < n_curves ; track++ )
    {

      int tmpNbPoints = this->pointsPerTrack[ track ] ;
      file.write( reinterpret_cast<char*>( &( tmpNbPoints ) ),
                                                           sizeof( int32_t ) ) ;

      for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
      {

        float tmpPoint = this->matrixTracks[ 0 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;
        tmpPoint = this->matrixTracks[ 1 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;
        tmpPoint = this->matrixTracks[ 2 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;

      }

      offsetPoints += 3 * this->pointsPerTrack[ track ] ;

    }

    file.close() ;

  }


  // If .trk format
  if ( outIsTrk )
  {


    if ( !bundleInfo.isTrk )
    {

      const char* outMessage = "BundlesData::write : Problem writing .trk  "
                               "BundlesMinf isTrk attribute is false \n" ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    std::ofstream file ;
    file.open( bundlesFilenameStr, std::ios::binary | std::ios::out ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesData::write : Problem reading file "
                    << bundlesFilenameStr << std::endl ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }


    for ( char c : bundleInfo.id_string )
    {

      file.write( reinterpret_cast<char*>( &( c ) ), sizeof( c ) ) ;

    }


    for( short int i : bundleInfo.size )
    {

      file.write( reinterpret_cast<char*>( &( i ) ), sizeof( i ) ) ;

    }


    for ( float f : bundleInfo.resolution )
    {

      file.write( reinterpret_cast<char*>( &( f ) ), sizeof( f ) ) ;

    }


    for ( float f : bundleInfo.origin )
    {

      file.write( reinterpret_cast<char*>( &( f ) ), sizeof( f ) ) ;

    }


    short int n_scalars = bundleInfo.n_scalars ;
    file.write( reinterpret_cast<char*>( &( n_scalars ) ),
                                                         sizeof( n_scalars ) ) ;


    for ( std::vector<char> vc : bundleInfo.scalar_name )
    {

      for ( char c : vc )
      {

        file.write( reinterpret_cast<char*>( &( c ) ), sizeof( c ) ) ;

      }

    }


    short int n_properties = bundleInfo.n_properties ;
    file.write( reinterpret_cast<char*>( &( n_properties ) ),
                                                      sizeof( n_properties ) ) ;




    for ( std::vector<char> vc : bundleInfo.property_name )
    {

      for ( char c : vc )
      {

        file.write( reinterpret_cast<char*>( &( c ) ), sizeof( c ) ) ;

      }

    }



    for ( std::vector<float> vf : bundleInfo.vox_to_ras )
    {

      for ( float f : vf )
      {

        file.write( reinterpret_cast<char*>( &( f ) ), sizeof( f ) ) ;

      }

    }



    for ( char c : bundleInfo.reserved )
    {

      file.write( reinterpret_cast<char*>( &( c ) ), sizeof( c ) ) ;

    }



    for ( char c : bundleInfo.voxel_order )
    {

      file.write( reinterpret_cast<char*>( &( c ) ), sizeof( c ) ) ;

    }



    for ( char c : bundleInfo.pad2 )
    {

      file.write( reinterpret_cast<char*>( &( c ) ), sizeof( c ) ) ;

    }



    for ( float f : bundleInfo.image_orientation_patient )
    {

      file.write( reinterpret_cast<char*>( &( f ) ), sizeof( f ) ) ;

    }



    for ( char c : bundleInfo.pad1 )
    {

      file.write( reinterpret_cast<char*>( &( c ) ), sizeof( c ) ) ;

    }


    unsigned char invert_x = bundleInfo.invert_x ;
    file.write( reinterpret_cast<char*>( &( invert_x ) ), sizeof( invert_x ) ) ;

    unsigned char invert_y = bundleInfo.invert_y ;
    file.write( reinterpret_cast<char*>( &( invert_y ) ), sizeof( invert_y ) ) ;

    unsigned char invert_z = bundleInfo.invert_z ;
    file.write( reinterpret_cast<char*>( &( invert_z) ),  sizeof( invert_z ) ) ;

    unsigned char swap_xy = bundleInfo.swap_xy ;
    file.write( reinterpret_cast<char*>( &( swap_xy ) ), sizeof( swap_xy ) ) ;

    unsigned char swap_yz = bundleInfo.swap_yz ;
    file.write( reinterpret_cast<char*>( &( swap_yz ) ), sizeof( swap_yz ) ) ;

    unsigned char swap_zx = bundleInfo.swap_zx ;
    file.write( reinterpret_cast<char*>( &( swap_zx ) ), sizeof( swap_zx ) ) ;

    int curves_count = bundleInfo.curves_count ;
    file.write( reinterpret_cast<char*>( &( curves_count ) ),
                                                      sizeof( curves_count ) ) ;

    int version = bundleInfo.version ;
    file.write( reinterpret_cast<char*>( &( version ) ), sizeof( version ) ) ;

    int hdr_size = bundleInfo.hdr_size ;
    file.write( reinterpret_cast<char*>( &( hdr_size ) ), sizeof( hdr_size ) ) ;

    int offsetPoints = 0 ;
    for( int track = 0 ; track < this->curves_count ; track++ )
    {

      int nb_points = this->pointsPerTrack[ track ] ;

      file.write( reinterpret_cast<char*>( &nb_points ), sizeof( int ) ) ;

      for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
      {

        float tmpPoint = this->matrixTracks[ 0 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;
        tmpPoint = this->matrixTracks[ 1 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;
        tmpPoint = this->matrixTracks[ 2 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;

        if ( bundleInfo.n_scalars > 0 )
        {

          // Cannot use "&" here because it binds trackScalars to const
          // value_type
          std::vector<float> trackScalars =  this->tracksScalars[ track ] ;

          for ( int scalar = 0 ; scalar < bundleInfo.n_scalars ; scalar++ )
          {

            file.write( reinterpret_cast<char*>( &( trackScalars[ scalar +
                         bundleInfo.n_scalars * point ] ) ), sizeof( float ) ) ;

          }

        }

      }

      if ( bundleInfo.n_properties > 0 )
      {

        // Cannot use "&" here because it binds trackProperties to const
        // value_type
        std::vector<float> trackProperties = this->tracksProperties[ track ] ;

        for ( int property = 0 ; property < bundleInfo.n_properties ;
                                                                    property++ )
        {

          file.write( reinterpret_cast<char*>(
                                             &( trackProperties[ property ] ) ),
                                                             sizeof( float ) ) ;

        }

      }

      offsetPoints += 3 * this->pointsPerTrack[ track ] ;

    }


    file.close() ;

  }


  // If .tck format
  if ( outIsTck )
  {

    BundlesMinf newBundleInfo( bundleInfo ) ;

    newBundleInfo.computeHeaderTck( this->curves_count ) ;
    newBundleInfo.computeTckOffsetBinary() ;

    std::ofstream file ;
    file.open( bundlesFilenameStr, std::ios::binary | std::ios::out ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "BundlesData::write : Problem reading file "
                    << bundlesFilenameStr << std::endl ;

      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }


    int nbElementsHeader = newBundleInfo.tckHeaderInfo.size() ;
    int sizeBytesHeader = 0 ;
    for ( int i = 0 ; i < nbElementsHeader ; i++ )
    {

      std::string lineToWrite = newBundleInfo.tckHeaderInfo[ i ] + "\n" ;
      sizeBytesHeader += lineToWrite.size() ;

    }

    for ( int i = 0 ; i < nbElementsHeader ; i++ )
    {

      std::string line = newBundleInfo.tckHeaderInfo[ i ] ;
      std::vector<std::string> words ;
      std::string spaceDelimiter = " " ;
      std::size_t pos = 0 ;
      while ( ( pos = line.find( spaceDelimiter ) ) != std::string::npos )
      {

        words.push_back( line.substr( 0, pos ) ) ;
        line.erase( 0, pos + spaceDelimiter.length() ) ;

      }

      words.push_back( line ) ;

      if ( words[ 0 ] == "file:" )
      {

        file << words[ 0 ] << " " << words[ 1 ] << " " << sizeBytesHeader
                                                                  << std::endl ;

      }
      else
      {

        file << newBundleInfo.tckHeaderInfo[ i ] << std::endl ;

      }

    }

    ////////////////////////////////////////////////////////////////////////////
    int nbCurves = this->curves_count ;
    float nanValue = NAN ;
    int offsetPoints = 0 ;
    for ( int track = 0 ; track < nbCurves ; track++ )
    {

      int nbPoints = this->pointsPerTrack[ track ] ;

      for ( int point = 0 ; point < nbPoints ; point++ )
      {

        float tmpPoint = this->matrixTracks[ 0 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;
        tmpPoint = this->matrixTracks[ 1 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;
        tmpPoint = this->matrixTracks[ 2 + 3 * point + offsetPoints ] ;
        file.write( reinterpret_cast<char*>( &( tmpPoint ) ),
                                                             sizeof( float ) ) ;

      }

      offsetPoints += 3 * nbPoints ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        file.write( reinterpret_cast<char*>( &nanValue ), sizeof( float ) ) ;

      }

    }

    float infinityValue = INFINITY ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      file.write( reinterpret_cast<char*>( &infinityValue ), sizeof( float ) ) ;

    }

    file.close() ;

  }

  if ( bundleInfo.haveMinf )
  {

    bundleInfo.write( bundlesFilenameInfoStr.c_str() ) ;

  }

}


////////////////////////////////////////////////////////////////////////////////
// For details on .tck http://trackvis.org/docs/?subsect=fileformat
// For details on .trk http://trackvis.org/docs/?subsect=fileformat
void BundlesData::write( const char* bundlesFilename,
                         const BundlesMinf& bundleInfo,
                         const char* referenceAnatomyPath ) const
{

  NiftiImage referenceAnatomyInfo( referenceAnatomyPath ) ;

  BundlesMinf newBundleInfo( bundleInfo ) ;
  newBundleInfo.resolution = referenceAnatomyInfo.resolution ;
  newBundleInfo.size = referenceAnatomyInfo.size ;

  std::string bundlesFilenameStr = bundlesFilename ;
  if ( endswith( bundlesFilenameStr, ".bundles" ) )
  {

    newBundleInfo.isBundles = true ;
    newBundleInfo.isTrk = false ;
    newBundleInfo.isTck = false ;

  }
  else if ( endswith( bundlesFilenameStr, ".trk" ) )
  {

    newBundleInfo.isBundles = false ;
    newBundleInfo.isTrk = true ;
    newBundleInfo.isTck = false ;

  }
  else if ( endswith( bundlesFilenameStr, ".tck" ) )
  {

    // Not necessary for .tck because the header info is filled automatically
    newBundleInfo.isBundles = false ;
    newBundleInfo.isTrk = false ;

  }
  else
  {

    const char* outMessage = "BundlesData::write : the only supported formats "
                             "are .bundles, .trk, .tck \n" ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  this->write( bundlesFilename, newBundleInfo ) ;

}


////////////////////////////////////////////////////////////////////////////////
std::string BundlesData::getFormat() const
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

    std::string outMessage = "BundlesData::getFormat : attributes isBundles, " \
                             "isTrk and isTck are all set to false, not " \
                             "supported format found \n" ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

}


////////////////////////////////////////////////////////////////////////////////
bool BundlesData::matrixTracksEquals( std::vector<float>& matrixTracks ) const
{

  if ( matrixTracks.size() == this->matrixTracks.size() )
  {

    for ( int i = 0 ; i < matrixTracks.size() ; i++ )
    {

      if ( matrixTracks[ i ] != this->matrixTracks[ i ] )
      {

        std::cout << matrixTracks[ i ] << "\t|\t" << this->matrixTracks[ i ] << std::endl ;

        return( false ) ;

      }

    }

  }
  else
  {


    return( false ) ;

  }


}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::computeLengthFiber(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int nbPoints ) const
{

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  float lengthFiber = 0 ;
  for ( int point = 1 ; point < nbPoints ; point++ )
  {

    float tmpDistancePoints = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

       tmpDistancePoints += pow(
               tractogramFibers[ 3 * ( point - 1 ) + i + offsetTractogram ] -
               tractogramFibers[ 3 * point + i + offsetTractogram ], 2 ) ;

    }

    tmpDistancePoints = sqrt( tmpDistancePoints ) ;

    lengthFiber += tmpDistancePoints ;

  }

  return lengthFiber ;

}

//----------------------------------------------------------------------------//
float BundlesData::computeLengthFiber( const std::vector<float>& fiber ) const
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;

  float lengthFiber = computeLengthFiber( fiber, fiberIndex, nbPoints ) ;

  return lengthFiber ;

}

//----------------------------------------------------------------------------//
float BundlesData::computeLengthFiber( int fiberIndex ) const
{

  float lengthFiber = computeLengthFiber( this->matrixTracks, fiberIndex,
                                          this->pointsPerTrack[ fiberIndex ] ) ;

  return lengthFiber ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::resampleFiberEquidistant(
                                           const std::vector<float>& inputFiber,
                                           std::vector<float>& outputFiber,
                                           int nbPointsToResample ) const
{

  int nbPointsInputFiber = ( int )( inputFiber.size() / 3 ) ;

  float lengthFiber = computeLengthFiber( inputFiber ) ;

  float stepSize = lengthFiber / ( nbPointsToResample - 1 ) ;

  outputFiber.resize( 3 * nbPointsToResample, 0 ) ;

  std::vector<float> pointToTranslate( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    int point = 0 ;
    outputFiber[ 3 * point + i ] = inputFiber[ 3 * point + i ] ;
    pointToTranslate[ i ] = inputFiber[ 3 * point + i ] ;

  }


  int pointInInputFiber = 1 ;
  for ( int point = 1 ; point < nbPointsToResample ; point++ )
  {

    float translationDistance = 0 ;
    while ( translationDistance < stepSize )
    {

      std::vector<float> translationVector( 3, 0 ) ;
      std::vector<float> point1( 3, 0 ) ;
      std::vector<float> point2( 3, 0 ) ;

      for ( int i = 0 ; i < 3 ; i++)
      {

        point1[ i ] = pointToTranslate[ i ] ;
        point2[ i ] = inputFiber[ 3 * pointInInputFiber + i ] ;

      }
      vectorFromPoints( point2, point1, translationVector ) ;
      float distanceBetweenPoints = normVector( translationVector ) ;
      normalizeVector( translationVector ) ;

      float residualDistance = stepSize - translationDistance ;

      if ( distanceBetweenPoints < residualDistance )
      {

        if ( pointInInputFiber == nbPointsInputFiber - 1 )
        {

          translatePoint( pointToTranslate, translationVector,
                                          residualDistance, pointToTranslate ) ;
          translationDistance += residualDistance ;

        }
        else
        {

          translatePoint( pointToTranslate, translationVector,
                                     distanceBetweenPoints, pointToTranslate ) ;
          translationDistance += distanceBetweenPoints ;
          pointInInputFiber++ ;

        }

      }
      else
      {

        translatePoint( pointToTranslate, translationVector,
                                       residualDistance, pointToTranslate ) ;

        translationDistance += residualDistance ;

      }

    }

    for ( int i = 0 ; i < 3 ; i++ )
    {

      outputFiber[ 3 * point + i ] = pointToTranslate[ i ] ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::computeMedialPointFiberTractogram(
                          const std::vector<float>& tractogramFibers,
                          int fiberIndex,
                          int nbPoints,
                          std::vector<float>& medialPointFiberTractogram ) const
{

  int64_t offsetTractogram = nbPoints * 3 * fiberIndex ;

  if ( nbPoints % 2 != 0 )
  {

    int64_t tmpIndexTract = 3 * ( nbPoints - 1 ) / 2 + offsetTractogram ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointFiberTractogram[ i ] = tractogramFibers[ tmpIndexTract + i ] ;

    }

  }
  else
  {

    int64_t tmpIndexTract = 3 * nbPoints / 2 + offsetTractogram ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointFiberTractogram[ i ] = ( tractogramFibers[ tmpIndexTract + i ]
                        + tractogramFibers[ tmpIndexTract + i + 3 ] ) / 2 ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesData::computeMedialPointFiberTractogram(
                          int fiberIndex,
                          std::vector<float>& medialPointFiberTractogram ) const
{

  computeMedialPointFiberTractogram( this->matrixTracks, fiberIndex,
                                     this->pointsPerTrack[ fiberIndex ],
                                                  medialPointFiberTractogram ) ;

}

// -------------------------------------------------------------------------- //
void BundlesData::computeMedialPointFiberTractogram(
                          const std::vector<float>& fiber,
                          std::vector<float>& medialPointFiberTractogram ) const
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;
  computeMedialPointFiberTractogram( fiber, fiberIndex, nbPoints,
                                                  medialPointFiberTractogram ) ;


}

////////////////////////////////////////////////////////////////////////////////
int BundlesData::computeMedialPointFiberWithDistance(
                          const std::vector<float>& tractogramFibers,
                          int fiberIndex,
                          int nbPoints,
                          std::vector<float>& medialPointFiberTractogram ) const
{

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  float lengthFiber = computeLengthFiber( tractogramFibers, fiberIndex,
                                                                    nbPoints ) ;

  float dist = 0 ;
  int point = 0 ;
  while ( dist < lengthFiber / 2 )
  {

    point++ ;

    float tmpDist = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      tmpDist += pow( tractogramFibers[ 3 * point + i + offsetTractogram ] -
             tractogramFibers[ 3 * ( point - 1 ) + i + offsetTractogram ], 2 ) ;
    }

    tmpDist = sqrt( tmpDist ) ;

    dist += tmpDist ;

  }

  std::vector<float> translationVector( 3 , 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    int index1 = 3 * ( point - 1 ) + i + offsetTractogram ;
    int index2 = 3 * point + i + offsetTractogram ;

    // Vector going from index2 to index1 because we are going to search the
    // center going from dist > length / 2 until dist < length by changing
    // point2 where dist > length ( and not point1 where dist < length / 2 )
    translationVector[ i ] = tractogramFibers[ index1 ] -
                                                    tractogramFibers[ index2 ] ;

  }

  float normTranslationVector = scalarProduct( translationVector,
                                                           translationVector ) ;
  normTranslationVector = sqrt( normTranslationVector ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translationVector[ i ] /= normTranslationVector ; // Unitary vector

  }

  normTranslationVector = scalarProduct( translationVector,
                                                           translationVector ) ;
  normTranslationVector = sqrt( normTranslationVector ) ;

  float step = 0.1 ;
  int numberOfSteps = 0 ;
  while ( dist > lengthFiber / 2 )
  {

    dist -= step ;
    numberOfSteps++ ;

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    int index = 3 * point + i + offsetTractogram ;
    medialPointFiberTractogram[ i ] = tractogramFibers[ index ] + step *
                                        numberOfSteps * translationVector[ i ] ;
                                        // It's + step and not - step because
                                        // it's a vector translation

  }

  return ( point - 1 ) ; // Return point where dist < length / 2

}

//----------------------------------------------------------------------------//
int BundlesData::computeMedialPointFiberWithDistance(
                          int fiberIndex,
                          std::vector<float>& medialPointFiberTractogram ) const
{

  int point = computeMedialPointFiberWithDistance( this->matrixTracks,
                                 fiberIndex, this->pointsPerTrack[ fiberIndex ],
                                                  medialPointFiberTractogram ) ;

  return point ; // It is point and no point - 1 because the function use to
                 // compute point returns already point - 1

}

//----------------------------------------------------------------------------//
int BundlesData::computeMedialPointFiberWithDistance(
                          const std::vector<float>& fiber,
                          std::vector<float>& medialPointFiberTractogram ) const
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;
  int point = computeMedialPointFiberWithDistance( fiber, fiberIndex, nbPoints,
                                                  medialPointFiberTractogram ) ;

  return point ; // It is point and no point - 1 because the function use to
                 // compute point returns already point - 1


}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::computeGravityCenterFiber(
                                  const std::vector<float>& tractogramFibers,
                                  int fiberIndex,
                                  int nbPoints,
                                  std::vector<float>& gravityCenterFiber ) const
{

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterFiber[ i ] = 0 ;

  }

  int sampling = nbPoints ;

  float dis = computeLengthFiber( tractogramFibers, fiberIndex, nbPoints ) ;

  float step = dis / sampling ;

  std::vector<float> pointLeft( 3, 0 ) ;
  std::vector<float> pointRight( 3 , 0 ) ;
  std::vector<float> pointStep( 3, 0 ) ;
  std::vector<float> direction( 3, 0 ) ;

  int pointLeftIndex = 0 ;
  int pointRightIndex = 1 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    pointLeft[ i ] = tractogramFibers[ 3 * pointLeftIndex + i +
                                                            offsetTractogram ] ;
    pointRight[ i ] = tractogramFibers[ 3 * pointRightIndex + i +
                                                            offsetTractogram ] ;

  }

  vectorFromPoints( pointRight, pointLeft, direction ) ;
  float segmentLength = normVector( direction ) ;
  normalizeVector( direction ) ;

  translatePoint( pointLeft, direction, step, pointStep ) ;

  if ( segmentLength < step )
  {

    float translationDistance = 0 ;

    while ( segmentLength < step )
    {

      pointLeftIndex++ ;
      pointRightIndex++ ;

      if ( pointRightIndex > ( nbPoints - 1 ) )
      {

        dis = -1 ;
        break ;

      }

      pointLeft = pointRight ;
      for ( int i = 0 ; i < 3 ; i++ )
      {

        pointRight[ i ] = tractogramFibers[ 3 * pointRightIndex + i +
                                                            offsetTractogram ] ;

      }

      vectorFromPoints( pointRight, pointLeft, direction ) ;
      translationDistance = step - segmentLength ;
      segmentLength += normVector( direction ) ;

    }

    translatePoint( pointLeft, direction, translationDistance, pointStep ) ;


  }

  int count = 0 ;


  while ( dis > 0 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      gravityCenterFiber[ i ] += pointStep[ i ] ;

    }

    pointLeft = pointStep ;

    vectorFromPoints( pointRight, pointLeft, direction ) ;
    segmentLength = normVector( direction ) ;
    normalizeVector( direction ) ;

    translatePoint( pointLeft, direction, step, pointStep ) ;

    if ( segmentLength < step )
    {

      float translationDistance = 0 ;

      while ( segmentLength < step )
      {

        pointLeftIndex++ ;
        pointRightIndex++ ;

        if ( pointRightIndex > ( nbPoints - 1 ) )
        {

          dis = -1 ;
          break ;

        }

        pointLeft = pointRight ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          pointRight[ i ] = tractogramFibers[ 3 * pointRightIndex + i +
                                                            offsetTractogram ] ;

        }

        vectorFromPoints( pointRight, pointLeft, direction ) ;
        translationDistance = step - segmentLength ;
        segmentLength += normVector( direction ) ;

      }

      translatePoint( pointLeft, direction, translationDistance, pointStep ) ;

    }

    dis -= step ;

    count++ ;

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterFiber[ i ] /= count ;

  }

  // std::cout << "Sampling = " << sampling << "   |   Count = " << count << "\n" ;

}

// -------------------------------------------------------------------------- //
void BundlesData::computeGravityCenterFiber(
                                  int fiberIndex,
                                  std::vector<float>& gravityCenterFiber ) const
{

  computeGravityCenterFiber( this->matrixTracks,
                             fiberIndex,
                             this->pointsPerTrack[ fiberIndex ],
                             gravityCenterFiber ) ;

}
// -------------------------------------------------------------------------- //
void BundlesData::computeGravityCenterFiber(
                                  const std::vector<float>& fiber,
                                  std::vector<float>& gravityCenterFiber ) const
{

  int fiberIndex = 0 ;
  int nbPoints = ( int )( fiber.size() / 3 ) ;


  computeGravityCenterFiber( fiber,
                             fiberIndex,
                             nbPoints,
                             gravityCenterFiber ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::computeNormalVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& normalVector ) const
{


  // From https://gist.github.com/ialhashim/0a2554076a6cf32831ca
  // copy coordinates to  matrix in Eigen format
	Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic > coord( 3, nbPoints ) ;
  int offsetTractogram = 3 * nbPoints * fiberIndex ;
	for ( size_t point = 0 ; point < nbPoints ; ++point )
  {

    // coord.col( point ) = tractogramFibers[ 3 * point + 0 + offsetTractogram ] ;
    coord( 0, point ) = tractogramFibers[ 3 * point + 0 + offsetTractogram ] ;
    coord( 1, point ) = tractogramFibers[ 3 * point + 1 + offsetTractogram ] ;
    coord( 2, point ) = tractogramFibers[ 3 * point + 2 + offsetTractogram ] ;

  }

	// calculate centroid
	std::vector<float> centroid{ coord.row( 0 ).mean(),
                               coord.row( 1 ).mean(),
                               coord.row( 2 ).mean() } ;

	// subtract centroid
	coord.row( 0 ).array() -= centroid[ 0 ] ;
  coord.row( 1 ).array() -= centroid[ 1 ] ;
  coord.row( 2 ).array() -= centroid[ 2 ] ;

	// we only need the left-singular matrix here
	//  http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
	auto svd = coord.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeThinV ) ;
	// Eigen::Vector3 plane_normal = svd.matrixU().rightCols<1>() ;
	std::vector<float> plane_normal{ svd.matrixU().rightCols<1>()[ 0 ],
                                   svd.matrixU().rightCols<1>()[ 1 ],
                                   svd.matrixU().rightCols<1>()[ 2 ] } ;
  normalVector[ 0 ] = plane_normal[ 0 ] ;
  normalVector[ 1 ] = plane_normal[ 1 ] ;
  normalVector[ 2 ] = plane_normal[ 2 ] ;

  normalizeVector( normalVector ) ;


}
// void BundlesData::computeNormalVectorFiberTractogram(
//                            const std::vector<float>& tractogramFibers,
//                            const std::vector<float>& medialPointFiberTractogram,
//                            int fiberIndex,
//                            int nbPoints,
//                            std::vector<float>& normalVector ) const
// {
//
//   int64_t offset = 3 * nbPoints * fiberIndex ;
//
//   std::vector<float> vector1( 3, 0 ) ;
//   std::vector<float> vector2( 3, 0 ) ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     vector1[ i ] = tractogramFibers[ 3 * 0 + i + offset ] -
//                                                medialPointFiberTractogram[ i ] ;
//     vector2[ i ] = tractogramFibers[ 3 * ( nbPoints - 1 ) + i + offset ] -
//                                                medialPointFiberTractogram[ i ] ;
//
//   }
//
//   normalizeVector( vector1 ) ;
//   normalizeVector( vector2 ) ;
//
//   crossProduct( vector1, vector2, normalVector ) ;
//
//   normalizeVector( normalVector ) ;
//
//
// }

// -------------------------------------------------------------------------- //
void BundlesData::computeNormalVectorFiberTractogram(
                                        int fiberIndex,
                                        std::vector<float>& normalVector ) const
{

  int nbPoints = this->pointsPerTrack[ fiberIndex ];
  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // this->computeMedialPointFiberTractogram( fiberIndex,
  //                                          medialPointFiberTractogram ) ;
  int medialPointIndex = computeMedialPointFiberWithDistance(
                                                  fiberIndex,
                                                  medialPointFiberTractogram ) ;

  computeNormalVectorFiberTractogram( this->matrixTracks,
                                      medialPointFiberTractogram,
                                      fiberIndex,
                                      nbPoints,
                                      normalVector ) ;

}

// -------------------------------------------------------------------------- //
void BundlesData::computeNormalVectorFiberTractogram(
                                        const std::vector<float>& fiber,
                                        std::vector<float>& normalVector ) const
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;

  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber, medialPointFiberTractogram ) ;
  int medialPointIndex = computeMedialPointFiberWithDistance(
                                                  fiber,
                                                  medialPointFiberTractogram ) ;

  int fiberIndex = 0 ;
  computeNormalVectorFiberTractogram( fiber,
                                      medialPointFiberTractogram,
                                      fiberIndex,
                                      nbPoints,
                                      normalVector ) ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::computeDirectionVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           const std::vector<float>& normalVector,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& directionVector ) const
{


  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

  std::vector<float> point1( 3, 0 ) ;
  std::vector<float> point2( 3, 0 ) ;

  std::vector<float> vector1( 3, 0 ) ;
  std::vector<float> vector2( 3, 0 ) ;
  std::vector<float> vector3( 3, 0 ) ;


  for ( int i = 0 ; i < 3 ; i++ )
  {

    point1[ i ] = tractogramFibers[ 3 * 0 + i + offsetTractogram ] ;
    point2[ i ] = tractogramFibers[ 3 * ( nbPoints - 1 ) + i +
                                                          offsetTractogram ] ;
    vector1[ i ] = point1[ i ] - medialPointFiberTractogram[ i ] ;
    vector2[ i ] = point2[ i ] - medialPointFiberTractogram[ i ] ;



    // point1[ i ] = tractogramFibers[ 3 * 0 + i + offsetTractogram ] ;
    // point2[ i ] = tractogramFibers[ 3 * ( nbPoints - 1 ) + i +
    //                                                       offsetTractogram ] ;
    // vector1[ i ] = point2[ i ] - point1[ i ] ;

  }

  normalizeVector( vector1 ) ;
  normalizeVector( vector2 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    vector3[ i ] = vector1[ i ] + vector2[ i ] ;

  }


  // crossProduct( vector1, normalVector, vector3 ) ;

  directionVector = vector3 ;

  // projectVectorToPlane( normalVector, vector3, directionVector ) ;

  normalizeVector( directionVector ) ;

}



// -------------------------------------------------------------------------- //
void BundlesData::computeDirectionVectorFiberTractogram(
                                     int fiberIndex,
                                     const std::vector<float>& normalVector,
                                     std::vector<float>& directionVector ) const
{

  int nbPoints = this->pointsPerTrack[ 0 ] ;

  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // this->computeMedialPointFiberTractogram( fiberIndex,
  //                                          medialPointFiberTractogram ) ;
  this->computeMedialPointFiberWithDistance( fiberIndex,
                                             medialPointFiberTractogram ) ;

  computeDirectionVectorFiberTractogram( this->matrixTracks,
                                         medialPointFiberTractogram,
                                         normalVector,
                                         fiberIndex,
                                         nbPoints,
                                         directionVector ) ;

}


// -------------------------------------------------------------------------- //
void BundlesData::computeDirectionVectorFiberTractogram(
                                     const std::vector<float>& fiber,
                                     const std::vector<float>& normalVector,
                                     std::vector<float>& directionVector ) const
{

  int nbPoints = ( int )( fiber.size() / 3 ) ;
  int fiberIndex = 0 ;

  std::vector<float> medialPointFiberTractogram( 3, 0 ) ;
  // this->computeMedialPointFiberTractogram( fiber,
  //                                          medialPointFiberTractogram ) ;
  int medialPointIndex = computeMedialPointFiberWithDistance( fiber,
                                             medialPointFiberTractogram ) ;

  computeDirectionVectorFiberTractogram( fiber,
                                         medialPointFiberTractogram,
                                         normalVector,
                                         fiberIndex,
                                         nbPoints,
                                         directionVector ) ;

}


////////////////////////////////////////////////////////////////////////////////
float BundlesData::computeMDADBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) const
{

  std::vector<float> fiber1( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_1, fiberIndex_1, nbPoints, fiber1 ) ;
  std::vector<float> resampledFiber1( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber1, resampledFiber1, nbPoints ) ;

  std::vector<float> fiber2( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_2, fiberIndex_2, nbPoints, fiber2 ) ;
  std::vector<float> resampledFiber2( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber2, resampledFiber2, nbPoints ) ;

  float dMDA = 0 ;

  ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
  std::vector<float> translation( 3, 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointTractFiber_1[ i ] -
                                                  medialPointTractFiber_2[ i ] ;

  }

  float tmpDirectDistance = 0 ;
  float tmpFlippedDistance = 0 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    tmpDirectDistance += pow( resampledFiber1[ i + 3 * 0 ] -
                                ( resampledFiber2[ i + 3 * 0 ] +
                                                       translation[ i ] ), 2 ) ;

    tmpFlippedDistance += pow( resampledFiber1[
                                i + 3 * 0 ] -
                                ( resampledFiber2[ i + 3 * ( nbPoints - 1 ) ]
                                                     + translation[ i ] ), 2 ) ;

  }

  bool isDirectSens = true ;
  if ( tmpFlippedDistance < tmpDirectDistance )
  {

    isDirectSens = false ;

  }

  // Computing MDA when distance
  dMDA = 0 ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    int k = point ;
    if ( !isDirectSens )
    {

      k = nbPoints - point - 1 ;

    }

    float tmpMDA = 0 ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      tmpMDA += pow( resampledFiber1[ i + 3 * point ] - ( resampledFiber2[
                                         i + 3 * k ] + translation[ i ] ), 2 ) ;

    }

    tmpMDA = sqrt( tmpMDA ) ;

    if ( tmpMDA > dMDA )
    {

      dMDA = tmpMDA ;

    }

  }

  return dMDA ;

}
// float BundlesData::computeMDADBetweenTwoFibers(
//                               const std::vector<float>& tractogramFibers_1,
//                               const std::vector<float>& tractogramFibers_2,
//                               const std::vector<float>& medialPointTractFiber_1,
//                               const std::vector<float>& medialPointTractFiber_2,
//                               int fiberIndex_1,
//                               int fiberIndex_2,
//                               int nbPoints )
// {
//
//   float dMDA = 0 ;
//
//   int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;
//   int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;
//
//   ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
//   std::vector<float> translation( 3, 0 ) ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     translation[ i ] = medialPointTractFiber_1[ i ] -
//                                                   medialPointTractFiber_2[ i ] ;
//
//   }
//
//   float tmpDirectDistance = 0 ;
//   float tmpFlippedDistance = 0 ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     tmpDirectDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * 0 +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//     tmpFlippedDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * ( nbPoints - 1 ) +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//   }
//
//   bool isDirectSens = true ;
//   if ( tmpFlippedDistance < tmpDirectDistance )
//   {
//
//     isDirectSens = false ;
//
//   }
//
//   // Computing MDA when distance
//   dMDA = 0 ;
//
//   for ( int point = 0 ; point < nbPoints ; point++ )
//   {
//
//     int k = point ;
//     if ( !isDirectSens )
//     {
//
//       k = nbPoints - point - 1 ;
//
//     }
//
//     float tmpMDA = 0 ;
//
//     for ( int i = 0 ; i < 3 ; i++ )
//     {
//
//       tmpMDA += pow( tractogramFibers_1[ i + 3 * point +
//                      offsetTractogram_1 ] - ( tractogramFibers_2[
//                      i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
//                                                                            2 ) ;
//
//     }
//
//     tmpMDA = sqrt( tmpMDA ) ;
//
//     if ( tmpMDA > dMDA )
//     {
//
//       dMDA = tmpMDA ;
//
//     }
//
//   }
//
//   return dMDA ;
//
// }

//----------------------------------------------------------------------------//
float BundlesData::computeMDADBetweenTwoFibers(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int fiberIndex_1,
                                   int fiberIndex_2,
                                   int nbPoints ) const
{

  ///////////// Searching the medial point of tractogram 1st fiber /////////////
  int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_1( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_1,
                                       fiberIndex_1,
                                       nbPoints,
                                       medialPointTractFiber_1 ) ;


  ///////////// Searching the medial point of tractogram 2nd fiber /////////////
  int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_2( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_2,
                                       fiberIndex_2,
                                       nbPoints,
                                       medialPointTractFiber_2 ) ;

  float dMDA = computeMDADBetweenTwoFibers( tractogramFibers_1,
                                            tractogramFibers_2,
                                            medialPointTractFiber_1,
                                            medialPointTractFiber_2,
                                            fiberIndex_1,
                                            fiberIndex_2,
                                            nbPoints ) ;

  return dMDA ;


}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::computeMDFBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) const
{


  std::vector<float> fiber1( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_1, fiberIndex_1, nbPoints, fiber1 ) ;
  std::vector<float> resampledFiber1( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber1, resampledFiber1, nbPoints ) ;

  std::vector<float> fiber2( 3 * nbPoints, 0 ) ;
  getFiberFromTractogram( tractogramFibers_2, fiberIndex_2, nbPoints, fiber2 ) ;
  std::vector<float> resampledFiber2( 3 * nbPoints, 0 ) ;
  resampleFiberEquidistant( fiber2, resampledFiber2, nbPoints ) ;

  float dMDF = 0 ;

  ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
  std::vector<float> translation( 3, 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointTractFiber_1[ i ] -
                                                  medialPointTractFiber_2[ i ] ;

  }

  float tmpDirectDistance = 0 ;
  float tmpFlippedDistance = 0 ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    tmpDirectDistance += pow( resampledFiber1[ i + 3 * 0 ] -
                                                ( resampledFiber2[ i + 3 * 0 ] +
                                                       translation[ i ] ), 2 ) ;

    tmpFlippedDistance += pow( resampledFiber1[ i + 3 * 0 ] -
                                 ( resampledFiber2[ i + 3 * ( nbPoints - 1 ) ] +
                                                       translation[ i ] ), 2 ) ;

  }

  bool isDirectSens = true ;
  if ( tmpFlippedDistance < tmpDirectDistance )
  {

    isDirectSens = false ;

  }

  // Computing MDA when distance
  dMDF = 0 ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    int k = point ;
    if ( !isDirectSens )
    {

      k = nbPoints - point - 1 ;

    }

    float tmpMDF = 0 ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      // tmpMDF += pow( resampledFiber1[ i + 3 * point ] - ( resampledFiber2[
      //                                 i + 3 * k ] + translation[ i ] ), 2 ) ;
      tmpMDF += pow( resampledFiber1[ i + 3 * point ] -  resampledFiber2[
                                                              i + 3 * k ], 2 ) ;

    }

    tmpMDF = sqrt( tmpMDF ) ;

    dMDF += tmpMDF ;

  }

  dMDF /= nbPoints ;

  return dMDF ;


}

// float BundlesData::computeMDFBetweenTwoFibers(
//                               const std::vector<float>& tractogramFibers_1,
//                               const std::vector<float>& tractogramFibers_2,
//                               const std::vector<float>& medialPointTractFiber_1,
//                               const std::vector<float>& medialPointTractFiber_2,
//                               int fiberIndex_1,
//                               int fiberIndex_2,
//                               int nbPoints )
// {
//
//   float dMDF = 0 ;
//
//   int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;
//   int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;
//
//   ///////// Computing translation from tract 1 fiber to tract 2 fiber //////////
//   std::vector<float> translation( 3, 0 ) ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     translation[ i ] = medialPointTractFiber_1[ i ] -
//                                                   medialPointTractFiber_2[ i ] ;
//
//   }
//
//   float tmpDirectDistance = 0 ;
//   float tmpFlippedDistance = 0 ;
//
//   for ( int i = 0 ; i < 3 ; i++ )
//   {
//
//     tmpDirectDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * 0 +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//     tmpFlippedDistance += pow( tractogramFibers_1[
//                                 i + 3 * 0 + offsetTractogram_1 ] -
//                                 ( tractogramFibers_2[ i + 3 * ( nbPoints - 1 ) +
//                                 offsetTractogram_2 ] + translation[ i ] ), 2 ) ;
//
//   }
//
//   bool isDirectSens = true ;
//   if ( tmpFlippedDistance < tmpDirectDistance )
//   {
//
//     isDirectSens = false ;
//
//   }
//
//   // Computing MDA when distance
//   dMDF = 0 ;
//
//   for ( int point = 0 ; point < nbPoints ; point++ )
//   {
//
//     int k = point ;
//     if ( !isDirectSens )
//     {
//
//       k = nbPoints - point - 1 ;
//
//     }
//
//     float tmpMDF = 0 ;
//
//     for ( int i = 0 ; i < 3 ; i++ )
//     {
//
//       // tmpMDF += pow( tractogramFibers_1[ i + 3 * point +
//       //                offsetTractogram_1 ] - ( tractogramFibers_2[
//       //                i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
//       //                                                                      2 ) ;
//       tmpMDF += pow( tractogramFibers_1[ i + 3 * point +
//                      offsetTractogram_1 ] -  tractogramFibers_2[
//                                          i + 3 * k + offsetTractogram_2 ], 2 ) ;
//
//     }
//
//     tmpMDF = sqrt( tmpMDF ) ;
//
//     dMDF += tmpMDF ;
//
//   }
//
//   dMDF /= nbPoints ;
//
//   return dMDF ;
//
//
// }

//----------------------------------------------------------------------------//
float BundlesData::computeMDFBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) const
{

  ///////////// Searching the medial point of tractogram 1st fiber /////////////
  int64_t offsetTractogram_1 = fiberIndex_1 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_1( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_1,
                                       fiberIndex_1,
                                       nbPoints,
                                       medialPointTractFiber_1 ) ;


  ///////////// Searching the medial point of tractogram 2nd fiber /////////////
  int64_t offsetTractogram_2 = fiberIndex_2 * nbPoints * 3 ;

  std::vector<float> medialPointTractFiber_2( 3, 0 ) ;
  computeMedialPointFiberWithDistance( tractogramFibers_2,
                                       fiberIndex_2,
                                       nbPoints,
                                       medialPointTractFiber_2 ) ;

  float dMDF = computeMDFBetweenTwoFibers( tractogramFibers_1,
                                           tractogramFibers_2,
                                           medialPointTractFiber_1,
                                           medialPointTractFiber_2,
                                           fiberIndex_1,
                                           fiberIndex_2,
                                           nbPoints ) ;

  return dMDF ;


}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::compareDisimilarityBundles(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int nbFibersTract_1,
                                   int nbFibersTract_2,
                                   int nbPoints1,
                                   int nbPoints2 ) const
{

  if ( nbPoints1 != nbPoints2 )
  {

    std::cout << "ERROR : in compareDisimilarityBundles( const float* "
              << "tractogramFibers_1, const float* tractogramFibers_2, int "
              << "nbFibersTract_1, int nbFibersTract_2, int nbPoints1, int "
              << "nbPoints2 ) the number of points of the fibers must be the "
              << "same " << std::endl ;
    exit( 1 ) ;

  }
  int nbPoints = nbPoints1 ;

  float disimilarity = 0 ;

  #pragma omp parallel for reduction( + : disimilarity )
  for ( int fiberIndex_1 = 0 ; fiberIndex_1 < nbFibersTract_1 ; fiberIndex_1++ )
  {

    float minDMDA = 1000 ;

    for ( int fiberIndex_2 = 0 ; fiberIndex_2 < nbFibersTract_2 ; fiberIndex_2++ )
    {

      float dMDA = computeMDADBetweenTwoFibers( tractogramFibers_1,
                                                tractogramFibers_2,
                                                fiberIndex_1,
                                                fiberIndex_2,
                                                nbPoints ) ;

      if ( minDMDA > dMDA )
      {

        minDMDA = dMDA ;

      }

    }

    disimilarity += minDMDA ;

  }

  disimilarity /= nbFibersTract_1 ;

  return disimilarity ;

}

//----------------------------------------------------------------------------//
float BundlesData::compareDisimilarityBundles(
                                     const std::vector<float>& tractogramFibers,
                                     int nbFibersTract,
                                     int nbPoints ) const
{

  int nbFibersTract_1 = this->curves_count ;
  int nbFibersTract_2 = nbFibersTract ;
  int nbPoints1 = this->pointsPerTrack[ 0 ] ;
  for ( int fiber = 1 ; fiber < nbFibersTract_2 ; fiber++ )
  {

    if ( this->pointsPerTrack[ fiber ] != nbPoints1 )
    {

      std::cout << "ERROR : in compareDisimilarityBundles( const float* "
                << "tractogramFibers, int nbFibersTract, int nbPoints ) the "
                << "number of points of all the fibers in tractogramFibers "
                << "has to be the same " << std::endl ;
      exit( 1 ) ;
    }

  }

  if ( nbPoints != nbPoints1 )
  {

    std::cout << "ERROR : in compareDisimilarityBundles( const float* "
              << "tractogramFibers, int nbFibersTract, int nbPoints ) the "
              << "number of points of the fibers must be the same "
              << std::endl ;
    exit( 1 ) ;

  }


  float disimilarity = compareDisimilarityBundles( this->matrixTracks,
                                                   tractogramFibers,
                                                   nbFibersTract_1,
                                                   nbFibersTract_2,
                                                   nbPoints1,
                                                   nbPoints ) ;

  return disimilarity ;

}


////////////////////////////////////////////////////////////////////////////////
double BundlesData::distanceBetweenBundles( const std::vector<float>& bundle1,
                                            const std::vector<float>& bundle2,
                                            int nbFibersBundle1,
                                            int nbFibersBundle2,
                                            int nbPointsTract_1,
                                            int nbPointsTract_2 ) const
{

  double disimilarity_1 = compareDisimilarityBundles( bundle1,
                                                      bundle2,
                                                      nbFibersBundle1,
                                                      nbFibersBundle2,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity_2 = compareDisimilarityBundles( bundle2,
                                                      bundle1,
                                                      nbFibersBundle2,
                                                      nbFibersBundle1,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity = 0 ;

  if ( disimilarity_1 > disimilarity_2 )
  {

    disimilarity = disimilarity_1 ;

  }
  else
  {

    disimilarity = disimilarity_2 ;

  }

  return disimilarity ;

}


//----------------------------------------------------------------------------//


double BundlesData::distanceBetweenBundles( const std::vector<float>& bundle,
                                            int nbFibersBundle,
                                            int nbPoints ) const
{

  const std::vector<float>& bundle1 = bundle ;
  int nbFibersBundle1 = nbFibersBundle ;
  int nbPointsTract_1 = nbPoints ;

  const std::vector<float>& bundle2 = this->matrixTracks ;
  int nbFibersBundle2 = this->curves_count ;
  int nbPointsTract_2 = this->pointsPerTrack[ 0 ] ;



  double disimilarity_1 = compareDisimilarityBundles( bundle1,
                                                      bundle2,
                                                      nbFibersBundle1,
                                                      nbFibersBundle2,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity_2 = compareDisimilarityBundles( bundle2,
                                                      bundle1,
                                                      nbFibersBundle2,
                                                      nbFibersBundle1,
                                                      nbPointsTract_1,
                                                      nbPointsTract_2 ) ;

  double disimilarity = 0 ;

  if ( disimilarity_1 > disimilarity_2 )
  {

    disimilarity = disimilarity_1 ;

  }
  else
  {

    disimilarity = disimilarity_2 ;

  }

  return disimilarity ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::computeNumberAdjacentFibersBundle1toBundle2(
                                const BundlesData& bundle1,
                                const BundlesData& bundle2,
                                float threshold,
                                std::vector<int>& nbAdjacentFibersBundle ) const
{

  if ( threshold <= 0 )
  {
    std::cout << "Warning : the threshold value in bundlesAdjacency method "
              << "must be greater than 0 but got " << threshold
              << "... Using default value of 10 " << std::endl ;

    threshold = 10 ;

  }

  const std::vector<float>& bundlesData1 = bundle1.matrixTracks ;
  int nbFibersBundle1 = bundle1.curves_count ;
  int nbPoints1 = bundle1.pointsPerTrack[ 0 ] ;

  const std::vector<float>& bundleData2 = bundle2.matrixTracks ;
  int nbFibersBundle2 = bundle2.curves_count ;
  int nbPoints2 = bundle2.pointsPerTrack[ 0 ] ;

  if ( nbPoints2 != nbPoints1 )
  {

    std::cout << "ERROR : Problem in compareDisimilarityBundles, not the same "
              << "number of points per track in each bundle " << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPoints1 ;

  #pragma omp parallel for
  for ( int fiberIndex1 = 0 ; fiberIndex1 < nbFibersBundle1 ; fiberIndex1++ )
  {

    int64_t offsetTractogram_1 = fiberIndex1 * nbPoints * 3 ;

    // Searching the medial point of tractogram 1 fiber
    std::vector<float> medialPointTractFiber1( 3, 0 ) ;
    bundle1.computeMedialPointFiberWithDistance( fiberIndex1,
                                                 medialPointTractFiber1 ) ;

    float minDMDA = 1000 ;
    float dMDA = 0 ;
    nbAdjacentFibersBundle[ fiberIndex1 ] = 0 ;

    for ( int fiberIndex2 = 0 ; fiberIndex2 < nbFibersBundle2 ; fiberIndex2++ )
    {


      int64_t offsetTractogram_2 = fiberIndex2 * nbPoints * 3 ;

      // Searching the medial point of tractogram 2 fiber
      std::vector<float> medialPointTractFiber2( 3, 0 ) ;
      bundle2.computeMedialPointFiberWithDistance( fiberIndex2,
                                                   medialPointTractFiber2 ) ;


      // Computing translation from tract 1 fiber to tract 2 fiber
      std::vector<float> translation( 3, 0 ) ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        translation[ i ] = medialPointTractFiber1[ i ] -
                                                   medialPointTractFiber2[ i ] ;

      }

      float tmpDirectDistance = 0 ;
      float tmpFlippedDistance = 0 ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        tmpDirectDistance += pow(
                              bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                              ( bundleData2[ i + 3 * 0 + offsetTractogram_2 ] +
                                                       translation[ i ] ), 2 ) ;
        tmpFlippedDistance += pow(
                                bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                                ( bundleData2[ i + 3 * ( nbPoints - 1 ) +
                                offsetTractogram_2 ] + translation[ i ] ), 2 ) ;

      }

      bool isDirectSens = true ;
      if ( tmpFlippedDistance < tmpDirectDistance )
      {

        isDirectSens = false ;

      }

      // Computing MDA when distance between middle points is below threshold
      dMDA = 0 ;

      for ( int point = 0 ; point < nbPoints ; point++ )
      {

        int k = point ;
        if ( !isDirectSens )
        {

          k = nbPoints - point - 1 ;

        }

        float tmpMDA = 0 ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          tmpMDA += pow( bundlesData1[ i + 3 * point +
                         offsetTractogram_1 ] - ( bundleData2[
                         i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
                                                                           2 ) ;

        }

        tmpMDA = sqrt( tmpMDA ) ;

        if ( tmpMDA > dMDA )
        {

          dMDA = tmpMDA ;

        }

      }

      if ( dMDA < threshold )
      {

        nbAdjacentFibersBundle[ fiberIndex1 ] += 1 ;

      }

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesData::computeNumberAdjacentFibersBundle1toBundle2(
                                const std::vector<float>& bundle1,
                                const std::vector<float>& bundle2,
                                int nbFibersBundle1,
                                int nbFibersBundle2,
                                int nbPoints,
                                float threshold,
                                std::vector<int>& nbAdjacentFibersBundle ) const
{

  if ( threshold <= 0 )
  {
    std::cout << "Warning : the threshold value in bundlesAdjacency method "
              << "must be greater than 0 but got " << threshold
              << "... Using default value of 10 " << std::endl ;

    threshold = 10 ;

  }

  if ( nbAdjacentFibersBundle.size() < nbFibersBundle1 )
  {

    nbAdjacentFibersBundle.resize( nbFibersBundle1, 0 ) ;

  }

  const std::vector<float>& bundlesData1 = bundle1 ;

  const std::vector<float>& bundleData2 = bundle2 ;

  #pragma omp parallel for
  for ( int fiberIndex1 = 0 ; fiberIndex1 < nbFibersBundle1 ; fiberIndex1++ )
  {

    int64_t offsetTractogram_1 = fiberIndex1 * nbPoints * 3 ;

    // Searching the medial point of tractogram 1 fiber
    std::vector<float> medialPointTractFiber1( 3, 0 ) ;
    this->computeMedialPointFiberWithDistance( bundle1, fiberIndex1, nbPoints,
                                                 medialPointTractFiber1 ) ;

    float minDMDA = 1000 ;
    float dMDA = 0 ;
    nbAdjacentFibersBundle[ fiberIndex1 ] = 0 ;

    for ( int fiberIndex2 = 0 ; fiberIndex2 < nbFibersBundle2 ; fiberIndex2++ )
    {


      int64_t offsetTractogram_2 = fiberIndex2 * nbPoints * 3 ;

      // Searching the medial point of tractogram 2 fiber
      std::vector<float> medialPointTractFiber2( 3, 0 ) ;
      this->computeMedialPointFiberWithDistance( bundle2, fiberIndex2, nbPoints,
                                                   medialPointTractFiber2 ) ;


      // Computing translation from tract 1 fiber to tract 2 fiber
      std::vector<float> translation( 3, 0 ) ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        translation[ i ] = medialPointTractFiber1[ i ] -
                                                   medialPointTractFiber2[ i ] ;

      }

      float tmpDirectDistance = 0 ;
      float tmpFlippedDistance = 0 ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        tmpDirectDistance += pow(
                              bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                              ( bundleData2[ i + 3 * 0 + offsetTractogram_2 ] +
                                                       translation[ i ] ), 2 ) ;
        tmpFlippedDistance += pow(
                                bundlesData1[ i + 3 * 0 + offsetTractogram_1 ] -
                                ( bundleData2[ i + 3 * ( nbPoints - 1 ) +
                                offsetTractogram_2 ] + translation[ i ] ), 2 ) ;

      }

      bool isDirectSens = true ;
      if ( tmpFlippedDistance < tmpDirectDistance )
      {

        isDirectSens = false ;

      }

      // Computing MDA when distance between middle points is below threshold
      dMDA = 0 ;

      for ( int point = 0 ; point < nbPoints ; point++ )
      {

        int k = point ;
        if ( !isDirectSens )
        {

          k = nbPoints - point - 1 ;

        }

        float tmpMDA = 0 ;

        for ( int i = 0 ; i < 3 ; i++ )
        {

          tmpMDA += pow( bundlesData1[ i + 3 * point +
                         offsetTractogram_1 ] - ( bundleData2[
                         i + 3 * k + offsetTractogram_2 ] + translation[ i ] ),
                                                                           2 ) ;

        }

        tmpMDA = sqrt( tmpMDA ) ;

        if ( tmpMDA > dMDA )
        {

          dMDA = tmpMDA ;

        }

      }

      if ( dMDA < threshold )
      {

        nbAdjacentFibersBundle[ fiberIndex1 ] += 1 ;

      }

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::coverageBundle1toBundle2( const BundlesData& bundle1,
                                             const BundlesData& bundle2,
                                             float threshold ) const
{

  int nbFibersBundle = bundle1.curves_count ;

  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersBundle1toBundle2( bundle1,
                                                     bundle2,
                                                     threshold,
                                                     nbAdjacentFibersBundle ) ;

  float coverage = 0 ;
  // #pragma omp parallel for reduction( + : coverage )
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersBundle[ i ] > 0 )
    {

      coverage += 1.0 ;

    }

  }

  coverage /= nbFibersBundle ;

  return coverage ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::overlapBundle1toBundle2( const BundlesData& bundle1,
                                            const BundlesData& bundle2,
                                            float threshold ) const
{

  int nbFibersBundle = bundle1.curves_count ;

  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersBundle1toBundle2( bundle1,
                                                     bundle2,
                                                     threshold,
                                                     nbAdjacentFibersBundle ) ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
  // #pragma omp parallel for reduction( + : nbAdjacentFibersRecognizedBundle, nbAdjacentFibersAtlasBundle )
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersBundle[ i ] > 0 )
    {

      nbAdjacentFibersRecognizedBundle += 1 ;
      nbAdjacentFibersAtlasBundle += nbAdjacentFibersBundle[ i ] ;

    }

  }

  if ( nbAdjacentFibersAtlasBundle == 0 )
  {

    std::cout << "Warning : there are not adjacent fibers of bundle 1 "
              << "to atlas, overlap is not defined \n" ;

    return 0 ;

  }

  // Overlap is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::bundlesAdjacency( const BundlesData& bundle1,
                                     const BundlesData& bundle2,
                                     float threshold ) const
{

  float coverageBundle1toBundle2 =
                           this->coverageBundle1toBundle2( bundle1,
                                                           bundle2,
                                                           threshold ) ;

  float coverageBundle2toBundle1 =
                           this->coverageBundle1toBundle2( bundle2,
                                                           bundle1,
                                                           threshold ) ;

  float bundlesAdjacencyMeasure = ( coverageBundle1toBundle2 +
                                    coverageBundle2toBundle1 ) / 2 ;

  return bundlesAdjacencyMeasure ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::computeAngleBetweenVectors(
                                       const std::vector<float>& vector1,
                                       const std::vector<float>& vector2 ) const
{

  float dotProduct = scalarProduct( vector1, vector2 ) ;
  float normVector1 = normVector( vector1 ) ;
  float normVector2 = normVector( vector2 ) ;

  float angle = acos( dotProduct / ( normVector1 * normVector2 ) ) *
                                                                  180 / M_PI ;

  return angle ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::computeAngleBetweenPlanes(
                                       const std::vector<float>& vector1,
                                       const std::vector<float>& vector2 ) const
{

  float vector1DotVector2 = scalarProduct( vector1, vector2 ) ;

  std::vector<float> correctedVector2( 3, 0 ) ;
  if ( vector1DotVector2 < 0 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      correctedVector2[ i ] = -vector2[ i ] ;

    }

  }
  else
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      correctedVector2[ i ] = vector2[ i ] ;

    }

  }

  float angle = computeAngleBetweenVectors( vector1, correctedVector2 ) ;

  if ( angle > 90 )
  {

    angle = 180 - angle ;

  }

  // If is angle = nan then angle = 0
  if ( angle != angle )
  {

    angle = 0 ;

  }

  return angle ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::computeAngleBetweenDirections(
                                       const std::vector<float>& vector1,
                                       const std::vector<float>& vector2 ) const
{

  float angle = computeAngleBetweenVectors( vector1, vector2 ) ;

  float vector1DotVector2 = scalarProduct( vector1, vector2 ) ;


  if ( vector1DotVector2 < 0 && angle < 90 )
  {

    angle = 180 - angle ;

  }

  // If is angle = nan then angle = 0
  if ( angle != angle )
  {

    angle = 0 ;

  }

  return angle ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::projectVectorToPlane(
                                     const std::vector<float>& normalToPlane,
                                     const std::vector<float>& inputVector,
                                     std::vector<float>& projectedVector ) const
{

  // Normalisation of normal vector
  float normalVectorNorm = sqrt( pow( normalToPlane[ 0 ], 2 ) +
                                 pow( normalToPlane[ 1 ], 2 ) +
                                 pow( normalToPlane[ 2 ], 2 ) ) ;

  std::vector<float> unitNormalVector( 3, 0 ) ;
  if ( normalVectorNorm == 1 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitNormalVector[ i ] = normalToPlane[ i ] ;

    }

  }
  else
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitNormalVector[ i ] = normalToPlane[ i ] / normalVectorNorm ;

    }

  }

  // Scalar product of normalToPlane by inputVector
  float unitNormalVectorDotInputVector = scalarProduct(
                                               unitNormalVector, inputVector ) ;

  // Projection
  for ( int i = 0 ; i < 3 ; i++ )
  {

    projectedVector[ i ] = inputVector[ i ] - unitNormalVectorDotInputVector *
                                                         unitNormalVector[ i ] ;

  }

}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::computeRotationMatrixFromVectorAndAngle(
                                      const std::vector<float>& vector,
                                      float angle, // in rad
                                      std::vector<float>& rotationMatrix ) const
{

  /////////////////////// Normalisation of normal vector ///////////////////////
  float vectorNorm = sqrt( pow( vector[ 0 ], 2 ) + pow( vector[ 1 ], 2 ) +
                                                       pow( vector[ 2 ], 2 ) ) ;

  std::vector<float> unitVector( 3, 0 ) ;
  if ( vectorNorm == 1 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitVector[ i ] = vector[ i ] ;

    }
  }
  else
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      unitVector[ i ] = vector[ i ] / vectorNorm ;

    }

  }

  float cosAngle = cos( angle ) ;
  float sinAngle = sin( angle ) ;

  float m00 = cosAngle + pow( unitVector[ 0 ], 2 ) * ( 1 - cosAngle ) ;
  float m01 = unitVector[ 0 ] * unitVector[ 1 ] * ( 1 - cosAngle ) -
                                                    unitVector[ 2 ] * sinAngle ;
  float m02 = unitVector[ 0 ] * unitVector[ 2 ] * ( 1 - cosAngle ) +
                                                    unitVector[ 1 ] * sinAngle ;
  float m10 = unitVector[ 0 ] * unitVector[ 1 ] * ( 1 - cosAngle ) +
                                                    unitVector[ 2 ] * sinAngle ;
  float m11 = cosAngle + pow( unitVector[ 1 ], 2 ) * ( 1 - cosAngle ) ;
  float m12 = unitVector[ 1 ] * unitVector[ 2 ] * ( 1 - cosAngle ) -
                                                    unitVector[ 0 ] * sinAngle ;
  float m20 = unitVector[ 0 ] * unitVector[ 2 ] * ( 1 - cosAngle ) -
                                                    unitVector[ 1 ] * sinAngle ;
  float m21 = unitVector[ 1 ] * unitVector[ 2 ] * ( 1 - cosAngle ) +
                                                    unitVector[ 0 ] * sinAngle ;
  float m22 = cosAngle + pow( unitVector[ 2 ], 2 ) * ( 1 - cosAngle ) ;



  //////////////////// Computing elements of rotationMatrix ////////////////////
  // ( 0, 0 )
  rotationMatrix[ 0 ] = m00 ;
  // ( 0, 1 )
  rotationMatrix[ 1 ] = m01 ;
  // ( 0, 2 )
  rotationMatrix[ 2 ] = m02 ;
  // ( 1, 0 )
  rotationMatrix[ 3 ] = m10 ;
  // ( 1, 1 )
  rotationMatrix[ 4 ] = m11 ;
  // ( 1, 2 )
  rotationMatrix[ 5 ] = m12 ;
  // ( 2, 0 )
  rotationMatrix[ 6 ] = m20 ;
  // ( 2, 1 )
  rotationMatrix[ 7 ] = m21 ;
  // ( 2, 2 )
  rotationMatrix[ 8 ] = m22 ;


}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::applyRotationMatrixToVector(
                                       const std::vector<float>& vector,
                                       const std::vector<float>& rotationMatrix,
                                       std::vector<float>& rotatedVector ) const
{

  for ( int i = 0 ; i < 3 ; i++ )
  {

    rotatedVector[ i ] = rotationMatrix[ 3 * i + 0 ] * vector[ 0 ] +
                         rotationMatrix[ 3 * i + 1 ] * vector[ 1 ] +
                         rotationMatrix[ 3 * i + 2 ] * vector[ 2 ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::applyRotationMatrixToFiber(
                                     const std::vector<float>& fiber,
                                     const std::vector<float>& rotationMatrix,
                                     const std::vector<float>& medialPointFiber,
                                     int nbPoints,
                                     std::vector<float>& rotatedFiber ) const
{

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = fiber[ 3 * point + i ] - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesData::applyRotationMatrixToFiber(
                                       const BundlesData& inputBundlesData,
                                       const std::vector<float>& rotationMatrix,
                                       int fiberIndex,
                                       int nbPoints,
                                       std::vector<float>& rotatedFiber ) const
{

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // inputBundlesData.computeMedialPointFiberTractogram(
  //                                               fiberIndex, medialPointFiber ) ;
  inputBundlesData.computeMedialPointFiberWithDistance(
                                                fiberIndex, medialPointFiber ) ;

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = inputBundlesData[ offsetTractogram + 3 * point + i ]
                                                       - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesData::applyRotationMatrixToFiber(
                                     const std::vector<float>& tractogramFibers,
                                     const std::vector<float>& rotationMatrix,
                                     int fiberIndex,
                                     int nbPoints,
                                     std::vector<float>& rotatedFiber ) const
{

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( tractogramFibers,
  //                                    fiberIndex,
  //                                    nbPoints,
  //                                    medialPointFiber ) ;
  computeMedialPointFiberWithDistance( tractogramFibers,
                                       fiberIndex,
                                       nbPoints,
                                       medialPointFiber ) ;

  int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = tractogramFibers[ offsetTractogram + 3 * point + i ]
                                                       - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

// -------------------------------------------------------------------------- //
void BundlesData::applyRotationMatrixToFiber(
                                       const std::vector<float>& rotationMatrix,
                                       int fiberIndex,
                                       int nbPoints,
                                       std::vector<float>& rotatedFiber ) const
{

  applyRotationMatrixToFiber( this->matrixTracks, rotationMatrix, fiberIndex,
                                                      nbPoints, rotatedFiber ) ;


}

// -------------------------------------------------------------------------- //
void BundlesData::applyRotationMatrixToFiber(
                                       const std::vector<float>& fiber,
                                       const std::vector<float>& rotationMatrix,
                                       int nbPoints,
                                       std::vector<float>& rotatedFiber ) const
{

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber, medialPointFiber ) ;
  computeMedialPointFiberWithDistance( fiber, medialPointFiber ) ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    std::vector<float> pointFiber( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      pointFiber[ i ] = fiber[ 3 * point + i ] - medialPointFiber[ i ] ;

    }

    std::vector<float> rotatedPointFiber( 3, 0 ) ;
    applyRotationMatrixToVector( pointFiber,
                                 rotationMatrix,
                                 rotatedPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      rotatedFiber[ 3 * point + i ] = rotatedPointFiber[ i ] +
                                                         medialPointFiber[ i ] ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::putFibersInSamePlane(
                               const std::vector<float>& normalVector1,
                               const std::vector<float>& normalVector2,
                               const std::vector<float>& tractogramFibers2,
                               int fiberIndex2,
                               int nbPoints,
                               std::vector<float>& fiber2ToPlane1,
                               std::vector<float>& newNormalVectorFiber2 ) const
{

  float vector1DotVector2 = scalarProduct( normalVector1, normalVector2 ) ;

  std::vector<float> correctedNormalVector2( 3, 0 ) ;
  if ( vector1DotVector2 < 0 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      correctedNormalVector2[ i ] = -normalVector2[ i ] ;

    }

  }
  else
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      correctedNormalVector2[ i ] = normalVector2[ i ] ;

    }

  }

  float angleNormaleVector = computeAngleBetweenPlanes(
                                                      normalVector1,
                                                      correctedNormalVector2 ) ;

  if ( angleNormaleVector < 1 )
  {

    int64_t offsetTractogram = 3 * nbPoints * fiberIndex2 ;

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        fiber2ToPlane1[ 3 * point + i ] = tractogramFibers2[ 3 * point + i +
                                                            offsetTractogram ] ;

      }

    }

    for ( int i = 0 ; i < 3 ; i++ )
    {

      newNormalVectorFiber2[ i ] = correctedNormalVector2[ i ] ;

    }

    return ;

  }

  angleNormaleVector *= M_PI / 180 ;

  std::vector<float> rotationAxisSamePlane( 3, 0 );
  crossProduct( normalVector1, correctedNormalVector2, rotationAxisSamePlane ) ;

  std::vector<float> rotationMatrixSamePlane( 9, 0 ) ;
  computeRotationMatrixFromVectorAndAngle( rotationAxisSamePlane,
                                           angleNormaleVector,
                                           rotationMatrixSamePlane ) ;

  applyRotationMatrixToFiber( tractogramFibers2,
                              rotationMatrixSamePlane,
                              fiberIndex2,
                              nbPoints,
                              fiber2ToPlane1 ) ;

  applyRotationMatrixToVector( correctedNormalVector2,
                               rotationMatrixSamePlane,
                               newNormalVectorFiber2 ) ;

  float residualAnglePlanes = computeAngleBetweenPlanes(
                                                       normalVector1,
                                                       newNormalVectorFiber2 ) ;

  if ( residualAnglePlanes > 1 )
  {

    std::vector<float> rotationMatrixSamePlaneTmp( 9, 0 );
    computeRotationMatrixFromVectorAndAngle( rotationAxisSamePlane,
                                             - angleNormaleVector,
                                             rotationMatrixSamePlaneTmp ) ;
    std::vector<float> fiber2ToPlane1Tmp( 3 * nbPoints, 0 ) ;
    applyRotationMatrixToFiber( tractogramFibers2,
                                rotationMatrixSamePlaneTmp,
                                fiberIndex2,
                                nbPoints,
                                fiber2ToPlane1Tmp ) ;

    std::vector<float> newNormalVectorTmp( 3, 0 ) ;

    applyRotationMatrixToVector( correctedNormalVector2,
                                 rotationMatrixSamePlaneTmp,
                                 newNormalVectorTmp ) ;

    float residualAnglePlanesTmp = computeAngleBetweenPlanes(
                                                          normalVector1,
                                                          newNormalVectorTmp ) ;

    // The comparison f != f is only true if f = nan
    if ( residualAnglePlanesTmp < residualAnglePlanes ||
                              residualAnglePlanesTmp != residualAnglePlanesTmp )
    {

      fiber2ToPlane1 = fiber2ToPlane1Tmp ;
      newNormalVectorFiber2 = newNormalVectorTmp ;
      residualAnglePlanes = residualAnglePlanesTmp ;

    }

  }


  if ( residualAnglePlanes > 5 )
  {

    std::cout << "\nERROR : in BundlesData::putFibersInSamePlane "
              << "could not put fibers in same plan, got a final angle of "
              << residualAnglePlanes << std::endl ;

    exit( 1 ) ;

  }

}

// -------------------------------------------------------------------------- //
void BundlesData::putFibersInSamePlane(
                               const std::vector<float>& normalVector1,
                               const std::vector<float>& normalVector2,
                               const std::vector<float>& fiber2,
                               int nbPoints,
                               std::vector<float>& fiber2ToPlane1,
                               std::vector<float>& newNormalVectorFiber2 ) const
{

  int fiberIndex2 = 0 ;

  putFibersInSamePlane( normalVector1,
                        normalVector2,
                        fiber2,
                        fiberIndex2,
                        nbPoints,
                        fiber2ToPlane1,
                        newNormalVectorFiber2 ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::putFiberInPlaneXY(
                                const std::vector<float>& normalVector,
                                const std::vector<float>& tractogramFibers,
                                int fiberIndex,
                                int nbPoints,
                                std::vector<float>& fiberToPlaneXY,
                                std::vector<float>& newNormalVectorFiber ) const
{

  std::vector<float> zAxisVector{ 0, 0, 1 } ;

  putFibersInSamePlane( zAxisVector,
                        normalVector,
                        tractogramFibers,
                        fiberIndex,
                        nbPoints,
                        fiberToPlaneXY,
                        newNormalVectorFiber ) ;

}

//----------------------------------------------------------------------------//
void BundlesData::putFiberInPlaneXY(
                                const std::vector<float>& normalVector,
                                const std::vector<float>& fiber,
                                int nbPoints,
                                std::vector<float>& fiberToPlaneXY,
                                std::vector<float>& newNormalVectorFiber ) const
{

  std::vector<float> zAxisVector{ 0, 0, 1 } ;

  putFibersInSamePlane( zAxisVector,
                        normalVector,
                        fiber,
                        nbPoints,
                        fiberToPlaneXY,
                        newNormalVectorFiber ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::putFibersInSameDirection(
                                  const std::vector<float>& normalVector1,
                                  const std::vector<float>& normalVector2,
                                  const std::vector<float>& directionVector1,
                                  const std::vector<float>& fiber2,
                                  int nbPoints,
                                  std::vector<float>& fiber2ToDirection1 ) const
{

  std::vector<float> medialPointFiber2( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber2, medialPointFiber2 ) ;
  computeMedialPointFiberWithDistance( fiber2, medialPointFiber2 ) ;

  float angleBetweenPlanes = computeAngleBetweenPlanes( normalVector1,
                                                        normalVector2 ) ;
  // if ( angleBetweenPlanes > 5 )
  // {
  //
  //   std::cout << "WARNING : in putFibersInSameDirection, fiber 1 and fiber 2 "
  //             << "MUST be in the same plane but got an angle between their "
  //             << "normal vectors of " << angleBetweenPlanes << ". The "
  //             << "computations will not be accurate." << std::endl ;
  //
  // }

  std::vector<float> projectedDirectionVectorFiber2( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber2,
                                         normalVector2,
                                         projectedDirectionVectorFiber2 ) ;

  float angleDirectionVectors = computeAngleBetweenDirections(
                                              directionVector1,
                                              projectedDirectionVectorFiber2 ) ;

  if ( angleDirectionVectors < 1 )
  {

    fiber2ToDirection1 = fiber2 ;
    return ;

  }

  angleDirectionVectors *= M_PI / 180 ;

  std::vector<float> rotationAxisSameDirection( 3, 0 ) ;
  crossProduct( directionVector1, projectedDirectionVectorFiber2,
                                                   rotationAxisSameDirection ) ;

  std::vector<float> rotationMatrixSameDirection( 9, 0 ) ;
  computeRotationMatrixFromVectorAndAngle( rotationAxisSameDirection,
                                           angleDirectionVectors,
                                           rotationMatrixSameDirection) ;

  applyRotationMatrixToFiber( fiber2,
                              rotationMatrixSameDirection,
                              nbPoints,
                              fiber2ToDirection1 ) ;


  std::vector<float> newDirectionVector( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber2ToDirection1, normalVector1,
                                                          newDirectionVector ) ;

  float residualAngleDirections = computeAngleBetweenDirections(
                                                          directionVector1,
                                                          newDirectionVector ) ;

  if ( residualAngleDirections > 3 )
  {

    std::vector<float> rotationMatrixSameDirectionTmp( 9, 0 ) ;
    computeRotationMatrixFromVectorAndAngle( rotationAxisSameDirection,
                                             - angleDirectionVectors,
                                             rotationMatrixSameDirectionTmp ) ;
    std::vector<float> fiber2ToDirection1Tmp( 3 * nbPoints, 0 ) ;
    applyRotationMatrixToFiber( fiber2,
                                rotationMatrixSameDirectionTmp,
                                nbPoints,
                                fiber2ToDirection1Tmp ) ;

    std::vector<float> newDirectionVectorTmp( 3, 0 ) ;
    computeDirectionVectorFiberTractogram( fiber2ToDirection1Tmp,
                                        normalVector1, newDirectionVectorTmp ) ;

    float residualAngleDirectionsTmp = computeAngleBetweenDirections(
                                                       directionVector1,
                                                       newDirectionVectorTmp ) ;

    // The comparison f != f is only true if f = nan
    if ( residualAngleDirectionsTmp < residualAngleDirections ||
                      residualAngleDirectionsTmp != residualAngleDirectionsTmp )
    {

      fiber2ToDirection1 = fiber2ToDirection1Tmp ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::registerFiber(
                               const std::vector<float>& tractogramFibers2,
                               const std::vector<float>& normalVector1,
                               const std::vector<float>& normalVector2,
                               const std::vector<float>& directionVector1,
                               const std::vector<float>& medialPointFiber1,
                               const std::vector<float>& medialPointFiber2,
                               int fiberIndexTractogram2,
                               int nbPoints,
                               std::vector<float>& fiber2Tofiber1,
                               std::vector<float>& newNormalVectorFiber2 ) const
{

  // Compute translation
  std::vector<float> translation( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointFiber1[ i ] - medialPointFiber2[ i ] ;

  }

  // Create translated fiber 2
  std::vector<float> translatedFiber2( 3 * nbPoints, 0 ) ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ;  i < 3 ; i++ )
    {

      int indexTmp = fiberIndexTractogram2 * nbPoints * 3 + 3 * point + i ;
      translatedFiber2[ 3 * point + i ] = tractogramFibers2[ indexTmp ] +
                                                              translation[ i ] ;

    }

  }

  // Computing rotation matrix to put fiber 2 in same direction as fiber 1
  std::vector<float> sameDirectionTranslatedFiber2( 3 * nbPoints, 0 ) ;
  putFibersInSameDirection( normalVector1,
                            normalVector2,
                            directionVector1,
                            translatedFiber2,
                            nbPoints,
                            sameDirectionTranslatedFiber2 ) ;


  // Computing roation matrix to put fiber 2 in same plane as fiber 1
  std::vector<float> normalVectorSameDirectionTranslated2( 3, 0 ) ;
  computeNormalVectorFiberTractogram( sameDirectionTranslatedFiber2,
                                      normalVectorSameDirectionTranslated2 ) ;
  putFibersInSamePlane( normalVector1,
                        normalVectorSameDirectionTranslated2,
                        sameDirectionTranslatedFiber2,
                        nbPoints,
                        fiber2Tofiber1,
                        newNormalVectorFiber2 ) ;

  // // Computing residual angles
  // std::vector<float> directionMovedFiber2( 3, 0 ) ;
  // computeDirectionVectorFiberTractogram( fiber2Tofiber1,
  //                                        newNormalVectorFiber2,
  //                                        directionMovedFiber2 ) ;
  //
  // float angleNormaleVectorMoved = computeAngleBetweenPlanes(
  //                                                    normalVector1,
  //                                                    newNormalVectorFiber2 ) ;
  //
  // float angleDirectionVectorsMoved = computeAngleBetweenDirections(
  //                                                     directionVector1,
  //                                                     directionMovedFiber2 ) ;
  //
  // std::cout << "Residual angles : \n   Angle between planes : "
  //           << angleNormaleVectorMoved << "\n   Angle between direction : "
  //           << angleDirectionVectorsMoved << std::endl ;


}

// -------------------------------------------------------------------------- //
void BundlesData::registerFiber(
                               const std::vector<float>& tractogramFibers1,
                               const std::vector<float>& tractogramFibers2,
                               int fiberIndexTractogram1,
                               int fiberIndexTractogram2,
                               int nbPoints,
                               std::vector<float>& fiber2Tofiber1,
                               std::vector<float>& newNormalVectorFiber2 ) const
{

  std::vector<float> medialPointFiber1( 3, 0 ) ;
  // computeMedialPointFiberTractogram( tractogramFibers1, fiberIndexTractogram1,
  //                                                nbPoints, medialPointFiber1 ) ;
  int medialPoint1Index = computeMedialPointFiberWithDistance( tractogramFibers1,
                          fiberIndexTractogram1, nbPoints, medialPointFiber1 ) ;
  std::vector<float> normalVector1( 3, 0 ) ;
  computeNormalVectorFiberTractogram( tractogramFibers1, medialPointFiber1,
                              fiberIndexTractogram1, nbPoints, normalVector1 ) ;
  // computeNormalVectorFiberTractogram( tractogramFibers1, medialPointFiber1,
  //          medialPoint1Index, fiberIndexTractogram1, nbPoints, normalVector1 ) ;
  std::vector<float> directionVector1( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( tractogramFibers1, medialPointFiber1,
                                         normalVector1, fiberIndexTractogram1,
                                                  nbPoints, directionVector1 ) ;

  std::vector<float> medialPointFiber2( 3, 0 ) ;
  // computeMedialPointFiberTractogram( tractogramFibers2, fiberIndexTractogram2,
  //                                                nbPoints, medialPointFiber2 ) ;
  int medialPoint2Index = computeMedialPointFiberWithDistance( tractogramFibers2,
                          fiberIndexTractogram2, nbPoints, medialPointFiber2 ) ;
  std::vector<float> normalVector2( 3, 0 ) ;
  computeNormalVectorFiberTractogram( tractogramFibers2, medialPointFiber2,
                              fiberIndexTractogram2, nbPoints, normalVector2 ) ;
  // computeNormalVectorFiberTractogram( tractogramFibers2, medialPointFiber2,
  //          medialPoint2Index, fiberIndexTractogram2, nbPoints, normalVector2 ) ;


  registerFiber( tractogramFibers2,
                 normalVector1,
                 normalVector2,
                 directionVector1,
                 medialPointFiber1,
                 medialPointFiber2,
                 fiberIndexTractogram2,
                 nbPoints,
                 fiber2Tofiber1,
                 newNormalVectorFiber2 ) ;


}
// -------------------------------------------------------------------------- //
void BundlesData::registerFiber(
                               const std::vector<float>& fiber2,
                               const std::vector<float>& normalVector1,
                               const std::vector<float>& normalVector2,
                               const std::vector<float>& directionVector1,
                               const std::vector<float>& medialPointFiber1,
                               const std::vector<float>& medialPointFiber2,
                               int nbPoints,
                               std::vector<float>& fiber2Tofiber1,
                               std::vector<float>& newNormalVectorFiber2 ) const
{

  int fiberIndexTractogram2 = 0 ;

  registerFiber( fiber2,
                 normalVector1,
                 normalVector2,
                 directionVector1,
                 medialPointFiber1,
                 medialPointFiber2,
                 fiberIndexTractogram2,
                 nbPoints,
                 fiber2Tofiber1,
                 newNormalVectorFiber2 ) ;

}

// -------------------------------------------------------------------------- //
void BundlesData::registerFiber(
                               const std::vector<float>& fiber1,
                               const std::vector<float>& fiber2,
                               int nbPoints,
                               std::vector<float>& fiber2Tofiber1,
                               std::vector<float>& newNormalVectorFiber2 ) const
{

  std::vector<float> medialPointFiber1( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber1, medialPointFiber1 ) ;
  computeMedialPointFiberWithDistance( fiber1, medialPointFiber1 ) ;
  std::vector<float> normalVector1( 3, 0 ) ;
  computeNormalVectorFiberTractogram( fiber1, normalVector1 ) ;
  std::vector<float> directionVector1( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber1, normalVector1,
                                                            directionVector1 ) ;

  std::vector<float> medialPointFiber2( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber2, medialPointFiber2 ) ;
  computeMedialPointFiberWithDistance( fiber2, medialPointFiber2 ) ;
  std::vector<float> normalVector2( 3, 0 ) ;
  computeNormalVectorFiberTractogram( fiber2, normalVector2 ) ;

  registerFiber( fiber2,
                 normalVector1,
                 normalVector2,
                 directionVector1,
                 medialPointFiber1,
                 medialPointFiber2,
                 nbPoints,
                 fiber2Tofiber1,
                 newNormalVectorFiber2 ) ;

}


////////////////////////////////////////////////////////////////////////////////
void BundlesData::registerFiberToPlaneXYAndDirectionX(
                               const std::vector<float>& tractogramFibers,
                               const std::vector<float>& normalVector,
                               const std::vector<float>& medialPointFiber,
                               int fiberIndexTractogram,
                               int nbPoints,
                               std::vector<float>& fiberToPlaneXYAndDirectionX,
                               std::vector<float>& newNormalVectorFiber ) const
{

  std::vector<float> xAxisVector{ 1, 0, 0 } ;
  std::vector<float> zAxisVector{ 0, 0, 1 } ;
  std::vector<float> origin{ 0, 0, 0 } ;

  registerFiber( tractogramFibers,
                 zAxisVector,
                 normalVector,
                 xAxisVector,
                 origin,
                 medialPointFiber,
                 fiberIndexTractogram,
                 nbPoints,
                 fiberToPlaneXYAndDirectionX,
                 newNormalVectorFiber ) ;

}

// -------------------------------------------------------------------------- //
void BundlesData::registerFiberToPlaneXYAndDirectionX(
                               const std::vector<float>& fiber,
                               const std::vector<float>& normalVector,
                               const std::vector<float>& medialPointFiber,
                               int nbPoints,
                               std::vector<float>& fiberToPlaneXYAndDirectionX,
                               std::vector<float>& newNormalVectorFiber ) const
{

  int fiberIndexTractogram = 0 ;

  registerFiberToPlaneXYAndDirectionX( fiber,
                                       normalVector,
                                       medialPointFiber,
                                       fiberIndexTractogram,
                                       nbPoints,
                                       fiberToPlaneXYAndDirectionX,
                                       newNormalVectorFiber ) ;


}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::getFiberFromTractogram(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int nbPoints,
                                     std::vector<float>& fiber ) const
{

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      fiber[ 3 * point + i ] = tractogramFibers[ 3 * nbPoints *
                                               fiberIndex + 3 * point + i ] ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::getFiberWithVectors(
                                    const std::vector<float>& fiber,
                                    const std::vector<float>& referenceFiber,
                                    int nbPoints,
                                    std::vector<float>& fiberWithVectors ) const
{

  std::vector<float> medialPointReferenceFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( referenceFiber,
  //                                                  medialPointReferenceFiber ) ;
  computeMedialPointFiberWithDistance( referenceFiber,
                                                   medialPointReferenceFiber ) ;

  std::vector<float> medialPointFiber( 3, 0 ) ;
  // computeMedialPointFiberTractogram( fiber, medialPointFiber ) ;
  computeMedialPointFiberWithDistance( fiber, medialPointFiber ) ;

  std::vector<float> translation( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = medialPointReferenceFiber[ i ] - medialPointFiber[ i ] ;

  }

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      fiberWithVectors[ 3 * point + i ] = fiber[  3 * point + i ] +
                                                              translation[ i ] ;

    }

  }


  std::vector<float> normalVector( 3, 0 ) ;
  computeNormalVectorFiberTractogram( fiber, normalVector ) ;

  std::vector<float> directionVector( 3, 0 ) ;
  computeDirectionVectorFiberTractogram( fiber, normalVector,
                                                             directionVector ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {


    fiberWithVectors[ 3 * nbPoints + i ] = normalVector[ i ] +
                                      medialPointFiber[ i ] + translation[ i ] ;
    fiberWithVectors[ 3 * nbPoints + 3 + i ] = medialPointFiber[ i ] +
                                                              translation[ i ] ;

    fiberWithVectors[ 3 * nbPoints + 3 * 2 + i ] = directionVector[ i ] +
                                      medialPointFiber[ i ] + translation[ i ] ;
    fiberWithVectors[ 3 * nbPoints + 3 * 3 + i ] = medialPointFiber[ i ] +
                                                              translation[ i ] ;

  }


}

// -------------------------------------------------------------------------- //
void BundlesData::getFiberWithVectors(
                                    const std::vector<float>& fiber,
                                    int nbPoints,
                                    std::vector<float>& fiberWithVectors ) const
{

  getFiberWithVectors( fiber, fiber, nbPoints, fiberWithVectors ) ;

}





////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Optimization functions ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void BundlesData::getMatrixAndVectorForLeastSquares(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int nbPoints,
                                     Eigen::MatrixXf& AMatrix,
                                     Eigen::VectorXf& BVector ) const
{

  int offsetTractogram = 3 * nbPoints * fiberIndex ;

  int nbRows = nbPoints;
  int nbColumns = 3 ;
  AMatrix.resize( nbRows, nbColumns ) ;
  BVector.resize( nbRows ) ;

  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    AMatrix( point, 0 ) = tractogramFibers[ 3 * point + 0 + offsetTractogram ] ;
    AMatrix( point, 1 ) = tractogramFibers[ 3 * point + 1 + offsetTractogram ] ;
    AMatrix( point, 2 ) = 1 ;
    BVector( point ) = tractogramFibers[ 3 * point + 2 + offsetTractogram ] ;

  }


}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Algebra functions //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void BundlesData::vectorFromPoints( const std::vector<float>& point1,
                                    const std::vector<float>& point2,
                                    std::vector<float>& vector ) const
{

  for ( int i = 0 ; i < 3 ; i++ )
  {

    vector[ i ] = point1[ i ] - point2[ i ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::scalarProduct( const std::vector<float>& vector1,
                                  const std::vector<float>& vector2 ) const
{

  float result = 0 ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    result += vector1[ i ] * vector2[ i ] ;

  }

  return result ;

}

////////////////////////////////////////////////////////////////////////////////
float BundlesData::normVector( const std::vector<float>& vector ) const
{

  float norm = scalarProduct( vector, vector ) ;
  norm = sqrt( norm ) ;

  return norm ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::normalizeVector( std::vector<float>& vector ) const
{

  float norm = normVector( vector ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    vector[ i ] /= norm ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::translatePoint(
                             const std::vector<float>& point,
                             const std::vector<float>& unitaryTranslationVector,
                             float translationDistance,
                             std::vector<float>& translatedPoint ) const
{

  for ( int i = 0 ; i < 3 ; i++ )
  {

    translatedPoint[ i ] = point[ i ] + translationDistance *
                                                 unitaryTranslationVector[ i ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::translatePoint( const std::vector<float>& point,
                                  const std::vector<float>& translationVector,
                                  std::vector<float>& translatedPoint ) const
{

  translatePoint( point, translationVector, 1.0, translatedPoint ) ;

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::crossProduct( const std::vector<float>& vector1,
                                const std::vector<float>& vector2,
                                std::vector<float>& result ) const
{

  result[ 0 ] = vector1[ 1 ] * vector2[ 2 ] - vector1[ 2 ] * vector2[ 1 ] ;
  result[ 1 ] = vector1[ 2 ] * vector2[ 0 ] - vector1[ 0 ] * vector2[ 2 ] ;
  result[ 2 ] = vector1[ 0 ] * vector2[ 1 ] - vector1[ 1 ] * vector2[ 0 ] ;

}



////////////////////////////////////////////////////////////////////////////////
int BundlesData::checkIfFiberPointCanBeResampled(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int point,
                                     int nbPoints ) const
{

  // Returns :
  //           * 0 if the fiber cannot be resample
  //           * 1 if using two consecutive points after the point
  //           * 2 if using two consecutive points before the point
  //           * 3 if using one point before and one point after the point


  if ( nbPoints < 4 )
  {

    std::cout << "ERROR in BundlesData::checkIfFiberPointCanBeResampled :"
              << " imposible to resample fiber with NaN and less than 4 points "
              << std::endl ;
    return 0 ;

  }

  int64_t offset = 3 * nbPoints * fiberIndex ;

  // Check if the fiber has the points to do the interpolation :
  // - If NaN in first point : two consecutive points after the first point
  // - If NaN is last point : two consecutive points before the last point
  // - Else : one point beore and one point after or two consecutive points
  //          before or after (execpt for second and before last points which
  //          can only be resampled with before and after point or two
  //          consecutive points after the point for the second point, and two
  //          consecutive points before the point for the before last point


  // Test if method 1 is possible
  bool isMethod1 = true ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( point < ( nbPoints - 2 ) )
    {

      if ( isnan( tractogramFibers[ offset + 3 * ( point + 1 ) + i ] ) ||
           isnan( tractogramFibers[ offset + 3 * ( point + 2 ) + i ] ) )
      {

        isMethod1 = false ;
        break ;

      }

    }
    else
    {

      isMethod1 = false ;
      break ;

    }

  }

  if ( isMethod1 )
  {

    return 1 ;

  }
  else
  {

    if ( point == 0 )
    {

      return 0 ;

    }

  }

  // Test if method 2 is possible
  bool isMethod2 = true ;
  for ( int i = 0 ; i < 3 ; i++ )
  {
    if ( point > 2 )
    {

      if ( isnan( tractogramFibers[ offset + 3 * ( point - 1 ) + i ] ) ||
           isnan( tractogramFibers[ offset + 3 * ( point - 2 ) + i ] ) )
      {

        isMethod2 = false ;
        break ;

      }

    }
    else
    {

      isMethod2 = false ;
      break ;

    }

  }

  if ( isMethod2 )
  {

    return 2 ;

  }
  else
  {

    if ( point == ( nbPoints - 1 ) )
    {

      return 0 ;

    }

  }

  // Test if method 3 is possible
  bool isMethod3 = true ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( point > 0 && point < ( nbPoints - 1 ) )
    {

      if ( isnan( tractogramFibers[ offset + 3 * ( point - 1 ) + i ] ) ||
           isnan( tractogramFibers[ offset + 3 * ( point + 1 ) + i ] ) )
      {

        isMethod3 = false ;
        break ;

      }

    }
    else
    {

      isMethod3 = false ;
      break ;

    }

  }

  if ( isMethod3 )
  {

    return 3 ;

  }
  else
  {

    return 0 ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void BundlesData::resampleFiberWithNan( std::vector<float>& tractogramFibers,
                                        int fiberIndex,
                                        int nbPoints ) const
{

  if ( nbPoints < 4 )
  {

    std::cout << "ERROR in BundlesData::resampleFiberWithNan : "
              << "imposible to resample fiber with NaN and less than 4 points "
              << std::endl ;
    return ;

  }

  int64_t offset = 3 * nbPoints * fiberIndex ;

  float lengthFiberSegment = computeLengthFiber( tractogramFibers, fiberIndex,
                                                  nbPoints ) / (nbPoints - 1 ) ;

  // Check if the fiber has the points to do the interpolation :
  // - If NaN in first point : two consecutive points after the first point
  // - If NaN is last point : two consecutive points before the last point
  // - Else : one point before and one point after or two consecutive points
  //          before or after (execpt for second and before last points which
  //          can only be resampled with before and after point or two
  //          consecutive points after the point for the second point, and two
  //          consecutive points before the point for the before last point
  std::vector<int> nanPointIndices ;
  std::vector<int> isPointOkForInterpolation( nbPoints, 0 ) ;
  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    // Checking if point is NaN
    for ( int i = 0 ; i < 3 ; i++ )
    {

      if ( isnan( tractogramFibers[ offset + 3 * point + i ] ) )
      {

        nanPointIndices.push_back( point ) ;
        break ;

      }

    }

    // Checking if points can be resampled
    isPointOkForInterpolation[ point ] = checkIfFiberPointCanBeResampled(
                                                               tractogramFibers,
                                                               fiberIndex,
                                                               point,
                                                               nbPoints ) ;

  }

  int nbNanPoints = nanPointIndices.size() ;
  bool isResamplingPossible = false ;
  for ( int point = 0 ; point < nbNanPoints ; point++ )
  {

    if ( isPointOkForInterpolation[ nanPointIndices[ point ] ] )
    {

      isResamplingPossible = true ;
      break ;

    }

  }

  if ( !isResamplingPossible )
  {

    std::cout << "ERROR in BundlesData::resampleFiberWithNan : fiber "
              << fiberIndex << " with NaN values cannot be resampled "
              << "setting values to 0 \n" ;

    for ( int _tmpPointIndex = 0 ; _tmpPointIndex < nbPoints ;
                                                              _tmpPointIndex++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        tractogramFibers[ offset + 3 * _tmpPointIndex + i ] =
                                              tractogramFibers[ offset + 3 *
                                                  ( _tmpPointIndex + 1 ) + i ] ;

      }

    }

    return ;

  }



  // Resampling
  int nbResampledPoints = 0 ;
  int nbIterations = 0 ;
  std::vector<bool> nanPointResampled( nbNanPoints ,false ) ;
  while ( nbResampledPoints < nbNanPoints && nbIterations <= nbNanPoints )
  {

    nbIterations++ ;

    for ( int nanPoint = 0 ; nanPoint < nbNanPoints ; nanPoint++ )
    {

      int nanPointIndex = nanPointIndices[ nanPoint ] ;
      if ( !nanPointResampled[ nanPoint ] )
      {

        int isInterpolationPossible = isPointOkForInterpolation[
                                                               nanPointIndex ] ;

        if (  isInterpolationPossible == 1 )
        {

          for ( int i = 0 ; i < 3 ; i++ )
          {

            tractogramFibers[ offset + 3 * nanPointIndex + i ] =
                                                  tractogramFibers[ offset + 3 *
                                                   ( nanPointIndex + 1 ) + i ] ;

            tractogramFibers[ offset + 3 * ( nanPointIndex + 1 ) + i ] =
                          ( tractogramFibers[ offset + 3 * nanPointIndex + i ] +
                            tractogramFibers[ offset + 3 * ( nanPointIndex + 2 )
                                                                   + i ] ) / 2 ;


          }

          nanPointResampled[ nanPoint ] = true ;
          nbResampledPoints++ ;


        }

        if (  isInterpolationPossible == 2 )
        {

          for ( int i = 0 ; i < 3 ; i++ )
          {

            tractogramFibers[ offset + 3 * nanPointIndex + i ] =
                                                tractogramFibers[ offset + 3 * (
                                                     nanPointIndex - 1 ) + i ] ;

            tractogramFibers[ offset + 3 * ( nanPointIndex - 1 ) + i ] =
                          ( tractogramFibers[ offset + 3 * nanPointIndex + i ] +
                            tractogramFibers[ offset + 3 * ( nanPointIndex - 2 )
                                                                   + i ] ) / 2 ;


          }

          nanPointResampled[ nanPoint ] = true ;
          nbResampledPoints++ ;

        }

        if (  isInterpolationPossible == 3 )
        {

          for ( int i = 0 ; i < 3 ; i++ )
          {

            tractogramFibers[ offset + 3 * nanPointIndex + i ] =
                  ( tractogramFibers[ offset + 3 * ( nanPointIndex + 0 ) + i ] +
                    tractogramFibers[ offset + 3 * ( nanPointIndex - 1 ) + i ] )
                                                                           / 2 ;


          }

          nanPointResampled[ nanPoint ] = true ;
          nbResampledPoints++ ;

        }

      }

    }

    // Checking for new points that can be interpolated
    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      isPointOkForInterpolation[ point ] = checkIfFiberPointCanBeResampled(
                                                               tractogramFibers,
                                                               fiberIndex,
                                                               point,
                                                               nbPoints ) ;

    }

  }

  if ( nbResampledPoints < nbNanPoints )
  {

    std::cout << "WARNING : Not all NaN points in fiber " << fiberIndex
              << " could be resampled " << std::endl ;

  }

}
