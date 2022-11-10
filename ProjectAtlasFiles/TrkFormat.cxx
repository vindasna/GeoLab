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

#include "TrkFormat.h"


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
TrkFormat::TrkFormat()
{

  memcpy( this->id_string, "TRACK", 6 * sizeof( char ) ) ;
  memcpy( this->voxel_order, "LPI", 4 * sizeof( char ) ) ;

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

}

TrkFormat::TrkFormat( const char* trkFilename, int verbose )
{

  trkReading( trkFilename, verbose ) ;

}

TrkFormat::TrkFormat( char* id_string,
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
                      std::vector< std::vector < float > >  tracksProperties,
                      bool isOk )
{

  memcpy( this->id_string, id_string, 6 * sizeof( char ) ) ;
  memcpy( this->dim, dim, 3 * sizeof( short int ) ) ;
  memcpy( this->voxel_size, voxel_size, 3 * sizeof( float ) ) ;
  memcpy( this->origin, origin, 3 * sizeof( float ) ) ;
  this->n_scalars = n_scalars ;
  memcpy( this->scalar_name, scalar_name, 10 * 20 * sizeof( char ) ) ;
  this->n_properties = n_properties ;
  memcpy( this->property_name, property_name, 10 * 20 * sizeof( char ) ) ;
  memcpy( this->vox_to_ras, vox_to_ras, 4 * 4 * sizeof( float ) ) ;
  memcpy( this->reserved, reserved, 444 * sizeof( char ) ) ;
  memcpy( this->voxel_order, voxel_order, 4 * sizeof( char ) ) ;
  memcpy( this->pad2, pad2, 4 * sizeof( char ) ) ;
  memcpy( this->image_orientation_patient, image_orientation_patient,
                                                         6 * sizeof( float ) ) ;
  memcpy( this->pad1, pad1, 2 * sizeof( char ) ) ;
  this->invert_x = invert_x ;
  this->invert_y = invert_y ;
  this->invert_z = invert_z ;
  this->swap_xy = swap_xy ;
  this->swap_yz = swap_yz ;
  this->swap_zx = swap_zx ;
  this->curves_count = curves_count ;
  this->version = version ;
  this->hdr_size = hdr_size ;
  this->matrixTracks = matrixTracks ;
  this->pointsPerTrack = pointsPerTrack ;
  this->tracksScalars = tracksScalars ;
  this->tracksProperties = tracksProperties ;
  this->isOk = isOk ;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
TrkFormat::~TrkFormat() {}



////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void TrkFormat::trkReading( const char* trkFilename, int verbose )
{

  std::ifstream file ;
  file.open( trkFilename, std::ios::binary | std::ios::in ) ;
  if ( file.fail() )
  {

    std::cout << "Problem reading file : " << trkFilename << std::endl ;
    exit( 1 ) ;

  }

  file.read( reinterpret_cast<char*>( &( this->id_string ) ),
                                                   sizeof( this->id_string ) ) ;
  file.read( reinterpret_cast<char*>( &( this->dim ) ), sizeof( this->dim ) ) ;
  file.read( reinterpret_cast<char*>( &( this->voxel_size ) ),
                                                  sizeof( this->voxel_size ) ) ;
  file.read( reinterpret_cast<char*>( &( this->origin ) ),
                                                      sizeof( this->origin ) ) ;
  file.read( reinterpret_cast<char*>( &( this->n_scalars ) ),
                                                   sizeof( this->n_scalars ) ) ;
  file.read( reinterpret_cast<char*>( &( this->scalar_name ) ),
                                                 sizeof( this->scalar_name ) ) ;
  file.read( reinterpret_cast<char*>( &( this->n_properties ) ),
                                                sizeof( this->n_properties ) ) ;
  file.read( reinterpret_cast<char*>( &( this->property_name ) ),
                                               sizeof( this->property_name ) ) ;
  file.read( reinterpret_cast<char*>( &( this->vox_to_ras ) ),
                                                  sizeof( this->vox_to_ras ) ) ;
  file.read( reinterpret_cast<char*>( &( this->reserved ) ),
                                                    sizeof( this->reserved ) ) ;
  file.read( reinterpret_cast<char*>( &( this->voxel_order ) ),
                                                 sizeof( this->voxel_order ) ) ;
  file.read( reinterpret_cast<char*>( &( this->pad2 ) ),
                                                        sizeof( this->pad2 ) ) ;
  file.read( reinterpret_cast<char*>( &( this->image_orientation_patient ) ),
                                   sizeof( this->image_orientation_patient ) ) ;
  file.read( reinterpret_cast<char*>( &( this->pad1 ) ),
                                                        sizeof( this->pad1 ) ) ;
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
  file.read( reinterpret_cast<char*>( &( this->curves_count ) ),
                                                sizeof( this->curves_count ) ) ;
  file.read( reinterpret_cast<char*>( &( this->version ) ),
                                                     sizeof( this->version ) ) ;
  file.read( reinterpret_cast<char*>( &( this->hdr_size ) ),
                                                    sizeof( this->hdr_size ) ) ;

  int64_t sizeMatrixTracks = 0 ;

  if ( this->curves_count == 0 )
  {

   std::cout << "Problem reading file, the number of traks was not stored "
             << "in the header " << std::endl ;

   this->isOk = false ;

  }
  else
  {

   this->pointsPerTrack.resize( this->curves_count ) ;

   // Getting number of points per curve
   int64_t offsetBytes = 1000 ;
   for( int track = 0 ; track < this->curves_count ; track++ )
   {

     if ( verbose > 1 && ( track % 1000 == 0 ||
                              ( track + 1 ) == this->curves_count ) )
     {

       printf("\rProcessing tracks [ %10d / %10d ]", track + 1 ,
                                                this->curves_count ) ;
       std::cout << "" << std::flush ;

     }

     file.read( reinterpret_cast<char*>( &( this->pointsPerTrack[ track ] ) ),
                                                           sizeof( int32_t ) ) ;

     offsetBytes += sizeof( int32_t ) + this->pointsPerTrack[ track ] *
                                 ( 3 + this->n_scalars ) * sizeof( float )
                                    + this->n_properties * sizeof( float ) ;

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

     std::vector< float >& trackScalars =  this->tracksScalars[ track ] ;
     std::vector< float >& trackProperties =
                                        this->tracksProperties[ track ] ;

     trackScalars.resize( this->n_scalars *
                                          this->pointsPerTrack[ track ] ) ;
     trackProperties.resize( this->n_properties ) ;

     file.seekg( sizeof( int32_t ), std::ios_base::cur ) ;

     for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
     {

       file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 0 +
                             3 * point + offsetPoints ] ) ), sizeof( float ) ) ;
       file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 1 +
                             3 * point + offsetPoints ] ) ), sizeof( float ) ) ;
       file.read( reinterpret_cast<char*>( &( this->matrixTracks[ 2 +
                             3 * point + offsetPoints ] ) ), sizeof( float ) ) ;


       for ( int scalar = 0 ; scalar < this->n_scalars ; scalar++ )
       {

         file.read( reinterpret_cast<char*>( &( trackScalars[ scalar +
                              this->n_scalars * point ] ) ), sizeof( float ) ) ;

       }

     }

     for ( int property = 0 ; property < this->n_properties ; property++ )
     {

       file.read( reinterpret_cast<char*>( &( trackProperties[ property ] ) ),
                                                             sizeof( float ) ) ;

     }


     offsetPoints += 3 * this->pointsPerTrack[ track ] ;


   }

   this->isOk = true ;

  }

  file.close() ;

}



void TrkFormat::trkWriting( const char* trkFilename, int verbose )
{

  std::ofstream file ;
  file.open( trkFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << trkFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }


  file.write( reinterpret_cast<char*>( &(this->id_string) ),
                                                   sizeof( this->id_string ) ) ;
  file.write( reinterpret_cast<char*>( &(this->dim) ), sizeof( this->dim ) ) ;
  file.write( reinterpret_cast<char*>( &(this->voxel_size) ),
                                                  sizeof( this->voxel_size ) ) ;
  file.write( reinterpret_cast<char*>( &(this->origin) ),
                                                      sizeof( this->origin ) ) ;
  file.write( reinterpret_cast<char*>( &(this->n_scalars) ),
                                                   sizeof( this->n_scalars ) ) ;
  file.write( reinterpret_cast<char*>( &(this->scalar_name) ),
                                                 sizeof( this->scalar_name ) ) ;
  file.write( reinterpret_cast<char*>( &(this->n_properties) ),
                                                sizeof( this->n_properties ) ) ;
  file.write( reinterpret_cast<char*>( &(this->property_name) ),
                                               sizeof( this->property_name ) ) ;
  file.write( reinterpret_cast<char*>( &(this->vox_to_ras) ),
                                                  sizeof( this->vox_to_ras ) ) ;
  file.write( reinterpret_cast<char*>( &(this->reserved) ),
                                                    sizeof( this->reserved ) ) ;
  file.write( reinterpret_cast<char*>( &(this->voxel_order) ),
                                                 sizeof( this->voxel_order ) ) ;
  file.write( reinterpret_cast<char*>( &(this->pad2) ), sizeof( this->pad2 ) ) ;
  file.write( reinterpret_cast<char*>( &(this->image_orientation_patient) ),
                                   sizeof( this->image_orientation_patient ) ) ;
  file.write( reinterpret_cast<char*>( &(this->pad1) ), sizeof( this->pad1 ) ) ;
  file.write( reinterpret_cast<char*>( &(this->invert_x) ),
                                                    sizeof( this->invert_x ) ) ;
  file.write( reinterpret_cast<char*>( &(this->invert_y) ),
                                                    sizeof( this->invert_y ) ) ;
  file.write( reinterpret_cast<char*>( &(this->invert_z) ),
                                                    sizeof( this->invert_z ) ) ;
  file.write( reinterpret_cast<char*>( &(this->swap_xy) ),
                                                     sizeof( this->swap_xy ) ) ;
  file.write( reinterpret_cast<char*>( &(this->swap_yz) ),
                                                     sizeof( this->swap_yz ) ) ;
  file.write( reinterpret_cast<char*>( &(this->swap_zx) ),
                                                     sizeof( this->swap_zx ) ) ;
  file.write( reinterpret_cast<char*>( &(this->curves_count) ),
                                                sizeof( this->curves_count ) ) ;
  file.write( reinterpret_cast<char*>( &(this->version) ),
                                                     sizeof( this->version ) ) ;
  file.write( reinterpret_cast<char*>( &(this->hdr_size) ),
                                                    sizeof( this->hdr_size ) ) ;

  int offsetPoints = 0 ;
  for( int track = 0 ; track < ( int )this->curves_count ; track++ )
  {

    int nb_points = this->pointsPerTrack[ track ] ;
    std::vector< float >& trackScalars =  this->tracksScalars[ track ] ;
    std::vector< float >& trackProperties = this->tracksProperties[ track ] ;


    file.write( reinterpret_cast<char*>( &nb_points ), sizeof( int ) ) ;

    for ( int point = 0 ; point < this->pointsPerTrack[ track ] ; point++ )
    {

      file.write( reinterpret_cast<char*>( &( this->matrixTracks[ 0 + 3 * point
                                       + offsetPoints ] ) ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>( &( this->matrixTracks[ 1 + 3 * point
                                       + offsetPoints ] ) ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>( &( this->matrixTracks[ 2 + 3 * point
                                       + offsetPoints ] ) ), sizeof( float ) ) ;

      if ( this->n_scalars > 0 )
      {

        for ( int scalar = 0 ; scalar < this->n_scalars ; scalar++ )
        {

          file.write( reinterpret_cast<char*>( &( trackScalars[ scalar +
                              this->n_scalars * point ] ) ), sizeof( float ) ) ;

        }

      }

    }

    if ( this->n_properties > 0 )
    {

      for ( int property = 0 ; property < this->n_properties ; property++ )
      {

        file.write( reinterpret_cast<char*>( &( trackProperties[ property ] ) ),
                                                             sizeof( float ) ) ;

      }

    }

    offsetPoints += 3 * this->pointsPerTrack[ track ] ;

  }

  file.close() ;

}



////////////////////////////////////////////////////////////////////////////////
void TrkFormat::printTrkHeader( )
{

  std::cout << "\nId_string : " << (this->id_string) << std::endl ;

  std::cout << "Dim : " ;
  for (int i = 0 ; i < 3 ; i++)
  {

   std::cout << (this->dim[i]) << "  ";

  }
  std::cout << "\n" ;

  std::cout << "Voxel size : " ;
  for (int i = 0 ; i < 3 ; i++)
  {

   std::cout << (this->voxel_size[ i ]) << "  ";

  }
  std::cout << "\n" ;

  std::cout << "Origin : " ;
  for (int i = 0 ; i < 3 ; i++)
  {

   std::cout << (this->origin[ i ]) << "  ";

  }
  std::cout << "\n" ;

  std::cout << "Number of scalars : " << (this->n_scalars) << std::endl;

  std::cout << "Scalar names : " << std::endl ;
  if ( this->n_scalars > 0 )
  {

    for ( int i = 0 ; i < 10 ; i++ )
    {

      std::cout << this->scalar_name[ i ] << std::endl;

    }
    std::cout << "\n" ;

  }

  std::cout << "Number of properties : " << this->n_properties << std::endl;

  std::cout << "Properties names : " << std::endl ;
  if ( this->n_properties > 0 )
  {

    for ( int i = 0 ; i < 10 ; i ++)
    {

       std::cout << this->property_name[ i ] << std::endl ;

    }
    std::cout << "\n";

  }

  std::cout << "Vox to ras : " << std::endl ;
  for ( int i = 0 ; i < 4 ; i ++)
  {

    for ( int j = 0 ; j < 4 ; j ++)
    {

      std::cout << this->vox_to_ras[ i ][ j ] << " " ;

    }

    std::cout << "\n" ;

  }

  std::cout << "Reserved : " << this->reserved << std::endl ;

  // std::cout << "Voxel order : " << this->voxel_order << std::endl ;
  std::cout << "Voxel order : " ;
  for ( int i = 0 ; i < 4 ; i++ )
  {

    std::cout << this->voxel_order[ i ] ;

  }
  std::cout << "\n" ;

  std::cout << "Padding 2 : " << this->pad2 << std::endl ;

  std::cout << "Image orientation : " << std::endl ;
  for ( int i = 0 ; i < 6 ; i++ )
  {

    std::cout << this->image_orientation_patient[ i ] << "  " ;

  }
  std::cout << "\n" ;

  std::cout << "Padding 1 : " << this->pad1 << std::endl ;

  std::cout << "Invert x : " << this->invert_x << std::endl ;

  std::cout << "Invert y : " << this->invert_y << std::endl ;

  std::cout << "Invert z : " << this->invert_z << std::endl ;

  std::cout << "Swap xy : " << this->swap_xy << std::endl ;

  std::cout << "Swap yz : " << this->swap_yz << std::endl ;

  std::cout << "Swap zx : " << this->swap_zx << std::endl ;

  std::cout << "Number of tracks stored : " << this->curves_count << std::endl ;

  std::cout << "Version : " << this->version << std::endl ;

  std::cout << "Header size : " << this->hdr_size << std::endl ;

}


////////////////////////////////////////////////////////////////////////////////
void TrkFormat::toBundles( BundlesFormat& bundlesInfo,
                           BundlesDataFormat& bundlesData )
{

  bundlesData.matrixTracks = this->matrixTracks ;
  bundlesData.pointsPerTrack = this->pointsPerTrack ;
  bundlesData.curves_count = this->curves_count ;

  bundlesInfo.binary = ( int )1 ;
  bundlesInfo.bundles = "None" ;
  bundlesInfo.byte_order = "DCBA" ;
  bundlesInfo.curves_count = this->curves_count ;
  bundlesInfo.radio = 0 ;
  bundlesInfo.length = 0 ;
  bundlesInfo.nSubjects = 0 ;
  bundlesInfo.data_file_name = "*.bundlesdata" ;
  bundlesInfo.format = "bundles_1.0" ;
  bundlesInfo.io_mode = "None" ;
  bundlesInfo.item_count = 0 ;
  bundlesInfo.label_type = "None" ;
  bundlesInfo.labels = "None" ;
  bundlesInfo.object_type = "BundleMap" ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    bundlesInfo.resolution[ i ] = this->voxel_size[ i ] ;
    bundlesInfo.size[ i ] = this->dim[ i ] ;

  }

  bundlesInfo.space_dimension = 3 ;
  bundlesInfo.maxRadius = 0 ;
  bundlesInfo.averageAngle = 0 ;
  bundlesInfo.averageDirectionAngle = 0 ;
  bundlesInfo.averageShapeAngle = 0 ;
  bundlesInfo.averageLength = 0 ;
  bundlesInfo.averageDisimilarity = 0 ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    bundlesInfo.centerBundle[ i ] = 0  ;

  }
  bundlesInfo.density = 0 ;


}
