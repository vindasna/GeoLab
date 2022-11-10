///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------- Libraries ----------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <ncurses.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <zlib.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

#include "ConvertBundleFormat.h"


////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to check if file exists //////////////////////////
////////////////////////////////////////////////////////////////////////////////

inline bool is_file( const std::string& path )
{

  struct stat buffer ;
  return( stat ( path.c_str(), &buffer ) == 0 ) ;

}


///////////////////////////////////////////////////////////////////////////////
//////////////////////////// TRK Reading function /////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// For details see http://trackvis.org/docs/?subsect=fileformat
trkFormat trkReading( const char* trkFilename )
{

   std::streampos size ;
   trkFormat trkData ;

   std::ifstream file ;
   file.open( trkFilename, std::ios::binary | std::ios::in ) ;
   if ( file.fail() )
   {

     std::cout << "Problem reading file : " << trkFilename << std::endl ;
     exit( 1 ) ;

   }

   file.read( reinterpret_cast<char*>( &( trkData.id_string ) ),
                                                 sizeof( trkData.id_string ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.dim ) ),
                                                       sizeof( trkData.dim ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.voxel_size ) ),
                                                sizeof( trkData.voxel_size ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.origin ) ),
                                                    sizeof( trkData.origin ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.n_scalars ) ),
                                                 sizeof( trkData.n_scalars ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.scalar_name ) ),
                                               sizeof( trkData.scalar_name ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.n_properties ) ),
                                              sizeof( trkData.n_properties ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.property_name ) ),
                                             sizeof( trkData.property_name ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.vox_to_ras ) ),
                                                sizeof( trkData.vox_to_ras ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.reserved ) ),
                                                  sizeof( trkData.reserved ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.voxel_order ) ),
                                               sizeof( trkData.voxel_order ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.pad2 ) ),
                                                      sizeof( trkData.pad2 ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.image_orientation_patient ) ),
                                 sizeof( trkData.image_orientation_patient ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.pad1 ) ),
                                                      sizeof( trkData.pad1 ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.invert_x ) ),
                                                  sizeof( trkData.invert_x ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.invert_y ) ),
                                                  sizeof( trkData.invert_y ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.invert_z ) ),
                                                  sizeof( trkData.invert_z ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.swap_xy ) ),
                                                   sizeof( trkData.swap_xy ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.swap_yz ) ),
                                                   sizeof( trkData.swap_yz ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.swap_zx ) ),
                                                   sizeof( trkData.swap_zx ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.n_count ) ),
                                                   sizeof( trkData.n_count ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.version ) ),
                                                   sizeof( trkData.version ) ) ;
   file.read( reinterpret_cast<char*>( &( trkData.hdr_size ) ),
                                                  sizeof( trkData.hdr_size ) ) ;



   if ( trkData.n_count == 0 )
   {

     std::cout << "Problem reading file, the number of traks was not stored "
               << "in the header " << std::endl ;

     trkData.isOk = false ;
     exit( 1 ) ;

   }
   else
   {

     trkData.pointsPerTrack.resize( trkData.n_count ) ;

     // Getting number of points per curve
     int64_t sizeMatrixTracks = 0 ;
     int64_t sizeScalars = 0 ;
     int64_t offsetBytes = 1000 ;
     for( int track = 0 ; track < trkData.n_count ; track++ )
     {

       file.read( reinterpret_cast<char*>(
                       &( trkData.pointsPerTrack[ track ] ) ), sizeof( int ) ) ;

       sizeMatrixTracks += 3 * trkData.pointsPerTrack[ track ] ;
       sizeScalars += trkData.n_scalars * trkData.pointsPerTrack[ track ] ;

       offsetBytes += sizeof( int ) + ( 3 + trkData.n_scalars ) *
                             sizeof( float ) * trkData.pointsPerTrack[ track ] +
                                        sizeof( float ) * trkData.n_properties ;

       file.seekg( offsetBytes, std::ios_base::beg ) ;

     }
     int64_t sizeProperties = trkData.n_count * trkData.n_properties ;

     // Getting vector with fiber coordinates
     trkData.tracksProperties.resize( sizeProperties, 0 ) ;
     trkData.matrixTracks.resize( sizeMatrixTracks, 0 ) ;
     trkData.tracksScalars.resize( sizeScalars, 0 ) ;
     file.seekg( 1000, std::ios_base::beg ) ;
     int offsetProperties = 0 ;
     int offsetScalars = 0 ;
     int offsetPoints = 0 ;
     for( int track = 0 ; track < trkData.n_count ; track++ )
     {

       file.seekg( sizeof( int ), std::ios_base::cur ) ;

       for ( int point = 0 ; point < trkData.pointsPerTrack[ track ] ; point++ )
       {

         // Points
         file.read( reinterpret_cast<char*>(
                                        &( trkData.matrixTracks[ 0 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
         file.read( reinterpret_cast<char*>(
                                        &( trkData.matrixTracks[ 1 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
         file.read( reinterpret_cast<char*>(
                                        &( trkData.matrixTracks[ 2 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
         // Scalars
         if ( trkData.n_scalars > 0 )
         {

           for ( int scalar = 0 ; scalar < trkData.n_scalars ; scalar++ )
           {

             file.read( reinterpret_cast<char*>(
                                        &( trkData.tracksScalars[ scalar +
                                        offsetScalars ] ) ), sizeof( float ) ) ;

           }

         }

       }
       // Properties
       if ( trkData.n_properties > 0 )
       {

         for ( int property = 0 ; property < trkData.n_properties ; property++ )
         {

           file.read( reinterpret_cast<char*>(
                                     &( trkData.tracksProperties[ property +
                                     offsetProperties ] ) ), sizeof( float ) ) ;

         }

       }

       offsetScalars += trkData.n_scalars * trkData.pointsPerTrack[ track ] ;
       offsetPoints += 3 * trkData.pointsPerTrack[ track ] ;
       offsetProperties += trkData.n_properties ;

     }

     trkData.isOk = true ;

   }

   file.close() ;

   return trkData ;

}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// TCK Reading function /////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// For details see http://trackvis.org/docs/?subsect=fileformat

tckFormat tckReading( const char* tckFilename )
{


  tckFormat tckData ;

  uint64_t curvesCount = 0 ;
  std::string dataTypeBinary = "" ;

  std::fstream infile ;
  infile.open( tckFilename, std::ios::in ) ;

  if ( infile.is_open() )
  {

    std::string line ;
    bool notFinished = true ;
    std::string spaceDelimiter = " " ;
    while ( std::getline( infile, line ) && notFinished )
    {

      tckData.headerInfo.push_back( line ) ;

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

        tckData.offsetBinary = std::stoi( words[ words.size() - 1 ] ) ;

      }

      if ( words[ 0 ] == "datatype:" )
      {

        dataTypeBinary = words[ words.size() - 1 ] ;

      }

      if ( words[ 0 ] == "count:" )
      {

        curvesCount = std::stoi( words[ words.size() - 1 ] ) ;

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

    std::cout << "ERROR : Problem reading " << tckFilename << std::endl ;
    exit( 1 ) ;

  }


  if ( dataTypeBinary == "Float32LE" || dataTypeBinary == "Float32BE" )
  {

    tckData.sizeDataType = sizeof( float ) ;

  }
  else
  {

    tckData.sizeDataType = sizeof( double ) ;

  }


   /////////////////////////////////////////////////////////////////////////////
   if ( verbose)
   {

     std::cout << "Reading tractogram : " << tckFilename << std::endl ;

   }

   std::ifstream file ;
   file.open( tckFilename, std::ios::binary | std::ios::in ) ;
   if ( file.fail() )
   {

     std::cout << "Problem reading file : " << tckFilename << std::endl ;
     exit( 1 ) ;

   }

   file.seekg( tckData.offsetBinary, std::ios_base::beg ) ;

   bool notEndFile = true ;
   uint64_t currentFiber = 0 ;
   int nbPoints = 0 ;
   while ( notEndFile )
   {

     if ( ( ( currentFiber % 1000 == 0 ) || ( currentFiber + 1 == curvesCount )
         || ( currentFiber == 0 ) ) && curvesCount > 0 && verbose > 1 )
     {

       printf( "Reading fiber : [ %ld / %ld ] \r", currentFiber, curvesCount ) ;
       std::cout << "" << std::flush ;

     }

     float track[ 3 ] ;
     file.read( reinterpret_cast<char*>( &( track ) ), 3 * tckData.sizeDataType ) ;


     if ( std::isinf( track[ 0 ] ) && std::isinf( track[ 1 ] ) &&
                                                      std::isinf( track[ 2 ] ) )
     {

       notEndFile = false ;

     }
     else if ( std::isnan( track[ 0 ] ) && std::isnan( track[ 1 ] ) &&
                                                      std::isnan( track[ 2 ] ) )
     {

       tckData.pointsPerTrack.push_back( nbPoints ) ;
       nbPoints = 0 ;

       currentFiber++ ;

     }
     else
     {

       tckData.matrixTracks.push_back( track[ 0 ] ) ;
       tckData.matrixTracks.push_back( track[ 1 ] ) ;
       tckData.matrixTracks.push_back( track[ 2 ] ) ;

       nbPoints++ ;

     }

   }

   file.close() ;


   if ( curvesCount != 0 )
   {

     std::cout << "\n" ;

     // int tmpNbPoints = tckData.pointsPerTrack[ 0 ] ;
     // for ( int _tmp = 1 ; _tmp < tckData.pointsPerTrack.size() ; _tmp++ )
     // {
     //
     //   if ( tmpNbPoints != tckData.pointsPerTrack[ _tmp ] )
     //   {
     //
     //     std::cout << "Number of points : " << tmpNbPoints << "\t|\t" << tckData.pointsPerTrack[ _tmp ] << std::endl ;
     //   }
     //
     // }

     if ( curvesCount != tckData.pointsPerTrack.size() )
     {

       std::cout << "ERROR : Problem reading " << tckFilename << "   ->  "
                 << "the number of fibers stored in the header is different "
                 << "from the number of fibers in the binary data, got "
                 << curvesCount << " from header and "
                 << tckData.pointsPerTrack.size() << " from binary "
                 << "( size pointsPerTrack vector ) \n";
       exit( 1 ) ;

     }

   }

   tckData.curves_count = curvesCount ;


   // Check the the coordinates where read correctly
   int64_t sizeMatrixTracks = 0 ;
   for ( int i = 0 ; i < curvesCount ; i++ )
   {

     sizeMatrixTracks += 3 * tckData.pointsPerTrack[ i ] ;

   }

   if ( tckData.matrixTracks.size() != sizeMatrixTracks )
   {

     std::cout << "ERROR : Problem reading " << tckFilename << "   ->  "
               << "the size of matrixTracks according to the header is "
               << "different from the size of the matrix when reading the "
               << "the binary, got " << sizeMatrixTracks << " from header and "
               << tckData.matrixTracks.size() << " from binary \n" ;
     exit( 1 ) ;

   }

   return tckData ;

}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// TCK Reading function /////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void tckWriting( const char* tckFilename,
                 tckFormat& tckData )
{

  std::ofstream outFile ;
  outFile.open( tckFilename, std::ios::binary | std::ios::out ) ;
  if ( outFile.is_open() )
  {

    int nbElementsHeader = tckData.headerInfo.size() ;
    int sizeBytesHeader = 0 ;
    for ( int i = 0 ; i < nbElementsHeader ; i++ )
    {

      std::string lineToWrite = tckData.headerInfo[ i ] + "\n" ;
      sizeBytesHeader += lineToWrite.size() ;

    }

    for ( int i = 0 ; i < nbElementsHeader ; i++ )
    {

      std::string line = tckData.headerInfo[ i ] ;
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

        outFile << words[ 0 ] << " " << words[ 1 ] << " " << sizeBytesHeader << std::endl ;

      }
      else
      {

        outFile << tckData.headerInfo[ i ] << std::endl ;

      }


    }

    ////////////////////////////////////////////////////////////////////////////
    int nbCurves = tckData.curves_count ;
    float nanValue = NAN ;
    int offsetPoints = 0 ;
    for ( int track = 0 ; track < nbCurves ; track++ )
    {

      int nbPoints = tckData.pointsPerTrack[ track ] ;

      for ( int point = 0 ; point < nbPoints ; point++ )
      {

        outFile.write( reinterpret_cast<char*>(
                                        &tckData.matrixTracks[ 0 + 3 * point +
                                            offsetPoints ] ), sizeof( float ) ) ;
        outFile.write( reinterpret_cast<char*>(
                                        &tckData.matrixTracks[ 1 + 3 * point +
                                            offsetPoints ] ), sizeof( float ) ) ;
        outFile.write( reinterpret_cast<char*>(
                                        &tckData.matrixTracks[ 2 + 3 * point +
                                            offsetPoints ] ), sizeof( float ) ) ;

      }

      offsetPoints += 3 * nbPoints ;

      for ( int i = 0 ; i < 3 ; i++ )
      {

        outFile.write( reinterpret_cast<char*>( &nanValue ), sizeof( float ) ) ;

      }

      // if ( track < nbCurves )
      // {
      //
      //   for ( int i = 0 ; i < 3 ; i++ )
      //   {
      //
      //     outFile.write( reinterpret_cast<char*>( &nanValue ), sizeof( float ) ) ;
      //
      //   }
      //
      // }

    }

    float infinityValue = INFINITY ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      outFile.write( reinterpret_cast<char*>( &infinityValue ), sizeof( float ) ) ;

    }
    ////////////////////////////////////////////////////////////////////////////

    outFile.close() ;

  }
  else
  {

    std::cout << "ERROR : Problem reading " << tckFilename << std::endl ;
    exit( 1 ) ;

  }

}



///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Read .bundles /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
BundlesFormat bundlesReading( const char* bundlesFilename )
{

  BundlesFormat bundlesInfo ;

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

      bundlesInfo.binary = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'bundles\'" )
    {

      bundlesInfo.bundles = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'byte_order\'" )
    {

      bundlesInfo.byte_order = out[ 1 ] ;


    }

    if ( out[ 0 ] == "\'curves_count\'" )
    {

      bundlesInfo.curves_count = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'radio\'" )
    {

      bundlesInfo.radio = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'length\'" )
    {

      bundlesInfo.length = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'nSubjects\'" )
    {

      bundlesInfo.nSubjects = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'data_file_name\'" )
    {

      bundlesInfo.data_file_name = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'format\'" )
    {

      bundlesInfo.format = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'io_mode\'" )
    {

      bundlesInfo.io_mode = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'item_count\'" )
    {

      bundlesInfo.item_count = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'label_type\'" )
    {

      bundlesInfo.label_type = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'labels\'" )
    {

      bundlesInfo.labels = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'object_type\'" )
    {

      bundlesInfo.object_type = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'resolutionX\'" )
    {

      bundlesInfo.resolution[ 0 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'resolutionY\'" )
    {

      bundlesInfo.resolution[ 1 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'resolutionZ\'" )
    {

      bundlesInfo.resolution[ 2 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeX\'" )
    {

      bundlesInfo.size[ 0 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeY\'" )
    {

      bundlesInfo.size[ 1 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeZ\'" )
    {

      bundlesInfo.size[ 2 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'space_dimension\'" )
    {

      bundlesInfo.space_dimension = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageRadius\'" )
    {

      bundlesInfo.averageRadius = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minRadius\'" )
    {

      bundlesInfo.minRadius = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxRadius\'" )
    {

      bundlesInfo.maxRadius = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageAngle\'" )
    {

      bundlesInfo.averageAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minAngle\'" )
    {

      bundlesInfo.minAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxAngle\'" )
    {

      bundlesInfo.maxAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageDirectionAngle\'" )
    {

      bundlesInfo.averageDirectionAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minDirectionAngle\'" )
    {

      bundlesInfo.minDirectionAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxDirectionAngle\'" )
    {

      bundlesInfo.maxDirectionAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageShapeAngle\'" )
    {

      bundlesInfo.averageShapeAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minShapeAngle\'" )
    {

      bundlesInfo.minShapeAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxShapeAngle\'" )
    {

      bundlesInfo.maxShapeAngle = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageLength\'" )
    {

      bundlesInfo.averageLength = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minLength\'" )
    {

      bundlesInfo.minLength = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxLength\'" )
    {

      bundlesInfo.maxLength = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'averageDisimilarity\'" )
    {

      bundlesInfo.averageDisimilarity = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'minDisimilarity\'" )
    {

      bundlesInfo.minDisimilarity = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'maxDisimilarity\'" )
    {

      bundlesInfo.maxDisimilarity = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'centerBundleX\'" )
    {

      bundlesInfo.centerBundle[ 0 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'centerBundleY\'" )
    {

      bundlesInfo.centerBundle[ 1 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'centerBundleZ\'" )
    {

      bundlesInfo.centerBundle[ 2 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'density\'" )
    {

      bundlesInfo.density = std::stof( out[ 1 ] ) ;

    }

  }
  minfFile.close() ;

  return bundlesInfo ;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////// Read .bundlesdata ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////
BundlesDataFormat bundlesdataReading( const char* bundlesdataFilename,
                                      const char* bundlesFilename )
{

  BundlesFormat bundlesInfo = bundlesReading( bundlesFilename ) ;

  int n_curves = bundlesInfo.curves_count ;

  if ( n_curves == 0 )
  {

    std::cout << "Problem reading file, the number of traks was not stored "
              << "in the header " << std::endl ;

    std::exit( 1 ) ;

  }

  BundlesDataFormat bundlesdata ;

   bundlesdata.curves_count = n_curves ;

  int64_t sizeMatrixTracks = 0 ;

  bundlesdata.pointsPerTrack.resize( n_curves, 0 ) ;

  std::ifstream file ;
  file.open( bundlesdataFilename, std::ios::binary | std::ios::in ) ;
  if ( file.fail() )
  {

    std::cout << "Problem reading file : " << bundlesdataFilename << std::endl ;
    exit( 1 ) ;

  }

  // Getting number of points per curve
  int64_t offsetBytes = 0 ;
  for( int track = 0 ; track < bundlesdata.curves_count ; track++ )
  {

    if ( verbose > 1 && ( track % 1000 == 0 ||
                                  ( track + 1 ) == bundlesdata.curves_count ) )
    {

      printf("\rProcessing tracks [ %10d / %10d ]",
                                        track + 1 , bundlesdata.curves_count ) ;
      std::cout << "" << std::flush ;

    }

    file.read( reinterpret_cast<char*>(
                                     &( bundlesdata.pointsPerTrack[ track ] ) ),
                                                           sizeof( int32_t ) ) ;

    offsetBytes += sizeof( int32_t ) +
                           bundlesdata.pointsPerTrack[ track ] * 3 *
                                                               sizeof( float ) ;

    file.seekg( offsetBytes, std::ios_base::beg ) ;

    sizeMatrixTracks += 3 * bundlesdata.pointsPerTrack[ track ] ;

  }

  // Getting vector with fiber coordinates
  bundlesdata.matrixTracks.resize( sizeMatrixTracks, 0 ) ;
  file.seekg( 0, std::ios_base::beg ) ;
  int offsetPoints = 0 ;
  for ( int track = 0 ; track < bundlesdata.curves_count ; track++ )
  {

    file.seekg( sizeof( int32_t ), std::ios_base::cur ) ;
    for ( int point = 0 ; point < bundlesdata.pointsPerTrack[ track ] ;
                                                                       point++ )
    {

      file.read( reinterpret_cast<char*>(
                                    &( bundlesdata.matrixTracks[ 0 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
      file.read( reinterpret_cast<char*>(
                                    &( bundlesdata.matrixTracks[ 1 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
      file.read( reinterpret_cast<char*>(
                                    &( bundlesdata.matrixTracks[ 2 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;

    }

    offsetPoints += 3 * bundlesdata.pointsPerTrack[ track ] ;

  }

  file.close() ;

  if ( verbose > 1 )
  {

    std::cout << "\n" ;

  }

  return bundlesdata ;

}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Read .bundlesmap.minf /////////////////////////////
///////////////////////////////////////////////////////////////////////////////


bundlesMapMinf bundlesMapMinfReading( const char* bundlesMapMinfFilename )
{

  bundlesMapMinf bundlesMapMinfData ;

  const char delim = ':' ;
  std::string line ;
  std::ifstream minfFile ;
  minfFile.open( bundlesMapMinfFilename ) ;
  if ( minfFile.fail() )
  {

    std::cout << "Problem reading file : " << bundlesMapMinfFilename
                                                                << std::endl ;
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

    if ( out[ 0 ] == "\'byte_order\'" )
    {

      bundlesMapMinfData.byte_order = out[ 1 ] ;


    }

    if ( out[ 0 ] == "\'curve3d_counts\'" )
    {

      std::string tmp_string = out[ 1 ] ;
      tmp_string.erase(std::remove( tmp_string.begin(), tmp_string.end(),
                                                    '[' ), tmp_string.end() ) ;
      tmp_string.erase(std::remove( tmp_string.begin(), tmp_string.end(),
                                                    ']' ), tmp_string.end() ) ;
      if ( verbose )
      {

        std::cout << "curve3d_counts : " << tmp_string << std::endl ;
      }

      bundlesMapMinfData.curve3d_counts = std::stoi( tmp_string ) ;

    }

    if ( out[ 0 ] == "\'io_mode\'" )
    {

      bundlesMapMinfData.io_mode = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'item_count\'" )
    {

      bundlesMapMinfData.item_count = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'label_type\'" )
    {

      bundlesMapMinfData.label_type = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'labels\'" )
    {

      bundlesMapMinfData.labels = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'object_type\'" )
    {

      bundlesMapMinfData.object_type = out[ 1 ] ;

    }

    if ( out[ 0 ] == "\'resolutionX\'" )
    {

      bundlesMapMinfData.resolution[ 0 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'resolutionY\'" )
    {

      bundlesMapMinfData.resolution[ 1 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'resolutionZ\'" )
    {

      bundlesMapMinfData.resolution[ 2 ] = std::stof( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeX\'" )
    {

      bundlesMapMinfData.size[ 0 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeY\'" )
    {

      bundlesMapMinfData.size[ 1 ] = std::stoi( out[ 1 ] ) ;

    }

    if ( out[ 0 ] == "\'sizeZ\'" )
    {

      bundlesMapMinfData.size[ 2 ] = std::stoi( out[ 1 ] ) ;

    }

  }

  minfFile.close() ;
  return bundlesMapMinfData ;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////// Read .bundlesmap ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

bundlesMapFormat bundlesMapReading( const char* bundlesMapFilename )
{

  bundlesMapFormat bundlesMapData ;

  std::ifstream file ;
  file.open( bundlesMapFilename, std::ios::binary | std::ios::in ) ;
  if ( file.fail() )
  {

    std::cout << "Problem reading file : " << bundlesMapFilename << std::endl ;
    exit( 1 ) ;

  }

  file.read( reinterpret_cast<char*>( &( bundlesMapData.labels ) ),
                                                           sizeof( int16_t ) ) ;

  if ( verbose )
  {

    std::cout << "Labels : " << bundlesMapData.labels << std::endl ;

  }

  file.read( reinterpret_cast<char*>( &( bundlesMapData.curves_count ) ),
                                                           sizeof( int32_t ) ) ;

  if ( bundlesMapData.curves_count == 0 )
  {

    std::cout << "Problem reading file, the number of traks was not stored "
              << "in the header " << std::endl ;

    std::exit( 1 ) ;

  }

  if ( verbose )
  {

    std::cout << "Curves count : " << bundlesMapData.curves_count << std::endl ;

  }

  // Getting number of points per curve
  int64_t sizeMatrixTracks = 0 ;
  bundlesMapData.pointsPerTrack.resize( bundlesMapData.curves_count ) ;
  int64_t offsetBytes = sizeof( int16_t ) + sizeof( int32_t ) ;
  for( int track = 0 ; track <  bundlesMapData.curves_count ; track++ )
  {

    if ( verbose > 1 && ( track % 1000 == 0 ||
                               ( track + 1 ) ==  bundlesMapData.curves_count ) )
    {

      printf("\rProcessing tracks [ %10d / %10d ]", track + 1 ,
                                                 bundlesMapData.curves_count ) ;
      std::cout << "" << std::flush ;

    }

    file.read( reinterpret_cast<char*>(
                                  &( bundlesMapData.pointsPerTrack[ track ] ) ),
                                                           sizeof( int32_t ) ) ;

    offsetBytes += sizeof( int32_t ) + bundlesMapData.pointsPerTrack[ track ] *
                                                           3 * sizeof( float ) ;

    file.seekg( offsetBytes, std::ios_base::beg ) ;

    sizeMatrixTracks += 3 *  bundlesMapData.pointsPerTrack[ track ] ;

  }


  // Getting vector with fiber coordinates
  bundlesMapData.matrixTracks.resize( sizeMatrixTracks, 0 ) ;
  file.seekg( sizeof( int16_t ) + sizeof( int32_t ), std::ios_base::beg ) ;
  int offsetPoints = 0 ;
  for ( int track = 0 ; track < bundlesMapData.curves_count ; track++ )
  {

    file.seekg( sizeof( int32_t ), std::ios_base::cur ) ;
    for ( int point = 0 ; point < bundlesMapData.pointsPerTrack[ track ] ;
                                                                       point++ )
    {

      file.read( reinterpret_cast<char*>(
                                 &( bundlesMapData.matrixTracks[ 0 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
      file.read( reinterpret_cast<char*>(
                                 &( bundlesMapData.matrixTracks[ 1 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;
      file.read( reinterpret_cast<char*>(
                                 &( bundlesMapData.matrixTracks[ 2 + 3 * point +
                                         offsetPoints ] ) ), sizeof( float ) ) ;

    }

    offsetPoints += 3 * bundlesMapData.pointsPerTrack[ track ] ;

    file.seekg( sizeof( int8_t ), std::ios_base::cur ) ;

  }

  file.close() ;

  if ( verbose > 1 )
  {

    std::cout << "\n" ;

  }


  return bundlesMapData ;

}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Write .trk //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void trkWriting( const char* trkFilename, trkFormat& trkData )
{

  std::ofstream file ;
  file.open( trkFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << trkFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  int size_check = 0 ;

  file.write( reinterpret_cast<char*>( &( trkData.id_string ) ),
                                                 sizeof( trkData.id_string ) ) ;
  size_check += sizeof( trkData.id_string ) ;
  file.write( reinterpret_cast<char*>( &( trkData.dim ) ),
                                                       sizeof( trkData.dim ) ) ;
  size_check += sizeof( trkData.dim ) ;
  file.write( reinterpret_cast<char*>( &( trkData.voxel_size ) ),
                                                sizeof( trkData.voxel_size ) ) ;
  size_check += sizeof( trkData.voxel_size ) ;
  file.write( reinterpret_cast<char*>( &( trkData.origin ) ),
                                                    sizeof( trkData.origin ) ) ;
  size_check += sizeof( trkData.origin ) ;
  file.write( reinterpret_cast<char*>( &( trkData.n_scalars ) ),
                                                 sizeof( trkData.n_scalars ) ) ;
  size_check += sizeof( trkData.n_scalars ) ;
  file.write( reinterpret_cast<char*>( &( trkData.scalar_name ) ),
                                               sizeof( trkData.scalar_name ) ) ;
  size_check += sizeof( trkData.scalar_name ) ;
  file.write( reinterpret_cast<char*>( &( trkData.n_properties ) ),
                                              sizeof( trkData.n_properties ) ) ;
  size_check += sizeof( trkData.n_properties ) ;
  file.write( reinterpret_cast<char*>( &( trkData.property_name ) ),
                                             sizeof( trkData.property_name ) ) ;
  size_check += sizeof( trkData.property_name ) ;
  file.write( reinterpret_cast<char*>( &( trkData.vox_to_ras ) ),
                                                sizeof( trkData.vox_to_ras ) ) ;
  size_check += sizeof( trkData.vox_to_ras ) ;
  file.write( reinterpret_cast<char*>( &( trkData.reserved ) ),
                                                  sizeof( trkData.reserved ) ) ;
  size_check += sizeof( trkData.reserved ) ;
  file.write( reinterpret_cast<char*>( &( trkData.voxel_order ) ),
                                               sizeof( trkData.voxel_order ) ) ;
  size_check += sizeof( trkData.voxel_order ) ;
  file.write( reinterpret_cast<char*>( &( trkData.pad2 ) ),
                                                      sizeof( trkData.pad2 ) ) ;
  size_check += sizeof( trkData.pad2 ) ;
  file.write( reinterpret_cast<char*>( &( trkData.image_orientation_patient ) ),
                                 sizeof( trkData.image_orientation_patient ) ) ;
  size_check += sizeof( trkData.image_orientation_patient ) ;
  file.write( reinterpret_cast<char*>( &( trkData.pad1 ) ),
                                                      sizeof( trkData.pad1 ) ) ;
  size_check += sizeof( trkData.pad1 ) ;
  file.write( reinterpret_cast<char*>( &( trkData.invert_x ) ),
                                                  sizeof( trkData.invert_x ) ) ;
  size_check += sizeof( trkData.invert_x ) ;
  file.write( reinterpret_cast<char*>( &( trkData.invert_y ) ),
                                                  sizeof( trkData.invert_y ) ) ;
  size_check += sizeof( trkData.invert_y ) ;
  file.write( reinterpret_cast<char*>( &( trkData.invert_z ) ),
                                                  sizeof( trkData.invert_z ) ) ;
  size_check += sizeof( trkData.invert_z ) ;
  file.write( reinterpret_cast<char*>( &( trkData.swap_xy ) ),
                                                   sizeof( trkData.swap_xy ) ) ;
  size_check +=  sizeof( trkData.swap_xy ) ;
  file.write( reinterpret_cast<char*>( &( trkData.swap_yz ) ),
                                                   sizeof( trkData.swap_yz ) ) ;
  size_check +=  sizeof( trkData.swap_yz ) ;
  file.write( reinterpret_cast<char*>( &( trkData.swap_zx ) ),
                                                   sizeof( trkData.swap_zx ) ) ;
  size_check +=  sizeof( trkData.swap_zx ) ;
  file.write( reinterpret_cast<char*>( &( trkData.n_count ) ),
                                                   sizeof( trkData.n_count ) ) ;
  size_check +=  sizeof( trkData.n_count ) ;
  file.write( reinterpret_cast<char*>( &( trkData.version ) ),
                                                   sizeof( trkData.version ) ) ;
  size_check +=  sizeof( trkData.version ) ;
  file.write( reinterpret_cast<char*>( &( trkData.hdr_size ) ),
                                                  sizeof( trkData.hdr_size ) ) ;
  size_check +=  sizeof( trkData.hdr_size ) ;


  int offsetProperties = 0 ;
  int offsetScalars = 0 ;
  int offsetPoints = 0 ;
  for( int track = 0 ; track < trkData.n_count ; track++ )
  {

    int nb_points = trkData.pointsPerTrack[ track ] ;

    file.write( reinterpret_cast<char*>( &( trkData.pointsPerTrack[ track ] ) ),
                                                               sizeof( int ) ) ;

    for ( int point = 0 ; point < nb_points ; point++ )
    {

      file.write( reinterpret_cast<char*>( &trkData.matrixTracks[ 0 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>( &trkData.matrixTracks[ 1 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>( &trkData.matrixTracks[ 2 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;

      if ( trkData.n_scalars > 0 )
      {

        for ( int scalar = 0 ; scalar < trkData.n_scalars ; scalar++ )
        {

          file.write( reinterpret_cast<char*>( &trkData.tracksScalars[ scalar +
                                          offsetScalars ] ), sizeof( float ) ) ;

        }

      }

    }

    if ( trkData.n_properties > 0 )
    {

      for ( int property = 0 ; property < trkData.n_properties ; property++ )
      {

        file.write( reinterpret_cast<char*>( &trkData.tracksProperties[ property
                                     + offsetProperties ] ), sizeof( float ) ) ;

      }

    }

    offsetScalars += trkData.n_scalars * trkData.pointsPerTrack[ track ] ;
    offsetPoints += 3 * trkData.pointsPerTrack[ track ] ;
    offsetProperties += trkData.n_properties ;

  }


  file.close() ;

}


///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Write .bundles ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void bundlesWriting( const char* bundlesFilename, BundlesFormat& bundlesInfo )
{

  if ( verbose > 1 )
  {

    std::cout << "\nWriting " << bundlesFilename << std::endl ;

  }

  std::ofstream file( bundlesFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path " ;
    exit( 1 );

  }

  file << "attributes = {\n" ;

  file << "    'binary' : "
       << bundlesInfo.binary
       << ",\n" ;

  file << "    'bundles' : "
       << bundlesInfo.bundles
       << ",\n" ;

  file << "    'byte_order' : "
       << bundlesInfo.byte_order
       << ",\n" ;

  file << "    'curves_count' : "
       << bundlesInfo.curves_count
       << ",\n" ;

  file << "    'radio' : "
       << bundlesInfo.radio
       << ",\n" ;

  file << "    'length' : "
       << bundlesInfo.length
       << ",\n" ;

  file << "    'nSubjects' : "
       << bundlesInfo.nSubjects
       << ",\n" ;

  file << "    'data_file_name' : "
       << bundlesInfo.data_file_name
       << ",\n" ;

  file << "    'format' : "
       << bundlesInfo.format
       << ",\n" ;

  file << "    'io_mode' : "
       << bundlesInfo.io_mode
       << ",\n" ;

  file << "    'item_count' : "
       << bundlesInfo.item_count
       << ",\n" ;

  file << "    'label_type' : "
       << bundlesInfo.label_type
       << ",\n" ;

  file << "    'labels' : "
       << bundlesInfo.labels
       << ",\n" ;

  file << "    'object_type' : "
       << bundlesInfo.object_type
       << ",\n" ;

  file << "    'averageRadius' : "
       << bundlesInfo.averageRadius
       << ",\n" ;

  file << "    'minRadius' : "
       << bundlesInfo.minRadius
       << ",\n" ;

  file << "    'maxRadius' : "
       << bundlesInfo.maxRadius
       << ",\n" ;

  file << "    'averageAngle' : "
       << bundlesInfo.averageAngle
       << ",\n" ;

  file << "    'minAngle' : "
       << bundlesInfo.minAngle
       << ",\n" ;

  file << "    'maxAngle' : "
       << bundlesInfo.maxAngle
       << ",\n" ;

  file << "    'averageDirectionAngle' : "
       << bundlesInfo.averageDirectionAngle
       << ",\n" ;

  file << "    'minDirectionAngle' : "
       << bundlesInfo.minDirectionAngle
       << ",\n" ;

  file << "    'maxDirectionAngle' : "
       << bundlesInfo.maxDirectionAngle
       << ",\n" ;

  file << "    'averageShapeAngle' : "
       << bundlesInfo.averageShapeAngle
       << ",\n" ;

  file << "    'minShapeAngle' : "
       << bundlesInfo.minShapeAngle
       << ",\n" ;

  file << "    'maxShapeAngle' : "
       << bundlesInfo.maxShapeAngle
       << ",\n" ;

  file << "    'averageLength' : "
       << bundlesInfo.averageLength
       << ",\n" ;

  file << "    'minLength' : "
       << bundlesInfo.minLength
       << ",\n" ;

  file << "    'maxLength' : "
       << bundlesInfo.maxLength
       << ",\n" ;

  file << "    'averageDisimilarity' : "
       << bundlesInfo.averageDisimilarity
       << ",\n" ;

  file << "    'minDisimilarity' : "
       << bundlesInfo.minDisimilarity
       << ",\n" ;

  file << "    'maxDisimilarity' : "
       << bundlesInfo.maxDisimilarity
       << ",\n" ;

  file << "    'density' : "
       << bundlesInfo.density
       << ",\n" ;

  file << "    'centerBundleX' : "
       << bundlesInfo.centerBundle[ 0 ]
       << ",\n" ;

  file << "    'centerBundleY' : "
       << bundlesInfo.centerBundle[ 1 ]
       << ",\n" ;

  file << "    'centerBundleZ' : "
       << bundlesInfo.centerBundle[ 2 ]
       << ",\n" ;

  file << "    'resolutionX' : "
       << bundlesInfo.resolution[ 0 ]
       << ",\n" ;

  file << "    'resolutionY' : "
       << bundlesInfo.resolution[ 1 ]
       << ",\n" ;

  file << "    'resolutionZ' : "
       << bundlesInfo.resolution[ 2 ]
       << ",\n" ;

  file << "    'sizeX' : "
       << bundlesInfo.size[ 0 ]
       << ",\n" ;

  file << "    'sizeY' : "
       << bundlesInfo.size[ 1 ]
       << ",\n" ;

  file << "    'sizeZ' : "
       << bundlesInfo.size[ 2 ]
       << ",\n" ;

  file << "    'space_dimension' : "
       << bundlesInfo.space_dimension
       << "\n    }" ;

  file.close() ;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////// Write .bundlesdata //////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void bundlesdataWriting( const char* bundlesdataFilename,
                                                BundlesDataFormat& bundlesdata )
{
  if ( verbose > 1 )
  {

    std::cout << "Writing " << bundlesdataFilename << std::endl ;

  }

  std::ofstream file ;
  file.open( bundlesdataFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << bundlesdataFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  int n_curves = bundlesdata.curves_count ;


  int offsetPoints = 0 ;


  for ( int track = 0 ; track < n_curves ; track++ )
  {

    file.write( reinterpret_cast<char*>( &bundlesdata.pointsPerTrack[ track ] ),
                                                           sizeof( int32_t ) ) ;

    for ( int point = 0 ; point < bundlesdata.pointsPerTrack[ track ] ;
                                                                       point++ )
    {

      file.write( reinterpret_cast<char*>(
                                      &bundlesdata.matrixTracks[ 0 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>(
                                      &bundlesdata.matrixTracks[ 1 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>(
                                      &bundlesdata.matrixTracks[ 2 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;

    }

    offsetPoints += 3 * bundlesdata.pointsPerTrack[ track ] ;

  }

  file.close() ;

}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Write .bundlesmap.minf ////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void bundlesMapMinfWriting( const char* bundlesMapMinfFilename,
                            bundlesMapMinf& bundlesMapMinfData )
{

  std::cout << "Writing " << bundlesMapMinfFilename << std::endl ;


  std::ofstream file( bundlesMapMinfFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path " ;
    exit( 1 );

  }

  file << "attributes = {\n" ;

  file << "    'byte_order' : "
       << bundlesMapMinfData.byte_order
       << ",\n" ;

  file << "    'curve3d_counts' : "
       << bundlesMapMinfData.curve3d_counts
       << ",\n" ;

  file << "    'io_mode' : "
       << bundlesMapMinfData.io_mode
       << ",\n" ;

  file << "    'item_count' : "
       << bundlesMapMinfData.item_count
       << ",\n" ;

  file << "    'label_type' : "
       << bundlesMapMinfData.label_type
       << ",\n" ;

  file << "    'labels' : "
       << bundlesMapMinfData.labels
       << ",\n" ;

  file << "    'object_type' : "
       << bundlesMapMinfData.object_type
       << ",\n" ;

  file << "    'resolutionX' : "
       << bundlesMapMinfData.resolution[ 0 ]
       << ",\n" ;

  file << "    'resolutionY' : "
       << bundlesMapMinfData.resolution[ 1 ]
       << ",\n" ;

  file << "    'resolutionZ' : "
       << bundlesMapMinfData.resolution[ 2 ]
       << ",\n" ;

  file << "    'sizeX' : "
       << bundlesMapMinfData.size[ 0 ]
       << ",\n" ;

  file << "    'sizeY' : "
       << bundlesMapMinfData.size[ 1 ]
       << ",\n" ;

  file << "    'sizeZ' : "
       << bundlesMapMinfData.size[ 2 ]
       << ",\n" ;

  file.close() ;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////// Write .bundlesmap //////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void bundlesMapWriting( const char* bundlesMapFilename,
                        bundlesMapFormat& bundlesMapData )
{

  if ( verbose > 1 )
  {

    std::cout << "Writing " << bundlesMapFilename << std::endl ;

  }

  std::ofstream file ;
  file.open( bundlesMapFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << bundlesMapFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }


  file.write( reinterpret_cast<char*>( &bundlesMapData.labels ),
                                                           sizeof( int16_t ) ) ;

  file.write( reinterpret_cast<char*>( &bundlesMapData.curves_count ),
                                                           sizeof( int32_t ) ) ;

  int n_curves = bundlesMapData.curves_count ;

  int offsetPoints = 0 ;
  for ( int track = 0 ; track < n_curves ; track++ )
  {

    file.write( reinterpret_cast<char*>(
                                      &bundlesMapData.pointsPerTrack[ track ] ),
                                                           sizeof( int32_t ) ) ;

    for ( int point = 0 ; point < bundlesMapData.pointsPerTrack[ track ] ;
                                                                       point++ )
    {

      file.write( reinterpret_cast<char*>(
                                   &bundlesMapData.matrixTracks[ 0 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>(
                                   &bundlesMapData.matrixTracks[ 1 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;
      file.write( reinterpret_cast<char*>(
                                   &bundlesMapData.matrixTracks[ 2 + 3 * point +
                                           offsetPoints ] ), sizeof( float ) ) ;

    }

    offsetPoints += 3 * bundlesMapData.pointsPerTrack[ track ] ;

  }

  file.close() ;


}

///////////////////////////////////////////////////////////////////////////////
//////////////// .bundlesdata to .trk with reference image  ///////////////////
///////////////////////////////////////////////////////////////////////////////
void bundlesData2Trk( const char* inputBundlesdataFormat,
                      const char* inputBundlesFormat,
                      const char* outputFile,
                      const char* refFile,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{

  std::string refString( refFile ) ;
  std::vector<std::vector<float>> vox_to_ras ;
  vox_to_ras.resize( 4 ) ;

  short dim[ 8 ] ;
  float pixdim[ 8 ] ;

  if ( refString.find( ".nii" ) != std::string::npos )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << refFile << std::endl ;

    }

    std::ifstream file ;
    file.open( refFile, std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::cout << "Problem reading file : " << refFile << std::endl ;
      exit( 1 ) ;

    }

    file.seekg( 40, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &dim ), sizeof( dim ) ) ;

    file.seekg( 76, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &pixdim ), sizeof( pixdim ) ) ;

    file.seekg( 280, std::ios_base::beg ) ;
    for ( int i = 0 ; i < 3 ; i ++ )
    {

      vox_to_ras[ i ].resize( 4, 0 ) ;

      float buffer[ 4 ] ;
      file.read( reinterpret_cast<char*>( &buffer ), sizeof( buffer ) ) ;
      for ( int j = 0 ; j < 4 ; j++ )
      {

        vox_to_ras[ i ][ j ] = buffer[ j ] ;

      }

    }
    vox_to_ras[ 3 ].resize( 4, 0.0f ) ;
    vox_to_ras[ 3 ][ 3 ] = 1.0f ;

    file.close() ;

  }
  else
  {

    std::cout << "ERROR : Reference image must be .nii " << std::endl ;
    exit( 1 ) ;

  }

  std::cout << "\n" ;
  for ( int i = 0 ; i < 4 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++)
    {

      std::cout << vox_to_ras[ i ][ j ] << "   " ;

    }
     std::cout << "\n" ;

  }

  if ( verbose )
  {

    std::cout << "Reading : " << inputBundlesFormat << std::endl ;

  }

  BundlesFormat bundlesInfo = bundlesReading( inputBundlesFormat );


  if ( bundlesInfo.resolution[ 0 ] == 0 || bundlesInfo.resolution[ 1 ] == 0
                                        || bundlesInfo.resolution[ 2 ] == 0 )
  {

    if ( verbose )
    {

      std::cout << "No resolution information in input .bundles... "
                << "Using reference image header information " << std::endl ;

    }

    for ( int i = 0 ; i < 3 ; i++ )
    {

      bundlesInfo.resolution[ i ] = pixdim[ i + 1 ] ;

    }


  }

  if ( bundlesInfo.size[ 0 ] == 0 || bundlesInfo.size[ 1 ] == 0
                                  || bundlesInfo.size[ 2 ] == 0 )
  {

    if ( verbose )
    {

      std::cout << "No size information in input .bundles... "
                << "Using reference image header information " << std::endl ;

    }

    for ( int i = 0 ; i < 3 ; i++ )
    {

      bundlesInfo.size[ i ] = dim[ i + 1 ] ;

    }

  }


  if ( verbose )
  {

    std::cout << "Reading : " << inputBundlesdataFormat << std::endl ;

  }

  BundlesDataFormat bundlesdata = bundlesdataReading( inputBundlesdataFormat,
                                                      inputBundlesFormat  ) ;

  BundlesDataFormat bundlesdataFile ;
  // Applying vox_to_ras to fibers
  applyVoxToRasToTractogram( bundlesdata.matrixTracks, vox_to_ras,
                                                bundlesdataFile.matrixTracks ) ;

  bundlesdataFile.pointsPerTrack = bundlesdata.pointsPerTrack ;


  BundlesDataFormat bundlesdataFileFlipped ;
  flipTractogram( bundlesdataFile,
                  bundlesInfo,
                  bundlesdataFileFlipped,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;
  bundlesdataFileFlipped.pointsPerTrack = bundlesdataFile.pointsPerTrack ;

  trkFormat trkData ;

  memcpy( trkData.id_string, "TRACK", sizeof( trkData.id_string ) ) ;

  if ( bundlesInfo.size[ 0 ] > 1e-5 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      trkData.dim[ i ] = bundlesInfo.size[ i ] ;


    }

  }
  else
  {

    std::cout << "Not size information in input .bundles... "
              << "Using reference image header information" << std::endl ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      trkData.dim[ i ] = dim[ i + 1 ] ;

    }

  }

  if ( bundlesInfo.resolution[ 0 ] > 1e-5 )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      trkData.voxel_size[ i ] = bundlesInfo.resolution[ i ] ;
    }

  }
  else
  {

    std::cout << "Not resolution information in input .bundles... "
              << "Using reference image header information " << std::endl ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      trkData.voxel_size[ i ] = pixdim[ i + 1 ] ;

    }

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    trkData.origin[ i ] = 0 ;

  }

  trkData.n_scalars = ( short int )0;

  trkData.n_properties = ( short int )0 ;

  for ( int i = 0 ; i < 4 ; i ++ )
  {

    for ( int j = 0 ; j < 4 ; j ++ )
    {

      trkData.vox_to_ras[ i ][ j ] = vox_to_ras[ i ][ j ] ;

    }

  }

  // memcpy( trkData.voxel_order, "LPI", sizeof( trkData.voxel_order ) ) ;
  memcpy( trkData.voxel_order, "RAS", sizeof( trkData.voxel_order ) ) ;


  for ( int i = 0 ; i < 6 ; i++ )
  {

    trkData.image_orientation_patient[ i ] = 0 ;

  }


  trkData.n_count = bundlesInfo.curves_count ;
  trkData.version = 2 ;
  trkData.hdr_size = ( int )1000 ;
  trkData.matrixTracks = bundlesdataFileFlipped.matrixTracks ;
  trkData.pointsPerTrack = bundlesdataFileFlipped.pointsPerTrack ;
  // trkData.matrixTracks = bundlesdataFile.matrixTracks ;
  // trkData.pointsPerTrack = bundlesdataFile.pointsPerTrack ;

  if ( verbose )
  {

    std::cout << "TKR Header info : " << std::endl ;
    printTrkHeader( trkData ) ;

  }

  if ( verbose )
  {

    std::cout << "\nWriting file : " << outputFile << std::endl ;

  }
  trkWriting( outputFile, trkData ) ;


}


///////////////////////////////////////////////////////////////////////////////
/////////////// Fonction to print information in .trk header ///////////////////
///////////////////////////////////////////////////////////////////////////////

void printTrkHeader( trkFormat& trkData )
{

  std::cout << "\nId_string : " << (trkData.id_string) << std::endl ;

  std::cout << "Dim : " ;
  for (int i = 0 ; i < 3 ; i++)
  {

   std::cout << (trkData.dim[i]) << "  ";

  }
  std::cout << "\n" ;

  std::cout << "Voxel size : " ;
  for (int i = 0 ; i < 3 ; i++)
  {

   std::cout << (trkData.voxel_size[ i ]) << "  ";

  }
  std::cout << "\n" ;

  std::cout << "Origin : " ;
  for (int i = 0 ; i < 3 ; i++)
  {

   std::cout << (trkData.origin[ i ]) << "  ";

  }
  std::cout << "\n" ;

  std::cout << "Number of scalars : " << (trkData.n_scalars) << std::endl;

  std::cout << "Scalar names : " << std::endl ;
  if ( trkData.n_scalars > 0 )
  {

    for ( int i = 0 ; i < 10 ; i++ )
    {

      std::cout << trkData.scalar_name[ i ] << std::endl;

    }
    std::cout << "\n" ;

  }

  std::cout << "Number of properties : " << trkData.n_properties << std::endl;

  std::cout << "Properties names : " << std::endl ;
  if ( trkData.n_properties > 0 )
  {

    for ( int i = 0 ; i < 10 ; i ++)
    {

       std::cout << trkData.property_name[ i ] << std::endl ;

    }
    std::cout << "\n";

  }

  std::cout << "Vox to ras : " << std::endl ;
  for ( int i = 0 ; i < 4 ; i ++)
  {

    for ( int j = 0 ; j < 4 ; j ++)
    {

      std::cout << trkData.vox_to_ras[ i ][ j ] << " " ;

    }

    std::cout << "\n" ;

  }

  std::cout << "Reserved : " << trkData.reserved << std::endl ;

  // std::cout << "Voxel order : " << trkData.voxel_order << std::endl ;
  std::cout << "Voxel order : " ;
  for ( int i = 0 ; i < 4 ; i++ )
  {

    std::cout << trkData.voxel_order[ i ] ;

  }
  std::cout << "\n" ;

  std::cout << "Padding 2 : " << trkData.pad2 << std::endl ;

  std::cout << "Image orientation : " << std::endl ;
  for ( int i = 0 ; i < 6 ; i++ )
  {

    std::cout << trkData.image_orientation_patient[ i ] << "  " ;

  }
  std::cout << "\n" ;

  std::cout << "Padding 1 : " << trkData.pad1 << std::endl ;

  std::cout << "Invert x : " << trkData.invert_x << std::endl ;

  std::cout << "Invert y : " << trkData.invert_y << std::endl ;

  std::cout << "Invert z : " << trkData.invert_z << std::endl ;

  std::cout << "Swap xy : " << trkData.swap_xy << std::endl ;

  std::cout << "Swap yz : " << trkData.swap_yz << std::endl ;

  std::cout << "Swap zx : " << trkData.swap_zx << std::endl ;

  std::cout << "Number of tracks stored : " << trkData.n_count << std::endl ;

  std::cout << "Version : " << trkData.version << std::endl ;

  std::cout << "Header size : " << trkData.hdr_size << std::endl ;

}



///////////////////////////////////////////////////////////////////////////////
//////////////////////// Function trk -> bundlesdata //////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Trk2BundlesData( const char* inputFile,
                      const char* refFile,
                      const char* bundlesFilename,
                      const char* bundlesdataFilename,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{

  // ----------------------------- Reading .trk ----------------------------- //

  trkFormat trkFile = trkReading( inputFile ) ;

  // ----------------------- Reading reference image ------------------------ //
  std::string refString( refFile ) ;
  short dim[ 8 ] ;
  float pixdim[ 8 ] ;
  std::vector<std::vector<float>> vox_to_ras ;
  vox_to_ras.resize( 4 ) ;
  if ( refString.find( ".nii" ) != std::string::npos )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << refFile << std::endl ;

    }

    std::ifstream file ;
    file.open( refFile, std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::cout << "Problem reading file : " << refFile << std::endl ;
      exit( 1 ) ;

    }

    file.seekg( 40, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &dim ), sizeof( dim ) ) ;

    file.seekg( 76, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &pixdim ), sizeof( pixdim ) ) ;

    file.seekg( 280, std::ios_base::beg ) ;
    for ( int i = 0 ; i < 3 ; i ++ )
    {

      vox_to_ras[ i ].resize( 4, 0 ) ;

      float buffer[ 4 ] ;
      file.read( reinterpret_cast<char*>( &buffer ), sizeof( buffer ) ) ;
      for ( int j = 0 ; j < 4 ; j++ )
      {

        vox_to_ras[ i ][ j ] = buffer[ j ] ;

      }

    }
    vox_to_ras[ 3 ].resize( 4, 0.0f ) ;
    vox_to_ras[ 3 ][ 3 ] = 1.0f ;

    file.close() ;

  }
  else
  {

    std::cout << "ERROR : Reference image must be .nii " << std::endl ;
    exit( 1 ) ;

  }

  // ----------------------- Print header if verbose ------------------------ //

  if ( verbose )
  {

    printTrkHeader( trkFile ) ;

  }

  // -------------------------- Writing .bundles ---------------------------- //

  BundlesFormat bundlesFile ;

  bundlesFile.binary = 1 ;
  bundlesFile.bundles = "[ '255', 0 ]" ;
  bundlesFile.byte_order = "'DCBA'" ;
  bundlesFile.curves_count = ( int )trkFile.n_count ;
  bundlesFile.data_file_name = "'*.bundlesdata'" ;
  bundlesFile.format = "'bundles_1.0'" ;
  bundlesFile.io_mode = "'binary'" ;
  bundlesFile.item_count = ( int )1 ;
  bundlesFile.label_type = "'std_string'" ;
  bundlesFile.labels = "[ '255' ]" ;
  bundlesFile.object_type = "'BundleMap'" ;
  bundlesFile.resolution[ 0 ] = trkFile.voxel_size[ 0 ] ;
  bundlesFile.resolution[ 1 ] = trkFile.voxel_size[ 1 ] ;
  bundlesFile.resolution[ 2 ] = trkFile.voxel_size[ 2 ] ;
  bundlesFile.size[ 0 ] = trkFile.dim[ 0 ] ;
  bundlesFile.size[ 1 ] = trkFile.dim[ 1 ] ;
  bundlesFile.size[ 2 ] = trkFile.dim[ 2 ] ;
  bundlesFile.space_dimension = ( int )3 ;

  if ( verbose )
  {

    std::cout << "Writing bundles file : " << bundlesFilename << std::endl ;

  }


  bundlesWriting( bundlesFilename, bundlesFile ) ;


  // ------------------------ Writing .bundlesdata -------------------------- //
  std::cout << "\n" ;
  for ( int i = 0 ; i < 4 ; i++ )
  {

    for ( int j = 0 ; j < 4 ; j++)
    {

      std::cout << vox_to_ras[ i ][ j ] << "   " ;

    }
    std::cout << "\n" ;

  }


  BundlesDataFormat bundlesdataFileFlipped ;
  flipTractogram( trkFile.matrixTracks,
                  bundlesFile.resolution,
                  bundlesFile.size,
                  bundlesdataFileFlipped.matrixTracks,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;

  bundlesdataFileFlipped.pointsPerTrack = trkFile.pointsPerTrack ;
  bundlesdataFileFlipped.curves_count = trkFile.n_count ;

  std::vector<std::vector<float>> ras_to_vox ;
  computeInverseVoxToRas( vox_to_ras, ras_to_vox ) ;

  BundlesDataFormat bundlesdataFile ;
  // Applying vox_to_ras to tck fibers
  applyVoxToRasToTractogram( bundlesdataFileFlipped.matrixTracks, ras_to_vox,
                                                bundlesdataFile.matrixTracks ) ;

  bundlesdataFile.pointsPerTrack = bundlesdataFileFlipped.pointsPerTrack ;
  bundlesdataFile.curves_count = bundlesdataFileFlipped.curves_count ;
  // applyVoxToRasToTractogram( trkFile.matrixTracks, ras_to_vox,
  //                                               bundlesdataFile.matrixTracks ) ;
  //
  // bundlesdataFile.pointsPerTrack = trkFile.pointsPerTrack ;
  // bundlesdataFile.curves_count = trkFile.n_count ;

  if ( verbose )
  {

    std::cout << "Writing bundlesdata file : "
              << bundlesdataFilename << std::endl ;

  }

  bundlesdataWriting( bundlesdataFilename, bundlesdataFile ) ;

}

///////////////////////////////////////////////////////////////////////////////
//////////////////////// Function trk -> bundlesdata //////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Tck2BundlesData( const char* inputFile,
                      const char* refFile,
                      const char* bundlesMinfFilename,
                      const char* bundlesFilename,
                      const char* bundlesdataFilename,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{

  // ----------------------------- Reading .tck ----------------------------- //
  tckFormat tckFile = tckReading( inputFile ) ;

  // ----------------------- Reading reference image ------------------------ //
  std::string refString( refFile ) ;
  short dim[ 8 ] ;
  float pixdim[ 8 ] ;
  std::vector<std::vector<float>> vox_to_ras ;
  vox_to_ras.resize( 4 ) ;
  if ( refString.find( ".nii" ) != std::string::npos )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << refFile << std::endl ;

    }

    std::ifstream file ;
    file.open( refFile, std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::cout << "Problem reading file : " << refFile << std::endl ;
      exit( 1 ) ;

    }

    file.seekg( 40, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &dim ), sizeof( dim ) ) ;

    file.seekg( 76, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &pixdim ), sizeof( pixdim ) ) ;

    file.seekg( 280, std::ios_base::beg ) ;
    for ( int i = 0 ; i < 3 ; i ++ )
    {

      vox_to_ras[ i ].resize( 4, 0 ) ;

      float buffer[ 4 ] ;
      file.read( reinterpret_cast<char*>( &buffer ), sizeof( buffer ) ) ;
      for ( int j = 0 ; j < 4 ; j++ )
      {

        vox_to_ras[ i ][ j ] = buffer[ j ] ;

      }

    }
    vox_to_ras[ 3 ].resize( 4, 0.0f ) ;
    vox_to_ras[ 3 ][ 3 ] = 1.0f ;

    file.close() ;

  }
  else
  {

    std::cout << "ERROR : Reference image must be .nii " << std::endl ;
    exit( 1 ) ;

  }

  // ---------------------------- Writing .minf ----------------------------- //
  std::cout << "Writing " << bundlesMinfFilename << std::endl ;


  std::ofstream file( bundlesMinfFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path "
              << bundlesMinfFilename << std::endl ;
    exit( 1 );

  }

  for ( std::string line : tckFile.headerInfo )
  {

    file << line << std::endl ;

  }

  file.close() ;


  // --------------------------- Writing .bundles --------------------------- //

  BundlesFormat bundlesFile ;

  bundlesFile.binary = 1 ;
  bundlesFile.bundles = "[ '255', 0 ]" ;
  bundlesFile.byte_order = "'DCBA'" ;
  bundlesFile.curves_count = ( int )tckFile.curves_count ;
  bundlesFile.data_file_name = "'*.bundlesdata'" ;
  bundlesFile.format = "'bundles_1.0'" ;
  bundlesFile.io_mode = "'binary'" ;
  bundlesFile.item_count = ( int )1 ;
  bundlesFile.label_type = "'std_string'" ;
  bundlesFile.labels = "[ '255' ]" ;
  bundlesFile.object_type = "'BundleMap'" ;
  bundlesFile.resolution[ 0 ] = pixdim[ 0 + 1 ] ;
  bundlesFile.resolution[ 1 ] = pixdim[ 1 + 1 ] ;
  bundlesFile.resolution[ 2 ] = pixdim[ 2 + 1 ] ;
  bundlesFile.size[ 0 ] = dim[ 0 + 1 ] ;
  bundlesFile.size[ 1 ] = dim[ 1 + 1 ] ;
  bundlesFile.size[ 2 ] = dim[ 2 + 1 ] ;
  bundlesFile.space_dimension = ( int )3 ;

  bundlesWriting( bundlesFilename, bundlesFile ) ;


  // ------------------------- Writing .bundlesdata ------------------------- //
  std::vector<std::vector<float>> ras_to_vox ;
  computeInverseVoxToRas( vox_to_ras, ras_to_vox ) ;

  BundlesDataFormat bundlesdataFile ;
  // Applying vox_to_ras to tck fibers
  applyVoxToRasToTractogram( tckFile.matrixTracks, ras_to_vox,
                                                bundlesdataFile.matrixTracks ) ;

  // bundlesdataFile.matrixTracks = tckFile.matrixTracks ;
  bundlesdataFile.pointsPerTrack = tckFile.pointsPerTrack ;
  bundlesdataFile.curves_count = tckFile.curves_count ;

  BundlesDataFormat bundlesdataFileFlipped ;
  flipTractogram( bundlesdataFile,
                  bundlesFile,
                  bundlesdataFileFlipped,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;

  bundlesdataFileFlipped.pointsPerTrack = bundlesdataFile.pointsPerTrack ;
  bundlesdataFileFlipped.curves_count = bundlesdataFile.curves_count ;

  bundlesdataWriting( bundlesdataFilename, bundlesdataFileFlipped ) ;


}

///////////////////////////////////////////////////////////////////////////////
//////////////////////// Function bundlesdata -> tck //////////////////////////
///////////////////////////////////////////////////////////////////////////////
void BundlesData2Tck( const char* inputBundlesdataFormat,
                      const char* inputBundlesFormat,
                      const char* inputBundlesMinf,
                      const char* refFile,
                      const char* outputTckFilename,
                      int flip_x,
                      int flip_y,
                      int flip_z )
{

  tckFormat tckData ;

  // --------------------------- Reading .bundles --------------------------- //
  BundlesFormat inputBundles( inputBundlesFormat, 1 ) ;

  // ------------------------- Reading .bundlesdata ------------------------- //
  BundlesDataFormat inputBundlesdata( inputBundlesdataFormat,
                                      inputBundlesFormat,
                                      1 ) ;

  // ----------------------------- Reading .minf ---------------------------- //
  std::cout << "Reading .minf" << std::endl ;
  std::string inputBundlesMinfString = inputBundlesMinf ;
  if ( is_file( inputBundlesMinfString ) )
  {

    std::string dataTypeBinary = "" ;
    std::fstream infile ;
    infile.open( inputBundlesMinf, std::ios::in ) ;
    if ( infile.is_open() )
    {

      std::string line ;
      bool notFinished = true ;
      std::string spaceDelimiter = " " ;
      while ( std::getline( infile, line ) && notFinished )
      {

        tckData.headerInfo.push_back( line ) ;

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

          tckData.offsetBinary = std::stoi( words[ words.size() - 1 ] ) ;

        }

        if ( words[ 0 ] == "datatype:" )
        {

          dataTypeBinary = words[ words.size() - 1 ] ;

        }

        if ( words[ 0 ] == "count:" )
        {

          tckData.curves_count = std::stoi( words[ words.size() - 1 ] ) ;

        }

        // if ( ( tckData.offsetBinary > 0 && dataTypeBinary != "" &&
        //        tckData.curves_count > 0 ) || ( words[ 0 ] == "END" ) )
        if (  words[ 0 ] == "END" )
        {

          notFinished = false ;

        }

      }

      infile.close() ;

    }
    else
    {

      std::cout << "ERROR : Problem reading " << inputBundlesMinf << std::endl ;
      exit( 1 ) ;

    }


  }
  else
  {

    std::cout << "Creating .minf information" << std::endl ;

    std::string matrixTracksString = "mrtrix tracks    " ;

    std::string timeStampString = "timestamp: 1661165504.3768622875" ;

    std::string datatypeString = "datatype: Float32LE" ;

    std::ostringstream countStringOss ;
    countStringOss << "count: " << inputBundlesdata.curves_count ;
    std::string countString = countStringOss.str() ;

    std::ostringstream totalCountStringOss ;
    totalCountStringOss << "total_count: " << inputBundlesdata.curves_count ;
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

    tckData.headerInfo.push_back( matrixTracksString ) ;
    tckData.headerInfo.push_back( timeStampString ) ;
    tckData.headerInfo.push_back( datatypeString ) ;
    tckData.headerInfo.push_back( fileString ) ;
    tckData.headerInfo.push_back( countString ) ;
    tckData.headerInfo.push_back( totalCountString ) ;
    tckData.headerInfo.push_back( endString ) ;

    tckData.offsetBinary = tmpSizeBytes + 1 ;

  }



  // ----------------------- Reading reference image ------------------------ //
  std::string refString( refFile ) ;
  short dim[ 8 ] ;
  float pixdim[ 8 ] ;
  std::vector<std::vector<float>> vox_to_ras ;
  vox_to_ras.resize( 4 ) ;
  if ( refString.find( ".nii" ) != std::string::npos )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << refFile << std::endl ;

    }

    std::ifstream file ;
    file.open( refFile, std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::cout << "Problem reading file : " << refFile << std::endl ;
      exit( 1 ) ;

    }

    file.seekg( 40, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &dim ), sizeof( dim ) ) ;

    file.seekg( 76, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &pixdim ), sizeof( pixdim ) ) ;

    file.seekg( 280, std::ios_base::beg ) ;
    for ( int i = 0 ; i < 3 ; i ++ )
    {

      vox_to_ras[ i ].resize( 4, 0 ) ;

      float buffer[ 4 ] ;
      file.read( reinterpret_cast<char*>( &buffer ), sizeof( buffer ) ) ;
      for ( int j = 0 ; j < 4 ; j++ )
      {

        vox_to_ras[ i ][ j ] = buffer[ j ] ;

      }

    }
    vox_to_ras[ 3 ].resize( 4, 0.0f ) ;
    vox_to_ras[ 3 ][ 3 ] = 1.0f ;

    file.close() ;

  }
  else
  {

    std::cout << "ERROR : Reference image must be .nii " << std::endl ;
    exit( 1 ) ;

  }

  // Correcting resolution of .bundles if not present
  bool isResolution = false ;
  for ( float _value : inputBundles.resolution )
  {

    if ( _value !=0 )
    {

      isResolution = true ;

    }

  }
  if ( !isResolution )
  {
    int counter = 1 ;  // Resolution starts at index 1 in .nii
    for ( float& _value : inputBundles.resolution )
    {
      _value = pixdim[ counter ] ;
      counter++ ;

    }

  }

  // Correcting size of .bundles if not present
  bool isSize = false ;
  for ( short int _value : inputBundles.size )
  {

    if ( _value !=0 )
    {

      isSize = true ;

    }

  }
  if ( !isSize )
  {
    int counter = 1 ;  // Size starts at index 1 in .nii
    for ( short int& _value : inputBundles.size )
    {
      _value = dim[ counter ] ;
      counter++ ;

    }

  }

  // ---------------------------- Writing .tck ----------------------------- //
  std::vector<float> flippedMatrixTracks ;
  flipTractogram( inputBundlesdata.matrixTracks,
                  inputBundles.resolution,
                  inputBundles.size,
                  flippedMatrixTracks,
                  - flip_x,
                  - flip_y,
                  - flip_z ) ;


  // Applying vox_to_ras to tck fibers
  applyVoxToRasToTractogram( flippedMatrixTracks, vox_to_ras,
                                                        tckData.matrixTracks ) ;

  tckData.pointsPerTrack = inputBundlesdata.pointsPerTrack ;
  tckData.curves_count = inputBundlesdata.curves_count ;
  tckData.sizeDataType = sizeof( float ) ;

  std::cout << "Writing .tck" << std::endl ;
  tckWriting( outputTckFilename, tckData ) ;


}


///////////////////////////////////////////////////////////////////////////////
///////////// Function to apply vox_to_ras to point matrices //////////////////
///////////////////////////////////////////////////////////////////////////////
void applyVoxToRasToTractogram(
                              const std::vector<float>& tractogram,
                              const std::vector<std::vector<float>>& vox_to_ras,
                              std::vector<float>& newTractogram )
{

  int64_t sizeTractogram = tractogram.size() ;
  newTractogram.resize( sizeTractogram ) ;

  int64_t nbPointsTractogram = ( int64_t )( sizeTractogram / 3 ) ;

  #pragma omp parallel for
  for ( int pointIndex = 0 ; pointIndex < nbPointsTractogram ; pointIndex++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      newTractogram[ 3 * pointIndex + i ] = 0 ;
      for ( int j = 0 ; j < 3 ; j++ )
      {

        newTractogram[ 3 * pointIndex + i ] += tractogram[ 3 * pointIndex + j ]
                                                        * vox_to_ras[ i ][ j ] ;
      }
      newTractogram[ 3 * pointIndex + i ] += vox_to_ras[ i ][ 3 ] ;

    }

  }

}


///////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute adjacent of 3x3 matrix ////////////////////////
///////////////////////////////////////////////////////////////////////////////
void computeAdjugateMatrix( const std::vector<std::vector<float>>& matrix,
                            std::vector<std::vector<float>>& adjugate )
{

  adjugate.resize( 3 ) ;
  adjugate[ 0 ].resize( 3 ) ;
  adjugate[ 1 ].resize( 3 ) ;
  adjugate[ 2 ].resize( 3 ) ;

  // Sanity checks
  if ( matrix.size() != 3 )
  {

    std::cout << "ERRROR in computeAdjugateMatrix : The input matrix must be "
              << "a 3x3 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 3 )
    {

      std::cout << "ERRROR in computeAdjugateMatrix : The input matrix must be "
                << "a 3x3 matrix" << std::endl ;
      exit( 1 ) ;

    }

  }

  // Transposing matrix
  std::vector<std::vector<float>> matrixTranspose ;
  matrixTranspose.resize( 3 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    matrixTranspose[ i ].resize( 3 ) ;
    for ( int j = 0 ; j < 3 ; j++ )
    {
      matrixTranspose[ i ][ j ] = matrix[ j ][ i ] ;

    }

  }

  // Computing adjugate matrix
  // #pragma omp parallel for
  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      std::vector<std::vector<float>> minorMatrix ;
      minorMatrix.resize( 2 ) ;
      minorMatrix[ 0 ].resize( 2 ) ;
      minorMatrix[ 1 ].resize( 2 ) ;

      int lineIndex = 0 ;
      int columnIndex = 0 ;
      bool stop = false ;
      for ( int n = 0 ; n < 3 ; n++ )
      {

        for ( int p = 0 ; p < 3 ; p++ )
        {

          if ( n != i && p != j )
          {

            minorMatrix[ lineIndex ][ columnIndex ] =
                                                     matrixTranspose[ n ][ p ] ;
            if ( columnIndex == 1 )
            {

              if ( lineIndex < 1 )
              {

                lineIndex++ ;
                columnIndex = 0 ;

              }
              else
              {

                stop = true ;
                break ;

              }

            }
            else
            {

              columnIndex++ ;

            }

          }

        }

        if ( stop )
        {

          break ;

        }

      }

      float sign = 0 ;
      if ( ( ( i + j ) % 2 == 0 ) )
      {

        sign = 1 ;

      }
      else
      {

        sign = -1 ;

      }

      adjugate[ i ][ j ] = minorMatrix[ 0 ][ 0 ] * minorMatrix[ 1 ][ 1 ] -
                           minorMatrix[ 1 ][ 0 ] * minorMatrix[ 0 ][ 1 ] ;
      adjugate[ i ][ j ] *= sign ;

    }

  }

}

///////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute determinant 3x3 matrix ////////////////////////
///////////////////////////////////////////////////////////////////////////////
float computeDeterminantMatrix( const std::vector<std::vector<float>>& matrix )
{

  // Sanity checks
  if ( matrix.size() != 3 )
  {

    std::cout << "ERRROR in computeDeterminantMatrix : The input matrix must "
              << "be a 3x3 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 3 )
    {

      std::cout << "ERRROR in computeDeterminantMatrix : The input matrix must "
                << "be a 3x3 matrix" << std::endl ;
      exit( 1 ) ;

    }

  }


  // Computing determinant
  float determinant = 0 ;
  int j = 0 ; // Select first row to compute determinant
  // #pragma omp parallel for
  for ( int i = 0 ; i < 3 ; i++ )
  {

    std::vector<std::vector<float>> minorMatrix ;
    minorMatrix.resize( 2 ) ;
    minorMatrix[ 0 ].resize( 2 ) ;
    minorMatrix[ 1 ].resize( 2 ) ;

    int lineIndex = 0 ;
    int columnIndex = 0 ;
    bool stop = false ;
    for ( int n = 0 ; n < 3 ; n++ )
    {

      for ( int p = 0 ; p < 3 ; p++ )
      {

        if ( n != i && p != j )
        {

          minorMatrix[ lineIndex ][ columnIndex ] = matrix[ n ][ p ] ;
          if ( columnIndex == 1 )
          {

            if ( lineIndex < 1 )
            {

              lineIndex++ ;
              columnIndex = 0 ;

            }
            else
            {

              stop = true ;
              break ;

            }

          }
          else
          {

            columnIndex++ ;

          }

        }

      }

      if ( stop )
      {

        break ;

      }

    }

    float sign = 0 ;
    if ( ( ( i + j ) % 2 == 0 ) )
    {

      sign = 1 ;

    }
    else
    {

      sign = -1 ;

    }

    determinant += sign * matrix[ i ][ j ] * (
                               minorMatrix[ 0 ][ 0 ] * minorMatrix[ 1 ][ 1 ] -
                               minorMatrix[ 1 ][ 0 ] * minorMatrix[ 0 ][ 1 ] ) ;

  }

  return determinant ;

}


///////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute inverse 3x3 matrix ////////////////////////
///////////////////////////////////////////////////////////////////////////////
void computeInverseMatrix( const std::vector<std::vector<float>>& matrix,
                           std::vector<std::vector<float>>& inverseMatrix )
{

  // Sanity checks
  if ( matrix.size() != 3 )
  {

    std::cout << "ERRROR in computeInverseMatrix : The input matrix must be a "
              << "3x3 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( std::vector<float> line : matrix )
  {

    if ( line.size() != 3 )
    {

      std::cout << "ERRROR in computeInverseMatrix : The input matrix must be "
                << "a 3x3 matrix" << std::endl ;
      exit( 1 ) ;

    }

  }

  float determinant = computeDeterminantMatrix( matrix ) ;
  if ( abs( determinant ) < 1e-30 )
  {

    std::cout << "ERROR in computeInverseMatrix : determinant is 0, singular "
              << "matrix are not invertible" << std::endl ;
    exit( 1 ) ;

  }

  std::vector<std::vector<float>> adjugateMatrix ;
  computeAdjugateMatrix( matrix, adjugateMatrix ) ;

  inverseMatrix.resize( 3 ) ;
  inverseMatrix[ 0 ].resize( 3, 0 ) ;
  inverseMatrix[ 1 ].resize( 3, 0 ) ;
  inverseMatrix[ 2 ].resize( 3, 0 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      inverseMatrix[ i ][ j ] = 1 / determinant * adjugateMatrix[ i ][ j ] ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
/////////////////////// Compute inverse vox_to_ras /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void computeInverseVoxToRas( const std::vector<std::vector<float>>& vox_to_ras,
                             std::vector<std::vector<float>>& ras_to_vox )
{

  // Sanity checks
  if ( vox_to_ras.size() != 4 )
  {

    std::cout << "ERRROR in computeInverseVoxToRas : The input matrix must be "
              << "a 4x4 matrix" << std::endl ;
    exit( 1 ) ;

  }
  for ( int i = 0 ; i < 4 ; i++ )
  {

    if ( vox_to_ras[ i ].size() != 4 )
    {

      std::cout << "ERRROR in computeInverseVoxToRas : The input matrix must "
                << "be a 4x4 matrix" << std::endl ;
      exit( 1 ) ;

    }

    if ( i == 3 )
    {

      if ( vox_to_ras[ i ][ 0 ] != 0.0f || vox_to_ras[ i ][ 1 ] != 0.0f ||
           vox_to_ras[ i ][ 2 ] != 0.0f || vox_to_ras[ i ][ 3 ] != 1.0f )
        {

          std::cout << "ERRROR in computeInverseVoxToRas : The input matrix is "
          << "not a 4x4 vox_to_ras matrix" << std::endl ;
          exit( 1 ) ;

        }

    }

  }

  // Resizeing ras_to_vox
  ras_to_vox.resize( 4 ) ;
  ras_to_vox[ 0 ].resize( 4, 0 ) ;
  ras_to_vox[ 1 ].resize( 4, 0 ) ;
  ras_to_vox[ 2 ].resize( 4, 0 ) ;
  ras_to_vox[ 3 ].resize( 4, 0 ) ;

  // Getting 3x3 subTransform
  std::vector<std::vector<float>> subTransform ;
  subTransform.resize( 3 ) ;
  subTransform[ 0 ].resize( 3 ) ;
  subTransform[ 1 ].resize( 3 ) ;
  subTransform[ 2 ].resize( 3 ) ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      subTransform[ i ][ j ] = vox_to_ras[ i ][ j ] ;

    }

  }

  // Computing inverse 3x3 subTransform
  std::vector<std::vector<float>> inverseSubTransform ;
  computeInverseMatrix( subTransform, inverseSubTransform ) ;

  // Getting translation
  std::vector<float> translation( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    translation[ i ] = vox_to_ras[ i ][ 3 ] ;

  }

  // Compute inverseSubTransform.dot( - traslation )
  std::vector<float> translation2( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      translation2[ i ] += - inverseSubTransform[ i ][ j ] * translation[ j ] ;

    }

  }


  // Filling ras_to_vox
  for ( int i = 0 ; i < 3 ; i++ )
  {

    for ( int j = 0 ; j < 3 ; j++ )
    {

      ras_to_vox[ i ][ j ] = inverseSubTransform[ i ][ j ] ;

    }

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    ras_to_vox[ i ][ 3 ] = translation2[ i ] ;

  }

  ras_to_vox[ 3 ][ 3 ] = 1.0f ;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////// Function bundlesmap -> bundlesdata //////////////////////
///////////////////////////////////////////////////////////////////////////////

void BundlesMap2BundlesData( const char* inputBundlesMapFilename,
                             const char* inputBundlesMapMinfFilename,
                             const char* bundlesFilename,
                             const char* bundlesdataFilename )
{

  if ( verbose )
  {

    std::cout << "Reading : " << inputBundlesMapFilename << "..." << std::endl ;

  }

  bundlesMapFormat bundlesMapFormatData = bundlesMapReading(
                                           inputBundlesMapFilename ) ;


  if ( verbose )
  {

    std::cout << "Reading : " << inputBundlesMapMinfFilename << "..."
                                                             << std::endl ;

  }

  bundlesMapMinf bundlesMapMinfData = bundlesMapMinfReading(
                                                inputBundlesMapMinfFilename ) ;


  BundlesDataFormat bundlesdataFile ;
  bundlesdataFile.matrixTracks = bundlesMapFormatData.matrixTracks ;
  bundlesdataFile.pointsPerTrack = bundlesMapFormatData.pointsPerTrack ;


  if ( verbose )
  {

    std::cout << "Writing bundlesdata file : "
              << bundlesdataFilename << std::endl ;

  }

  bundlesdataWriting( bundlesdataFilename, bundlesdataFile ) ;


  BundlesFormat bundlesFile ;

  bundlesFile.binary = 1 ;
  bundlesFile.bundles = "[ '255', 0 ]" ;
  bundlesFile.byte_order = "'DCBA'" ;
  bundlesFile.curves_count = ( int )bundlesMapMinfData.curve3d_counts ;
  bundlesFile.data_file_name = "'*.bundlesdata'" ;
  bundlesFile.format = "'bundles_1.0'" ;
  bundlesFile.io_mode = "'binary'" ;
  bundlesFile.item_count = ( int )1 ;
  bundlesFile.label_type = "'std_string'" ;
  bundlesFile.labels = "[ '255' ]" ;
  bundlesFile.object_type = "'BundleMap'" ;
  bundlesFile.resolution[ 0 ] = bundlesMapMinfData.resolution[ 0 ] ;
  bundlesFile.resolution[ 1 ] = bundlesMapMinfData.resolution[ 1 ] ;
  bundlesFile.resolution[ 2 ] = bundlesMapMinfData.resolution[ 2 ] ;
  bundlesFile.size[ 0 ] = bundlesMapMinfData.size[ 0 ] ;
  bundlesFile.size[ 1 ] = bundlesMapMinfData.size[ 1 ] ;
  bundlesFile.size[ 2 ] = bundlesMapMinfData.size[ 2 ] ;
  bundlesFile.space_dimension = 3 ;

  if ( verbose )
  {

    std::cout << "Writing bundles file : " << bundlesFilename << std::endl ;

  }

  bundlesWriting( bundlesFilename, bundlesFile ) ;


  // std::cout << "Conversion .bundlesmap -> .bundlesdata (.bunldes) not "
  //           << "implemented." << std::endl ;
  // exit( 1 ) ;

}


///////////////////////////////////////////////////////////////////////////////
////////////////////////////// Flip tractogrtam ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void flipTractogram( BundlesDataFormat& inputTractogram,
                     BundlesFormat& inputTractogramInfo,
                     BundlesDataFormat& outputTractogram,
                     int flip_x,
                     int flip_y,
                     int flip_z )
{

  int nbFibers = inputTractogramInfo.curves_count ;
  int64_t nbPointsTractogram = ( int64_t )( inputTractogram.matrixTracks.size()
                                                                         / 3 ) ;

  if ( nbFibers != inputTractogram.pointsPerTrack.size() )
  {

    std::cout << "ERROR : the number of fibers in the .bundles file is "
             << "different from the number of fibers in .bundlesdata file \n" ;
    exit( 1 ) ;

  }

  for ( int fiber = 0 ; fiber < nbFibers ; fiber++ )
  {

    int nbPoints = inputTractogram.pointsPerTrack[ fiber ] ;
    if ( nbPoints < 0 )
    {

      std::cout << "ERROR : the number of points for each fiber has to be "
                << "strictly positive " << std::endl ;
      exit( 1 ) ;

    }

  }


  outputTractogram.pointsPerTrack = inputTractogram.pointsPerTrack ;
  outputTractogram.curves_count = inputTractogram.curves_count ;
  outputTractogram.matrixTracks.resize( inputTractogram.matrixTracks.size() ) ;

  float translation[ 3 ] = { 0 } ;
  int flip[ 3 ] = { 0 } ;
  flip[ 0 ] = flip_x ;
  flip[ 1 ] = flip_y ;
  flip[ 2 ] = flip_z ;


  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( flip[ i ] == -1 )
    {

      int tmpMiddle = ceilf( inputTractogramInfo.size[ i ] / 2 ) + 1 ;

      if ( inputTractogramInfo.size[ i ] % 2 != 0 )
      {

        translation[ i ] = tmpMiddle * inputTractogramInfo.resolution[ i ] ;

      }
      else
      {

        translation[ i ] = ( ( 2 * tmpMiddle - 1 ) / 2 ) *
                                          inputTractogramInfo.resolution[ i ] ;

      }

    }

  }

  std::vector< float >& tractogramFibers = inputTractogram.matrixTracks ;
  for ( int point = 0 ; point < nbPointsTractogram ; point++ )
  {

    outputTractogram.matrixTracks[ 3 * point + 0 ]  =
                                      flip_x * tractogramFibers[ 3 * point + 0 ]
                                                        + 2 * translation[ 0 ] ;

    outputTractogram.matrixTracks[ 3 * point + 1 ]  =
                                      flip_y * tractogramFibers[ 3 * point + 1 ]
                                                        + 2 * translation[ 1 ] ;

    outputTractogram.matrixTracks[ 3 * point + 2 ]  =
                                      flip_z * tractogramFibers[ 3 * point + 2 ]
                                                        + 2 * translation[ 2 ] ;


  }

}



// -------------------------------------------------------------------------- //

void flipTractogram( const std::vector<float>& inputTractogram,
                     const std::vector<float>& resolution,
                     const std::vector<short int>& size,
                     std::vector<float>& outputTractogram,
                     int flip_x,
                     int flip_y,
                     int flip_z )
{

  int64_t nbPointsTractogram = ( int64_t )( inputTractogram.size() / 3 ) ;

  outputTractogram.resize( inputTractogram.size() ) ;

  float translation[ 3 ] = { 0 } ;
  int flip[ 3 ] = { 0 } ;
  flip[ 0 ] = flip_x ;
  flip[ 1 ] = flip_y ;
  flip[ 2 ] = flip_z ;


  for ( int i = 0 ; i < 3 ; i++ )
  {

    if ( flip[ i ] == -1 )
    {

      int tmpMiddle = ceilf( size[ i ] / 2 ) + 1 ;

      if ( size[ i ] % 2 != 0 )
      {

        translation[ i ] = tmpMiddle * resolution[ i ] ;

      }
      else
      {

        translation[ i ] = ( ( 2 * tmpMiddle - 1 ) / 2 ) * resolution[ i ] ;

      }

    }

  }


  for ( int point = 0 ; point < nbPointsTractogram ; point++ )
  {

    outputTractogram[ 3 * point + 0 ]  =
                                      flip_x * inputTractogram[ 3 * point + 0 ]
                                                        + 2 * translation[ 0 ] ;

    outputTractogram[ 3 * point + 1 ]  =
                                      flip_y * inputTractogram[ 3 * point + 1 ]
                                                        + 2 * translation[ 1 ] ;

    outputTractogram[ 3 * point + 2 ]  =
                                      flip_z * inputTractogram[ 3 * point + 2 ]
                                                        + 2 * translation[ 2 ] ;


  }

}


///////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ///////////////
///////////////////////////////////////////////////////////////////////////////

int getFlagPosition( int argc, char* argv[], const std::string& flag )
{

  std::string cmd ;
  for ( int i = 0 ; i < argc ; i++ )
  {

    std::string arg = argv[ i ] ;
    if ( 0 == arg.find( flag ) )
    {

      return i ;

    }

  }
  return 0 ;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  int index_input, index_output, index_reference, index_x, index_y, index_z,
                                                    index_verbose, index_help ;
  index_input =   getFlagPosition( argc, argv, "-i" ) ;
  index_output =   getFlagPosition( argc, argv, "-o" ) ;
  index_reference =   getFlagPosition( argc, argv, "-r" ) ;
  index_x =   getFlagPosition( argc, argv, "-x" ) ;
  index_y =   getFlagPosition( argc, argv, "-y" ) ;
  index_z =   getFlagPosition( argc, argv, "-z" ) ;
  index_verbose =   getFlagPosition( argc, argv, "-v" ) ;
  index_help =   getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : " << std::endl
              << "-i : Input file name " << std::endl
              << "-o : Out file name " << std::endl
              << "[-r] : Reference image of the bundles (only when "
              << "converting to .trk) " << std::endl
              << "[-x] : Flip around x brain axis " << std::endl
              << "[-y] : Flip around y brain axis " << std::endl
              << "[-z] : Flip around z brain axis " << std::endl
              << "[-v] : set verbosity level at 1 " << std::endl
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_reference )
  {

    std::cout << "-r argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_x )
  {

    flip_x = -1 ;

  }
  else
  {

    flip_x = 1 ;

  }

  if ( index_y )
  {

    flip_y = -1 ;

  }
  else
  {

    flip_y = 1 ;

  }

  if ( index_z )
  {

    flip_z = -1 ;

  }
  else
  {

    flip_z = 1 ;

  }

  if ( index_verbose )
  {

    verbose = 1 ;

  }

  // Checking format input reference
  char* refFilename = argv[ index_reference + 1 ] ;
  std::string refString( refFilename ) ;

  if ( refString.find( ".nii" ) == std::string::npos )
  {

    std::cout << "ERROR: Only supported format for reference image "
    << "is .nii " << std::endl ;
    exit( 1 ) ;

  }

  if ( refString.find( ".nii.gz" ) != std::string::npos )
  {

    std::cout << "ERROR: Only supported format for reference image "
    << "is .nii " << std::endl ;
    exit( 1 ) ;

  }



  // Different scenarios depending on the inputs
  std::string inputFilename ( argv[ index_input + 1 ] ) ;
  std::string outputFilename( argv[ index_output + 1 ] ) ;

  std::string bundlesdataFilename ;
  std::string bundlesMinfFilename ;
  std::string bundlesFilename ;


  if ( outputFilename.find( ".trk" ) != std::string::npos )
  {


    if ( inputFilename.find( ".bundlesdata" ) != std::string::npos )
    {

      bundlesdataFilename = inputFilename ;

      bundlesFilename = inputFilename ;
      size_t index = inputFilename.find( ".bundlesdata" ) ;
      bundlesFilename.replace( index, 12, ".bundles" ) ;

    }
    else if ( inputFilename.find( ".bundles" ) != std::string::npos )
    {

      bundlesFilename = inputFilename ;

      bundlesdataFilename = inputFilename ;
      size_t index = inputFilename.find( ".bundles" ) ;
      bundlesdataFilename.replace( index, 8, ".bundlesdata" ) ;

    }
    else
    {

      std::cout << "ERROR: Only conversions supported .trk -> .bundles "
      << "(.bundlesdata) and .bundles (.bundlesdata) -> .trk \n" ;
      exit( 1 ) ;

    }

    std::cout << "Conversion .bundlesdata (.bundles) -> .trk ... \n" ;

    bundlesData2Trk( bundlesdataFilename.c_str(),
                     bundlesFilename.c_str(),
                     outputFilename.c_str(),
                     refFilename,
                     flip_x,
                     flip_y,
                     flip_z ) ;

    return 0 ;


  }
  else if ( outputFilename.find( ".bundles" ) != std::string::npos ||
                 outputFilename.find( ".bundlesdata" )  != std::string::npos )
  {

    if ( outputFilename.find( ".bundlesdata" ) != std::string::npos )
    {

      bundlesdataFilename = outputFilename ;

      bundlesFilename = outputFilename ;
      size_t index = outputFilename.find( ".bundlesdata" ) ;
      bundlesFilename.replace( index, 12, ".bundles" ) ;

    }
    else if ( outputFilename.find( ".bundles" ) != std::string::npos )
    {

      bundlesFilename = outputFilename ;

      bundlesdataFilename = outputFilename ;
      size_t index = outputFilename.find( ".bundles" ) ;
      bundlesdataFilename.replace( index, 8, ".bundlesdata" ) ;

    }

    if ( inputFilename.find( ".trk" ) != std::string::npos )
    {

      std::cout << "Conversion .trk -> .bundlesdata (.bundles) ... \n" ;

      size_t index = bundlesFilename.find( ".bundles" ) ;
      std::string minfFilename = bundlesFilename ;
      minfFilename.replace( index, 8, ".minf" ) ;


      Trk2BundlesData( inputFilename.c_str(),
                       refFilename,
                       bundlesFilename.c_str(),
                       bundlesdataFilename.c_str(),
                       flip_x,
                       flip_y,
                       flip_z ) ;

    }
    else if ( inputFilename.find( ".tck" ) != std::string::npos )
    {

      std::cout << "Conversion .tck -> .bundlesdata (.bundles) ... \n" ;

      size_t index = bundlesFilename.find( ".bundles" ) ;
      std::string minfFilename = bundlesFilename ;
      minfFilename.replace( index, 8, ".minf" ) ;


      Tck2BundlesData( inputFilename.c_str(),
                       refFilename,
                       minfFilename.c_str(),
                       bundlesFilename.c_str(),
                       bundlesdataFilename.c_str(),
                       flip_x,
                       flip_y,
                       flip_z ) ;


    }
    else if ( inputFilename.find( ".bundlemap" ) != std::string::npos ||
                inputFilename.find( ".bundlemap.minf" )  != std::string::npos )
    {

      std::cout << "Conversion .bundlemap -> .bundlesdata (.bundles) ... \n" ;

      std::string bundlesMapFilename ;
      std::string bundlesMapMinfFilename ;

      if ( inputFilename.find( ".bundlemap.minf" ) != std::string::npos )
      {

        bundlesMapMinfFilename = inputFilename ;

        bundlesMapFilename = inputFilename ;
        size_t index = inputFilename.find( ".bundlemap.minf" ) ;
        bundlesMapFilename.replace( index, 15, ".bundlemap" ) ;

      }
      else if ( inputFilename.find( ".bundlemap" ) != std::string::npos )
      {

        bundlesMapFilename = inputFilename ;

        bundlesMapMinfFilename = inputFilename ;
        size_t index = inputFilename.find( ".bundlemap" ) ;
        bundlesMapMinfFilename.replace( index, 10, ".bundlemap.minf" ) ;

      }

      BundlesMap2BundlesData( bundlesMapFilename.c_str(),
                              bundlesMapMinfFilename.c_str(),
                              bundlesFilename.c_str(),
                              bundlesdataFilename.c_str() ) ;

    }

    return 0 ;

  }
  else if ( outputFilename.find( ".tck" ) != std::string::npos )
  {

    if ( inputFilename.find( ".bundlesdata" ) != std::string::npos )
    {

      bundlesdataFilename = inputFilename ;

      bundlesFilename = inputFilename ;
      bundlesMinfFilename = inputFilename ;
      size_t index = inputFilename.find( ".bundlesdata" ) ;
      bundlesFilename.replace( index, 12, ".bundles" ) ;
      bundlesMinfFilename.replace( index, 12, ".minf" ) ;

    }
    else if ( inputFilename.find( ".bundles" ) != std::string::npos )
    {

      bundlesFilename = inputFilename ;

      bundlesdataFilename = inputFilename ;
      bundlesMinfFilename = inputFilename ;
      size_t index = inputFilename.find( ".bundles" ) ;
      bundlesdataFilename.replace( index, 8, ".bundlesdata" ) ;
      bundlesMinfFilename.replace( index, 8, ".minf" ) ;

    }
    else
    {

      std::cout << "ERROR: Only conversion supported .tck -> .bundles \n " ;
      exit( 1 ) ;

    }

    std::cout << "Conversion .bundlesdata (.bundles) -> .tck ... \n" ;

    BundlesData2Tck( bundlesdataFilename.c_str(),
                     bundlesFilename.c_str(),
                     bundlesMinfFilename.c_str(),
                     refFilename,
                     outputFilename.c_str(),
                     flip_x,
                     flip_y,
                     flip_z ) ;

    return 0 ;

  }
  else
  {

    std::cout << "Only conversions supported : \n "
              << ".trk <-> .bundles/.bundlesdata \n "
              << ".bundlemap <-> .bundles/.bundlesdata"
              << ".tck -> .bundles/.bundlesdata \n " << std::endl ;
    exit( 1 ) ;

  }

}
