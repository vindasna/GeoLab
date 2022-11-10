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

#include "AtlasBundles.h"


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
AtlasBundles::AtlasBundles(){}

AtlasBundles::AtlasBundles( const char* atlasDirectory,
                            bool isBundlesFormat,
                            bool isTRKFormat,
                            int verbose )
{

  atlasReading( atlasDirectory, isBundlesFormat, isTRKFormat, verbose ) ;

}

AtlasBundles::AtlasBundles( std::vector< BundlesFormat > bundles,
                            std::vector< BundlesDataFormat > bundlesData,
                            std::vector< std::string > bundlesNames )
{

  this->bundles = bundles ;
  this->bundlesData = bundlesData ;
  this->bundlesNames = bundlesNames ;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
AtlasBundles::~AtlasBundles(){}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::atlasReading( const char* atlasDirectory,
                                 bool isBundlesFormat,
                                 bool isTRKFormat,
                                 int verbose )
{

  int sizeBytesFullAtlas = 0 ;
  std::vector< std::string > atlasTRKFilenames ;
  std::vector< std::string > atlasBundlesFilenames ;
  std::vector< std::string > atlasBundlesDataFilenames ;

  if ( isBundlesFormat )
  {

    for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                              atlasDirectory ) )
    {

      std::string tmpBundlesDataFilename ;
      std::string tmpBundlesFilename ;

      tmpBundlesDataFilename = file.path() ;
      tmpBundlesFilename = tmpBundlesDataFilename ;
      std::string key (".bundlesdata") ;

      if ( tmpBundlesDataFilename.find( ".bundlesdata" ) != std::string::npos )
      {

        std::size_t found = tmpBundlesFilename.rfind( key ) ;
        tmpBundlesFilename.replace( found, key.length(), ".bundles" ) ;

        atlasBundlesDataFilenames.push_back( tmpBundlesDataFilename ) ;
        atlasBundlesFilenames.push_back( tmpBundlesFilename ) ;

      }
      else if ( tmpBundlesDataFilename.find( ".trk" ) != std::string::npos )
      {

        std::cout << "Format error : Input tractogram is .bundles but .trk "
                  << "files found in the atlas directory. Projection between "
                  << "different formats is not supported yet \n" ;
        exit( 1 ) ;

      }

    }

    int nbBundlesDataFiles = atlasBundlesDataFilenames.size() ;
    int nbBundlesFiles = atlasBundlesFilenames.size() ;

    if ( nbBundlesDataFiles != nbBundlesFiles )
    {

      std::cout << "There has to be a .bundles file for each .bundlesdata : \n"
                << "Number of .bundles files : " << nbBundlesFiles
                << "\nNumber of .bundlesdata files : " << nbBundlesDataFiles
                << std::endl ;
      exit( 1 ) ;

    }

    this->bundles.resize( nbBundlesFiles ) ;
    this->bundlesData.resize( nbBundlesFiles ) ;
    this->bundlesNames.resize( nbBundlesFiles ) ;

    if ( verbose )
    {

      std::cout << "Reading atlas : " << atlasDirectory << std::endl ;

    }

    for ( int k = 0 ; k < nbBundlesFiles ; k++ )
    {

      if ( verbose )
      {

        printf("\rReading bundle [ %d / %d ]", k + 1 , nbBundlesFiles ) ;
        std::cout << "" << std::flush ;

      }

      std::string bundlesFilename = atlasBundlesFilenames[ k ] ;
      std::string bundlesdataFilename = atlasBundlesDataFilenames[ k ] ;

      if ( verbose > 1 )
      {

        std::cout << "\nReading : " << bundlesFilename << std::endl ;

      }


      BundlesFormat bundleInfo( bundlesFilename.c_str(), verbose ) ;

      // this->bundles.push_back( bundleInfo ) ;
      this->bundles[ k ] = bundleInfo ;

      if ( verbose > 1 )
      {

        std::cout << "Reading : " << bundlesdataFilename << std::endl ;

      }


      BundlesDataFormat bundleData( bundlesdataFilename.c_str(),
                                    bundlesFilename.c_str(),
                                    verbose - 1 ) ;


      // this->bundlesData.push_back( bundleData ) ;
      this->bundlesData[ k ] = bundleData  ;

      for ( int fiber = 0 ; fiber < this->bundlesData[ k ].curves_count ;
                                                                       fiber++ )
      {

        sizeBytesFullAtlas +=
                         3 * this->bundlesData[ k ].pointsPerTrack[ fiber ] *
                                                               sizeof( float ) ;

      }

      this->bundlesNames[ k ] = getFilenameNoExtension( bundlesFilename ) ;

    }

    if ( verbose == 1 )
    {

      std::cout << "\n" ;

    }

  }
  else if ( isTRKFormat )
  {

    std::string tmpTRKFilename ;

    for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                              atlasDirectory ) )
    {

      tmpTRKFilename = file.path() ;

      if ( tmpTRKFilename.find( ".trk" ) != std::string::npos )
      {

        atlasTRKFilenames.push_back( tmpTRKFilename ) ;

      }

    }


    int nbTRKFiles = atlasTRKFilenames.size() ;

    if ( nbTRKFiles == 0 )
    {

      std::cout << "ERROR : problem reading TRK atlas, no .trk files found \n" ;
      exit( 1 ) ;

    }


    this->bundles.resize( nbTRKFiles ) ;
    this->bundlesData.resize( nbTRKFiles ) ;
    this->bundlesNames.resize( nbTRKFiles ) ;

    if ( verbose )
    {

      std::cout << "Reading atlas : " << atlasDirectory << std::endl ;
    }


    for ( int k = 0 ; k < nbTRKFiles ; k++ )
    {

      if ( verbose )
      {

        printf("\rReading bundle [ %d / %d ]", k + 1 , nbTRKFiles ) ;
        std::cout << "" << std::flush ;

      }

      std::string trkFilename = atlasTRKFilenames[ k ] ;

      if ( verbose > 1 )
      {

        std::cout << "\nReading : " << trkFilename << std::endl ;

      }


      TrkFormat trkData( trkFilename.c_str(), verbose ) ;

      trkData.toBundles( this->bundles[ k ], this->bundlesData[ k ] ) ;

      for ( int fiber = 0 ; fiber < this->bundlesData[ k ].curves_count ;
                                                                      fiber++ )
      {

        sizeBytesFullAtlas +=
                        3 * this->bundlesData[ k ].pointsPerTrack[ fiber ] *
                                                               sizeof( float ) ;

      }

      this->bundlesNames[ k ] = getFilenameNoExtension( trkFilename ) ;

    }

    if ( verbose == 1 )
    {

      std::cout << "\n" ;

    }

  }

}


////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::atlasWriting( const char* outDirectory,
                                 bool isBundlesFormat,
                                 bool isTRKFormat,
                                 int verbose )
{

  char lastChar = outDirectory[
                                 ( ( std::string )outDirectory ).size()  - 1 ] ;

  int sizeAtlas = bundles.size() ;

  std::cout << "Saving atlas in " << outDirectory << std::endl ;

  for ( int i = 0 ; i < sizeAtlas ; i++ )
  {

    if ( isBundlesFormat )
    {

      printf( "\rSaving bundles [ %d / %d ]", i + 1, sizeAtlas ) ;

      std::string bundlesOutFile = outDirectory ;
      std::string bundlesDataOutFile = outDirectory ;

      if ( lastChar != '/' )
      {

        bundlesOutFile += "/" ;
        bundlesDataOutFile += "/" ;

      }

      bundlesOutFile += bundlesNames[ i ] + ".bundles" ;
      bundlesDataOutFile += bundlesNames[ i ] + ".bundlesdata" ;

      bundles[ i ].bundlesWriting( bundlesOutFile.c_str(), verbose ) ;
      bundlesData[ i ].bundlesdataWriting( bundlesDataOutFile.c_str(), verbose ) ;

    }
    else if ( isTRKFormat )
    {

      // TO DO
      std::cout << "Writing .trk atlas is not supported yet" << std::endl ;

    }

  }
  std::cout << "\nDone" << std::endl ;

}

////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::analysisWriting(
                         const char* analysisFilename,
                         std::vector<float>& distancesToCenterAtlasBundle,
                         std::vector<float>& lengthsAtlasBundle,
                         std::vector<float>& disimilaritiesAtlasBundle,
                         std::vector<float>& anglesAtlasBundle,
                         std::vector<float>& directionAnglesAtlasBundle,
                         std::vector<float>& shapeAnglesAtlasBundle,
                         std::vector<float>& distancesBetweenMedialPointsBundle,
                         int nbFibersAtlasBundle,
                         int numberDisimilaritiesAtlasBundle,
                         int numberDistancesBetweenMedialPoints,
                         int verbose )
{

  if ( verbose > 1 )
  {

    std::cout << "Writing " << analysisFilename << std::endl ;

  }

  std::ofstream file ;
  file.open( analysisFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << analysisFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  file.write( reinterpret_cast<char*>( &nbFibersAtlasBundle ), sizeof( int ) ) ;

  file.write( reinterpret_cast<char*>( &numberDisimilaritiesAtlasBundle ),
                                                               sizeof( int ) ) ;

  file.write( reinterpret_cast<char*>( &numberDistancesBetweenMedialPoints ),
                                                               sizeof( int ) ) ;

    for ( int i = 0 ; i < nbFibersAtlasBundle ; i++ )
    {


      file.write( reinterpret_cast<char*>(
                   &( distancesToCenterAtlasBundle[ i ] ) ), sizeof( float ) ) ;

    }

    for ( int i = 0 ; i < nbFibersAtlasBundle ; i++ )
    {


      file.write( reinterpret_cast<char*>( &( lengthsAtlasBundle[ i ] ) ),
                                                             sizeof( float ) ) ;

    }

    for ( int i = 0 ; i < numberDisimilaritiesAtlasBundle ; i++ )
    {


      file.write( reinterpret_cast<char*>(
                      &( disimilaritiesAtlasBundle[ i ] ) ), sizeof( float ) ) ;

    }

    for ( int i = 0 ; i < numberDisimilaritiesAtlasBundle ; i++ )
    {

      file.write( reinterpret_cast<char*>( &( anglesAtlasBundle[ i ] ) ),
                                                             sizeof( float ) ) ;

    }

    for ( int i = 0 ; i < numberDisimilaritiesAtlasBundle ; i++ )
    {

      file.write( reinterpret_cast<char*>(
                     &( directionAnglesAtlasBundle[ i ] ) ), sizeof( float ) ) ;

    }

    for ( int i = 0 ; i < nbFibersAtlasBundle ; i++ )
    {

      if ( verbose > 1 )
      {

        std::cout << "Shape angle of fiber " << i << " : "
                  << shapeAnglesAtlasBundle[ i ] << std::endl ;

      }

      file.write( reinterpret_cast<char*>( &( shapeAnglesAtlasBundle[ i ] ) ),
                                                             sizeof( float ) ) ;

    }

    for ( int i = 0 ; i < numberDistancesBetweenMedialPoints ; i++ )
    {


      file.write( reinterpret_cast<char*>(
                                 &( distancesBetweenMedialPointsBundle[ i ] ) ),
                                                             sizeof( float ) ) ;

    }

    file.close() ;

}

////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::buildFullAtlas( BundlesFormat& bundlesFullAtlas,
                                   BundlesDataFormat& bundlesDataFullAtlas )
{

  int64_t sizeFullAtlas = 0 ;
  int sizeAtlas = this->bundles.size() ;

  bundlesDataFullAtlas.curves_count = 0 ;

  for ( int bundleIndex = 0 ; bundleIndex < sizeAtlas ; bundleIndex++ )
  {

    int nbFibersBundle = this->bundlesData[ bundleIndex ].curves_count ;
    bundlesDataFullAtlas.curves_count += nbFibersBundle ;
    for ( int fiber = 0 ; fiber < nbFibersBundle ; fiber++ )
    {

      sizeFullAtlas +=
                  3 * this->bundlesData[ bundleIndex ].pointsPerTrack[ fiber ] ;

    }

  }

  bundlesDataFullAtlas.matrixTracks.resize( sizeFullAtlas, 0 ) ;

  bundlesDataFullAtlas.pointsPerTrack.resize(
                                        bundlesDataFullAtlas.curves_count, 0 ) ;

  int fiberIndexFullAtlas = 0 ;
  int64_t offsetFullAtlas = 0 ;
  for ( int bundleIndex = 0 ; bundleIndex < sizeAtlas ; bundleIndex++ )
  {

    int tmpCurveCount = this->bundlesData[ bundleIndex ].curves_count ;

    for ( int fiber = 0 ; fiber < tmpCurveCount ; fiber++ )
    {

      bundlesDataFullAtlas.pointsPerTrack[ fiberIndexFullAtlas ] =
                      this->bundlesData[ bundleIndex ].pointsPerTrack[ fiber ] ;

      fiberIndexFullAtlas += 1 ;

    }

    std::copy( this->bundlesData[ bundleIndex ].matrixTracks.begin(),
               this->bundlesData[ bundleIndex ].matrixTracks.end(),
               bundlesDataFullAtlas.matrixTracks.begin() + offsetFullAtlas ) ;

    offsetFullAtlas += ( int )( 3 *
        this->bundlesData[ bundleIndex ].pointsPerTrack[ 0 ] * tmpCurveCount ) ;

  }

  bundlesFullAtlas.binary = this->bundles[ 0 ].binary ;
  bundlesFullAtlas.byte_order = this->bundles[ 0 ].byte_order ;
  bundlesFullAtlas.data_file_name = this->bundles[ 0 ].data_file_name ;
  bundlesFullAtlas.format = this->bundles[ 0 ].format ;
  bundlesFullAtlas.object_type = this->bundles[ 0 ].object_type ;
  bundlesFullAtlas.space_dimension = this->bundles[ 0 ].space_dimension ;
  bundlesFullAtlas.bundles = ( std::string )"[ 'Full_Atlas', 0 ]" ;
  bundlesFullAtlas.curves_count = bundlesDataFullAtlas.curves_count ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    bundlesFullAtlas.resolution[ i ] = this->bundles[ 0 ].resolution[ i ] ;
    bundlesFullAtlas.size[ i ] = this->bundles[ 0 ].size[ i ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::computeLengthsAtlasBundleFibers(
                                  BundlesDataFormat& atlasBundleData,
                                  int nbPoints,
                                  std::vector<float>& lengthsAtlasBundleFibers )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {


    float lengthFiber = atlasBundleData.computeLengthFiber(
                                                   atlasBundleData.matrixTracks,
                                                   atlasBundleFiberIndex,
                                                   nbPoints ) ;

    lengthsAtlasBundleFibers[ atlasBundleFiberIndex ] = lengthFiber ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::findCenterBundle(
                              BundlesDataFormat& atlasBundleData,
                              int nbPoints,
                              std::vector<float>& medialPointAtlasBundle,
                              std::vector<float>& medialPointAtlasBundleFibers )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    // Searching the medial point of atlas fiber
    std::vector<float> medialPointFiber( 3, 0 ) ;
    atlasBundleData.computeMedialPointFiberWithDistance( atlasBundleFiberIndex,
                                                         medialPointFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointAtlasBundleFibers[ 3 * atlasBundleFiberIndex + i ] =
                                                         medialPointFiber[ i ] ;
      medialPointAtlasBundle[ i ] += medialPointFiber[ i ] ;

    }

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    medialPointAtlasBundle[ i ] /= nbFibersAtlasBundle ;

  }

}

// -------------------------------------------------------------------------- //
void AtlasBundles::findCenterBundle(
                              int bundleIndex,
                              std::vector<float>& medialPointAtlasBundle,
                              std::vector<float>& medialPointAtlasBundleFibers )
{


    int nbFibersAtlasBundle = this->bundlesData[ bundleIndex ].curves_count ;
    BundlesDataFormat& atlasBundleData = this->bundlesData[ bundleIndex ] ;
    int nbPoints = this->bundlesData[ bundleIndex ].pointsPerTrack[ 0 ] ;

    findCenterBundle( atlasBundleData,
                      nbPoints,
                      medialPointAtlasBundle,
                      medialPointAtlasBundleFibers ) ;

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::computeNormalVectorFibersAtlasBundle(
                         BundlesDataFormat& atlasBundleData,
                         int nbPoints,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         std::vector<float>& normalVectors )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    std::vector<float> medialPointFiberAtlas( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointFiberAtlas[ i ] = medialPointAtlasBundleFibers[
                                               3 * atlasBundleFiberIndex + i ] ;

    }

    std::vector<float> normalVector( 3, 0 ) ;
    atlasBundleData.computeNormalVectorFiberTractogram(
                                                  atlasBundleData.matrixTracks,
                                                  medialPointFiberAtlas,
                                                  atlasBundleFiberIndex,
                                                  nbPoints,
                                                  normalVector ) ;

    normalVectors[ 3 * atlasBundleFiberIndex + 0 ] = normalVector[ 0 ] ;
    normalVectors[ 3 * atlasBundleFiberIndex + 1 ] = normalVector[ 1 ] ;
    normalVectors[ 3 * atlasBundleFiberIndex + 2 ] = normalVector[ 2 ] ;

  }

}

////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::computeDirectionVectorFibersAtlasBundle(
                        BundlesDataFormat& atlasBundleData,
                        const std::vector<float>& normalVectorsAtlasBundle,
                        int nbPoints,
                        const std::vector<float>& medialPointAtlasBundleFibers,
                        std::vector<float>& directionVectors )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                       atlasBundleFiberIndex++ )
  {

    std::vector<float> medialPointFiberAtlas( 3, 0 ) ;
    std::vector<float> normalVectorFiberAtlas( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointFiberAtlas[ i ] = medialPointAtlasBundleFibers[
                                               3 * atlasBundleFiberIndex + i ] ;
      normalVectorFiberAtlas[ i ] = normalVectorsAtlasBundle[
                                               3 * atlasBundleFiberIndex + i ] ;

    }

    std::vector<float> directionVector( 3, 0 ) ;
    atlasBundleData.computeDirectionVectorFiberTractogram(
                                                   atlasBundleData.matrixTracks,
                                                   medialPointFiberAtlas,
                                                   normalVectorFiberAtlas,
                                                   atlasBundleFiberIndex,
                                                   nbPoints,
                                                   directionVector ) ;

    directionVectors[ 3 * atlasBundleFiberIndex + 0 ] = directionVector[ 0 ] ;
    directionVectors[ 3 * atlasBundleFiberIndex + 1 ] = directionVector[ 1 ] ;
    directionVectors[ 3 * atlasBundleFiberIndex + 2 ] = directionVector[ 2 ] ;

  }

}
////////////////////////////////////////////////////////////////////////////////
double AtlasBundles::compareDisimilarityBundles(
                         BundlesDataFormat& bundle1,
                         BundlesDataFormat& bundle2,
                         const std::vector<float>& medialPointAtlasBundleFibers,
                         int nbPoints )
{

  int nbFibersBundle1 = bundle1.curves_count ;
  int nbFibersBundle2 = bundle2.curves_count ;

  double disimilarity = 0 ;

  #pragma omp parallel for reduction( + : disimilarity )
  for ( int fiberIndex1 = 0 ; fiberIndex1 < nbFibersBundle1 ; fiberIndex1++ )
  {

    // Searching the medial point of tractogram 1 fiber
    std::vector<float> medialPointTractFiber1( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointTractFiber1[ i ] =  medialPointAtlasBundleFibers[
                                                         3 * fiberIndex1 + i ] ;

    }


    float minDMDA = 1000 ;

    for ( int fiberIndex2 = 0 ; fiberIndex2 < nbFibersBundle2 ; fiberIndex2++ )
    {

      // Registering fibers
      int verbose = 0 ;
      std::vector<float> fiber2Tofiber1( 3 * nbPoints, 0 ) ;
      std::vector<float> newNormalVectorFiber2( 3, 0 ) ;
      bundle1.registerFiber( bundle1.matrixTracks,
                             bundle2.matrixTracks,
                             fiberIndex1,
                             fiberIndex2,
                             nbPoints,
                             fiber2Tofiber1,
                             newNormalVectorFiber2,
                             verbose ) ;


      // Computing MDA when distance between middle points is below threshold
      // float dMDA = bundle1.computeMDADBetweenTwoFibers( bundle1.matrixTracks,
      float dMDA = bundle1.computeMDFBetweenTwoFibers( bundle1.matrixTracks,
                                                       fiber2Tofiber1,
                                                       medialPointTractFiber1,
                                                       medialPointTractFiber1,
                                                       fiberIndex1,
                                                       0,
                                                       nbPoints ) ;

      if ( minDMDA > dMDA )
      {

        minDMDA = dMDA ;

      }

    }

    disimilarity += minDMDA ;

  }

  disimilarity /= nbFibersBundle1 ;

  return disimilarity ;

}

///                           ###################                            ///
double AtlasBundles::compareDisimilarityBundles( BundlesDataFormat& bundle1,
                                                 BundlesDataFormat& bundle2,
                                                 int nbPoints )
{

  int nbFibersBundle1 = bundle1.curves_count ;
  int nbFibersBundle2 = bundle2.curves_count ;

  double disimilarity = 0 ;

  #pragma omp parallel for reduction( + : disimilarity )
  for ( int fiberIndex1 = 0 ; fiberIndex1 < nbFibersBundle1 ; fiberIndex1++ )
  {

    // Searching the medial point of tractogram 1 fiber
    std::vector<float> medialPointTractFiber1( 3, 0 ) ;
    bundle1.computeMedialPointFiberWithDistance( fiberIndex1,
                                                 medialPointTractFiber1 ) ;


    float minDMDA = 1000 ;


    for ( int fiberIndex2 = 0 ; fiberIndex2 < nbFibersBundle2 ; fiberIndex2++ )
    {

      // Registering fibers
      int verbose = 0 ;
      std::vector<float> fiber2Tofiber1( 3 * nbPoints, 0 ) ;
      std::vector<float> newNormalVectorFiber2( 3, 0 ) ;
      bundle1.registerFiber( bundle1.matrixTracks,
                             bundle2.matrixTracks,
                             fiberIndex1,
                             fiberIndex2,
                             nbPoints,
                             fiber2Tofiber1,
                             newNormalVectorFiber2,
                             verbose ) ;

      // Computing MDA when distance between middle points is below threshold
      // float dMDA = bundle1.computeMDADBetweenTwoFibers( bundle1.matrixTracks,
      float dMDA = bundle1.computeMDFBetweenTwoFibers( bundle1.matrixTracks,
                                                       fiber2Tofiber1,
                                                       medialPointTractFiber1,
                                                       medialPointTractFiber1,
                                                       fiberIndex1,
                                                       0,
                                                       nbPoints ) ;

      if ( minDMDA > dMDA )
      {

        minDMDA = dMDA ;

      }

    }

    disimilarity += minDMDA ;

  }

  disimilarity /= nbFibersBundle1 ;

  return disimilarity ;

}

////////////////////////////////////////////////////////////////////////////////
double AtlasBundles::distanceBetweenBundles( BundlesDataFormat& bundle1,
                                             BundlesDataFormat& bundle2,
                                             int nbFibersBundle1,
                                             int nbFibersBundle2,
                                             int nbPoints )
{

  double disimilarity_1 = compareDisimilarityBundles( bundle1,
                                                      bundle2,
                                                      nbPoints ) ;

  double disimilarity_2 = compareDisimilarityBundles( bundle2,
                                                      bundle1,
                                                      nbPoints ) ;

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

///                           ###################                            ///
double AtlasBundles::distanceBetweenBundles( int bundleIndex,
                                             BundlesDataFormat& bundle,
                                             int nbFibersBundle,
                                             int nbPoints )
{

  BundlesDataFormat& bundle1 = this->bundlesData[ bundleIndex ] ;
  int nbFibersBundle1 = this->bundlesData[ bundleIndex ].curves_count ;
  int nbPoints1 = this->bundlesData[ bundleIndex ].pointsPerTrack[ 0 ] ;

  if ( nbPoints != nbPoints1 )
  {

    std::cout << "ERROR : Problem in compareDisimilarityBundles, not the same "
              << "number of points per track in each bundle " << std::endl ;
    exit( 1 ) ;

  }

  BundlesDataFormat& bundle2 = bundle ;
  int nbFibersBundle2 = nbFibersBundle ;

  double disimilarity = distanceBetweenBundles( bundle1,
                                                bundle2,
                                                nbFibersBundle1,
                                                nbFibersBundle2,
                                                nbPoints ) ;

  return disimilarity ;

}


////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::computeNumberAdjacentFibersBundle1ToBundle2(
                                      const std::vector<float>& bundle1,
                                      const std::vector<float>& bundle2,
                                      int nbFibersBundle1,
                                      int nbFibersBundle2,
                                      int nbPoints, // Same for 2 bundles
                                      float threshold,
                                      std::vector<int>& nbAdjacentFibersBundle )
{

  if ( threshold <= 0 )
  {
    std::cout << "Warning : the threshold value in bundlesAdjacency method "
              << "must be greater than 0 but got " << threshold
              << "... Using default value of 10 " << std::endl ;

    threshold = 10 ;

  }

  // #pragma omp parallel for reduction( + : nbAdjacentFibersBundle )
  // No need for reduction because there is no race conditions on
  // nbAdjacentFibersBundle[ fiberIndex1 ] givent that it is modified at
  // the interior of other loops and the value of
  // nbAdjacentFibersBundle[ fiberIndex1 ] is not modified for a fiber
  // greater or less to fiberIndex1
  #pragma omp parallel for shared(nbAdjacentFibersBundle)
  for ( int fiberIndex1 = 0 ; fiberIndex1 < nbFibersBundle1 ; fiberIndex1++ )
  {

    // Getting fiber
    std::vector<float> fiber1( 3 * nbPoints, 0 ) ;
    this->bundlesData[ 0 ].getFiberFromTractogram( bundle1,
                                                   fiberIndex1,
                                                   nbPoints,
                                                   fiber1 ) ;

    // Searching the medial point of tractogram 1 fiber
    std::vector<float> medialPointTractFiber1( 3, 0 ) ;
    this->bundlesData[ 0 ].computeMedialPointFiberWithDistance(
                                                      fiber1,
                                                      medialPointTractFiber1 ) ;

    float minDMDA = 1000 ;
    nbAdjacentFibersBundle[ fiberIndex1 ] = 0 ;

    for ( int fiberIndex2 = 0 ; fiberIndex2 < nbFibersBundle2 ; fiberIndex2++ )
    {

      // Getting fiber
      std::vector<float> fiber2( 3 * nbPoints, 0 ) ;
      this->bundlesData[ 0 ].getFiberFromTractogram( bundle2,
                                                     fiberIndex2,
                                                     nbPoints,
                                                     fiber2 ) ;

      // Searching the medial point of tractogram 2 fiber
      std::vector<float> medialPointTractFiber2( 3, 0 ) ;
      this->bundlesData[ 0 ].computeMedialPointFiberWithDistance(
                                                      fiber2,
                                                      medialPointTractFiber2 ) ;


      // Computing MDA when distance between middle points is below threshold
      // float dMDA = this->bundlesData[ 0 ].computeMDADBetweenTwoFibers(
      float dMDA = this->bundlesData[ 0 ].computeMDFBetweenTwoFibers(
                                                         fiber1,
                                                         fiber2,
                                                         medialPointTractFiber1,
                                                         medialPointTractFiber2,
                                                         0,
                                                         0,
                                                         nbPoints ) ;

      if ( dMDA < threshold )
      {

        nbAdjacentFibersBundle[ fiberIndex1 ] += 1 ;

      }

    }

  }

}



////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::computeNumberAdjacentFibersRecognizedToAtlasBundles(
                                      BundlesDataFormat& bundle,
                                      std::string bundleName,
                                      float threshold,
                                      std::vector<int>& nbAdjacentFibersBundle )
{

  if ( threshold <= 0 )
  {
    std::cout << "Warning : the threshold value in bundlesAdjacency method "
              << "must be greater than 0 but got " << threshold
              << "... Using default value of 10 " << std::endl ;

    threshold = 10 ;

  }



  int indexBundleAtlas = this->findBundleIndexByName( bundleName ) ;

  if ( indexBundleAtlas == -1 )
  {

    std::cout << "Warning : Bundle with name " << bundleName << " not found "
              << "in atlas, cannot compute number of adjacent fibers \n" ;
    return ;

  }

  BundlesDataFormat& bundle1 = bundle ;
  int nbFibersBundle1 = bundle.curves_count ;
  int nbPoints1 = bundle.pointsPerTrack[ 0 ] ;

  BundlesDataFormat& bundle2 = this->bundlesData[ indexBundleAtlas ] ;
  int nbFibersBundle2 = this->bundlesData[ indexBundleAtlas ].curves_count ;
  int nbPoints2 = this->bundlesData[ indexBundleAtlas ].pointsPerTrack[ 0 ] ;

  if ( nbPoints2 != nbPoints1 )
  {

    std::cout << "ERROR : Problem in compareDisimilarityBundles, not the same "
              << "number of points per track in each bundle " << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPoints1 ;

  computeNumberAdjacentFibersBundle1ToBundle2( bundle1.matrixTracks,
                                               bundle2.matrixTracks,
                                               nbFibersBundle1,
                                               nbFibersBundle2,
                                               nbPoints,
                                               threshold,
                                               nbAdjacentFibersBundle ) ;

}

////////////////////////////////////////////////////////////////////////////////
void AtlasBundles::computeNumberAdjacentFibersAtlasToRecognizedBundles(
                                      BundlesDataFormat& bundle,
                                      std::string bundleName,
                                      float threshold,
                                      std::vector<int>& nbAdjacentFibersBundle )
{

  if ( threshold <= 0 )
  {
    std::cout << "Warning : the threshold value in bundlesAdjacency method "
              << "must be greater than 0 but got " << threshold
              << "... Using default value of 10 " << std::endl ;

    threshold = 10 ;

  }

  int indexBundleAtlas = this->findBundleIndexByName( bundleName ) ;

  if ( indexBundleAtlas == -1 )
  {

    std::cout << "Warning : Bundle with name " << bundleName << " not found "
              << "in atlas, cannot compute number of adjacent fibers \n" ;
    return ;

  }

  BundlesDataFormat& bundle1 =
                            this->bundlesData[ indexBundleAtlas ] ;
  int nbFibersBundle1 = this->bundlesData[ indexBundleAtlas ].curves_count ;
  int nbPoints1 = this->bundlesData[ indexBundleAtlas ].pointsPerTrack[ 0 ] ;

  BundlesDataFormat& bundle2 = bundle ;
  int nbFibersBundle2 = bundle.curves_count ;
  int nbPoints2 = bundle.pointsPerTrack[ 0 ] ;

  if ( nbPoints2 != nbPoints1 )
  {

    std::cout << "ERROR : Problem in compareDisimilarityBundles, not the same "
              << "number of points per track in each bundle " << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPoints1 ;

  computeNumberAdjacentFibersBundle1ToBundle2( bundle1.matrixTracks,
                                               bundle2.matrixTracks,
                                               nbFibersBundle1,
                                               nbFibersBundle2,
                                               nbPoints,
                                               threshold,
                                               nbAdjacentFibersBundle ) ;

}


////////////////////////////////////////////////////////////////////////////////
float AtlasBundles::coverageBundle1ToBundle2(
                                             const std::vector<float>& bundle1,
                                             const std::vector<float>& bundle2,
                                             int nbFibersBundle1,
                                             int nbFibersBundle2,
                                             int nbPoints, // Same for 2 bundles
                                             float threshold,
                                             int verbose )
{


  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle1, 0 ) ;

  this->computeNumberAdjacentFibersBundle1ToBundle2( bundle1,
                                                     bundle2,
                                                     nbFibersBundle1,
                                                     nbFibersBundle2,
                                                     nbPoints,
                                                     threshold,
                                                     nbAdjacentFibersBundle ) ;

  float coverage = 0 ;
  for ( int i = 0 ; i < nbFibersBundle1 ; i++ )
  {

    if ( nbAdjacentFibersBundle[ i ] > 0 )
    {

      coverage += 1.0 ;

    }

  }

  coverage /= nbFibersBundle1 ;

  return coverage ;

}

//----------------------------------------------------------------------------//
float AtlasBundles::coverageBundle1ToBundle2(
                      const std::vector<int>& nbAdjacentFibersBundle1ToBundle2 )
{

  int nbFibersBundle1 = nbAdjacentFibersBundle1ToBundle2.size() ;

  float coverage = 0 ;
  for ( int i = 0 ; i < nbFibersBundle1 ; i++ )
  {

    if ( nbAdjacentFibersBundle1ToBundle2[ i ] > 0 )
    {

      coverage += 1.0 ;

    }

  }

  coverage /= nbFibersBundle1 ;

  return coverage ;

}

////////////////////////////////////////////////////////////////////////////////
float AtlasBundles::coverageRecognizedToAtlasBundles( BundlesDataFormat& bundle,
                                                      std::string bundleName,
                                                      float threshold,
                                                      int verbose )
{

  int nbFibersBundle = bundle.curves_count ;

  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersRecognizedToAtlasBundles(
                                                      bundle,
                                                      bundleName,
                                                      threshold,
                                                      nbAdjacentFibersBundle ) ;

  float coverage = 0 ;
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

//----------------------------------------------------------------------------//

float AtlasBundles::coverageRecognizedToAtlasBundles(
              const std::vector<int>& nbAdjacentFibersRecognizedToAtlasBundles )
{

  int nbFibersBundle = nbAdjacentFibersRecognizedToAtlasBundles.size() ;

  float coverage = 0 ;
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersRecognizedToAtlasBundles[ i ] > 0 )
    {

      coverage += 1.0 ;

    }

  }

  coverage /= nbFibersBundle ;

  return coverage ;

}


////////////////////////////////////////////////////////////////////////////////
float AtlasBundles::coverageAtlasToRecognizedBundles( BundlesDataFormat& bundle,
                                                      std::string bundleName,
                                                      float threshold,
                                                      int verbose )
{

  int indexBundleAtlas = this->findBundleIndexByName( bundleName ) ;
  int nbFibersBundle = this->bundlesData[ indexBundleAtlas ].curves_count ;

  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersAtlasToRecognizedBundles(
                                                      bundle,
                                                      bundleName,
                                                      threshold,
                                                      nbAdjacentFibersBundle ) ;

  float coverage = 0 ;
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

//----------------------------------------------------------------------------//

float AtlasBundles::coverageAtlasToRecognizedBundles(
              const std::vector<int>& nbAdjacentFibersAtlasToRecognizedBundles )
{

  int nbFibersBundle = nbAdjacentFibersAtlasToRecognizedBundles.size() ;

  float coverage = 0 ;
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersAtlasToRecognizedBundles[ i ] > 0 )
    {

      coverage += 1.0 ;

    }

  }

  coverage /= nbFibersBundle ;

  return coverage ;

}

////////////////////////////////////////////////////////////////////////////////
float AtlasBundles::overlapBundle1ToBundle2(
                                             const std::vector<float>& bundle1,
                                             const std::vector<float>& bundle2,
                                             int nbFibersBundle1,
                                             int nbFibersBundle2,
                                             int nbPoints, // Same for 2 bundles
                                             float threshold,
                                             int verbose )
{

  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle1, 0 ) ;

  this->computeNumberAdjacentFibersBundle1ToBundle2( bundle1,
                                                     bundle2,
                                                     nbFibersBundle1,
                                                     nbFibersBundle2,
                                                     nbPoints,
                                                     threshold,
                                                     nbAdjacentFibersBundle ) ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
  for ( int i = 0 ; i < nbFibersBundle1 ; i++ )
  {

    if ( nbAdjacentFibersBundle[ i ] > 0 )
    {

      nbAdjacentFibersRecognizedBundle += 1 ;
      nbAdjacentFibersAtlasBundle += nbAdjacentFibersBundle[ i ] ;

    }

  }

  if ( nbAdjacentFibersAtlasBundle == 0 )
  {

    if ( verbose > 1 )
    {

      std::cout << "Warning : there are not adjacent fibers of bundle1 "
                << "to bundle2, overlap is not defined \n" ;

    }

    return 0 ;

  }

  // Overlap is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;


}

//----------------------------------------------------------------------------//

float AtlasBundles::overlapBundle1ToBundle2(
                      const std::vector<int>& nbAdjacentFibersBundle1ToBundle2 )
{

  int nbFibersBundle1 = nbAdjacentFibersBundle1ToBundle2.size() ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
  for ( int i = 0 ; i < nbFibersBundle1 ; i++ )
  {

    if ( nbAdjacentFibersBundle1ToBundle2[ i ] > 0 )
    {

      nbAdjacentFibersRecognizedBundle += 1 ;
      nbAdjacentFibersAtlasBundle += nbAdjacentFibersBundle1ToBundle2[ i ] ;

    }

  }

  if ( nbAdjacentFibersAtlasBundle == 0 )
  {


    std::cout << "Warning : there are not adjacent fibers of bundle1 "
              << "to bundle2, overlap is not defined \n" ;


    return 0 ;

  }

  // Overlap is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;


}


////////////////////////////////////////////////////////////////////////////////
float AtlasBundles::overlapRecognizedToAtlasBundles( BundlesDataFormat& bundle,
                                                     std::string bundleName,
                                                     float threshold,
                                                     int verbose )
{

  int nbFibersBundle = bundle.curves_count ;
  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersRecognizedToAtlasBundles(
                                                      bundle,
                                                      bundleName,
                                                      threshold,
                                                      nbAdjacentFibersBundle ) ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
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

    if ( verbose > 1 )
    {

      std::cout << "Warning : there are not adjacent fibers of recognized "
                << bundleName << " to atlas, overlap is not defined \n" ;

    }

    return 0 ;

  }

  // Overlap is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;


}


//----------------------------------------------------------------------------//


float AtlasBundles::overlapRecognizedToAtlasBundles(
              const std::vector<int>& nbAdjacentFibersRecognizedToAtlasBundles )
{

  int nbFibersBundle = nbAdjacentFibersRecognizedToAtlasBundles.size() ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersRecognizedToAtlasBundles[ i ] > 0 )
    {

      nbAdjacentFibersRecognizedBundle += 1 ;
      nbAdjacentFibersAtlasBundle +=
                                 nbAdjacentFibersRecognizedToAtlasBundles[ i ] ;

    }

  }

  if ( nbAdjacentFibersAtlasBundle == 0 )
  {


    std::cout << "Warning : there are not adjacent fibers of recognized to "
              << "atlas, overlap is not defined \n" ;


    return 0 ;

  }

  // Overlap is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;


}


////////////////////////////////////////////////////////////////////////////////
float AtlasBundles::overlapAtlasToRecognizedBundles( BundlesDataFormat& bundle,
                                                     std::string bundleName,
                                                     float threshold,
                                                     int verbose )
{

  int indexBundleAtlas = this->findBundleIndexByName( bundleName ) ;
  int nbFibersBundle = this->bundlesData[ indexBundleAtlas ].curves_count ;
  std::vector<int> nbAdjacentFibersBundle( nbFibersBundle, 0 ) ;

  this->computeNumberAdjacentFibersAtlasToRecognizedBundles(
                                                       bundle,
                                                       bundleName,
                                                       threshold,
                                                       nbAdjacentFibersBundle ) ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
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

    if ( verbose > 1 )
    {

      std::cout << "Warning : there are not adjacent fibers of recognized "
                << bundleName << " to atlas, overlap is not defined \n" ;

    }

    return 0 ;

  }

  // Overla is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;


}


//----------------------------------------------------------------------------//


float AtlasBundles::overlapAtlasToRecognizedBundles(
              const std::vector<int>& nbAdjacentFibersAtlasToRecognizedBundles )
{

  int nbFibersBundle = nbAdjacentFibersAtlasToRecognizedBundles.size() ;

  int nbAdjacentFibersRecognizedBundle = 0 ;
  int nbAdjacentFibersAtlasBundle = 0 ;
  for ( int i = 0 ; i < nbFibersBundle ; i++ )
  {

    if ( nbAdjacentFibersAtlasToRecognizedBundles[ i ] > 0 )
    {

      nbAdjacentFibersRecognizedBundle += 1 ;
      nbAdjacentFibersAtlasBundle +=
                                 nbAdjacentFibersAtlasToRecognizedBundles[ i ] ;

    }

  }

  if ( nbAdjacentFibersAtlasBundle == 0 )
  {


    std::cout << "Warning : there are not adjacent fibers of recognized to "
              << "atlas, overlap is not defined \n" ;


    return 0 ;

  }

  // Overla is greater or equal to 1
  float overlap = ( float )nbAdjacentFibersAtlasBundle /
                                    ( float )nbAdjacentFibersRecognizedBundle  ;

  return overlap ;


}

////////////////////////////////////////////////////////////////////////////////
float AtlasBundles::bundlesAdjacency( const std::vector<float>& bundle1,
                                      const std::vector<float>& bundle2,
                                      int nbFibersBundle1,
                                      int nbFibersBundle2,
                                      int nbPoints, // Same for 2 bundles
                                      float threshold,
                                      int verbose )
{

  float coverageRecognizedToAtlasBundlesMeasure =
                           this->coverageBundle1ToBundle2( bundle1,
                                                           bundle2,
                                                           nbFibersBundle1,
                                                           nbFibersBundle2,
                                                           nbPoints,
                                                           threshold,
                                                           verbose ) ;

  float coverageAtlasToRecognizedBundlesMeasure =
                           this->coverageBundle1ToBundle2( bundle2,
                                                           bundle1,
                                                           nbFibersBundle2,
                                                           nbFibersBundle1,
                                                           nbPoints,
                                                           threshold,
                                                           verbose ) ;

  float bundlesAdjacencyMeasure = ( coverageRecognizedToAtlasBundlesMeasure +
                                 coverageAtlasToRecognizedBundlesMeasure ) / 2 ;

  return bundlesAdjacencyMeasure ;

}


//----------------------------------------------------------------------------//


float AtlasBundles::bundlesAdjacency(  BundlesDataFormat& bundle,
                                       std::string bundleName,
                                       float threshold,
                                       int verbose )
{

  float coverageRecognizedToAtlasBundlesMeasure =
                           this->coverageRecognizedToAtlasBundles( bundle,
                                                                   bundleName,
                                                                   threshold,
                                                                   verbose ) ;

  float coverageAtlasToRecognizedBundlesMeasure =
                           this->coverageAtlasToRecognizedBundles( bundle,
                                                                   bundleName,
                                                                   threshold,
                                                                   verbose ) ;

  float bundlesAdjacencyMeasure = ( coverageRecognizedToAtlasBundlesMeasure +
                                 coverageAtlasToRecognizedBundlesMeasure ) / 2 ;

  return bundlesAdjacencyMeasure ;

}


//----------------------------------------------------------------------------//


float AtlasBundles::bundlesAdjacency(
                                 float coverageRecognizedToAtlasBundlesMeasure,
                                 float coverageAtlasToRecognizedBundlesMeasure )
{

  float bundlesAdjacencyMeasure = ( coverageRecognizedToAtlasBundlesMeasure +
                                 coverageAtlasToRecognizedBundlesMeasure ) / 2 ;

  return bundlesAdjacencyMeasure ;

}


////////////////////////////////////////////////////////////////////////////////
int AtlasBundles::findBundleIndexByName( std::string bundleName )
{

  int nbBundlesAtlas = this->bundles.size() ;

  for ( int i = 0 ; i < nbBundlesAtlas ; i++ )
  {

    if ( this->bundlesNames[ i ] == bundleName )
    {

      return i ;

    }

  }

  return -1 ;

}


////////////////////////////////////////////////////////////////////////////////
////////////////// Function to get filename without extension  /////////////////
////////////////////////////////////////////////////////////////////////////////
std::string getFilenameNoExtension( std::string path )
{

  size_t sep = path.find_last_of( "\\/" ) ;
  std::string filename ;
  if ( sep != std::string::npos )
  {

    filename = path.substr( sep + 1, path.size() - sep - 1 ) ;

  }

  size_t dot = filename.find_last_of( "." ) ;

  if ( dot != std::string::npos )
  {

    filename = filename.substr( 0, dot ) ;

  }

  return filename ;

}
