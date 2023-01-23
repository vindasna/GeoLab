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

#include "processAtlasInformation.h"
#include "ioWrapper.h"
#include "niftiImage.h"

///////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ///////////////
///////////////////////////////////////////////////////////////////////////////

int getFlagPosition( int argc, char* argv[], const std::string& flag )
{

  for ( int i = 0 ; i < argc ; i++ )
  {

    std::string arg = argv[ i ] ;
    // if ( 0 == arg.find( flag ) )
    if ( arg == flag )
    {

      return i ;

    }

  }
  return 0 ;

}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------//
//-------------------------- Find radius of bundle ---------------------------//
//----------------------------------------------------------------------------//
void computeCenterAtlasBundleFibers(
                             BundlesData& atlasBundleData,
                             std::vector<float>& medialPointsAtlasBundleFibers )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  #pragma omp parallel for
  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    // Searching the medial point of atlas fiber
    std::vector<float> medialPointAtlasBundleFiber( 3, 0 ) ;
    atlasBundleData.computeMedialPointFiberWithDistance(
                                                 atlasBundleFiberIndex,
                                                 medialPointAtlasBundleFiber ) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {

      medialPointsAtlasBundleFibers[ 3 * atlasBundleFiberIndex + i ] =
                                              medialPointAtlasBundleFiber[ i ] ;

    }

  }

}

//----------------------------------------------------------------------------//
//----------------------- Compute average fiber bundle -----------------------//
//----------------------------------------------------------------------------//
void computeAverageFiberBundle(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::vector<float>& medialPointAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;


  std::vector<float> referenceFiber( 3 * nbPoints, 0 ) ;
  int fiberIndex = 0 ;
  atlasBundleData.getFiberFromTractogram( atlasBundleData.matrixTracks,
                                          fiberIndex,
                                          nbPoints,
                                          referenceFiber ) ;



  std::vector<float> centerReferenceFiber( 3, 0 ) ;
  for ( int i = 0 ; i < 3 ; i++ )
  {

    centerReferenceFiber[ i ] = medialPointsAtlasBundleFibers[ 3 * 0 + i ] ;

  }


  // #pragma omp parallel for reduction( +:averageFiber[:nbElementsAverageFiber ] )
  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    int offsetAtlas = 3 * nbPoints * atlasBundleFiberIndex ;

    std::vector<float> translation( 3, 0 ) ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      translation[ i ] = centerReferenceFiber[ i ] -
               medialPointsAtlasBundleFibers[ 3 * atlasBundleFiberIndex + i ] ;

    }

    float directDistance = 0 ;
    float indirectDistance = 0 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      directDistance += pow( referenceFiber[ 3 * 0 + i ] - ( atlasBundleData[ 3
                             * 0 + offsetAtlas + i ] + translation[ i ] ), 2 ) ;
      indirectDistance += pow( referenceFiber[ 3 * 0 + i ] - ( atlasBundleData[
            3 * ( nbPoints - 1 ) + offsetAtlas + i ] + translation[ i ] ), 2 ) ;

    }

    bool isDirectSens = true ;
    if ( directDistance > indirectDistance )
    {

      isDirectSens = false ;

    }

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      int k = point ;
      if ( !isDirectSens )
      {

        k = nbPoints - point - 1 ;

      }

      for ( int i = 0 ; i < 3 ; i++ )
      {

        averageFiber[ 3 * point + i ] +=
                                    atlasBundleData[ 3 * k + offsetAtlas + i ] ;

      }

    }

  }


  for ( int point = 0 ; point < nbPoints ; point++ )
  {

    for ( int i = 0 ; i < 3 ; i++ )
    {

      averageFiber[ 3 * point + i ] /= nbFibersAtlasBundle ;

    }

  }



  // Searching the medial point of average atlas bundle fiber
  atlasBundleData.computeMedialPointFiberWithDistance(
                                                      averageFiber,
                                                      medialPointAtlasBundle ) ;

}


//----------------------------------------------------------------------------//
//--------------------- Compute gravity center of bundle ---------------------//
//----------------------------------------------------------------------------//
void computeGravityCenterAtlasBundle(
                                  BundlesData& atlasBundleData,
                                  int nbPoints,
                                  std::vector<float>& gravityCenterAtlasBundle )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterAtlasBundle[ i ] = 0 ;

  }

  for ( int atlasBundleFiberIndex = 0 ;
                             atlasBundleFiberIndex < nbFibersAtlasBundle ;
                                                      atlasBundleFiberIndex++ )
  {

    int offsetAtlas = 3 * nbPoints * atlasBundleFiberIndex ;

    for ( int point = 0 ; point < nbPoints ; point++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        gravityCenterAtlasBundle[ i ] +=
                                    atlasBundleData[ 3 * point + offsetAtlas + i ] ;

      }

    }

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    gravityCenterAtlasBundle[ i ] /= nbFibersAtlasBundle * nbPoints ;

  }

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  int index_atlas, index_format, index_reference, index_averageFibers,
                                       index_useMDF, index_verbose, index_help ;
  index_atlas =   getFlagPosition( argc, argv, "-a") ;
  index_format =   getFlagPosition( argc, argv, "-f") ;
  index_reference =   getFlagPosition( argc, argv, "-r") ;
  index_averageFibers = getFlagPosition( argc, argv, "-af") ;
  index_useMDF = getFlagPosition( argc, argv, "-useMDF" ) ;
  index_verbose =   getFlagPosition( argc, argv, "-v") ;
  index_help =   getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-a : Directory with the atlas (one file per bundle) \n"
              << "-f : Format of atlas bundles ( optios = [ .bundles, .trk, "
              << ".tck ] ) \n"
              << "[-r] : Reference .nii image where the atlas is \n"
              << "[-af] : Output directory where to save the average fibers of "
              << "the bundles \n"
              << "[-useMDF] : Use MDF distance as disimilarity measure\n"
              << "[-v] : set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  std::string atlasDirectory ;
  if ( !index_atlas )
  {

    std::cout << "-a argument required ..." << std::endl ;
    exit( 1 ) ;

  }
  else
  {

    atlasDirectory = argv[ index_atlas + 1 ] ;
    char lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      atlasDirectory = atlasDirectory + "/" ;

    }

  }

  std::string format ;
  if ( index_format )
  {

    format = argv[ index_format + 1 ] ;

    if ( format != ".bundles" && format != ".trk" && format != ".tck" )
    {

      std::cout << "The only supported formats for the atlas bundles are "
                << ".bundles, .trk and .tck" << std::endl ;
      exit( 1 ) ;

    }

  }
  else
  {

    std::cout << "-f argument required..." << std::endl ;
    exit( 1 ) ;

  }

  std::string niiFilename ;
  if ( index_reference )
  {

    niiFilename = argv[ index_reference + 1 ] ;
    char lastChar = niiFilename[ niiFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      niiFilename = niiFilename.substr( 0, niiFilename.size() - 1 ) ;

    }

  }

  if ( index_useMDF )
  {

    useMDFDistance = true ;

  }

  if ( index_verbose )
  {
    if ( argv[ index_verbose + 1 ] )
    {

      verbose = std::stoi( argv[ index_verbose + 1 ] ) ;

    }
    else
    {

      verbose = 1 ;

    }

  }

  /////////////////////////////////////////////////////////////////////////////

  // Cheking format reference image
  if ( index_reference )
  {

    if ( niiFilename.find( ".nii") == std::string::npos ||
                            niiFilename.find( ".nii.gz") != std::string::npos )
    {

      std::cout << "Error format : the only supported format for the reference "
                << "image is .nii " << std::endl ;
      exit( 1 ) ;

    }

  }

  std::string averageFibersDirectory ;
  if ( index_averageFibers )
  {

    averageFibersDirectory = argv[ index_averageFibers + 1 ] ;

    char lastChar = averageFibersDirectory[ averageFibersDirectory.size()
                                                                      - 1 ] ;
    if ( lastChar != '/' )
    {

      averageFibersDirectory = averageFibersDirectory + "/" ;

    }

  }



  //   xxxxxxxxxxxxxxxxxxxxxxxxxx Reading Atlas xxxxxxxxxxxxxxxxxxxxxxxxxx   //

  std::vector< std::string > atlasBundlesFilenames ;
  std::string tmpBundlesFilename ;

  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                    atlasDirectory.c_str() ) )
  {

    tmpBundlesFilename = file.path() ;

    if ( endswith( tmpBundlesFilename, format) )
    {
;
      atlasBundlesFilenames.push_back( tmpBundlesFilename ) ;

    }

  }

  bool isBundles = false ;
  bool isTrk = false ;
  bool isTck = false ;

  if ( format == ".bundles" )
  {

    isBundles = true ;

  }
  else if ( format == ".trk" )
  {

    isTrk = true ;

  }
  else
  {

    isTck = true ;

  }



  bool isBundlesFormat = true ;
  bool isTRKFormat = false ;
  AtlasBundles atlasData( atlasDirectory.c_str(),
                          isBundles,
                          isTrk,
                          isTck,
                          verbose ) ;



  //////////////////////////////// Sanity cheks ////////////////////////////////
  int nbPoints = atlasData.bundlesData[ 0 ].pointsPerTrack[ 0 ] ;
  int nbBundlesAtlas = atlasData.bundlesData.size() ;


  for ( int bundleIndex = 0 ; bundleIndex < nbBundlesAtlas ; bundleIndex++ )
  {

    BundlesData& bundleFibers = atlasData.bundlesData[ bundleIndex ] ;

    int nbFibersBundle = bundleFibers.curves_count ;

    for ( int fiber = 1 ; fiber < nbFibersBundle ; fiber++ )
    {

      if ( bundleFibers.pointsPerTrack[ fiber ] != nbPoints )
      {

        std::cout << "Error atlas : the number of points in each fiber "
                  << "of the atlas has to be the same, got "
                  << bundleFibers.pointsPerTrack[ fiber ]
                  << " and " << nbPoints << " for fiber " << fiber
                  << " in bundle " << bundleIndex << " and for fiber 0 in "
                  << "bundle 0" << std::endl ;
        exit( 1 ) ;

      }

    }

  }



  /////////////////////// Computing center Atlas bundles ///////////////////////

  // Looping over bundles
  for ( int atlasBundleIndex = 0 ; atlasBundleIndex < nbBundlesAtlas ;
                                                      atlasBundleIndex++ )
  {

    BundlesMinf& atlasBundleInfo = atlasData.bundlesMinf[ atlasBundleIndex ] ;
    BundlesData& atlasBundleData = atlasData.bundlesData[ atlasBundleIndex ] ;

    if ( verbose )
    {

      printf( "\rProcessing atlas bundles : [ %d  /  %d ]",
                                        atlasBundleIndex + 1, nbBundlesAtlas ) ;
      fflush( stdout ) ;

    }

    int nbFibersAtlasBundle = atlasBundleData.curves_count ;

    // Compute center points of atlas bundle fibers
    std::vector<float> medialPointsAtlasBundleFibers(
                                                  3 * nbFibersAtlasBundle, 0 ) ;

    computeCenterAtlasBundleFibers( atlasBundleData,
                                    medialPointsAtlasBundleFibers ) ;

    std::vector<float> gravityCenterAtlasBundle( 3, 0 ) ;
    computeGravityCenterAtlasBundle( atlasBundleData,
                                     nbPoints,
                                     gravityCenterAtlasBundle ) ;

    atlasBundleInfo.centerBundle[ 0 ] = gravityCenterAtlasBundle[ 0 ] ;
    atlasBundleInfo.centerBundle[ 1 ] = gravityCenterAtlasBundle[ 1 ] ;
    atlasBundleInfo.centerBundle[ 2 ] = gravityCenterAtlasBundle[ 2 ] ;

  }


  /////////////////////////////// Saving results ///////////////////////////////
  if ( verbose )
  {

    std::cout << "\n" ;

  }


  if ( verbose )
  {

    std::cout << "Checking other information in .bundles " << std::endl ;

  }

  // Checking other information in .bundles files
  if ( index_reference )
  {

    NiftiImage niiHeader( niiFilename.c_str() ) ;
    for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
    {

      for ( int i = 0 ; i < 3 ; i++ )
      {

        if ( atlasData.bundlesMinf[ bundle ].resolution[ i ] == 0 )
        {

          atlasData.bundlesMinf[ bundle ].resolution[ i ] =
                                                     niiHeader.resolution[ i ] ;

        }

        if ( atlasData.bundlesMinf[ bundle ].size[ i ] == 0 )
        {

          atlasData.bundlesMinf[ bundle ].size[ i ] = niiHeader.size[ i ] ;

        }


      }

      if ( atlasData.bundlesMinf[ bundle ].space_dimension == 0 )
      {

        atlasData.bundlesMinf[ bundle ].space_dimension =  3 ;

      }

    }

  }

  // Saving .bundles files

  if ( verbose )
  {

    std::cout << "Saving .bundles " << std::endl ;

  }

  for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
  {

    atlasData.bundlesMinf[ bundle ].write(
                                     atlasBundlesFilenames[ bundle ].c_str() ) ;

  }

  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }

}
