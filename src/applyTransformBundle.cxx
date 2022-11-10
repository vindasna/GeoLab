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

#include <Eigen/Core>
#include <LBFGS.h>

#include "./applyTransformBundle.h"


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
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  int index_input_bundle, index_transform, index_output_bundle,
                                                     index_verbose, index_help ;

  index_input_bundle = getFlagPosition( argc, argv, "-i") ;
  index_transform = getFlagPosition( argc, argv, "-t") ;
  index_output_bundle = getFlagPosition( argc, argv, "-o") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-i : Moving bundle \n"
              << "-t : Transform \n"
              << "-o : Output bundle \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input_bundle )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_transform )
  {

    std::cout << "-t argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output_bundle )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

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

  std::string movingFilename( argv[ index_input_bundle + 1 ] ) ;
  char lastChar = movingFilename[ movingFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    movingFilename = movingFilename.substr( 0, movingFilename.size() - 1 ) ;

  }

  std::string transformFilename( argv[ index_transform + 1 ] ) ;
  lastChar = transformFilename[ transformFilename.size()  - 1 ] ;
  if ( lastChar == '/' )
  {

    transformFilename = transformFilename.substr( 0,
                                                transformFilename.size() - 1 ) ;

  }

  std::string outputFilename( argv[ index_output_bundle + 1 ] ) ;
  lastChar = outputFilename[ outputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

  }


  //   xxxxxxxxxxxxxxxxxxxxxxx Reading moving bundle xxxxxxxxxxxxxxxxxxxxxx   //

  std::string movingBundlesFilename ;
  std::string movingBundlesDataFilename ;
  bool isBundlesFormat = false ;

  if ( movingFilename.find( ".bundlesdata" ) != std::string::npos )
  {

    movingBundlesDataFilename = movingFilename ;

    movingBundlesFilename = movingFilename ;
    size_t index = movingFilename.find( ".bundlesdata" ) ;
    movingBundlesFilename.replace( index, 12, ".bundles") ;

    isBundlesFormat = true ;

  }
  else if ( movingFilename.find( ".bundles" ) != std::string::npos )
  {

    movingBundlesFilename = movingFilename ;

    movingBundlesDataFilename = movingFilename ;
    size_t index = movingFilename.find( ".bundles" ) ;
    movingBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

    isBundlesFormat = true ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesDataFormat movingBundlesData ;
  BundlesFormat movingBundleInfo ;

  if ( isBundlesFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << movingBundlesFilename << std::endl ;

    }

    movingBundleInfo.bundlesReading( movingBundlesFilename.c_str(),
                                        verbose ) ;

    if ( verbose )
    {

      std::cout << "Reading : " << movingBundlesDataFilename << std::endl ;

    }
    movingBundlesData.bundlesdataReading( movingBundlesDataFilename.c_str(),
                                        movingBundlesFilename.c_str(),
                                        verbose ) ;

  }

  //   xxxxxxxxxxxxxxxxxxxxxxxx Reading transform xxxxxxxxxxxxxxxxxxxxxxxx   //

  if ( transformFilename.find( ".mat" ) != std::string::npos )
  {

    std::fstream file ;
    file.open( transformFilename, std::ios::in ) ;
    if ( file.is_open() )
    {

      std::string line ;
      int i = 0 ;
      while ( std::getline( file, line ) )
      {

        std::stringstream ss( line ) ;

        ss >> affineTransform[ i ] ;
        i++ ;

      }

      file.close() ;

    }


    if ( verbose )
    {

      std::cout << "Transform : [ " ;
      for ( int i = 0 ; i < 12 ; i++ )
      {

        if ( i < 11 )
        {

          std::cout << affineTransform[ i ] << ", " ;

        }
        else
        {

          std::cout << affineTransform[ i ] << " ]" << std::endl ;

        }

      }

    }

  }
  else
  {

    std::cout << "ERROR in transformation : The only format supported is .mat"
              << std::endl ;
    exit( 1 ) ;

  }


  //------------------------------ Sanity checks -----------------------------//
  int nbPointsMoving = movingBundlesData.pointsPerTrack[ 0 ] ;
  int nbCurvesMoving = movingBundlesData.curves_count ;
  for( int fiber = 1 ; fiber < nbCurvesMoving ; fiber++ )
  {

    if ( nbPointsMoving != movingBundlesData.pointsPerTrack[ fiber ] )
    {

      std::cout << "ERROR : All the fibers in the moving bundle must have the "
                << "same number of points, got " << nbPointsMoving
                << " for point 0 and "
                << movingBundlesData.pointsPerTrack[ fiber ]
                << " for point " << fiber << std::endl ;
      exit( 1 ) ;

    }

  }


  int nbPoints = nbPointsMoving ;

  if ( verbose )
  {

    std::cout << "Sanity cheks : OK " << std::endl ;

  }

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

  // Creating function to minimize
  BMDDistance fun( movingBundlesData, movingBundlesData ) ;


  BundlesDataFormat movedBundlesData( movingBundlesData ) ;
  fun.applyAffineToBundle( movingBundlesData.matrixTracks,
                           affineTransform,
                           movingBundlesData.curves_count,
                           nbPoints,
                           movedBundlesData.matrixTracks ) ;

  BundlesFormat movedBundleInfo( movingBundleInfo ) ;

  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //

  // Saving results
  std::string outputBundlesFilename ;
  std::string outputBundlesDataFilename ;

  if ( outputFilename.find( ".bundlesdata" ) != std::string::npos )
  {

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = outputFilename ;
    size_t index = outputFilename.find( ".bundlesdata" ) ;
    outputBundlesFilename.replace( index, 12, ".bundles") ;

  }
  else if ( outputFilename.find( ".bundles" ) != std::string::npos )
  {

    outputBundlesFilename = outputFilename ;

    outputBundlesDataFilename = outputFilename ;
    size_t index = outputFilename.find( ".bundles" ) ;
    outputBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }

  movedBundlesData.bundlesdataWriting( outputBundlesDataFilename.c_str(),
                                                                     verbose ) ;
  movedBundleInfo.bundlesWriting( outputBundlesFilename.c_str(), verbose ) ;


  return 0 ;

}
