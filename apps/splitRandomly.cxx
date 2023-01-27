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

#include <random>

#include <boost/process.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "splitRandomly.h"
#include "ioWrapper.h"

////////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ////////////////
////////////////////////////////////////////////////////////////////////////////

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
/////////////////////// Generates true with X probability //////////////////////
////////////////////////////////////////////////////////////////////////////////
bool generateTrueWithXprobability( float probability )
{

  // boost::random::mt19937 gen ;
  //
  // // Generator of integers between [1, 10]
  // boost::random::uniform_int_distribution<> dist( 1, 10 ) ;
  // int randomInteger = dist( gen ) ;

  std::random_device rd ;
  std::mt19937 gen( rd() ) ;
  std::uniform_real_distribution<> dist( 0, 100 ) ;

  int randomInteger = dist( gen ) ;

  if ( randomInteger <= probability * 100 )
  {

    return( true ) ;

  }

  return( false ) ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  const auto start_time = std::chrono::system_clock::now() ;

  int index_input, index_output, index_p, index_force, index_verbose,
                                                                    index_help ;

  index_input = getFlagPosition( argc, argv, "-i") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_p = getFlagPosition( argc, argv, "-p") ;
  index_force = getFlagPosition( argc, argv, "-force") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-i : input .bundles/.trk/.tck file \n"
              << "-o : Path of the output file (same format as input) \n"
              << "-p : Percentage to keep \n"
              << "[-force] : Force to overwrite files (default = false) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_p )
  {

    std::cout << "-p argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Checking arguments /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////// Force ////////////////////////////////////
  if ( index_force )
  {

    std::string _tmpIndexForce( argv[ index_force + 1 ] ) ;
    if ( _tmpIndexForce == "true" )
    {

      force = true ;

    }
    else if ( _tmpIndexForce == "false" )
    {

      force = false ;

    }
    else
    {

      std::cout << "Argument of -force must be either \"true\" or \"false\" "
                << std::endl ;
      exit( 1 ) ;

    }

  }

  ///////////////////////////////// Input file /////////////////////////////////
  std::string inputFilename = argv[ index_input + 1 ] ;
  char lastChar = inputFilename[ inputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputFilename = inputFilename.substr( 0,
                                             inputFilename.size() - 1 ) ;

  }
  if ( !is_file( inputFilename ) )
  {

    std::cout << "ERROR : Input file " << inputFilename << " does not "
                                                     << "exists " << std::endl ;
    exit( 1 );

  }


  bool isBundlesFormat = false ;
  bool isTrk = false ;
  bool isTck = false ;

  if ( endswith( inputFilename, ".bundlesdata" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".bundles" ) ;

    isBundlesFormat = true ;
    format = ".bundles" ;

  }
  else if ( endswith( inputFilename, ".bundles" ) )
  {

    inputBundlesFilename = inputFilename ;

    inputBundlesDataFilename = replaceExtension( inputFilename,
                                                              ".bundlesdata" ) ;

    isBundlesFormat = true ;
    format = ".bundles" ;

  }
  else if ( endswith( inputFilename, ".trk" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".minf" ) ;

    isTrk = true ;
    format = ".trk" ;

  }
  else if ( endswith( inputFilename, ".tck" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".minf" ) ;

    isTck = true ;
    format = ".tck" ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }




  ///////////////////////////////// Output file ////////////////////////////////
  std::string outputFilename = argv[ index_output + 1 ]  ;
  lastChar = outputFilename[ outputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

  }
  if ( is_file( outputFilename ) && !force )
  {

    std::cout << "WARNING : output file " << outputFilename
              << " already exists, and -force flag was not use, skipping "
              << "computations" << std::endl ;
    exit( 1 ) ;

  }

  if ( endswith( outputFilename, ".bundlesdata" ) )
  {

    if ( format != ".bundles" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outBundlesDataFilename = outputFilename ;

    outBundlesFilename = replaceExtension( outputFilename, ".bundles" ) ;

  }
  else if ( endswith( outputFilename, ".bundles" ) )
  {

    if ( format != ".bundles" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outBundlesFilename = outputFilename ;

    outBundlesDataFilename = replaceExtension( outputFilename,
                                                              ".bundlesdata" ) ;

  }
  else if ( endswith( outputFilename, ".trk" ) )
  {

    if ( format != ".trk" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outBundlesDataFilename = outputFilename ;

    outBundlesFilename = replaceExtension( outputFilename, ".minf" ) ;

  }
  else if ( endswith( outputFilename, ".tck" ) )
  {

    if ( format != ".tck" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outBundlesDataFilename = outputFilename ;

    outBundlesFilename = replaceExtension( outputFilename, ".minf" ) ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }



  ///////////////////////////// Percentage split ///////////////////////////////
  if ( index_p )
  {

    percentageSplit = std::stof( argv[ index_p + 1 ] ) ;
    if ( percentageSplit < 0 || percentageSplit > 1 )
    {

      std::cout << "Invalid argument for -p : value must belong to [0,1] "
                << std::endl ;
      exit( 1 ) ;

    }

  }

  ////////////////////////////////// Verbose ///////////////////////////////////
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


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  BundlesMinf tractogramInfo( inputBundlesFilename.c_str() ) ;
  BundlesData tractogramData( inputBundlesDataFilename.c_str() ) ;


  int nbFibers = tractogramData.curves_count ;
  // Checking points per track
  nbPointsPerFiber = tractogramData.pointsPerTrack[ 0 ] ;
  for ( int fiberIndex = 1 ; fiberIndex < nbFibers ; fiberIndex++ )
  {

    if ( nbPointsPerFiber != tractogramData.pointsPerTrack[ fiberIndex ] )
    {

      std::cout << "ERROR : all the fibers must have the same number of points "
                << "per track" << std::endl ;
      exit( 1 ) ;

    }

  }



  int nbFibersSplit = round( nbFibers * percentageSplit ) ;
  bool _stop = false ;
  std::vector<int> indexFibersToKeep ;
  while ( !_stop )
  {

    for ( int fiberIndex = 0 ; fiberIndex < nbFibers ; fiberIndex++ )
    {

      if ( std::find( indexFibersToKeep.begin(),
                              indexFibersToKeep.end(),
                                 fiberIndex ) != indexFibersToKeep.end() )
      {

        continue ;

      }

      if ( generateTrueWithXprobability( percentageSplit ) )
      {

        indexFibersToKeep.push_back( fiberIndex ) ;
        if ( indexFibersToKeep.size() == nbFibersSplit )
        {
          _stop = true ;
          break ;

        }

      }

    }

  }

  std::vector<float> matrixTracks ;
  std::vector<int32_t> pointsPerTrack  ;
  int curves_count = nbFibersSplit ;
  for ( int fiberIndex = 0 ; fiberIndex < nbFibersSplit ; fiberIndex++ )
  {

    int selectedFiberIndex = indexFibersToKeep[ fiberIndex ] ;
    int64_t offset = 3 * nbPointsPerFiber * selectedFiberIndex ;
    pointsPerTrack.push_back(
                         tractogramData.pointsPerTrack[ selectedFiberIndex ] ) ;
    for ( int _point = 0 ; _point < nbPointsPerFiber ; _point++ )
    {

      for ( int coord = 0 ; coord < 3 ; coord++ )
      {

        matrixTracks.push_back( tractogramData.matrixTracks[
                                               3 * _point + coord + offset ] ) ;

      }

    }

  }

  BundlesMinf tractogramInfoOut( tractogramInfo ) ;
  tractogramInfoOut.curves_count = nbFibersSplit ;

  BundlesData tractogramDataOut( tractogramData ) ;
  tractogramDataOut.matrixTracks = matrixTracks ;
  tractogramDataOut.pointsPerTrack = pointsPerTrack ;
  tractogramDataOut.curves_count = curves_count ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    tractogramInfoOut.haveMinf = true ;

  }
  tractogramDataOut.write( outBundlesDataFilename.c_str(), tractogramInfoOut ) ;

  return( 0 ) ;

}
