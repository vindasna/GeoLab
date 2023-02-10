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

#include <boost/process.hpp>

#include "convertAtlasFormat.h"
#include "ioWrapper.h"


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

  const auto start_time = std::chrono::system_clock::now() ;

  int index_input, index_output, index_input_format, index_output_format,
                        index_nbPoints, index_force, index_verbose, index_help ;

  index_input = getFlagPosition( argc, argv, "-i") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_input_format = getFlagPosition( argc, argv, "-if") ;
  index_output_format = getFlagPosition( argc, argv, "-of") ;
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints") ;
  index_force = getFlagPosition( argc, argv, "-force") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-i : Path to the input directory containing the bundles \n"
              << "-o : Path to the output directory \n"
              << "-if : input bundles format among .bundles/.trk/.tck \n"
              << "-of : ouput bundles format among .bundles/.trk/.tck \n"
              << "[-nbPoints] : Number of points per fiber (same number for all "
              << "fibers) \n"
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

  if ( !index_input_format )
  {

    std::cout << "-if argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output_format )
  {

    std::cout << "-of argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Checking arguments /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  ////////////////////////////// Input Directory ///////////////////////////////
  inputDirectoryPath = argv[ index_input + 1 ] ;
  char lastChar = inputDirectoryPath[ inputDirectoryPath.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    inputDirectoryPath = inputDirectoryPath + "/" ;

  }

  if ( !is_dir( inputDirectoryPath ) )
  {

    std::cout << "ERROR : input directory path " << inputDirectoryPath
              << " does not exists" << std::endl ;
    exit( 1 ) ;

  }

  ////////////////////////////// Output directory //////////////////////////////
  outputDirectoryPath = argv[ index_output + 1 ] ;
  lastChar = outputDirectoryPath[ outputDirectoryPath.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectoryPath = outputDirectoryPath + "/" ;


  }

  if ( !is_dir( outputDirectoryPath ) )
  {

    mkdir( outputDirectoryPath ) ;

  }


  //////////////////////////// Getting input format ////////////////////////////
  inputFormat = argv[ index_input_format + 1 ] ;
  if ( inputFormat == ".bundles" )
  {

    isBundlesFormat = true ;

  }
  else if ( inputFormat == ".trk" )
  {

    isTrkFormat = true ;

  }
  else if ( inputFormat == ".tck" )
  {

    isTckFormat = true ;

  }
  else
  {

    std::cout << "Argument -if must be either .bundles/.trk/.tck "
              << std::endl ;
    exit( 1 ) ;

  }

  //////////////////////////// Getting output format ///////////////////////////
  outputFormat = argv[ index_output_format + 1 ] ;
  if ( outputFormat == ".bundles" )
  {

    isBundlesFormat = true ;

  }
  else if ( outputFormat == ".trk" )
  {

    isTrkFormat = true ;

  }
  else if ( outputFormat == ".tck" )
  {

    isTckFormat = true ;

  }
  else
  {

    std::cout << "Argument -if must be either .bundles/.trk/.tck "
              << std::endl ;
    exit( 1 ) ;

  }



  ///////////////////////// Number of points per fiber /////////////////////////
  if ( index_nbPoints )
  {

    nbPointsPerFiber = std::stoi( argv[ index_nbPoints + 1 ] ) ;

    if ( nbPointsPerFiber <= 0 )
    {

      std::cout << "ERROR : -nbPoints must be greater than 0 " << std::endl ;
      exit( 1 ) ;

    }

  }


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
  //////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> bundlesPathsAtlasDir =
                           getFilesInDirectoryWithExtension( inputDirectoryPath,
                                                                 inputFormat ) ;

  std::vector<std::string> bundlesPathsOutDir ;
  for ( std::string _s : bundlesPathsAtlasDir )
  {
    std::cout << _s << std::endl ;
    std::string bundleName = getFilenameNoExtension( _s ) ;
    std::stringstream outBundlePathOss ;
    outBundlePathOss << outputDirectoryPath << bundleName << outputFormat ;
    std::cout << outBundlePathOss.str() << "\n" << std::endl ;
    bundlesPathsOutDir.push_back( outBundlePathOss.str() ) ;


  }



}
