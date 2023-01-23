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

#include "correctBundlesNan.h"
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


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  int index_input, index_output, index_verbose, index_help ;

  std::vector<std::string> possibleFlags{ "-i", "-o", "-v", "-h" } ;

  std::vector<bool> possibleFlagsNeedArgument{ true, true, true, false } ;



  for ( int k = 1 ; k < argc ; k++ )
  {

    std::string arg = argv[ k ] ;
    auto it = std::find( possibleFlags.begin(), possibleFlags.end(), arg ) ;
    if( it == possibleFlags.end() )
    {

      if ( k == 1 )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for ProjectAtlas command" << std::endl ;
        exit( 1 ) ;

      }

      std::string arg2 = argv[ k - 1 ] ;
      it = std::find( possibleFlags.begin(), possibleFlags.end(), arg2 ) ;
      if( it == possibleFlags.end() )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for ProjectAtlas command" << std::endl ;
        exit( 1 ) ;

      }
      int indexFound = it - possibleFlags.begin() ;
      if ( !possibleFlagsNeedArgument[ indexFound ] )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for ProjectAtlas command" << std::endl ;
        exit( 1 ) ;

      }

    }
    else
    {

      int indexFound = it - possibleFlags.begin() ;
      if ( possibleFlagsNeedArgument[ indexFound ] )
      {

        if ( k == argc - 1 )
        {

          std::cout << "ERROR in arguments : " << arg << " needs an input but "
                    << "none was given ( see -h )" << std::endl ;
          exit( 1 ) ;

        }

        std::string arg2 = argv[ k + 1 ] ;
        it = std::find( possibleFlags.begin(), possibleFlags.end(), arg2 ) ;
        if( it != possibleFlags.end() )
        {

          std::cout << "ERROR in arguments : " << arg << " needs an input but "
                    << "none was given ( see -h )" << std::endl ;
          exit( 1 ) ;

        }

      }

    }

  }

  index_input = getFlagPosition( argc, argv, "-i" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-i : Input tractogram \n"
              << "-o : Output corrected tractogram \n"
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
  std::string inputFilename( argv[ index_input + 1 ] ) ;
  char lastChar = inputFilename[ inputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputFilename = inputFilename.substr( 0, inputFilename.size() - 1 ) ;

  }


  std::string outputFilename( argv[ index_output + 1 ] ) ;
  lastChar = outputFilename[ outputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

  }




  //   xxxxxxxxxxxxxxxxxxxxx Reading Input tractogram xxxxxxxxxxxxxxxxxxxxx   //
  std::string inputBundlesFilename ;
  std::string inputBundlesDataFilename ;
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


  BundlesData inputBundlesData( inputBundlesDataFilename.c_str() ) ;
  BundlesMinf inputBundleInfo( inputBundlesFilename.c_str() ) ;



  ///////////////////////// Saving corrected tractogram ////////////////////////

  std::string outputBundlesFilename ;
  std::string outputBundlesDataFilename ;

  if ( endswith( outputFilename, ".bundlesdata" ) )
  {

    if ( format != ".bundles" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = replaceExtension( outputFilename, ".bundles" ) ;

  }
  else if ( endswith( outputFilename, ".bundles" ) )
  {

    if ( format != ".bundles" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesFilename = outputFilename ;

    outputBundlesDataFilename = replaceExtension( outputFilename,
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

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = replaceExtension( outputFilename, ".minf" ) ;

  }
  else if ( endswith( outputFilename, ".tck" ) )
  {

    if ( format != ".tck" )
    {

      std::cout << "ERROR : input and output format must be the same "
                << std::endl ;
      exit( 1 ) ;

    }

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = replaceExtension( outputFilename, ".minf" ) ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }

  inputBundlesData.write( outputBundlesDataFilename.c_str(), inputBundleInfo ) ;


  return 0 ;

}
