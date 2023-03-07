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

#include "computeCentroids.h"
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
///////////////////////// Function to close dipy server ////////////////////////
////////////////////////////////////////////////////////////////////////////////
void closeDipyServer( int portDipyServer )
{

  std::cout << "Closing dipy service... " ;
  std::ostringstream closeDipyServiceOss ;
  if ( index_cds )
  {

    closeDipyServiceOss << "python3 " << closeDipyServerClientFile
                                                                      << " " ;

  }
  else
  {

    closeDipyServiceOss << closeDipyServerClientFile << " " ;

  }
  closeDipyServiceOss << "-p " << portDipyServer << " " ;
  std::string closeDipyService = closeDipyServiceOss.str() ;

  int isCloseServerFail = run_sh_process( closeDipyService ) ;
  if ( isCloseServerFail )
  {

    std::cout << "WARNING : could not close dipy service... " << std::endl ;

  }

}


////////////////////////////////////////////////////////////////////////////////
////////////////// Function to get port number of dipy server //////////////////
////////////////////////////////////////////////////////////////////////////////
int getPortNumberDipyService( std::string& logFilePath )
{

  sleep( 1 ) ;

	const char delim = ':' ;
	// const char delim = '\n' ;
  std::string line ;
  std::ifstream logFile ;
	int portDipyServer = -1 ;
  logFile.open( logFilePath.c_str() ) ;
  if ( logFile.fail() )
  {

    std::cout << "Problem reading file : " << logFilePath
                                           << std::endl ;


    exit( 1 ) ;

  }
  while ( std::getline( logFile, line ) )
  {

    // std::vector< std::string > out ;
    std::stringstream ss( line ) ;
    std::string s ;
		bool isPort = false ;
    while ( std::getline( ss, s, delim ) )
    {

      // s.erase( std::remove( s.begin(), s.end(), ' ' ), s.end() ) ;
			// std::cout << "Line : " << i << " : \n" << s << std::endl ;
      // out.push_back( s ) ;
			if ( s == "Server port " )
			{

				isPort = true ;

			}

			if ( isPort )
			{

				try
				{

					portDipyServer = std::stoi( s ) ;
					isPort = false ;
					break ;

				}
				catch ( std::invalid_argument ){}

			}

			if ( portDipyServer >= 5000 )
			{

				break ;

			}

    }

  }

	logFile.close() ;

	return( portDipyServer ) ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  const auto start_time = std::chrono::system_clock::now() ;

  index_input = getFlagPosition( argc, argv, "-i" ) ;
  index_reference = getFlagPosition( argc, argv, "-ref" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_cc = getFlagPosition( argc, argv, "-cc" ) ;
  index_ods = getFlagPosition( argc, argv, "-ods" ) ;
  index_cds = getFlagPosition( argc, argv, "-cds" ) ;
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints" ) ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_format = getFlagPosition( argc, argv, "-format" ) ;
  index_force = getFlagPosition( argc, argv, "-force" ) ;
  index_thr = getFlagPosition( argc, argv, "-thr" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Compute centroids of atlas using QuickBundles : \n"
              << "-i : Path to the input tractogram \n"
              << "-ref : Path to the reference image where the tractogram is \n"
              << "-o : Path to the output directory \n"
              << "-nbPoints : Number of points per fiber (same number for all "
              << "fibers) \n"
              << "-format : Format of bundles in directory (.tck/.trk/.bundles)"
                                                                    << std::endl
              << "[-cc] : Path to the clientComputeCentroids.py file \n"
              << "[-ods] : Path to the dipyServer.py file \n"
              << "[-cds] : Path to the clientCloseServer.py file \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads "
              << "(default : number of cores ) \n"
              << "[-thr] : threshold for QuickBundles (default : based on .minf"
              <<                                             "averageRadius) \n"
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

  if ( !index_reference )
  {

    std::cout << "-ref argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_nbPoints )
  {

    std::cout << "-nbPoints argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_format )
  {

    std::cout << "-format argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  //////////////////////////////////////////////////////////////////////////////
  int nbCores = omp_get_num_procs() ;
  std::cout << "Number of cores in the system : " << nbCores << std::endl ;

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Checking arguments /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  ////////////////////////////// Getting format ////////////////////////////////
  format = argv[ index_format + 1 ] ;
  if ( format != ".bundles" && format != ".bundlesdata" && format != ".tck"
                                                           && format != ".trk" )
  {

    std::string outMessage = "the only supported formats are .bundles/.trk/" \
                                                                     ".tck \n" ;
    throw( std::invalid_argument( outMessage ) ) ;


  }

  ////////////////////////////// Input directory ///////////////////////////////
  inputDirPath = argv[ index_input + 1 ] ;
  char lastChar = inputDirPath[ inputDirPath.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    inputDirPath = inputDirPath + "/" ;

  }

  if ( !is_dir( inputDirPath ) )
  {

    std::string outMessage = "inputDir does not exists \n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }

  if ( countFilesDirectory( inputDirPath ) == 0 )
  {

    std::string outMessage = "inputDir is empty \n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }

  atlasBundlesDataPaths = getFilesInDirectoryWithExtension( inputDirPath,
                                                                      format ) ;
  if ( atlasBundlesDataPaths.size() == 0 )
  {

    std::string outMessage = "not files of specified format in inputDir \n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }


  atlasBundlesInfoPaths.resize(atlasBundlesDataPaths.size() ) ;
  for ( int i = 0 ; i < atlasBundlesDataPaths.size() ; i++ )
  {

    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      atlasBundlesDataPaths[ i ] = replaceExtension( atlasBundlesDataPaths[ i ],
                                                              ".bundlesdata" ) ;

      atlasBundlesInfoPaths[ i ] = replaceExtension( atlasBundlesDataPaths[ i ],
                                                                  ".bundles" ) ;

    }
    else
    {

      atlasBundlesInfoPaths[ i ] = replaceExtension( atlasBundlesDataPaths[ i ],
                                                                     ".minf" ) ;

    }

    if ( !is_file( atlasBundlesDataPaths[ i ] ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "File " << atlasBundlesDataPaths[ i ] << " does not "
                                                                 << "exists\n" ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

    if ( !is_file( atlasBundlesInfoPaths[ i ] ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "File " << atlasBundlesInfoPaths[ i ] << " does not "
                                                                 << "exists\n" ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }


  }

  /////////////////////////////// Reference image //////////////////////////////
  referenceFilename = argv[ index_reference + 1 ] ;
  lastChar = referenceFilename[ referenceFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    referenceFilename = referenceFilename.substr( 0,
                                                referenceFilename.size() - 1 ) ;

  }
  if ( !endswith( referenceFilename, ".nii" ) )
  {

    std::string outMessage = "The only reference image format supported is "\
                             ".nii\n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }
  else
  {

    std::cout << "Reference image : OK " << std::endl ;

  }

  if ( !is_file( referenceFilename ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ERROR : Reference image file " << referenceFilename
                                            << " does not exists" << std::endl ;
    std::string outMessage = outMessageOss.str() ;
    throw( std::invalid_argument( outMessage ) ) ;

  }

  ////////////////////////////// Output directory //////////////////////////////
  outputDirectory = argv[ index_output + 1 ] ;
  lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectory = outputDirectory + "/" ;

  }

  if ( is_dir( outputDirectory ) && countFilesDirectory( outputDirectory ) > 0
                                                                     && !force )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ERROR : output directory " << outputDirectory
              << " is not empty and the -force flag was not used" << std::endl ;
    std::string outMessage = outMessageOss.str() ;
    throw( std::invalid_argument( outMessage ) ) ;

  }


  if ( !is_dir( outputDirectory ) )
  {

    mkdir( outputDirectory ) ;

  }



  ////////////////////////// Compute centorids command /////////////////////////
  if ( index_cc )
  {

    computeCentroidsClientFilename = argv[ index_cc + 1 ] ;
    lastChar = computeCentroidsClientFilename[
                                   computeCentroidsClientFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      computeCentroidsClientFilename = computeCentroidsClientFilename.substr( 0,
                                   computeCentroidsClientFilename.size() - 1 ) ;

    }
    if ( !is_file( computeCentroidsClientFilename ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : clientComputeCentroids.py file "
                    << computeCentroidsClientFilename << " does not exists "
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

      exit( 1 ) ;

    }

  }
  else
  {

    computeCentroidsClientFilename = "clientComputeCentroids.py" ;

  }



  //////////////////////////// Open server command /////////////////////////////
  if ( index_ods )
  {

    openDipyServerClientFile = argv[ index_ods + 1 ] ;
    lastChar = openDipyServerClientFile[ openDipyServerClientFile.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      openDipyServerClientFile = openDipyServerClientFile.substr( 0,
                                         openDipyServerClientFile.size() - 1 ) ;

    }
    if ( !is_file( openDipyServerClientFile ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : dipyServer.py file "
               << openDipyServerClientFile << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

      exit( 1 ) ;

    }

  }
  else
  {

    openDipyServerClientFile = "dipyServer.py" ;

  }

  //////////////////////////// Close server command ////////////////////////////
  if ( index_cds )
  {

    closeDipyServerClientFile = argv[ index_cds + 1 ] ;
    lastChar = closeDipyServerClientFile[ closeDipyServerClientFile.size()
                                                                         - 1 ] ;
    if ( lastChar == '/' )
    {

      closeDipyServerClientFile = closeDipyServerClientFile.substr( 0,
                                        closeDipyServerClientFile.size() - 1 ) ;

    }
    if ( !is_file( closeDipyServerClientFile ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : clientCloseServer.py file "
              << closeDipyServerClientFile << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

      exit( 1 ) ;

    }


  }
  else
  {

    closeDipyServerClientFile = "clientCloseServer.py" ;

  }

  ///////////////////////// Number of points per fiber /////////////////////////
  nbPointsPerFiber = std::stoi( argv[ index_nbPoints + 1 ] ) ;



  ///////////////////////// Number of points per fiber /////////////////////////
  if ( index_thr )
  {

    thrQb = std::stof( argv[ index_thr + 1 ] ) ;

  }
  else
  {

    thrQb = -1 ;

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


  ///////////////////////////////// nbThreads //////////////////////////////////
  if ( index_nbThreads )
  {


    nbThreads = std::stoi( argv[ index_nbThreads + 1 ] ) ;
    if ( nbThreads <= 0 )
    {

      std::cout << "Invalid argument for -nbThreads : you must give a postive "
                << "integer " << std::endl ;
      exit( 1 ) ;

    }

  }
  else
  {

    // nbThreads = nbCores ;
    nbThreads = -1 ;

  }

  omp_set_num_threads( nbThreads ) ;

  #pragma omp parallel
  {

    nbThreadsUsed = omp_get_num_threads() ;

  }
  std::cout << "Number of threads : " << nbThreadsUsed << std::endl ;



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


  /////////////////////////// Launching dipy service ///////////////////////////
  std::cout << "Launching dipy service... " ;

  std::ostringstream serverLogFilePathOss ;
  serverLogFilePathOss << outputDirectory << "dipyServiceLog.txt" ;
  std::string serverLogFilePath = serverLogFilePathOss.str() ;


  std::ostringstream launchDipyServiceOss ;
  if ( index_ods )
  {

    launchDipyServiceOss << "python3 " << openDipyServerClientFile
                                                                      << " " ;

  }
  else
  {

    launchDipyServiceOss << openDipyServerClientFile << " " ;

  }
  launchDipyServiceOss << "-lf " << serverLogFilePath << " " ;
  std::string launchDipyService = launchDipyServiceOss.str() ;


  boost::process::ipstream out ; // To not pipe output in main process
  boost::process::ipstream err ; // To not pipe error in main process

  boost::process::child c( launchDipyService.c_str(),
                                               boost::process::std_out > out,
                                               boost::process::std_err > err ) ;

  std::cout << "" << std::flush ;

  while ( !is_file( serverLogFilePath ) ){}

  int portDipyServer = getPortNumberDipyService( serverLogFilePath ) ;

  std::cout << "Using port " << portDipyServer << std::endl ;
  std::cout << "" << std::flush ;



  ///////////////////////////// Computing Centroids ////////////////////////////
  if ( verbose )
  {

    std::cout << "Computing Centroids... "  << std::endl ;

  }

  const auto start_time_tract_neigh = std::chrono::system_clock::now() ;


  int nbBundles = atlasBundlesDataPaths.size() ;

  std::ostringstream clientLogFilePathOss ;
  clientLogFilePathOss << outputDirectory << "clientDipyLog.txt" ;
  std::string clientLogFilePath = clientLogFilePathOss.str() ;

  int bundleCounter = 0 ;
  #pragma omp parallel for schedule(dynamic)
  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    #pragma omp critical
    {

      printf( "\rNumber of bundles processed : [ %d  /  %d ]",
                                                    bundleCounter, nbBundles ) ;
      std::cout << "" << std::flush ;
      bundleCounter++ ;

    }

    float tmpThr = -1 ;
    if ( index_thr )
    {

      tmpThr = thrQb ;

    }
    else
    {


      tmpThr = getAverageRadiusAtlasBundle( atlasBundlesInfoPaths[ i ] ) ;

    }

    std::string bundleName = getFilenameNoExtension(
                                                  atlasBundlesDataPaths[ i ] ) ;

    std::stringstream outputCentroidsPathOss ;
    outputCentroidsPathOss << outputDirectory << bundleName << "_centroids"
                                                                     << format ;
    std::string outputCentroidsPath = outputCentroidsPathOss.str() ;

    std::string inputBundlePath ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      inputBundlePath = atlasBundlesInfoPaths [ i ] ;

    }
    else
    {

      inputBundlePath = atlasBundlesDataPaths[ i ] ;

    }


    std::stringstream computeCentroidsCommandOss ;
    if ( index_cc )
    {

      computeCentroidsCommandOss << "python3 " << computeCentroidsClientFilename
                                                                        << " " ;

    }
    else
    {

      computeCentroidsCommandOss << computeCentroidsClientFilename << " " ;

    }

    computeCentroidsCommandOss << "-i " << inputBundlePath << " "
                               << "-o " << outputCentroidsPath << " "
                               << "-r " << referenceFilename << " "
                               << "-thr " << tmpThr << " "
                               << "-nbPoints " << nbPointsPerFiber << " "
                               << "-lf " << clientLogFilePath << " "
                               << "-p " << portDipyServer << " "
                               << "-v 1 " ;
    std::string computeCentroidsCommand = computeCentroidsCommandOss.str() ;

    int _tmpNbFibers = getNbFibers( atlasBundlesDataPaths[ i ] ) ;
    bool isFailCentroids = 0 ;
    if ( _tmpNbFibers > 500 )
    {

      isFailCentroids = run_sh_process( computeCentroidsCommand ) ;
      if ( is_file( outputCentroidsPath ) )
      {

        isFailCentroids = 0 ;
        int _tmpNbFibers2 = getNbFibers( outputCentroidsPath ) ;
        if ( _tmpNbFibers2 < 10 )
        {

          if ( format == ".bundles" || format == ".bundlesdata" )
          {

            std::string tmpDataPath = replaceExtension( outputCentroidsPath,
                                                              ".bundlesdata" ) ;
            std::string tmpInfoPath = replaceExtension( outputCentroidsPath,
                                                                  ".bundles" ) ;
            if ( is_file( tmpDataPath ) )
            {

              rmfile( tmpDataPath ) ;

            }

            if ( is_file( tmpInfoPath ) )
            {

              rmfile( tmpInfoPath ) ;

            }


          }
          else
          {

            std::string tmpDataPath = replaceExtension( outputCentroidsPath,
                                                                      format ) ;
            std::string tmpInfoPath = replaceExtension( outputCentroidsPath,
                                                                     ".minf" ) ;
            if ( is_file( tmpDataPath ) )
            {

              rmfile( tmpDataPath ) ;

            }

            if ( is_file( tmpInfoPath ) )
            {

              rmfile( tmpInfoPath ) ;

            }


          }
          rmfile( outputCentroidsPath ) ;

          std::stringstream tmpComputeCentroidsCommandOss ;
          if ( index_cc )
          {

            tmpComputeCentroidsCommandOss << "python3 "
                                      << computeCentroidsClientFilename << " " ;

          }
          else
          {

            tmpComputeCentroidsCommandOss << computeCentroidsClientFilename
                                                                        << " " ;

          }
          tmpComputeCentroidsCommandOss
                                      << "-i " << inputBundlePath << " "
                                      << "-o " << outputCentroidsPath << " "
                                      << "-r " << referenceFilename << " "
                                      << "-thr " << tmpThr / 2.0 << " "
                                      << "-nbPoints " << nbPointsPerFiber << " "
                                      << "-lf " << clientLogFilePath << " "
                                      << "-p " << portDipyServer << " "
                                      << "-v 1 " ;
          std::string tmpComputeCentroidsCommand =
                                           tmpComputeCentroidsCommandOss.str() ;

          isFailCentroids = run_sh_process( tmpComputeCentroidsCommand ) ;
          if ( is_file( outputCentroidsPath ) )
          {

            isFailCentroids = 0 ;

          }
          else
          {

            isFailCentroids = 1 ;

          }


        }

      }
      else
      {

        isFailCentroids = 1 ;

      }

    }
    else
    {

      std::string outputCentroidsInfoPath ;
      if ( format == ".bundles" || format == ".bundlesdata" )
      {

        outputCentroidsInfoPath = replaceExtension( outputCentroidsPath,
                                                                  ".bundles" ) ;
        outputCentroidsPath = replaceExtension( outputCentroidsPath,
                                                              ".bundlesdata" ) ;

      }
      else
      {

        outputCentroidsInfoPath = replaceExtension( outputCentroidsPath,
                                                                     ".minf" ) ;

      }

      copy( atlasBundlesDataPaths[ i ], outputCentroidsPath ) ;
      // copy( atlasBundlesInfoPaths[ i ], outputCentroidsInfoPath ) ;
      isFailCentroids = 0 ;

    }


  }


  const std::chrono::duration< double > duration_tract_neigh =
                     std::chrono::system_clock::now() - start_time_tract_neigh ;

  if ( verbose )
  {

    std::cout << "\nDuration : " << duration_tract_neigh.count() << std::endl ;
    std::cout << "Done" << std::endl ;


  }


  return( 0 ) ;


}
