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
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include <boost/process.hpp>

#include "localSbrToAtlas.h"
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
/////////////////////////// Function to apply GeoLab ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
void applyLocalSBR( const std::string& movedTractogramNeighborhood,
                    const std::string& atlasBundleDirectory,
                    const std::string& atlasNeighborhoodFile,
                    const std::string& atlasNeighborhoodCentroidsFile,
                    const std::string& outputDirectory,
                    const std::string& referenceImage,
                    const std::string& format,
                    int nbPointsPerFiber,
                    int portDipyServer,
                    float& time_out,
                    int verbose )
{

  std::string _tmpAtlasBundlePath = atlasNeighborhoodFile ;

  char lastChar = _tmpAtlasBundlePath[ _tmpAtlasBundlePath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    _tmpAtlasBundlePath = _tmpAtlasBundlePath.substr( 0,
                                            _tmpAtlasBundlePath.size() - 1 ) ;

  }

  std::experimental::filesystem::path tmpPath(
                                 replaceExtension( _tmpAtlasBundlePath, "" ) ) ;
  std::string bundleName = tmpPath.stem() ;

  std::ostringstream clientsLogFilePathOss ;
  clientsLogFilePathOss << outputDirectory << "clientsDipyLog.txt" ;
  std::string clientsLogFilePath = clientsLogFilePathOss.str() ;


  std::string atlasNeighborhood ;
  if ( !( endswith( atlasNeighborhoodFile, format ) ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ERROR : in applyLocalSBR, atlasNeighborhoodFile must "
                  << "be in .bundles/.trk/.tck format, got "
                  << atlasNeighborhoodFile << std::endl ;
    std::string outMessage = outMessageOss.str() ;
    closeDipyServer( portDipyServer ) ;
    throw( std::invalid_argument( outMessage ) ) ;


  }
  else
  {

    atlasNeighborhood = atlasNeighborhoodFile ;

  }

  if ( !is_file( atlasNeighborhood ) )
  {

    std::cout << "ERROR : in applyLocalSBR, atlas neighborhood file "
              << atlasNeighborhood << " does not exists " << std::endl ;

    closeDipyServer( portDipyServer ) ;
    exit( 1 ) ;

  }


  // Computing centroids neighborhood atlas
  std::ostringstream atlasBundleFileOss ;
  if ( is_dir( atlasBundleDirectory ) )
  {

    atlasBundleFileOss << atlasBundleDirectory << bundleName << format ;

  }
  else
  {

      atlasBundleFileOss << atlasBundleDirectory ;

  }
  std::string atlasBundleFile = atlasBundleFileOss.str() ;
  float averageRadius = getAverageRadiusAtlasBundle( atlasBundleFile ) ;
  std::string atlasNeighborhoodCentroids ;
  if ( !isAtlasNeighborhoodCentroids )
  {

    atlasNeighborhoodCentroids = replaceExtension( atlasNeighborhood,
                                                   "_centroids.bundles" ) ;
    if ( format == ".trk" || format == ".tck" )
    {

      atlasNeighborhoodCentroids = replaceExtension( atlasNeighborhoodCentroids,
                                                                      format ) ;



    }

    std::ostringstream computeCentroidsCommandClientOss ;
    computeCentroidsCommandClientOss << computeCentroidsClientFilename << " " ;
    computeCentroidsCommandClientOss
                                   << "-i " << atlasNeighborhood << " "
                                   << "-o " << atlasNeighborhoodCentroids << " "
                                   << "-r " << referenceImage << " "
                                   << "-thr " << averageRadius << " "
                                   << "-nbPoints " << nbPointsPerFiber << " "
                                   << "-lf " << clientsLogFilePath << " "
                                   << "-p " << portDipyServer << " "
                                   << "-timeOut " << time_out << " "
                                   << "-v 1 " ;
    std::string computeCentroidsCommandClient =
                                        computeCentroidsCommandClientOss.str() ;
    int isFailCentroids = 0 ;
    if ( countFilesDirectory( atlasNeighborhoodCentroids ) > 5 && !force )
    {

      if ( verbose > 1 )
      {

        std::cout << "WARNING : output directory of "
                  << computeCentroidsClientFilename << " : "
                  << atlasNeighborhoodCentroids << " already exists and "
                  << "with more than 5 files and th -force flag was not used, "
                  << "trying computations with existing directory"
                  << std::endl ;

      }

    }
    else
    {

      int _tmpNbFibers = getNbFibers( atlasNeighborhood ) ;
      if ( _tmpNbFibers > 500 )
      {

        isFailCentroids = run_sh_process( computeCentroidsCommandClient ) ;
        if ( is_file( atlasNeighborhoodCentroids ) )
        {

          isFailCentroids = 0 ;

        }
        else
        {

          isFailCentroids = 1 ;

        }

      }
      else
      {

        atlasNeighborhoodCentroids = atlasNeighborhood ;
        isFailCentroids = 0 ;

      }


    }
    if ( isFailCentroids )
    {

      std::cout << "ERROR : could not compute centroids for "
                << atlasNeighborhood << ", got exit_code " << isFailCentroids
                                                                  << std::endl ;

      std::cout << "ERROR WITH COMMAND : \n " << computeCentroidsCommandClient
                                              << std::endl ;


      closeDipyServer( portDipyServer ) ;
      exit( 1 ) ;

    }

  }
  else
  {

    atlasNeighborhoodCentroids = atlasNeighborhoodCentroidsFile ;

  }

  // Compute centroids neighborhood moved tractogram
  // int nbClustersAtlasNeighborhood = getNbFibers( atlasNeighborhoodCentroids ) ;
  std::string movedTractogramNeighborhoodCentroids = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "_centroids.bundles" ) ;

  if ( ( format == ".trk" || format == ".tck" ) )
  {

    movedTractogramNeighborhoodCentroids = replaceExtension(
                                           movedTractogramNeighborhoodCentroids,
                                           format ) ;

  }

  std::ostringstream computeCentroidsCommandClient2Oss ;
  computeCentroidsCommandClient2Oss << computeCentroidsClientFilename << " " ;
  computeCentroidsCommandClient2Oss
                         << "-i " << movedTractogramNeighborhood << " "
                         << "-o " << movedTractogramNeighborhoodCentroids << " "
                         << "-r " << referenceImage << " "
                         << "-thr " << averageRadius << " "
                         << "-nbPoints " << nbPointsPerFiber << " "
                         // << "-mnc " << nbClustersAtlasNeighborhood << " "
                         << "-lf " << clientsLogFilePath << " "
                         << "-p " << portDipyServer << " "
                         << "-timeOut " << time_out << " "
                         << "-v 1 " ;
  std::string computeCentroidsCommandClient2 =
                                       computeCentroidsCommandClient2Oss.str() ;
  int isFailCentroids = 0 ;
  if ( countFilesDirectory( movedTractogramNeighborhoodCentroids ) > 5 &&
                                                                        !force )
  {

    if ( verbose > 1 )
    {

      std::cout << "WARNING : output directory of "
                << computeCentroidsClientFilename << " : "
                << movedTractogramNeighborhoodCentroids << " already exists "
                << "with more than 5 files and -force flag was not used,"
                << " trying computations  with existing directory"
                << std::endl ;

    }

  }
  else
  {

    int _tmpNbFibers = getNbFibers( movedTractogramNeighborhood ) ;
    if ( ( _tmpNbFibers > 500 &&
                              !is_file( movedTractogramNeighborhoodCentroids ) )
                                            || ( _tmpNbFibers > 500 && force ) )
    {

      isFailCentroids = run_sh_process( computeCentroidsCommandClient2 ) ;
      if ( is_file( movedTractogramNeighborhoodCentroids ) )
      {

        isFailCentroids = 0 ;

      }
      else
      {

        isFailCentroids = 1 ;

      }

    }
    else
    {

      movedTractogramNeighborhoodCentroids = movedTractogramNeighborhood ;
      isFailCentroids = 0 ;

    }


  }
  if ( isFailCentroids )
  {

    std::cout << "ERROR : could not compute centroids for "
              << movedTractogramNeighborhood << ", got exit_code "
                                               << isFailCentroids << std::endl ;

    std::cout << "ERROR WITH COMMAND : \n " << computeCentroidsCommandClient2
                                            << std::endl ;

    closeDipyServer( portDipyServer ) ;
    exit( 1 ) ;

  }


  // Registering centroids
  std::string neighborhoodRegistered = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "_moved.bundles" ) ;


  if ( ( format == ".trk" || format == ".tck" ) )
  {

    neighborhoodRegistered = replaceExtension( neighborhoodRegistered, format) ;

  }

  std::ostringstream registerBundlesClientCommadOss ;
  registerBundlesClientCommadOss << registerBundlesClientFile << " " ;
  registerBundlesClientCommadOss
                         << "-s " << replaceExtension(
			             atlasNeighborhoodCentroids, format ) << " "
                         << "-m " << replaceExtension(
                           movedTractogramNeighborhoodCentroids, format ) << " "
                         << "-ra " << referenceImage << " "
                         << "-b " << replaceExtension(
                                    movedTractogramNeighborhood, format ) << " "
                         << "-o " << neighborhoodRegistered << " "
                         << "-n " << nbPointsPerFiber << " "
                         << "-xfm " << "rigid" << " "
                         << "-lf " << clientsLogFilePath << " "
                         << "-p " << portDipyServer << " "
                         << "-timeOut " << time_out << " "
                         << "-v 1" ;
  std::string registerBundlesClientCommad =
                                          registerBundlesClientCommadOss.str() ;


  // int timeout = 100 ; // In s
  int timeout = 500 ; // In s
  std::string _tmpError1 ;
  if ( is_file( neighborhoodRegistered ) && !force )
  {

    if ( verbose > 1 )
    {

      std::cout << "WARNING : output file of " << registerBundlesClientFile
                << " : " << neighborhoodRegistered << " already exists and "
                << "the -force flag was not used, trying computations with "
                << "existing file" << std::endl ;

    }

  }
  else
  {

    _tmpError1 = run_sh_process_timeout( registerBundlesClientCommad,
                                                                     timeout ) ;
    if ( is_file( neighborhoodRegistered ) )
    {

      _tmpError1 = " OK" ;

    }
    else
    {

      _tmpError1 = "PROBLEM" ;

    }

  }

  int isFailRegister = 0 ;
  if ( !is_file( neighborhoodRegistered ) )
  {

    isFailRegister = 1 ;

  }
  if ( isFailRegister )
  {

    if ( verbose > 1 )
    {

      std::cout << "\nERROR : Fail with command :\n"
                << registerBundlesClientCommad << "\nFail to register "
                << movedTractogramNeighborhoodCentroids << " to "
                << atlasNeighborhoodCentroids << ", got : \"" << _tmpError1
                                                          << "\"" << std::endl ;

    }


    return ;

  }

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////// Function to close dipy server ////////////////////////
////////////////////////////////////////////////////////////////////////////////
void closeDipyServer( int portDipyServer )
{

  std::cout << "Closing dipy service... " ;
  std::ostringstream closeDipyServiceOss ;
  closeDipyServiceOss << closeDipyServerClientFile << " " ;
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
  index_atlas = getFlagPosition( argc, argv, "-a" ) ;
  index_reference = getFlagPosition( argc, argv, "-ref" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_cc = getFlagPosition( argc, argv, "-cc" ) ;
  index_rb = getFlagPosition( argc, argv, "-rb" ) ;
  index_ods = getFlagPosition( argc, argv, "-ods" ) ;
  index_cds = getFlagPosition( argc, argv, "-cds" ) ;
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints" ) ;
  index_fa = getFlagPosition( argc, argv, "-fa" ) ;
  index_an = getFlagPosition( argc, argv, "-an" ) ;
  index_anc = getFlagPosition( argc, argv, "-anc" ) ;
  index_tolThrCN = getFlagPosition( argc, argv, "-tolThrCN" ) ;
  index_cn = getFlagPosition( argc, argv, "-cn" ) ;
  index_slr = getFlagPosition( argc, argv, "-slr" ) ;
  index_force = getFlagPosition( argc, argv, "-force" ) ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_nbThreadsCN = getFlagPosition( argc, argv, "-nbThreadsCN" ) ;
  index_keep_tmp = getFlagPosition( argc, argv, "-ktf" ) ;
  index_time_out = getFlagPosition( argc, argv, "-timeOut" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using and SWM adapted version of"
              << " RecoBundles method : \n"
              << "-i : Path to the input tractogram \n"
              << "-o : Path to the output directory \n"
              << "-nbPoints : Number of points per fiber (same number for all "
              << "fibers) \n"
              << "[-a] : Path to the directory with the atlas, must be given if"
              << " not using the ESBA atlas \n"
              << "[-ref] : Path to the reference image where the tractogram is, "
              << "must be given if other than MNI152 \n"
              << "[-cc] : Path to the clientComputeCentroids.py file \n"
              << "[-rb] : Path to the clientRegisterBundles.py file \n"
              << "[-ods] : Path to the dipyServer.py file \n"
              << "[-cds] : Path to the clientCloseServer.py file \n"
              << "[-fa] : Path to full atlas tractogram (mandatory for global "
              << "SLR or if -anc is not given for other than ESBA atlas) \n"
              << "[-an] : Path to the atlas neighborhoods (must be given is -fa"
              << " is not given if other than ESBA atlas)\n"
              << "[-anc] : Path to the atlas neighborhoods centroids (must be "
              << "given is -fa is not given for other than ESBA atlas)\n"
              << "[-tolThrCN] : tolerance for computeNeighborhood threshold "
              << "(default = 2.0) \n"
              << "[-cn] : Path to the computeNeighborhood file \n"
              << "[-slr] : Do global SLR step (default : false)\n"
              << "[-force] : Force to overwrite files (default = false) \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads for all "
              << "computations except neighborhoods computations if "
              << "-nbThreadsCN is given (default : number of cores ) \n"
              << "[-nbThreadsCN] : Sets the value of omp_set_num_threads for the "
              << "computation of neighborhoods (default : value of -nbThreads )"
              << " \n"
              << "[-ktf] : Keep temp files (default = false) \n"
              << "[-timeOut] : Time out (in s) for computation of centroids and "
              << "local SBR registration (default = 50 s) \n"
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

  if ( !index_nbPoints )
  {

    std::cout << "-nbPoints argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  //////////////////////////////////////////////////////////////////////////////
  int nbCores = omp_get_num_procs() ;
  std::cout << "Number of cores in the system : " << nbCores << std::endl ;

  //////////////////////////////////////////////////////////////////////////////
  char lastChar ;
  char *tmpChar = std::getenv( "ESBA_DIR" ) ;
  if ( !tmpChar )
  {

    std::cout << "WARNING : environmental variable ESBA_DIR is not set"
              << std::endl ;

    default_ESBA_DIR = "" ;

  }
  else
  {

    default_ESBA_DIR = tmpChar ;

  }
  if ( default_ESBA_DIR.size() > 0 )
  {

    lastChar = default_ESBA_DIR[ default_ESBA_DIR.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      default_ESBA_DIR = default_ESBA_DIR + "/" ;

    }

  }

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Checking arguments /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  ////////////////////////////// Input tractogram //////////////////////////////
  std::string inputBundlesFilename( argv[ index_input + 1 ] ) ;

  if ( !endswith( inputBundlesFilename, ".bundles" ) &&
       !endswith( inputBundlesFilename, ".bundlesdata" ) &&
       !endswith( inputBundlesFilename, ".trk" ) &&
       !endswith( inputBundlesFilename, ".tck" ) )
  {

    std::string outMessage = "ProjectAtlasGeoLab : Only supported " \
                             "formats for input are .bundles/.bundlesdata, " \
                             ".trk, .tck \n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }


  lastChar = inputBundlesFilename[ inputBundlesFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputBundlesFilename = inputBundlesFilename.substr( 0,
                                             inputBundlesFilename.size() - 1 ) ;

  }


  if ( !is_file( inputBundlesFilename ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ProjectAtlasGeoLab : input file "
                  << inputBundlesFilename << "does not exists" << std::endl ;
    std::string outMessage = outMessageOss.str() ;
    throw( std::invalid_argument( outMessage ) ) ;


  }


  /// Getting format
  std::string format ;
  std::string formatUpper ;
  std::string inputBundlesMinfPath ;
  if ( endswith( inputBundlesFilename, ".bundles" ) ||
                              endswith( inputBundlesFilename, ".bundlesdata" ) )
  {

    format = ".bundles" ;
    formatUpper = "Bundles" ;

    inputBundlesMinfPath = replaceExtension( inputBundlesFilename,
                                                                  ".bundles" ) ;

  }
  else if ( endswith( inputBundlesFilename, ".trk" ) )
  {

    format = ".trk" ;
    formatUpper = "Trk" ;

    inputBundlesMinfPath = replaceExtension( inputBundlesFilename, ".minf" ) ;

  }
  else if ( endswith( inputBundlesFilename, ".tck" ) )
  {

    format = ".tck" ;
    formatUpper = "Tck" ;

    inputBundlesMinfPath = replaceExtension( inputBundlesFilename, ".minf" ) ;

  }
  // Already checked at the begginig of main if extension is one of those format

  if ( is_file( inputBundlesMinfPath ) )
  {

    haveMinf = true ;

  }


  ////////////////////////////// Atlas directory ///////////////////////////////
  std::string atlasDirectory ;
  if ( index_atlas )
  {

    isEsbaAtlas = false ;

    atlasDirectory = argv[ index_atlas + 1 ] ;
    lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      atlasDirectory = atlasDirectory + "/" ;

    }
    if ( !is_dir( atlasDirectory ) )
    {

      std::stringstream outMessageOss ;

      outMessageOss << "ERROR : Atlas directory " << atlasDirectory << " does not"
                                                      << " exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    isEsbaAtlas = true ;

    std::stringstream atlasDirectoryOss  ;
    atlasDirectoryOss << default_ESBA_DIR << "Esba" << formatUpper ;
    atlasDirectory = atlasDirectoryOss.str() ;
    lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      atlasDirectory = atlasDirectory + "/" ;

    }
    if ( !is_dir( atlasDirectory ) )
    {

      std::stringstream outMessageOss ;

      outMessageOss << "ERROR : Atlas directory " << atlasDirectory << " does "
                    << "not exists, it seems the environment variable ESBA_DIR "
                    << "was not correctly set" << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }

  /////////////////////////////// Reference image //////////////////////////////
  std::string referenceFilename ;
  if ( index_reference )
  {

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

  }
  else
  {

    std::stringstream referenceFilenameOss  ;
    referenceFilenameOss << default_ESBA_DIR
                         << "mni_icbm152_t1_tal_nlin_asym_09c_brain.nii" ;
    referenceFilename = referenceFilenameOss.str() ;
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
                    << " does not exists it seems the environment variable "
                    << "ESBA_DIR was not correctly set" << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }


  }


  ////////////////////////////// Output directory //////////////////////////////
  std::string outputDirectory( argv[ index_output + 1 ] ) ;
  lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectory = outputDirectory + "/" ;

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

    // std::ostringstream computeCentroidsClientFilenameOss ;
    // computeCentroidsClientFilenameOss << "python3 "
    //                                          << computeCentroidsClientFilename ;
    // computeCentroidsClientFilename = computeCentroidsClientFilenameOss.str() ;

  }
  else
  {

    computeCentroidsClientFilename = "clientComputeCentroids.py" ;

  }

  ////////////////////////// Register bundles command //////////////////////////
  if ( index_rb )
  {

    registerBundlesClientFile = argv[ index_rb + 1 ] ;
    lastChar = registerBundlesClientFile[ registerBundlesClientFile.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      registerBundlesClientFile = registerBundlesClientFile.substr( 0,
                                          registerBundlesClientFile.size() - 1 ) ;

    }
    if ( !is_file( registerBundlesClientFile ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : clientRegisterBundles.py file "
                << registerBundlesClientFile << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

      exit( 1 ) ;

    }

    // std::ostringstream registerBundlesClientFileOss ;
    // registerBundlesClientFileOss << "python3 " << registerBundlesClientFile ;
    // registerBundlesClientFile = registerBundlesClientFileOss.str() ;

  }
  else
  {

    registerBundlesClientFile = "clientRegisterBundles.py" ;

  }



  //////////////////////////// Open server command ////////////////////////////
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

    // std::ostringstream openDipyServerClientFileOss ;
    // openDipyServerClientFileOss << "python3 " << openDipyServerClientFile ;
    // openDipyServerClientFile = openDipyServerClientFileOss.str() ;


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

    // std::ostringstream closeDipyServerClientFileOss ;
    // closeDipyServerClientFileOss << "python3 " << closeDipyServerClientFile ;
    // closeDipyServerClientFile = closeDipyServerClientFileOss.str() ;


  }
  else
  {

    closeDipyServerClientFile = "clientCloseServer.py" ;

  }

  ///////////////////////// Number of points per fiber /////////////////////////
  nbPointsPerFiber = std::stoi( argv[ index_nbPoints + 1 ] ) ;


  ///////////////////////////// Atlas neighborhood /////////////////////////////
  std::string atlasNeighborhoodDirectory ;
  if ( index_an )
  {

    isAtlasNeighborhood = true ;

    atlasNeighborhoodDirectory = argv[ index_an + 1 ] ;
    lastChar = atlasNeighborhoodDirectory[ atlasNeighborhoodDirectory.size()
                                                                         - 1 ] ;
    if ( lastChar != '/' )
    {

      atlasNeighborhoodDirectory = atlasNeighborhoodDirectory + "/" ;

    }
    if ( !is_dir( atlasNeighborhoodDirectory ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : Atlas neighborhood directory "
                    << atlasNeighborhoodDirectory << " does not exists "
                    << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    if ( isEsbaAtlas )
    {

      std::stringstream atlasNeighborhoodDirectoryOss ;
      atlasNeighborhoodDirectoryOss << default_ESBA_DIR << "Neighborhood"
                                                                << formatUpper ;
      atlasNeighborhoodDirectory = atlasNeighborhoodDirectoryOss.str() ;

      lastChar = atlasNeighborhoodDirectory[ atlasNeighborhoodDirectory.size()
                                                                         - 1 ] ;
      if ( lastChar != '/' )
      {

        atlasNeighborhoodDirectory = atlasNeighborhoodDirectory + "/" ;

      }
      if ( !is_dir( atlasNeighborhoodDirectory ) )
      {

        std::stringstream outMessageOss ;
        outMessageOss << "WARNING : Default atlas neighborhood directory "
                      << atlasNeighborhoodDirectory << " does not exists, you "
                      << "might not have set the ESBA_DIR environment variable "
                      << "correctly, the atlas neighborhood will be computed on"
                      << " the fly" << std::endl ;
        std::string outMessage = outMessageOss.str() ;
        std::cout << outMessage ;

      }
      else
      {

        isAtlasNeighborhood = true ;

      }

    }

  }

  //////////////////////// Atlas neighborhood centroids ////////////////////////
  std::string atlasNeighborhoodCentroidsDirectory ;
  if ( index_anc )
  {

    isAtlasNeighborhoodCentroids = true ;

    atlasNeighborhoodCentroidsDirectory = argv[ index_anc + 1 ] ;
    lastChar = atlasNeighborhoodCentroidsDirectory[
                              atlasNeighborhoodCentroidsDirectory.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      atlasNeighborhoodCentroidsDirectory = atlasNeighborhoodCentroidsDirectory
                                                                         + "/" ;

    }
    if ( !is_dir( atlasNeighborhoodCentroidsDirectory ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : Atlas neighborhood centroids directory "
                   << atlasNeighborhoodCentroidsDirectory << " does not exists "
                   << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    if ( isEsbaAtlas )
    {

      std::stringstream atlasNeighborhoodCentroidsDirectoryOss ;
      atlasNeighborhoodCentroidsDirectoryOss << default_ESBA_DIR << "Centroids"
                                                                << formatUpper ;
      atlasNeighborhoodCentroidsDirectory =
                                  atlasNeighborhoodCentroidsDirectoryOss.str() ;
      lastChar = atlasNeighborhoodCentroidsDirectory[
                              atlasNeighborhoodCentroidsDirectory.size() - 1 ] ;
      if ( lastChar != '/' )
      {

        atlasNeighborhoodCentroidsDirectory = atlasNeighborhoodCentroidsDirectory
                                                                         + "/" ;

      }
      if ( !is_dir( atlasNeighborhoodCentroidsDirectory ) )
      {

        std::stringstream outMessageOss ;
        outMessageOss << "WARNING : Default atlas neighborhood centroids "
                      << "directory "<< atlasNeighborhoodCentroidsDirectory
                      << "does not exists, you might not have set the ESBA_DIR "
                      << "environment variable correctly, the atlas "
                      << "neighborhood centroids will be computed on the fly"
                      << std::endl ;
        std::string outMessage = outMessageOss.str() ;
        std::cout << outMessage ;

      }
      else
      {

        isAtlasNeighborhoodCentroids = true ;

      }

    }

  }

  ///////////////////////////////// Full atlas /////////////////////////////////
  // Must be after Checking -an and -anc
  std::string fullAtlasFilename ;
  if ( index_fa )
  {

    std::string tmpFullAtlasPath( argv[ index_fa + 1 ] ) ;
    isFullAtlas = true ;

    if ( !endswith( tmpFullAtlasPath, ".bundles" ) &&
         !endswith( tmpFullAtlasPath, ".bundlesdata" ) &&
         !endswith( tmpFullAtlasPath, ".trk" ) &&
         !endswith( tmpFullAtlasPath, ".tck" ) )
    {

      std::string outMessage = "ProjectAtlasGeoLab: Only supported formats " \
                               "for atlas are .bundles/.bundlesdata, .trk, " \
                               ".tck \n" ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

    char lastChar = tmpFullAtlasPath[ tmpFullAtlasPath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      tmpFullAtlasPath = tmpFullAtlasPath.substr( 0,
                                                tmpFullAtlasPath.size() - 1 ) ;

    }

    fullAtlasFilename = tmpFullAtlasPath ;


    if ( !is_file( fullAtlasFilename ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : Full atlas file " << fullAtlasFilename
                                            << " does not exists" << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


    }
    else
    {

      std::cout << "Full atlas tractogram : OK " << std::endl ;

    }

  }
  else if ( !index_fa && !isAtlasNeighborhood && !isAtlasNeighborhoodCentroids )
  {

    std::cout << "ERROR : input argument -fa must be given if -anc and -an "
              << "are not give" << std::endl ;
    exit( 1 ) ;


  }


  
  /////////////////////////////// Tolerances ///////////////////////////////////

  if ( index_tolThrCN )
  {

    toleranceThrComputeNeighborhood = std::stof( argv[ index_tolThrCN + 1 ] ) ;
    // if ( toleranceThrComputeNeighborhood <= 0 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "ERROR : argument -tolThrCN must be greater than 0"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    // }

  }


  

  ////////////////////////// computeNeighborhood file //////////////////////////
  if ( index_cn )
  {

    std::string tmpFile( argv[ index_cn + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : computeNeighborhood file " << tmpFile
                                           << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }
      exit( 1 ) ;

    }

    computeNeighborhoodFile = tmpFile ;

  }

  /////////////////////////////////// Do SLR ///////////////////////////////////
  if ( index_slr )
  {

    std::string _tmpIndexSLR( argv[ index_slr + 1 ] ) ;
    if ( _tmpIndexSLR == "true" )
    {

      doSLR = true ;

    }
    else if ( _tmpIndexSLR == "false" )
    {

      doSLR = false ;

    }
    else
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Argument of -slr must be either \"true\" or \"false\" "
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


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

    nbThreads = nbCores ;
    // nbThreads = -1 ;

  }

  omp_set_num_threads( nbThreads ) ;

  #pragma omp parallel
  {

    nbThreadsUsed = omp_get_num_threads() ;

  }
  std::cout << "Number of threads : " << nbThreadsUsed << std::endl ;

  omp_set_nested( 1 ) ;


  //////////////////////////////// nbThreadsCN /////////////////////////////////
  if ( index_nbThreadsCN )
  {


    nbThreadsCN = std::stoi( argv[ index_nbThreadsCN + 1 ] ) ;
    if ( nbThreadsCN <= 0 )
    {

      std::cout << "Invalid argument for -nbThreadsCN : you must give a postive"
                << " integer " << std::endl ;
      exit( 1 ) ;

    }

  }
  else
  {

    nbThreadsCN = nbThreads ;

  }


  ////////////////////////////// Keep temp files ///////////////////////////////
  if ( index_keep_tmp )
  {

    std::string _tmpIndexKeepTmpFiles( argv[ index_keep_tmp + 1 ] ) ;
    if ( _tmpIndexKeepTmpFiles == "true" )
    {

      keepTmpFiles = true ;

    }
    else if ( _tmpIndexKeepTmpFiles == "false" )
    {

      keepTmpFiles = false ;

    }
    else
    {

      std::cout << "Argument of -ktf must be either \"true\" or \"false\" "
                << std::endl ;
      exit( 1 ) ;

    }

  }

  ///////////////// Time out for centroids computations and SBR //////////////////
  if ( index_time_out )
  {

    time_out = std::stof( argv[ index_time_out + 1 ] ) ;
    if ( time_out <= 0 )
    {

      std::cout << "ERROR : -timeOut must be greater than 0, got "
                << time_out << std::endl ;
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
  if ( verbose )
  {
    
    std::cout << "\ttoleranceThrComputeNeighborhood : "
                               << toleranceThrComputeNeighborhood << std::endl ;

  }



  //////////////// Preparing atlas bundles paths (tmpAtlasDir) /////////////////
  if ( verbose )
  {

    std::cout << "#########################################################\n" ;
    std::cout << "############# Preparing atlas bundles paths #############"
                                                                  << std::endl ;
    std::cout << "#########################################################\n" ;

  }

  checkAtlasDirectory( atlasDirectory, format ) ;
  // std::vector<std::string> atlasBundleDirectories ;
  // std::string tmpAtlasDir = getAtlasBunldesPaths( outputDirectory,
  //                                                 atlasDirectory,
  //                                                 format,
  //                                                 atlasBundleDirectories ) ;
  std::vector<std::string> atlasBundleDirectories =
                    getFilesInDirectoryWithExtension( atlasDirectory, format ) ;

  if ( verbose )
  {

    std::cout << "Done " << std::endl ;

  }


  //////////////// Getting atlas neighborhood if input is given ////////////////
  std::vector<std::string> atlasNeighborhoodPaths ;
  if ( isAtlasNeighborhood )
  {

    getAtlasNeighborhoodCentroids( atlasNeighborhoodDirectory,
                                   atlasBundleDirectories,
                                   format,
                                   atlasNeighborhoodPaths ) ;

  }
  // else
  // {
  //
  //   std::cout << "ERROR : atlas neighborhood must be given using -an"
  //                                                                 << std::endl ;
  //   exit( 1 ) ;
  //
  // }

  ///////////////// Getting atlas centroids if input is given //////////////////
  std::vector<std::string> atlasNeighborhoodCentroidsPaths ;
  if ( isAtlasNeighborhoodCentroids )
  {

    getAtlasNeighborhoodCentroids( atlasNeighborhoodCentroidsDirectory,
                                   atlasBundleDirectories,
                                   format,
                                   atlasNeighborhoodCentroidsPaths ) ;

  }
  // else
  // {
  //
  //   std::cout << "ERROR : atlas neighborhood centroids must be given using -anc"
  //                                                                 << std::endl ;
  //   exit( 1 ) ;
  //
  // }


  ///////////////////////////////// Global SLR /////////////////////////////////
  std::string movedTractogram = inputBundlesFilename ;
  if ( doSLR )
  {

    if ( verbose )
    {

      std::cout << "#########################################################"
                << std::endl ;
      std::cout << "Doing global SLR..." << std::endl ;
      std::cout << "#########################################################"
                << std::endl ;
    }

    if ( format == ".bundles" && format == ".bundlesdata" )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ProjectAtlasGeoLab : global SLR not compatible "
                    << "with .bundles/.bundlesdata format\n" ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


    }

    std::ostringstream tmpSLRdirOss ;
    tmpSLRdirOss << outputDirectory << "tmpSLR/" ;
    std::string tmpSLRdir = tmpSLRdirOss.str() ;
    if ( !is_dir( tmpSLRdir ) )
    {

      mkdir( tmpSLRdir ) ;

    }


    std::ostringstream movedTractogramTrkOss ;
    movedTractogramTrkOss << tmpSLRdir << "moved" << format ;
    std::string movedTractogramTrk = movedTractogramTrkOss.str() ;

    std::ostringstream globalSLROss ;
    // globalSLROss << "dipy_slr "
    //              << fullAtlasFilename << " "
    //              << inputBundlesFilename << " "
    //              << "--greater_than 10 "
    //              << "--less_than 200 "
    //              << "--out_dir " << tmpSLRdir << " "
    //              << "--out_moved " << movedTractogramTrk << " "
    //              << "--force " ;
    globalSLROss << "dipy_slr "
                 << fullAtlasFilename << " "
                 << inputBundlesFilename << " "
                 << "--out_dir " << tmpSLRdir << " "
                 << "--out_moved " << movedTractogramTrk << " "
                 << "--force " ;
      std::string globalSLR = globalSLROss.str() ;

      int globaleSLRfail = 0 ;
      if ( is_file( movedTractogramTrk ) && !force )
      {

        if ( verbose > 1 )
        {

          std::cout << "WARNING : output file of dipy_slr : "
                    << movedTractogramTrk << " already exists and -force flag "
                    << "was not used, trying computations with existing file"
                    << std::endl ;

        }

      }
      else
      {

        globaleSLRfail = run_sh_process( globalSLR ) ;

      }
      if ( is_file( movedTractogramTrk ) )
      {

        globaleSLRfail = 0 ;


      }
      if ( globaleSLRfail )
      {

        std::cout << "ERROR : could not compute global SLR, got exit code "
                  << globaleSLRfail << std::endl ;
        // if ( is_dir( outputDirectory ) )
        // {
        //
        //   rmdir( outputDirectory ) ;
        //
        // }
        exit( 1 ) ;

      }

      std::ostringstream movedTractogramOss ;
      movedTractogramOss << tmpSLRdir << "moved" << format ;
      std::string movedTractogram = movedTractogramOss.str() ;

      if ( verbose )
      {

        std::cout << "Done" << std::endl ;

      }


  }


  ////////////////////// Computing neighborhood tractogram /////////////////////
  if ( verbose )
  {

    std::cout << "#########################################################\n" ;
    std::cout << "########### Computing tractogram neighborhood ###########"
                                                                  << std::endl ;
    std::cout << "#########################################################\n" ;

  }

  const auto start_time_tract_neigh = std::chrono::system_clock::now() ;

  std::ostringstream tmpNeighborhoodDirOss ;
  tmpNeighborhoodDirOss << outputDirectory << "tmp_neighborhood_dir/" ;
  std::string tmpNeighborhoodDir = tmpNeighborhoodDirOss.str() ;

  if ( !is_dir( tmpNeighborhoodDir ) )
  {

    mkdir( tmpNeighborhoodDir ) ;

  }

  std::ostringstream computeNeighborhoodCommandOss ;
  computeNeighborhoodCommandOss << computeNeighborhoodFile << " "
                                << "-i " << movedTractogram << " "
                                << "-a " << atlasDirectory << " "
                                << "-o " << tmpNeighborhoodDir << " "
				                        << "-minLen " << 10 << " "
				                        // << "-maxLen " << 100 << " "
				                        << "-maxLen " << 200 << " "
                                << "-tolThr " << toleranceThrComputeNeighborhood
                                                                          << " "
                                // << "-nbThreads " << nbCores << " "
                                << "-nbThreads " << nbThreadsCN << " "
                                << "-v " ;
  std::string computeNeighborhoodCommand = computeNeighborhoodCommandOss.str() ;
  int isNeighborhoodFail = 0 ;
  if ( countFilesDirectory( tmpNeighborhoodDir ) > 5 && !force )
  {

    if ( verbose > 1 )
    {

      std::cout << "WARNING : output directory of " << computeNeighborhoodFile
                << " : " <<  tmpNeighborhoodDir << " exists with more than 5 "
                << "files and  the -force flag was not used, trying "
                << "computations with existing directory" << std::endl ;

    }
    

  }
  else
  {

    isNeighborhoodFail = run_sh_process( computeNeighborhoodCommand ) ;

  }

  if ( isNeighborhoodFail )
  {

    std::cout << "ERROR : could not compute neighborhood of atlas bundles "
              << "in tractogram, got exit code "
              << isNeighborhoodFail << std::endl ;

    // if ( is_dir( outputDirectory ) )
    // {
    //
    //   rmdir( outputDirectory ) ;
    //
    // }
    exit( 1 ) ;

  }


  std::vector<std::string> neighborhoodFilenames ;
  getNeighborhoodFilenames( tmpNeighborhoodDir,
                            atlasBundleDirectories,
                            format,
                            neighborhoodFilenames ) ;


  const std::chrono::duration< double > duration_tract_neigh =
                     std::chrono::system_clock::now() - start_time_tract_neigh ;

  if ( verbose )
  {

    std::cout << "Duration : " << duration_tract_neigh.count() << std::endl ;
    std::cout << "Done" << std::endl ;


  }

  ///////////////////////// Computing neighborhood atlas ///////////////////////
  std::vector<std::string> neighborhoodAtlasFilenames ;
  std::string tmpNeighborhoodAtlasDir ;
  if ( ( !isAtlasNeighborhoodCentroids || !isAtlasNeighborhood ) &&
                                                                   isFullAtlas )
  {

    if ( verbose )
    {

      std::cout << "#########################################################"
                << std::endl ;
      std::cout << "Computing atlas neighborhood " << std::endl ;
      std::cout << "#########################################################"
                << std::endl ;

    }

    const auto start_time_atlas_neigh = std::chrono::system_clock::now() ;

    std::ostringstream tmpNeighborhoodAtlasDirOss ;
    tmpNeighborhoodAtlasDirOss << outputDirectory
                               << "tmp_neighborhood_atlas_dir/" ;
    tmpNeighborhoodAtlasDir = tmpNeighborhoodAtlasDirOss.str() ;

    if ( !is_dir( tmpNeighborhoodAtlasDir ) )
    {

      mkdir( tmpNeighborhoodAtlasDir ) ;

    }

    std::ostringstream computeAtlasNeighborhoodCommandOss ;
    computeAtlasNeighborhoodCommandOss << computeNeighborhoodFile << " "
                        << "-i " << fullAtlasFilename << " "
                        << "-a " << atlasDirectory << " "
                        << "-o " << tmpNeighborhoodAtlasDir << " "
			<< "-minLen " << 10 << " "
                        << "-maxLen " << 400 << " "
                        << "-tolThr " << toleranceThrComputeNeighborhood << " "
                        // << "-nbThreads " << nbCores << " "
                        << "-nbThreads " << nbThreadsCN << " "
                        << "-v " ;
    std::string computeAtlasNeighborhoodCommand =
                                      computeAtlasNeighborhoodCommandOss.str() ;
    int isAtlasNeighborhoodFail = 0 ;
    if ( countFilesDirectory( tmpNeighborhoodAtlasDir ) > 5 && !force )
    {

      if ( verbose > 1 )
      {

        std::cout << "WARNING : output directory of " << computeNeighborhoodFile
                  << " : " <<  tmpNeighborhoodAtlasDir << " already exists with"
                  << " more than 5 files and the -force flag was not used, "
                  << "trying computations with existing directory"
                  << std::endl ;

      }

    }
    else
    {

      isAtlasNeighborhoodFail = run_sh_process(
                                              computeAtlasNeighborhoodCommand ) ;

    }

    if ( isAtlasNeighborhoodFail )
    {
      std::cout << "ERROR : could not compute neighborhood of atlas bundles "
                << "in full atlas, got exit code " << isAtlasNeighborhoodFail
                                                                  << std::endl ;

      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }
      exit( 1 ) ;

    }

    getNeighborhoodFilenames( tmpNeighborhoodAtlasDir,
                              atlasBundleDirectories,
                              format,
                              neighborhoodAtlasFilenames ) ;


   const std::chrono::duration< double > duration_atlas_neigh =
                     std::chrono::system_clock::now() - start_time_atlas_neigh ;

    if ( verbose )
    {

      std::cout << "Duration : " << duration_atlas_neigh.count() << std::endl ;
      std::cout << "Done" << std::endl ;

    }

  }


  




  std::cout << "#########################################################\n" ;
  std::cout << "####################### Local SBR #######################"
                                                                  << std::endl ;
  std::cout << "#########################################################\n" ;
  const auto start_time_sbr = std::chrono::system_clock::now() ;

  // Launching dipy service
  std::cout << "Launching dipy service... " ;

  std::ostringstream serverLogFilePathOss ;
  serverLogFilePathOss << outputDirectory << "dipyServiceLog.txt" ;
  std::string serverLogFilePath = serverLogFilePathOss.str() ;


  std::ostringstream launchDipyServiceOss ;
  launchDipyServiceOss << openDipyServerClientFile << " " ;
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


  int nbBundles = 0 ;
  if ( !isAtlasNeighborhoodCentroids )
  {

    nbBundles = neighborhoodAtlasFilenames.size() ;

  }
  else
  {

    nbBundles = atlasNeighborhoodCentroidsPaths.size() ;

  }


  std::vector<std::vector<int>> indecesInNeighborhoodsReconized( nbBundles ) ;

  int nbFibersTractogram = getNbFibers( movedTractogram ) ;
  int nbBundlesProcessed = 1 ;
  while ( c.running() )
  {

    omp_set_num_threads( nbThreads ) ;
    #pragma omp parallel for schedule(dynamic)
    for ( int i = 0 ; i < nbBundles ; i++ )
    {

      std::string atlasBundleDirectory = atlasBundleDirectories[ i ] ;
      std::string movedTractogramNeighborhood = neighborhoodFilenames[ i ] ;
      std::string atlasNeighborhoodFile ;
      if ( isAtlasNeighborhood )
      {

        atlasNeighborhoodFile = atlasNeighborhoodPaths[ i ] ;


      }
      else
      {

        atlasNeighborhoodFile = neighborhoodAtlasFilenames[ i ] ;

      }

      std::string atlasNeighborhoodCentroidsFile ;
      if ( isAtlasNeighborhoodCentroids )
      {

        atlasNeighborhoodCentroidsFile = atlasNeighborhoodCentroidsPaths[ i ] ;


      }
      else
      {

        atlasNeighborhoodCentroidsFile = "" ;

      }

      if ( endswith( atlasNeighborhoodFile, ".minf" ) )
      {

        atlasNeighborhoodFile = replaceExtension( atlasNeighborhoodFile,
                                                                      format ) ;

      }


      applyLocalSBR( movedTractogramNeighborhood,
                     atlasBundleDirectory,
                     atlasNeighborhoodFile,
                     atlasNeighborhoodCentroidsFile,
                     outputDirectory,
                     referenceFilename,
                     format,
                     nbPointsPerFiber,
                     portDipyServer,
                     time_out,
                     verbose ) ;

      #pragma omp critical
      {

        printf( "\rNumber of bundles processed : [ %d  /  %d ]",
                                               nbBundlesProcessed, nbBundles ) ;
        std::cout << "" << std::flush ;

        nbBundlesProcessed++ ;

      }

    }

    // Closing dipy service
    std::cout << "\n" ;
    closeDipyServer( portDipyServer ) ;
    std::cout << "Done" << std::endl ;

    break ;

  }

  const std::chrono::duration< double > duration_sbr =
                          std::chrono::system_clock::now() - start_time_sbr ;

  if ( verbose )
  {

    std::cout << "Duration : " << duration_sbr.count() << std::endl ;

    std::cout << "Done" << std::endl ;

  }


  /*
  std::cout << "#########################################################\n" ;
  std::cout << "################## Cleaning temp files ##################"
                                                                  << std::endl ;
  std::cout << "#########################################################\n" ;


  if ( is_dir( tmpNeighborhoodDir) )
  {

    rmdir( tmpNeighborhoodDir ) ;

  }


  if ( !isAtlasNeighborhoodCentroids && isFullAtlas )
  {

    if ( is_dir( tmpNeighborhoodAtlasDir ) )
    {

      rmdir( tmpNeighborhoodAtlasDir ) ;

    }

  }
  */

  const std::chrono::duration< double > duration =
                              std::chrono::system_clock::now() - start_time ;

  if ( verbose )
  {

    std::cout << "Duration : " << duration.count() << std::endl ;

  }

  std::cout << "Done" << std::endl ;

  return( 0 ) ;

}
