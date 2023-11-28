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

#include "ProjectAtlasGeoLab.h"
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
void applyGeoLab( const std::string& movedTractogramNeighborhood,
                  const std::string& atlasBundleDirectory,
                  const std::string& atlasNeighborhoodFile,
                  const std::string& atlasNeighborhoodCentroidsFile,
                  const std::string& outputDirectory,
                  const std::string& referenceImage,
                  const std::string& format,
                  const int minimumNumberFibers,
                  std::vector<int>& indexInNeighborhoodRecognized,
                  float adjacency_classic,
                  int nbFibersClassic,
                  int nbPointsPerFiber,
                  int portDipyServer,
                  bool& keepClassic,
                  float& coverageGeoLab,
                  float& adjacencyGeoLab,
                  float& overlapGeoLab,
                  float& disimilarityGeoLab,
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

  coverageGeoLab = -1 ;
  adjacencyGeoLab = -1 ;
  overlapGeoLab = -1 ;
  disimilarityGeoLab = -1 ;

  if ( adjacency_classic >= 0.70 && nbFibersClassic >= minimumNumberFibers )
  {

    if ( verbose > 1 )
    {

      std::cout << "Keeping classic projection for bundle " << bundleName
                << std::endl ;

    }

    keepClassic = true ;

    return ;

  }
  else
  {

    if ( verbose > 1 )
    {

      std::cout << "Testing GeoLab projection for bundle "
                << bundleName << std::endl ;

    }

  }


  std::string atlasNeighborhood ;
  if ( !( endswith( atlasNeighborhoodFile, format ) ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ERROR : in applyGeoLab, atlasNeighborhoodFile must "
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

    std::cout << "ERROR : in applyGeoLab, atlas neighborhood file "
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

      // std::string atlasNeighborhoodCentroidsMinf = replaceExtension(
      //                                    atlasNeighborhoodCentroids, ".minf" ) ;
      //
      // std::string atlasNeighborhoodMinf = replaceExtension(
      //                                             atlasNeighborhood, ".minf" ) ;
      //
      // copy( atlasNeighborhoodMinf, atlasNeighborhoodCentroidsMinf ) ;



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


  // Projection
  std::ostringstream atlasInfoPathOss ;
  if ( is_dir( atlasBundleDirectory ) )
  {

    atlasInfoPathOss << atlasBundleDirectory << bundleName ;

    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      atlasInfoPathOss << ".bundles" ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      atlasInfoPathOss << ".minf" ;

    }

  }
  else
  {

    atlasInfoPathOss << atlasBundleDirectory ;

  }
  std::string atlasInfoPath = atlasInfoPathOss.str() ;
  float thrDistanceBetweenMedialPointsBundle = thrDistanceBetweenMedialPoints ;
  if ( !index_thrDBMP )
  {

    thrDistanceBetweenMedialPointsBundle =
                  getMaximumDistanceBetweenMedialPoints( atlasInfoPath ) ;
    thrDistanceBetweenMedialPointsBundle *= ( 1 +
                                        toleranceDistanceBetweenMedialPoints ) ;

  }


  RecognizedBundles recognized( 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                true, useMDF, true, useSimple,
                                toleranceP,
                                toleranceThr,
                                toleranceMaxAngle,
                                toleranceMaxDirectionAngle,
                                toleranceMinShapeAngle,
                                toleranceMaxShapeAngle,
                                toleranceLenght,
                                0.0,
                                thrPercentageSimilarity,
                                thrDistanceBetweenMedialPointsBundle,
                                minimumNumberFibers,
                                adjacencyForCompareBundles,
                                false, false, false, false, false, false, true,
                                false, false, saveBundlesSeparetly, false,
			                          false,
                                nbThreads ) ;
  
  if ( useDefaultThr )
  {

    recognized.useDefaultThr = true ;
    recognized.thrDistance = thrDistance ;


  }

  BundlesData atlasBundleData( atlasBundleDirectory.c_str() ) ;
  BundlesMinf atlasBundleInfo( atlasBundleDirectory.c_str() ) ;
  BundlesData subjectBundlesData( neighborhoodRegistered.c_str() ) ;

  recognized.projectBundle( atlasBundleData,
                            atlasBundleInfo,
                            subjectBundlesData,
                            indexInNeighborhoodRecognized,
                            coverageGeoLab,
                            adjacencyGeoLab,
                            overlapGeoLab,
                            disimilarityGeoLab,
                            useMeanForMDAD,
                            true ) ;


  if ( adjacency_classic > adjacencyGeoLab &&
                                         nbFibersClassic > minimumNumberFibers )
  {

    keepClassic = true ;

  }
  else
  {

    keepClassic = false ;

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
  index_thr = getFlagPosition( argc, argv, "-thr" ) ;
  index_thrCov = getFlagPosition( argc, argv, "-thrCov" ) ;
  index_thrAdj = getFlagPosition( argc, argv, "-thrAdj" ) ;
  index_thrDBMP = getFlagPosition( argc, argv, "-thrDBMP" )  ;
  index_tolP = getFlagPosition( argc, argv, "-tolP" )  ;
  index_tolThr = getFlagPosition( argc, argv, "-tolThr" )  ;
  index_tolMaxAngle = getFlagPosition( argc, argv, "-tolMaxAng" )  ;
  index_tolMaxDirectionAngle = getFlagPosition( argc, argv, "-tolMaxDirAng" )  ;
  index_tolMinShapeAngle = getFlagPosition( argc, argv, "-tolMinShapeAng" )  ;
  index_tolMaxShapeAngle = getFlagPosition( argc, argv, "-tolMaxShapeAng" )  ;
  index_tolLenght = getFlagPosition( argc, argv, "-tolLenght" )  ;
  index_tolThrCN = getFlagPosition( argc, argv, "-tolThrCN" ) ;
  index_tolDistBetMedPts = getFlagPosition( argc, argv, "-tolDBMP" )  ;
  index_minNbFibers = getFlagPosition( argc, argv, "-minNbFibers" ) ;
  index_thrSim = getFlagPosition( argc, argv, "-thrSim" ) ;
  index_adjCB = getFlagPosition( argc, argv, "-adjCB" ) ;
  index_pa = getFlagPosition( argc, argv, "-pa" ) ;
  index_cv = getFlagPosition( argc, argv, "-cv" ) ;
  index_cn = getFlagPosition( argc, argv, "-cn" ) ;
  index_slr = getFlagPosition( argc, argv, "-slr" ) ;
  index_cp = getFlagPosition( argc, argv, "-cp" )  ;
  index_mdf = getFlagPosition( argc, argv, "-mdf" )  ;
  index_simple = getFlagPosition( argc, argv, "-simple" )  ;
  index_sp = getFlagPosition( argc, argv, "-sp" ) ;
  index_force = getFlagPosition( argc, argv, "-force" ) ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_nbThreadsCN = getFlagPosition( argc, argv, "-nbThreadsCN" ) ;
  index_keep_tmp = getFlagPosition( argc, argv, "-ktf" ) ;
  index_time_out = getFlagPosition( argc, argv, "-timeOut" ) ;
  index_useMeanForMDAD = getFlagPosition( argc, argv, "-useMeanMDAD" ) ;
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
              << "[-thr] : Threshold (in mm) of the MDA distance to decide if "
              << "the fiber belong to a bundle in the atlas (Default = 10 mm) \n"
              << "[-thrCov] : Threshold to keep bundles where coverage is "
              << "greater than thrCov (default : 0 -> keep all bundles ) \n"
              << "[-thrDBMP] : Threshold for maximum distance between medial "
              << "points (default : based on bundles distribution) \n"
              << "[-tolP] : Tolerance for parameter p (for advanced users, "
              << "default = 0) \n"
              << "[-tolThr] : Tolerance for parameter thr (for advanced users, "
              << "default = 0) \n"
              << "[-tolMaxAng] : Tolerance for parameter max angle (for "
              << "advanced users,  default = 0.0) \n"
              << "[-tolMaxDirAng] : Tolerance for parameter max direction angle"
              << " (for advanced users,  default = 0.0) \n"
              << "[-tolMinShapeAng] : Tolerance for parameter min shape angle"
              << " (for advanced users,  default = 0.0) \n"
              << "[-tolMaxShapeAng] : Tolerance for parameter max shape angle"
              << " (for advanced users,  default = 0.0) \n"
              << "[-tolLenght] : Tolerance for parameter lenght (for advanced "
              << "users, default = 0.0) \n"
              << "[-tolThrCN] : tolerance for computeNeighborhood threshold "
              << "(default = 0.0) \n"
              << "[-tolDBMP] : Tolerance for distance between medial points "
              << "(for advanced users, default = -0.9) \n"
              << "[-thrAdj] : keep bundle with adjacency greater than given "
              << " value (default : 0 -> keep all bundles ) \n"
              << "[-minNbFibers] : Minimum number of fiber to consider a bundle"
              << " recognized ( default : 20 )\n"
              << "[-thrSim] : Threshold for percentage of similarity in "
              << "projection i.e nbAtlasBundleFibersSimilar / "
              << "nbFibersAtlasBundle (Range [ 0 ; 1 ], default = 0.00001) \n"
              << "[-adjCB] : adjacency for computing comparison with atlas "
              << "(default = 5 mm) \n"
              << "[-pa] : Path to the ProjectAtlas file \n"
              << "[-cv] : Path to the ConvertBundleFormat file \n"
              << "[-cn] : Path to the computeNeighborhood file \n"
              << "[-slr] : Do global SLR step (default : false)\n"
              << "[-cp] : Do first a classical projection without SBR (default "
	            << ": true)\n"
	            << "[-mdf]: Use MDF distance\n"
	            << "[-simple]: Do projection only using distance and length\n"
              << "[-sp] : Save recognized bundles separetly (default : true)\n"
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
              << "[-useMeanMDAD] : Use mean instead of max for MDAD (default : true ) \n"
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

  ////////////////////// Threshold on the MDAD/MDF distance //////////////////////
  if ( index_thr )
  {

    thrDistance = std::stof( argv[ index_thr + 1 ] ) ;
    if ( thrDistance <= 0 )
    {

      std::cout << "Error argument: thrDistance must be positive ( "
                << "thrDistance > 0 ), got thrDistance = " << thrDistance
                << std::endl ;
      exit( 1 ) ;

    }
    useDefaultThr = true ;


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


  ///////////////////////////// Coverage threshold /////////////////////////////
  if ( index_thrCov )
  {

    coverageThreshold = std::stof( argv[ index_thrCov + 1 ] ) ;

  }

  ///////////////////////////// Adjacency threshold ////////////////////////////
  if ( index_thrAdj )
  {

    adjacencyThreshold = std::stof( argv[ index_thrAdj + 1 ] ) ;

  }

  ////////////////// Threshold distance between medial points //////////////////
  if ( index_thrDBMP )
  {

    thrDistanceBetweenMedialPoints = std::stof( argv[ index_thrDBMP + 1 ] ) ;

    if ( thrDistanceBetweenMedialPoints < 0 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : thrDBMP must be positive "
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }


  }

  /////////////////////////////// Tolerances ///////////////////////////////////

  if ( index_tolP )
  {

    toleranceP = std::stof( argv[ index_tolP + 1 ] ) ;

    // if ( toleranceP < -1 || toleranceP > 1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolP must be in [-1;1]" << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    //
    // }

  }

  if ( index_tolThr )
  {

    toleranceThr = std::stof( argv[ index_tolThr + 1 ] ) ;

    // if ( toleranceThr < -1 || toleranceThr > 1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolThr must be in [-1;1]"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    // }

  }

  if ( index_tolMaxAngle )
  {

    toleranceMaxAngle = std::stof( argv[ index_tolMaxAngle + 1 ] ) ;

    // if ( toleranceMaxAngle < -1 || toleranceMaxAngle > 1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolMaxAng must be in [-1;1]"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    // }

  }

  if ( index_tolMaxDirectionAngle )
  {

    toleranceMaxDirectionAngle = std::stof( argv[
                                            index_tolMaxDirectionAngle + 1 ] ) ;

    // if ( toleranceMaxDirectionAngle < -1 || toleranceMaxDirectionAngle > 1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolMaxDirAng must be in [-1;1]"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    // }

  }

  if ( index_tolMinShapeAngle )
  {

    toleranceMinShapeAngle = std::stof( argv[ index_tolMinShapeAngle + 1 ] ) ;

    // if ( toleranceMinShapeAngle < -1 || toleranceMinShapeAngle > 1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolMinShapeAng must be in [-1;1]"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    //
    //
    // }

  }

  if ( index_tolMaxShapeAngle )
  {

    toleranceMaxShapeAngle = std::stof( argv[ index_tolMaxShapeAngle + 1 ] ) ;

    // if ( toleranceMaxShapeAngle < -1 || toleranceMaxShapeAngle > 1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolManShapeAng must be in [-1;1]"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    // }

  }

  if ( index_tolLenght )
  {

    toleranceLenght = std::stof( argv[ index_tolLenght + 1 ] ) ;

    // if ( toleranceLenght < -1 || toleranceLenght > 1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolLenght must be in [-1;1]"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    // }

  }


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


  if ( index_tolDistBetMedPts )
  {

    toleranceDistanceBetweenMedialPoints = std::stof(
                                          argv[ index_tolDistBetMedPts + 1 ] ) ;

    // if ( toleranceDistanceBetweenMedialPoints < -1 )
    // {
    //
    //   std::stringstream outMessageOss ;
    //   outMessageOss << "Error argument : -tolDBMP must be greater than -1"
    //                                                               << std::endl ;
    //   std::string outMessage = outMessageOss.str() ;
    //   throw( std::invalid_argument( outMessage ) ) ;
    //
    // }

  }


  /////////////////////////// Minimum number of fibers /////////////////////////
  if ( index_minNbFibers )
  {

    minimumNumberFibers = std::stoi( argv[ index_minNbFibers + 1 ] ) ;
    if ( minimumNumberFibers < 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : argument -minNbFibers must be > 0 "
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }

  ///////////////////////////// Threshold similarity ///////////////////////////
  if ( index_thrSim )
  {

    thrPercentageSimilarity = std::stof( argv[ index_thrSim + 1 ] ) ;
    if ( thrPercentageSimilarity <= 0 || thrPercentageSimilarity > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : argument -thrSim must be greater than 0 and "
                    << "lower or equal to 1 " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    thrPercentageSimilarity = 0.00001 ;

  }

  ////////////////////// Adjacency threshold for comparison ////////////////////
  if ( index_adjCB )
  {

    adjacencyForCompareBundles = std::stof( argv[ index_adjCB + 1 ] ) ;
    if ( adjacencyForCompareBundles <= 0 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : argument -adjCB must be greater than 0 "
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    adjacencyForCompareBundles = 5.0 ;

  }



  ///////////////////////////// Project atlas file /////////////////////////////
  if ( index_pa )
  {

    std::string tmpFile( argv[ index_pa + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : projecAtlas file " << tmpFile
                                           << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }


    }

    projectAtlasFile = tmpFile ;

  }

  ////////////////////////// ConvertBundleFormat file //////////////////////////
  if ( index_cv )
  {

    std::string tmpFile( argv[ index_cv + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : ConvertBundleFormat file " << tmpFile
                    << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }

    }

    convertBundleFormatsFile = tmpFile ;

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

  //////////////////////////// registerBunldes file ////////////////////////////
  if ( index_rb )
  {

    std::string tmpFile( argv[ index_rb + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : registerBunldes file " << tmpFile
                                           << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }

    }

    registerBundlesClientFile = tmpFile ;

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

  //////////////////// Do first a classical projection without SBR /////////////
  if ( index_cp )
  {

    std::string _tmpIndexCP( argv[ index_cp + 1 ] ) ;
    if ( _tmpIndexCP == "true" )
    {

      doClassical = true ;

    }
    else if ( _tmpIndexCP == "false" )
    {

      doClassical = false ;

    }
    else
    {

      std::cout << "Argument of -cp must be either \"true\" or \"false\" "
                << std::endl ;
      exit( 1 ) ;

    }


  }

  ///////////////////////////////// Use MDF distance //////////////////////////
  if ( index_mdf )
  {

    std::string _tmpIndexMDF( argv[ index_mdf + 1 ] ) ;
    if ( _tmpIndexMDF == "true" )
    {

      useMDF = true ;

    }
    else if ( _tmpIndexMDF == "false" )
    {

      useMDF = false ;

    }
    else
    {

      std::cout << "Argument of -mdf must be either \"true\" or \"false\" "
                << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( useMDF )
  {

    std::cout << "WARNING : Using MDF distance" << std::endl ;

  }

  ////////////////////////////// Use simple projection ////////////////////////
  if ( index_simple )
  {

    std::string _tmpIndexSimple( argv[ index_simple + 1 ] ) ;
    if ( _tmpIndexSimple == "true" )
    {

      useSimple = true ;

    }
    else if ( _tmpIndexSimple == "false" )
    {

      useSimple = false ;

    }
    else
    {

      std::cout << "Argument of -simple must be either \"true\" or \"false\" "
                << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( useSimple )
  {

    std::cout << "WARNING : Using simple projection" << std::endl ;

  }

  //////////////////////// Separate recognized bundles /////////////////////////
  if ( index_sp )
  {

    std::string _tmpIndexSP( argv[ index_sp + 1 ] ) ;
    if ( _tmpIndexSP == "true" )
    {

      saveBundlesSeparetly = true ;

    }
    else if ( _tmpIndexSP == "false" )
    {

      saveBundlesSeparetly = false ;

    }
    else
    {

      std::cout << "Argument of -sp must be either \"true\" or \"false\" "
                << std::endl ;
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


  /////////////////////// Use mean instead of max for MDAD ////////////////////////

  if ( index_useMeanForMDAD )
  {

    std::string _tmpUseMeanForMDAD( argv[ index_useMeanForMDAD + 1 ] ) ;
    if ( _tmpUseMeanForMDAD == "true" )
    {

      useMeanForMDAD = true ;

    }
    else if ( _tmpUseMeanForMDAD == "false" )
    {

      useMeanForMDAD = false ;

    }
    else
    {

      std::cout << "Argument of -useMeanMDAD must be either \"true\" or \"false\" "
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
  if ( verbose )
  {
    std::cout << "Do first a classical projection without SBR : " << doClassical
	                                                          << std::endl ;
    std::cout << "Tolerance on parameters : " << std::endl ;
    std::cout << "\ttoleranceP : " << toleranceP << std::endl ;
    std::cout << "\ttoleranceThr : " << toleranceThr << std::endl ;
    std::cout << "\ttoleranceMaxAngle : " << toleranceMaxAngle << std::endl ;
    std::cout << "\ttoleranceMaxDirectionAngle : "
                                    << toleranceMaxDirectionAngle << std::endl ;
    std::cout << "\ttoleranceMinShapeAngle : " << toleranceMinShapeAngle
                                                                  << std::endl ;
    std::cout << "\ttoleranceMaxShapeAngle : " << toleranceMaxShapeAngle
                                                                  << std::endl ;
    std::cout << "\ttoleranceLenght : " << toleranceLenght << std::endl ;
    std::cout << "\ttoleranceDistanceBetweenMedialPoints : "
                          << toleranceDistanceBetweenMedialPoints << std::endl ;
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
                        << "-maxLen " << 200 << " "
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


  //////////////////////// Projecting atlas without SBR ////////////////////////
  const auto start_time_no_sbr = std::chrono::system_clock::now() ;
  if ( verbose && doClassical )
  {

    std::cout << "#########################################################\n" ;
    std::cout << "############## Projecting atlas without SBR #############"
                                                                  << std::endl ;
    std::cout << "#########################################################\n" ;

  }

  int nbBundles = 0 ;
  if ( !isAtlasNeighborhoodCentroids )
  {

    nbBundles = neighborhoodAtlasFilenames.size() ;

  }
  else
  {

    nbBundles = atlasNeighborhoodCentroidsPaths.size() ;

  }

  if ( nbBundles == 0 && doClassical )
  {

    std::cout << "ERROR : no valid bundles in atlas directory" << std::endl ;
    exit( 1 ) ;

  }
  else
  {

    if ( doClassical )
    {

      std::cout << "Number of bundles in atlas : " << nbBundles << std::endl ;

    }


  }

  std::vector<std::vector<int>> indecesInNeighborhoodsReconizedClassic(
                                                                   nbBundles ) ;
  std::vector<float> coveragesClassic( nbBundles, -1 ) ;
  std::vector<float> adjacenciesClassic( nbBundles, -1 ) ;
  std::vector<float> overlapsClassic( nbBundles, -1 ) ;
  std::vector<float> disimilaritiesClassic( nbBundles, -1 ) ;
  std::vector<std::string> bundlesNames( nbBundles ) ;


  #pragma omp parallel for schedule(dynamic)
  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    std::string atlasBundleDirectory = atlasBundleDirectories[ i ] ;
    bundlesNames[ i ] = basenameNoExtension( atlasBundleDirectory ) ;
    std::string movedTractogramNeighborhood = neighborhoodFilenames[ i ] ;

    std::string _tmpAtlasBundlePath = atlasBundleDirectory ;

    char lastChar = _tmpAtlasBundlePath[ _tmpAtlasBundlePath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      _tmpAtlasBundlePath = _tmpAtlasBundlePath.substr( 0,
                                              _tmpAtlasBundlePath.size() - 1 ) ;

    }

    std::experimental::filesystem::path tmpPath( _tmpAtlasBundlePath ) ;
    std::string bundleName = tmpPath.stem() ;


    std::ostringstream atlasInfoPathOss ;
    if ( is_dir( atlasBundleDirectory ) )
    {

      atlasInfoPathOss << atlasBundleDirectory << bundleName ;

      if ( format == ".bundles" || format == ".bundlesdata" )
      {

        atlasInfoPathOss << ".bundles" ;

      }
      else if ( format == ".trk" || format == ".tck" )
      {

        atlasInfoPathOss << ".minf" ;

      }

    }
    else
    {

      atlasInfoPathOss << atlasBundleDirectory ;

    }
    std::string atlasInfoPath = atlasInfoPathOss.str() ;
    float thrDistanceBetweenMedialPointsBundle =
                                                thrDistanceBetweenMedialPoints ;

    if ( !index_thrDBMP )
    {

      thrDistanceBetweenMedialPointsBundle =
                  getMaximumDistanceBetweenMedialPoints( atlasInfoPath ) ;
      thrDistanceBetweenMedialPointsBundle *= ( 1 +
                                        toleranceDistanceBetweenMedialPoints ) ;

    }



    RecognizedBundles recognized( 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  true, useMDF, true, useSimple,
                                  toleranceP,
                                  toleranceThr,
                                  toleranceMaxAngle,
                                  toleranceMaxDirectionAngle,
                                  toleranceMinShapeAngle,
                                  toleranceMaxShapeAngle,
                                  toleranceLenght,
                                  0.0,
                                  thrPercentageSimilarity,
                                  thrDistanceBetweenMedialPointsBundle,
                                  minimumNumberFibers,
                                  adjacencyForCompareBundles,
                                  false, false, false, false, false, false,
                                  true, false, false, saveBundlesSeparetly,
  			                          false, false,
                                  nbThreads ) ;
    
    if ( useDefaultThr )
    {

      recognized.useDefaultThr = true ;
      recognized.thrDistance = thrDistance ;


    }

    BundlesData atlasBundleData( atlasBundleDirectory.c_str() ) ;
    BundlesMinf atlasBundleInfo( atlasBundleDirectory.c_str() ) ;
    BundlesData subjectBundlesData( movedTractogramNeighborhood.c_str() ) ;
    float coverageClassic = -1 ;
    float adjacencyClassic = -1 ;
    float overlapClassic = -1 ;
    float disimilarityClassic = -1 ;
    std::vector<int> indexInNeighborhoodRecognized ;
    if ( doClassical )
    {

      recognized.projectBundle( atlasBundleData,
                                atlasBundleInfo,
                                subjectBundlesData,
                                indexInNeighborhoodRecognized,
                                coverageClassic,
                                adjacencyClassic,
                                overlapClassic,
                                disimilarityClassic,
                                useMeanForMDAD,
                                true ) ;

    }
    #pragma omp critical
    {

      indecesInNeighborhoodsReconizedClassic[ i ] =
                                                 indexInNeighborhoodRecognized ;
      coveragesClassic[ i ] = coverageClassic ;
      adjacenciesClassic[ i ] = adjacencyClassic ;
      overlapsClassic[ i ] = overlapClassic ;
      disimilaritiesClassic[ i ] = disimilarityClassic ;

    }

  }


  const std::chrono::duration< double > duration_no_sbr =
                          std::chrono::system_clock::now() - start_time_no_sbr ;

  if ( verbose && doClassical )
  {

    std::cout << "Duration projection without SBR : " << duration_no_sbr.count()
                                                                  << std::endl ;

    std::cout << "Done" << std::endl ;

  }


  


  std::cout << "#########################################################\n" ;
  std::cout << "########### Projecting atlas bundles with SBR ###########"
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


  std::vector<std::vector<int>> indecesInNeighborhoodsReconized( nbBundles ) ;
  std::vector<bool> keepClassicProjection( nbBundles ) ;
  std::vector<float> coveragesGeoLab( nbBundles, -1 ) ;
  std::vector<float> adjacenciesGeoLab( nbBundles, -1 )  ;
  std::vector<float> overlapsGeoLab( nbBundles, -1 ) ;
  std::vector<float> disimilaritiesGeoLab( nbBundles, -1 ) ;

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

      std::vector<int> indexInNeighborhoodRecognized ;

      float adjacency_classic = adjacenciesClassic[ i ] ;
      int nbFibersClassic = indecesInNeighborhoodsReconizedClassic[ i ].size() ;
      bool keepClassic = true ;
      float coverageGeoLab = -1 ;
      float adjacencyGeoLab = -1 ;
      float overlapGeoLab = -1 ;
      float disimilarityGeoLab = -1 ;

      applyGeoLab( movedTractogramNeighborhood,
                   atlasBundleDirectory,
                   atlasNeighborhoodFile,
                   atlasNeighborhoodCentroidsFile,
                   outputDirectory,
                   referenceFilename,
                   format,
                   minimumNumberFibers,
                   indexInNeighborhoodRecognized,
                   adjacency_classic,
                   nbFibersClassic,
                   nbPointsPerFiber,
                   portDipyServer,
                   keepClassic,
                   coverageGeoLab,
                   adjacencyGeoLab,
                   overlapGeoLab,
                   disimilarityGeoLab,
                   time_out,
                   verbose ) ;

      #pragma omp critical
      {

        printf( "\rNumber of bundles processed : [ %d  /  %d ]",
                                               nbBundlesProcessed, nbBundles ) ;
        std::cout << "" << std::flush ;

        indecesInNeighborhoodsReconized[ i ] = indexInNeighborhoodRecognized ;
        keepClassicProjection[ i ] = keepClassic ;
        coveragesGeoLab[ i ] = coverageGeoLab ;
        adjacenciesGeoLab[ i ] = adjacencyGeoLab ;
        overlapsGeoLab[ i ] = overlapGeoLab ;
        disimilaritiesGeoLab[ i ] = disimilarityGeoLab ;

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

    std::cout << "Duration projection with SBR : " << duration_sbr.count()
                                                                  << std::endl ;

    std::cout << "Done" << std::endl ;

  }

  // Getting labels
  std::cout << "#########################################################\n" ;
  std::cout << "############# Saving labels in subject space ############\n"  ;
  std::cout << "#########################################################\n" ;
  std::vector<std::vector<int16_t>> labelsSubjectSpace( nbFibersTractogram ) ;
  std::vector<float> finalCoverages ;
  std::vector<float> finalAdjacencies ;
  std::vector<float> finalOverlaps ;
  std::vector<float> finalDisimilarities ;
  std::vector<int> finalNbFibersPerBundle ;
  std::vector<std::string> finalBundlesNames ;
  std::vector<std::string> finalBundlesNeigborhoods ;
  std::vector<bool> finalIsClassic ;
  std::vector<std::vector<int>> finalIndecesInNeighborhoodsReconized ;

  // std::vector<float> tmpAdjacenciesClassic ;
  // std::vector<float> tmpAdjacenciesGeoLab ;
  int tmpCpunter = 0 ;
  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    float finaleCoverage = -1 ;
    float finaleAdjacency = -1 ;
    float finaleOverlap = -1 ;
    float finaleDisimilarity = -1 ;
    std::vector<int> finalIndexInNeighborhoodRecognized ;
    bool finalKeepClassic = false ;

    if ( keepClassicProjection[ i ] )
    {

      finaleCoverage = coveragesClassic[ i ] ;
      finaleAdjacency = adjacenciesClassic[ i ] ;
      finaleOverlap = overlapsClassic[ i ] ;
      finaleDisimilarity = disimilaritiesClassic[ i ] ;
      finalIndexInNeighborhoodRecognized =
                                   indecesInNeighborhoodsReconizedClassic[ i ] ;
      finalKeepClassic = true ;

    }
    else
    {

      finaleCoverage = coveragesGeoLab[ i ] ;
      finaleAdjacency = adjacenciesGeoLab[ i ] ;
      finaleOverlap = overlapsGeoLab[ i ] ;
      finaleDisimilarity = disimilaritiesGeoLab[ i ] ;
      finalIndexInNeighborhoodRecognized =
                                          indecesInNeighborhoodsReconized[ i ] ;
      finalKeepClassic = false ;

    }


    if ( finalIndexInNeighborhoodRecognized.size() < minimumNumberFibers ||
                                          finaleAdjacency < adjacencyThreshold )
    {

      continue ;

    }

    // tmpAdjacenciesClassic.push_back( adjacenciesClassic[ i ] ) ;
    // tmpAdjacenciesGeoLab.push_back( adjacenciesGeoLab[ i ] ) ;

    finalCoverages.push_back( finaleCoverage ) ;
    finalAdjacencies.push_back( finaleAdjacency ) ;
    finalOverlaps.push_back( finaleOverlap ) ;
    finalDisimilarities.push_back( finaleDisimilarity ) ;
    finalBundlesNames.push_back( bundlesNames[ i ] ) ;
    finalBundlesNeigborhoods.push_back( neighborhoodFilenames[ i ] ) ;
    finalNbFibersPerBundle.push_back(
                                   finalIndexInNeighborhoodRecognized.size() ) ;
    finalIsClassic.push_back( finalKeepClassic ) ;
    finalIndecesInNeighborhoodsReconized.push_back(
                                          finalIndexInNeighborhoodRecognized ) ;


    std::string movedTractogramNeighborhood = neighborhoodFilenames[ i ] ;

    std::string neighborhoodTractogramBinPath = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "Index.bin" ) ;

    int nbFibersMovedTracotgramNeighborhood = getNbFibers(
                                                 movedTractogramNeighborhood ) ;

    std::vector<int64_t> indexInTractogram ;
    readIndexInTractogram( neighborhoodTractogramBinPath.c_str(),
                           indexInTractogram,
                           nbFibersMovedTracotgramNeighborhood ) ;


    for ( int fiberIndex : finalIndexInNeighborhoodRecognized )
    {

      int64_t tmpIndexInTractogram = indexInTractogram[ fiberIndex ] ;

      labelsSubjectSpace[ tmpIndexInTractogram ].push_back( tmpCpunter ) ;

    }

    tmpCpunter++ ;

  }

  std::cout << "Number of bundles recognized : " << finalIsClassic.size()
                                                                  << std::endl ;

  // for ( int i = 0 ; i < finalIsClassic.size() ; i++ )
  // {
  //
  //   std::cout << finalBundlesNames[ i ] << " : " << finalIsClassic[ i ]
  //             << " -> Classic : " << tmpAdjacenciesClassic[ i ]
  //             << "\t|\tGeoLab : " << tmpAdjacenciesGeoLab[ i ] << std::endl ;
  //
  // }

  // Free memory
  coveragesGeoLab = std::vector<float>() ;
  adjacenciesGeoLab = std::vector<float>() ;
  overlapsGeoLab = std::vector<float>() ;
  disimilaritiesGeoLab = std::vector<float>() ;
  indecesInNeighborhoodsReconized = std::vector<std::vector<int>>() ;

  coveragesClassic = std::vector<float>() ;
  adjacenciesClassic = std::vector<float>() ;
  overlapsClassic = std::vector<float>() ;
  disimilaritiesClassic = std::vector<float>() ;
  indecesInNeighborhoodsReconizedClassic = std::vector<std::vector<int>>() ;

  // Saving labels
  std::ostringstream labelsDictPathOss ;
  labelsDictPathOss << outputDirectory << "labels.dict" ;
  std::string labelsDictPath = labelsDictPathOss.str() ;
  saveLabelsDict( labelsDictPath.c_str(),
                  finalBundlesNames ) ;


  for ( int i = 0 ; i < nbFibersTractogram ; i++ )
  {

    if ( labelsSubjectSpace[ i ].empty() )
    {

      std::vector<int16_t> _tmpVectorNoLabel = { -1 } ;
      labelsSubjectSpace[ i ] = _tmpVectorNoLabel ;

    }

  }
  std::ostringstream labelsTxtPathOss ;
  labelsTxtPathOss << outputDirectory << "labels.txt" ;
  std::string labelsTxtPath = labelsTxtPathOss.str() ;
  saveLabels( labelsTxtPath.c_str(), labelsSubjectSpace ) ;

  std::cout << "Done" << std::endl ;


  //////////////////////// Saving labels local SBR space ///////////////////////
  std::cout << "#########################################################\n" ;
  std::cout << "########### Saving labels in local SBR sapce ############"
                                                                  << std::endl ;
  std::cout << "#########################################################\n" ;
  std::ostringstream labelsRecognizedSBRPathOss ;
  labelsRecognizedSBRPathOss << outputDirectory << "labelsSBR.txt" ;
  std::string labelsRecognizedSBRPath = labelsRecognizedSBRPathOss.str() ;

  std::ostringstream labelsDictRecognizedSBRPathOss ;
  labelsDictRecognizedSBRPathOss << outputDirectory << "labelsSBR.dict" ;
  std::string labelsDictRecognizedSBRPath =
                                          labelsDictRecognizedSBRPathOss.str() ;


  //
  int nbRecogniedBundles = finalBundlesNames.size() ;
  std::vector<std::vector<int16_t>> labelsSBR ;
  std::vector<float> matrixTracksSBR ;
  std::vector<int32_t> pointsPerTrackSBR ;
  std::vector<int> offsetFirsPointBundleSBR( nbRecogniedBundles ) ;
  std::vector<int> offsetFirstFiberBundleSBR( nbRecogniedBundles ) ;
  int curves_countSBR = 0 ;
  int _bundleLabel = 0 ;
  const auto start_time_tmp = std::chrono::system_clock::now() ;
  for ( int _bundleIndex = 0 ; _bundleIndex < nbRecogniedBundles ;
                                                                _bundleIndex++ )
  {


    std::vector<int> tmpIndexInNeighborhoodRecognized =
                          finalIndecesInNeighborhoodsReconized[ _bundleIndex ] ;

    offsetFirsPointBundleSBR[ _bundleIndex ] =
                                        3 * curves_countSBR * nbPointsPerFiber ;

    offsetFirstFiberBundleSBR[ _bundleIndex ] = curves_countSBR ;

    curves_countSBR += tmpIndexInNeighborhoodRecognized.size() ;

  }

  matrixTracksSBR.resize( 3 * curves_countSBR * nbPointsPerFiber ) ;
  pointsPerTrackSBR.resize( curves_countSBR, nbPointsPerFiber ) ;
  labelsSBR.resize( curves_countSBR ) ;
  #pragma omp parallel for schedule(dynamic)
  for ( int _bundleIndex = 0 ; _bundleIndex < nbRecogniedBundles ;
                                                                _bundleIndex++ )
  {

    std::string tmpBundleName = finalBundlesNames[ _bundleIndex ] ;

    std::string movedTractogramNeighborhood =
                                      finalBundlesNeigborhoods[ _bundleIndex ] ;

    std::string tmpLocalNeighborhood ;
    if ( finalIsClassic[ _bundleIndex ] )
    {

      tmpLocalNeighborhood = movedTractogramNeighborhood ;

    }
    else
    {

      tmpLocalNeighborhood = replaceExtension( movedTractogramNeighborhood,
                                               "_moved.bundles" ) ;
      tmpLocalNeighborhood = replaceExtension( tmpLocalNeighborhood, format ) ;

      if ( !is_file( tmpLocalNeighborhood ) )
      {

        std::cout << "ERROR : finalIsClassic for bundle " << tmpBundleName
                  << " is false but file SBR file " << tmpLocalNeighborhood
                  << " does not exists" << std::endl ;
      }

    }

    BundlesData tmpLocalNeighborhoodBundlesData(
                                                tmpLocalNeighborhood.c_str() ) ;

    std::vector<float>& tmpMatrixTracks =
                                  tmpLocalNeighborhoodBundlesData.matrixTracks ;
    std::vector<int> tmpIndexInNeighborhoodRecognized =
                          finalIndecesInNeighborhoodsReconized[ _bundleIndex ] ;

    for ( int i = 0 ; i < tmpIndexInNeighborhoodRecognized.size() ; i++ )
    {

      int fiber = tmpIndexInNeighborhoodRecognized[ i ] ;

      int64_t offset = 3 * fiber * nbPointsPerFiber ;

      int64_t offset_sbr = offsetFirsPointBundleSBR[ _bundleIndex ] +
                                                      3 * i * nbPointsPerFiber ;

      for ( int point = 0 ; point < nbPointsPerFiber ; point++ )
      {

        for ( int coord = 0 ; coord < 3 ; coord++ )
        {

          matrixTracksSBR[ offset_sbr + 3 * point + coord ] =
                                 tmpMatrixTracks[ 3 * point + coord + offset ] ;

        }

      }

      int64_t offsetLabelsSBR = offsetFirstFiberBundleSBR[ _bundleIndex ] ;

      std::vector<int16_t> _tmpVectorLabel = { (int16_t)_bundleIndex } ;
      labelsSBR[ offsetLabelsSBR + i ] = _tmpVectorLabel ;

    }

  }
  // for ( int _bundleIndex = 0 ; _bundleIndex < nbRecogniedBundles ;
  //                                                               _bundleIndex++ )
  // {
  //
  //   std::string tmpBundleName = finalBundlesNames[ _bundleIndex ] ;
  //
  //   std::string movedTractogramNeighborhood =
  //                                     finalBundlesNeigborhoods[ _bundleIndex ] ;
  //
  //   std::string tmpLocalNeighborhood ;
  //   if ( finalIsClassic[ _bundleIndex ] )
  //   {
  //
  //     tmpLocalNeighborhood = movedTractogramNeighborhood ;
  //
  //   }
  //   else
  //   {
  //
  //     tmpLocalNeighborhood = replaceExtension( movedTractogramNeighborhood,
  //                                              "_moved.bundles" ) ;
  //     tmpLocalNeighborhood = replaceExtension( tmpLocalNeighborhood, format ) ;
  //
  //     if ( !is_file( tmpLocalNeighborhood ) )
  //     {
  //
  //       std::cout << "ERROR : finalIsClassic for bundle " << tmpBundleName
  //                 << " is false but file SBR file " << tmpLocalNeighborhood
  //                 << " does not exists" << std::endl ;
  //     }
  //
  //   }
  //
  //   BundlesData tmpLocalNeighborhoodBundlesData(
  //                                               tmpLocalNeighborhood.c_str() ) ;
  //
  //   std::vector<float>& tmpMatrixTracks =
  //                                 tmpLocalNeighborhoodBundlesData.matrixTracks ;
  //   std::vector<int32_t>& tmpPointsPerTrack =
  //                               tmpLocalNeighborhoodBundlesData.pointsPerTrack ;
  //   int tmpCurves_count = tmpLocalNeighborhoodBundlesData.curves_count ;
  //   std::vector<int> tmpIndexInNeighborhoodRecognized =
  //                         finalIndecesInNeighborhoodsReconized[ _bundleIndex ] ;
  //
  //   for ( int fiber = 0 ; fiber < tmpCurves_count ; fiber++ )
  //   {
  //
  //     int64_t offset = 3 * fiber * nbPointsPerFiber ;
  //
  //     if ( std::find( tmpIndexInNeighborhoodRecognized.begin(),
  //                     tmpIndexInNeighborhoodRecognized.end(), fiber ) ==
  //                                       tmpIndexInNeighborhoodRecognized.end() )
  //     {
  //
  //       continue ;
  //
  //     }
  //
  //     pointsPerTrackSBR.push_back( nbPointsPerFiber ) ;
  //     for ( int point = 0 ; point < nbPointsPerFiber ; point++ )
  //     {
  //
  //       for ( int coord = 0 ; coord < 3 ; coord++ )
  //       {
  //
  //         matrixTracksSBR.push_back( tmpMatrixTracks[ 3 * point + coord +
  //                                                                   offset ] ) ;
  //
  //       }
  //
  //     }
  //
  //     std::vector<int16_t> _tmpVectorLabel = { (int16_t)_bundleIndex } ;
  //     labelsSBR.push_back( _tmpVectorLabel ) ;
  //
  //     curves_countSBR++ ;
  //
  //   }
  //
  // }

  const std::chrono::duration< double > duration_tmp =
                             std::chrono::system_clock::now() - start_time_tmp ;

  std::cout << "Duration : " << duration_tmp.count() << std::endl ;


  //
  BundlesMinf regroupedRecognizedBundles( movedTractogram.c_str() ) ;
  std::ostringstream regroupedRecognizedBundleOss ;
  regroupedRecognizedBundleOss << outputDirectory << "regroupedRecognized" ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    regroupedRecognizedBundleOss << ".bundles" ;

  }
  else if ( format == ".trk" || format == ".tck" )
  {

    regroupedRecognizedBundleOss << ".minf" ;

  }
  regroupedRecognizedBundles.curves_count = curves_countSBR ;


  std::vector<int64_t> fibersWithNans ;
  std::vector<std::vector<float>> tracksScalars ;
  std::vector<std::vector<float>> tracksProperties ;


  bool isBundles = false ;
  bool isTrk = false ;
  bool isTck = false ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    isBundles = true ;

  }
  else if ( format == ".trk" )
  {

    isTrk = true ;

  }
  else if ( format == ".tck" )
  {

    isTck = true ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ERROR : the only format supported are .bundles/.trk/.tck"
                  << std::endl ;
    std::string outMessage = outMessageOss.str() ;
    throw( std::invalid_argument( outMessage ) ) ;

  }



  BundlesData regroupedRecognizedBundlesData( matrixTracksSBR,
                                              pointsPerTrackSBR,
                                              fibersWithNans,
                                              tracksScalars,
                                              tracksProperties,
                                              curves_countSBR,
                                              isBundles,
                                              isTrk,
                                              isTck ) ;
  std::ostringstream regroupedRecognizedBundledataOss ;
  regroupedRecognizedBundledataOss << outputDirectory << "regroupedRecognized" ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    regroupedRecognizedBundledataOss << ".bundlesdata" ;

  }
  else if ( format == ".trk" || format == ".tck" )
  {

    regroupedRecognizedBundledataOss << format ;

  }


  std::string regroupedRecognizedBundledataPath =
                                      regroupedRecognizedBundledataOss.str() ;

  regroupedRecognizedBundlesData.write(
                                    regroupedRecognizedBundledataPath.c_str(),
                                    regroupedRecognizedBundles ) ;

  saveLabelsDict( labelsDictRecognizedSBRPath.c_str(),
                  finalBundlesNames ) ;

  saveLabels( labelsRecognizedSBRPath.c_str(),
              labelsSBR ) ;

  std::cout << "Done" << std::endl ;

  //////////////////////// Saving comparison with atlas ////////////////////////
  std::cout << "#########################################################\n" ;
  std::cout << "############# Saving comparison with atlas ##############"
                                                                  << std::endl ;
  std::cout << "#########################################################\n" ;
  std::ostringstream comparisonWithAtlasRecognizedSBRPathOss ;
  comparisonWithAtlasRecognizedSBRPathOss << outputDirectory
                                                  << "comparisonWithAtlas.tsv" ;
  std::string comparisonWithAtlasRecognizedSBRPath =
                                 comparisonWithAtlasRecognizedSBRPathOss.str() ;
  saveComparisonMeasuresWithAtlas(
                                finalCoverages,
                                finalAdjacencies,
                                finalOverlaps,
                                finalDisimilarities,
                                finalNbFibersPerBundle,
                                finalBundlesNames,
                                comparisonWithAtlasRecognizedSBRPath.c_str() ) ;



  if ( saveBundlesSeparetly )
  {

    std::ostringstream outputDirectoryBundlesOss ;
    outputDirectoryBundlesOss << outputDirectory << "bundles" ;
    std::string outputDirectoryBundles = outputDirectoryBundlesOss.str() ;

    std::ostringstream applyPredictedLabelsCommandOss ;
    applyPredictedLabelsCommandOss << "applyPredictedLabels "
                           << "-i " << regroupedRecognizedBundledataPath << " "
                           << "-l " << labelsRecognizedSBRPath << " "
                           << "-d " << labelsDictRecognizedSBRPath << " "
                           << "-o " << outputDirectoryBundles << " "
                           << "-v 1 " ;
    std::string applyPredictedLabelsCommand =
                                          applyPredictedLabelsCommandOss.str() ;

    int isFailApplyPredicted = run_sh_process( applyPredictedLabelsCommand ) ;

    if ( isFailApplyPredicted )
    {

      std::cout << "ERROR : problem saving separate bundles files"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( keepTmpFiles )
  {

    const std::chrono::duration< double > duration =
                                std::chrono::system_clock::now() - start_time ;

    if ( verbose )
    {

      std::cout << "Duration : " << duration.count() << std::endl ;

    }

    return( 0 ) ;

  }

  std::cout << "Done" << std::endl ;


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

  const std::chrono::duration< double > duration =
                              std::chrono::system_clock::now() - start_time ;

  if ( verbose )
  {

    std::cout << "\n\nDuration : " << duration.count() << std::endl ;

  }

  std::cout << "Done" << std::endl ;

  return( 0 ) ;

}
