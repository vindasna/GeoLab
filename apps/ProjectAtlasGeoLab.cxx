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
//////////////////////// Function to apply Recobundles /////////////////////////
////////////////////////////////////////////////////////////////////////////////
void applyRecoBundles( const std::string& movedTractogramNeighborhood,
                       const std::string& atlasBundleDirectory,
                       const std::string& atlasNeighborhoodFile,
                       const std::string& atlasNeighborhoodCentroidsFile,
                       const std::string& outputDirectory,
                       const std::string& referenceImage,
                       const std::string& format,
                       int nbPointsPerFiber,
                       int portDipyServer,
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


  std::ostringstream recognizedBundleOss ;
  recognizedBundleOss << outputDirectory << bundleName ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    recognizedBundleOss << ".bundles" ;

  }
  else if ( format == ".tck" || format == ".trk" )
  {

    recognizedBundleOss << ".minf" ;

  }
  std::string recognizedBundle = recognizedBundleOss.str() ;

  std::ostringstream recognizedBundledataOss ;
  recognizedBundledataOss << outputDirectory << bundleName ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    recognizedBundledataOss << ".bundlesdata" ;

  }
  else if ( format == ".tck" || format == ".trk" )
  {

    recognizedBundledataOss << format ;

  }
  std::string recognizedBundledata = recognizedBundledataOss.str() ;

  std::ostringstream recognizedBundleClassicOss ;
  recognizedBundleClassicOss << outputDirectory << bundleName << "_classic" ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    recognizedBundleClassicOss << ".bundles" ;

  }
  else if ( format == ".tck" || format == ".trk" )
  {

    recognizedBundleClassicOss << ".minf" ;

  }
  std::string recognizedBundleClassic = recognizedBundleClassicOss.str() ;

  std::ostringstream recognizedBundledataClassicOss ;
  recognizedBundledataClassicOss << outputDirectory << bundleName
                                                                 << "_classic" ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    recognizedBundledataClassicOss << ".bundlesdata" ;

  }
  else if ( format == ".tck" || format == ".trk" )
  {

    recognizedBundledataClassicOss << format ;

  }
  std::string recognizedBundledataClassic =
                                          recognizedBundledataClassicOss.str() ;

  std::ostringstream labelsRecognizedSBRPathOss ;
  labelsRecognizedSBRPathOss << outputDirectory << "labelsSBR_" << bundleName
                                                                     << ".txt" ;
  std::string labelsRecognizedSBRPath = labelsRecognizedSBRPathOss.str() ;

  std::ostringstream labelsDictRecognizedSBRPathOss ;
  labelsDictRecognizedSBRPathOss << outputDirectory << "labelsSBR_"
                                                      << bundleName << ".dict" ;
  std::string labelsDictRecognizedSBRPath =
                                          labelsDictRecognizedSBRPathOss.str() ;


  std::ostringstream clientsLogFilePathOss ;
  clientsLogFilePathOss << outputDirectory << "clientsDipyLog.txt" ;
  std::string clientsLogFilePath = clientsLogFilePathOss.str() ;

  float coverage_classic = -1 ;
  float adjacency_classic = -1 ;
  // if ( false )
  if ( is_file( recognizedBundle ) )
  {


    coverage_classic = getCoverageWithAtlas( recognizedBundle ) ;
    adjacency_classic = getAdjacencyWithAtlas( recognizedBundle ) ;
    int nbFibersClassic = getNbFibers( recognizedBundle ) ;
    // if ( coverage_classic > 0.80 && adjacency_classic > 0.80 )
    if ( adjacency_classic > 0.70 && nbFibersClassic > 10 )
    {

      if ( verbose > 1 )
      {

        std::cout << "Keeping classic projection for bundle " << bundleName
                  << std::endl ;

      }

      return ;

    }
    else
    {

      if ( verbose > 1 )
      {

        std::cout << "Testing RecoBundles projection for bundle "
                  << bundleName << std::endl ;

      }

      bool isRename = rename( recognizedBundle, recognizedBundleClassic ) ;
      if ( !isRename )
      {

        std::cout << "ERROR : could not copy " << recognizedBundle << " to "
                  << recognizedBundleClassic << std::endl ;

        // if ( is_dir( outputDirectory ) )
        // {
        //
        //   rmdir( outputDirectory ) ;
        //
        // }
        closeDipyServer( portDipyServer ) ;
        exit( 1 ) ;

      }

      isRename = rename( recognizedBundledata, recognizedBundledataClassic ) ;
      if ( !isRename )
      {

        std::cout << "ERROR : could not copy " << recognizedBundledata << " to "
                  << recognizedBundledataClassic << std::endl ;

        // if ( is_dir( outputDirectory ) )
        // {
        //
        //   rmdir( outputDirectory ) ;
        //
        // }
        closeDipyServer( portDipyServer ) ;
        exit( 1 ) ;

      }

    }

  }
  else
  {

    if ( verbose > 1 )
    {

      std::cout << "Testing RecoBundles projection for bundle "
                << bundleName << std::endl ;

    }

  }


  std::string atlasNeighborhood ;
  if ( !( endswith( atlasNeighborhoodFile, format ) ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ERROR : in applyRecoBundles, atlasNeighborhoodFile must "
                  << "be in .bundles/.trk/.tck format, got "
                  << atlasNeighborhoodFile << std::endl ;
    std::string outMessage = outMessageOss.str() ;
    throw( std::invalid_argument( outMessage ) ) ;


  }
  else
  {

    atlasNeighborhood = atlasNeighborhoodFile ;

  }

  if ( !is_file( atlasNeighborhood ) )
  {

    std::cout << "ERROR : in applyRecoBundles, atlas neighborhood file "
              << atlasNeighborhood << " does not exists " << std::endl ;

    // if ( is_dir( outputDirectory ) )
    // {
    //
    //   rmdir( outputDirectory ) ;
    //
    // }
    closeDipyServer( portDipyServer ) ;
    exit( 1 ) ;

  }


  // Computing centroids neighborhood atlas
  std::ostringstream atlasBundleFileOss ;
  atlasBundleFileOss << atlasBundleDirectory << bundleName << format ;
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

      std::string atlasNeighborhoodCentroidsMinf = replaceExtension(
                                         atlasNeighborhoodCentroids, ".minf" ) ;

      std::string atlasNeighborhoodMinf = replaceExtension(
                                                  atlasNeighborhood, ".minf" ) ;

      copy( atlasNeighborhoodMinf, atlasNeighborhoodCentroidsMinf ) ;



    }

    std::ostringstream computeCentroidsCommandClientOss ;
    if ( index_cc )
    {

      computeCentroidsCommandClientOss << "python3 "
                                      << computeCentroidsClientFilename << " " ;

    }
    else
    {

      computeCentroidsCommandClientOss << computeCentroidsClientFilename
                                                                        << " " ;

    }
    computeCentroidsCommandClientOss
                                   << "-i " << atlasNeighborhood << " "
                                   << "-o " << atlasNeighborhoodCentroids << " "
                                   << "-r " << referenceImage << " "
                                   << "-thr " << averageRadius << " "
                                   << "-nbPoints " << nbPointsPerFiber << " "
                                   << "-lf " << clientsLogFilePath << " "
                                   << "-p " << portDipyServer << " "
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

      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }
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

  if ( format == ".trk" || format == ".tck" )
  {

    movedTractogramNeighborhoodCentroids = replaceExtension(
                                           movedTractogramNeighborhoodCentroids,
                                           format ) ;

    std::string movedTractogramNeighborhoodCentroidsMinf = replaceExtension(
                               movedTractogramNeighborhoodCentroids, ".minf" ) ;

    std::string movedTractogramNeighborhoodMinf = replaceExtension(
                                        movedTractogramNeighborhood, ".minf" ) ;

    copy( movedTractogramNeighborhoodMinf,
                                    movedTractogramNeighborhoodCentroidsMinf ) ;



  }

  std::ostringstream computeCentroidsCommandClient2Oss ;
  if ( index_cc )
  {

    computeCentroidsCommandClient2Oss << "python3 "
                                      << computeCentroidsClientFilename << " " ;

  }
  else
  {

    computeCentroidsCommandClient2Oss << computeCentroidsClientFilename << " " ;

  }
  computeCentroidsCommandClient2Oss
                         << "-i " << movedTractogramNeighborhood << " "
                         << "-o " << movedTractogramNeighborhoodCentroids << " "
                         << "-r " << referenceImage << " "
                         << "-thr " << averageRadius << " "
                         << "-nbPoints " << nbPointsPerFiber << " "
                         // << "-mnc " << nbClustersAtlasNeighborhood << " "
                         << "-lf " << clientsLogFilePath << " "
                         << "-p " << portDipyServer << " "
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
    if ( _tmpNbFibers > 500 )
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

    // if ( is_dir( outputDirectory ) )
    // {
    //
    //   rmdir( outputDirectory ) ;
    //
    // }
    closeDipyServer( portDipyServer ) ;
    exit( 1 ) ;

  }


  // Registering centroids
  std::string neighborhoodRegistered = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "_moved.bundles" ) ;


  if ( format == ".trk" || format == ".tck" )
  {

    neighborhoodRegistered = replaceExtension( neighborhoodRegistered, format) ;

    std::string neighborhoodRegisteredMinf = replaceExtension(
                                             neighborhoodRegistered, ".minf" ) ;

    std::string movedTractogramNeighborhoodMinf = replaceExtension(
                                        movedTractogramNeighborhood, ".minf" ) ;

    copy( movedTractogramNeighborhoodMinf, neighborhoodRegisteredMinf ) ;



  }

  std::ostringstream registerBundlesClientCommadOss ;
  if ( index_rb )
  {

    registerBundlesClientCommadOss << "python3 "
                                           << registerBundlesClientFile << " " ;

  }
  else
  {

    registerBundlesClientCommadOss << registerBundlesClientFile << " " ;

  }
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

    if ( is_file( recognizedBundleClassic ) )
    {

      bool _tmpIsRename = rename( recognizedBundleClassic, recognizedBundle ) ;
      if ( !_tmpIsRename )
      {

        std::cout << "ERROR : could not move file " << recognizedBundleClassic
                  << " to " << recognizedBundle << std::endl ;
        closeDipyServer( portDipyServer ) ;
        exit( 1 ) ;

      }

    }

    if ( is_file( recognizedBundledataClassic ) )
    {

      bool _tmpIsRename = rename( recognizedBundledataClassic,
                                                      recognizedBundledata ) ;
      if ( !_tmpIsRename )
      {

        std::cout << "ERROR : could not move file "
                  << recognizedBundledataClassic << " to "
                   << recognizedBundledata << std::endl ;
        closeDipyServer( portDipyServer ) ;
        exit( 1 ) ;

      }

    }

    // if ( coverage_classic < coverageThreshold &&
    //                                   adjacency_classic < adjacencyThreshold )
    if ( adjacency_classic < adjacencyThreshold )
    {

      if ( is_file( recognizedBundle ) )
      {

        rmfile( recognizedBundle ) ;

      }

      if ( is_file( recognizedBundledata ) )
      {

        rmfile( recognizedBundledata ) ;

      }

    }

    return ;

  }


  // Projection
  std::ostringstream atlasInfoPathOss ;
  atlasInfoPathOss << atlasBundleDirectory << bundleName ;
  if ( format == ".bundles" || format == ".bundlesdata" )
  {

    atlasInfoPathOss << ".bundles" ;

  }
  else if ( format == ".trk" || format == ".tck" )
  {

    atlasInfoPathOss << ".minf" ;

  }
  std::string atlasInfoPath = atlasInfoPathOss.str() ;
  float thrDistanceBetweenMedialPointsBundle = thrDistanceBetweenMedialPoints ;
  if ( !index_thrDBMP )
  {

    thrDistanceBetweenMedialPointsBundle =
                  getAverageDistanceBetweenMedialPoints( atlasInfoPath ) ;
    thrDistanceBetweenMedialPointsBundle *= ( 1 +
                                      toleranceDistanceBetweenMedialPoints ) ;

  }

  std::ostringstream projectAtlasCommandOss ;
  projectAtlasCommandOss << projectAtlasFile << " "
                         << "-i " << neighborhoodRegistered << " "
                         << "-a " << atlasBundleDirectory << " "
                         << "-o " << outputDirectory << " "
                         << "-l " << "labelsSBR_" << bundleName << " "
                         << "-minNbFibers " << minimumNumberFibers << " "
                         << "-thrSim " << thrPercentageSimilarity << " "
                         << "-thrDBMP " << thrDistanceBetweenMedialPointsBundle
                                                                          << " "
                         << "-tolP " << toleranceP << " "
                         << "-tolThr " << toleranceThr << " "
                         << "-tolMaxAng " << toleranceMaxAngle << " "
                         << "-tolMaxDirAng " << toleranceMaxDirectionAngle
			                                                 << " "
                         << "-tolMinShapeAng " << toleranceMinShapeAngle << " "
                         << "-tolMaxShapeAng " << toleranceMaxShapeAngle << " "
                         << "-tolLenght " << toleranceLenght << " "
                         << "-thrAdj " << adjacencyForCompareBundles << " "
                         << "-cb "
                         << "-nbThreads " << nbThreads << " "
                         << "-v 1" ;

  std::string projectAtlasCommand = projectAtlasCommandOss.str() ;
  // Here we do not use force because we want to always do the RecoBundles
  // projection

  int isFailProjection = run_sh_process( projectAtlasCommand ) ;

  if ( isFailProjection )
  {

    std::cout << "ERROR : could not project " << neighborhoodRegistered
              << " to " << atlasBundleDirectory << ", got exit code "
              << isFailProjection << " with command " << projectAtlasCommand
                                                                  << std::endl ;

    // if ( is_dir( outputDirectory ) )
    // {
    //
    //   rmdir( outputDirectory ) ;
    //
    // }
    closeDipyServer( portDipyServer ) ;
    exit( 1 ) ;

  }

  // Getting labels
  if ( is_file( recognizedBundleClassic ) &&
                                        is_file( recognizedBundledataClassic ) )
  {

    if ( coverage_classic == -1 )
    {

      std::cout << "ERROR : invalid coverage_classic : -1 " << std::endl ;

      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }
      closeDipyServer( portDipyServer ) ;
      exit( 1 ) ;

    }

    if ( !is_file( recognizedBundle ) )
    {

      if ( is_file( labelsRecognizedSBRPath ) )
      {

        rmfile( labelsRecognizedSBRPath ) ;

      }

      if ( is_file( labelsDictRecognizedSBRPath ) )
      {

        rmfile( labelsDictRecognizedSBRPath ) ;

      }

      bool _tmpIsRename1 = rename( recognizedBundleClassic, recognizedBundle ) ;
      if ( !_tmpIsRename1 )
      {

        std::cout << "ERROR : could not move file " << recognizedBundleClassic
                  << " to " << recognizedBundle << std::endl ;
        closeDipyServer( portDipyServer ) ;
        exit( 1 ) ;

      }

      bool _tmpIsRename2 = rename( recognizedBundledataClassic,
                                                        recognizedBundledata ) ;
      if ( !_tmpIsRename2 )
      {

        std::cout << "ERROR : could not move file "
                  << recognizedBundledataClassic << " to "
                  << recognizedBundledata << std::endl ;
        closeDipyServer( portDipyServer ) ;
        exit( 1 ) ;

      }

      // if ( coverage_classic < coverageThreshold &&
      //                                   adjacency_classic < adjacencyThreshold )
      if ( adjacency_classic < adjacencyThreshold )
      {

        if ( is_file( recognizedBundle ) )
        {

          rmfile( recognizedBundle ) ;

        }

        if ( is_file( recognizedBundledata ) )
        {

          rmfile( recognizedBundledata ) ;

        }

      }

      return ;

    }

    float coverage_recobundles = getCoverageWithAtlas( recognizedBundle ) ;
    float adjacency_recobundles = getAdjacencyWithAtlas( recognizedBundle ) ;
    // if ( coverage_classic > coverage_recobundles )
    if ( adjacency_classic > adjacency_recobundles )
    {

      if ( is_file( labelsRecognizedSBRPath ) )
      {

        rmfile( labelsRecognizedSBRPath ) ;

      }

      if ( is_file( labelsDictRecognizedSBRPath ) )
      {

        rmfile( labelsDictRecognizedSBRPath ) ;

      }

      if ( is_file( recognizedBundle ) )
      {

        rmfile( recognizedBundle ) ;

      }

      if ( is_file( recognizedBundledata ) )
      {

        rmfile( recognizedBundledata ) ;

      }

      if ( is_file( recognizedBundleClassic ) )
      {


        bool _tmpIsRename = rename( recognizedBundleClassic,
                                                            recognizedBundle ) ;
        if ( !_tmpIsRename )
        {

          std::cout << "ERROR : could not move file " << recognizedBundleClassic
                    << " to " << recognizedBundle << std::endl ;
          closeDipyServer( portDipyServer ) ;
          exit( 1 ) ;

        }

      }

      if ( is_file( recognizedBundledataClassic ) )
      {

        bool _tmpIsRename = rename( recognizedBundledataClassic,
                                                        recognizedBundledata ) ;
        if ( !_tmpIsRename )
        {

          std::cout << "ERROR : could not move file "
                    << recognizedBundledataClassic << " to "
                     << recognizedBundledata << std::endl ;
          closeDipyServer( portDipyServer ) ;
          exit( 1 ) ;

        }

      }

      // if ( coverage_classic < coverageThreshold &&
      //                                   adjacency_classic < adjacencyThreshold )
      if ( adjacency_classic < adjacencyThreshold )
      {

        if ( is_file( recognizedBundle ) )
        {

          rmfile( recognizedBundle ) ;

        }

        if ( is_file( recognizedBundledata ) )
        {

          rmfile( recognizedBundledata ) ;

        }

      }


    }
    else
    {

      if ( !is_file( recognizedBundle ) )
      {

        if ( is_file( labelsRecognizedSBRPath ) )
        {

          rmfile( labelsRecognizedSBRPath ) ;

        }

        if ( is_file( labelsDictRecognizedSBRPath ) )
        {

          rmfile( labelsDictRecognizedSBRPath ) ;

        }

        if ( is_file( recognizedBundleClassic ) )
        {


          bool _tmpIsRename = rename( recognizedBundleClassic,
                                                            recognizedBundle ) ;
          if ( !_tmpIsRename )
          {

            std::cout << "ERROR : could not move " << recognizedBundleClassic
                      << " to " << recognizedBundle << std::endl ;
            closeDipyServer( portDipyServer ) ;
            exit( 1 ) ;

          }

        }

        if ( is_file( recognizedBundledataClassic ) )
        {

          if ( is_file( recognizedBundledata ) )
          {

            rmfile( recognizedBundledata ) ;

          }

          bool _tmpIsRename = rename( recognizedBundledataClassic,
                                                        recognizedBundledata ) ;
          if ( !_tmpIsRename )
          {

            std::cout << "ERROR : could not move file "
                      << recognizedBundledataClassic << " to "
                       << recognizedBundledata << std::endl ;
            closeDipyServer( portDipyServer ) ;
            exit( 1 ) ;

          }

        }

        // if ( coverage_classic < coverageThreshold &&
        //                                 adjacency_classic < adjacencyThreshold )
        if ( adjacency_classic < adjacencyThreshold )
        {

          if ( is_file( recognizedBundle ) )
          {

            rmfile( recognizedBundle ) ;

          }

          if ( is_file( recognizedBundledata ) )
          {

            rmfile( recognizedBundledata ) ;

          }

        }

        return ;

      }

      if ( is_file( recognizedBundleClassic ) )
      {

        rmfile( recognizedBundleClassic ) ;

      }

      if ( is_file( recognizedBundledataClassic ) )
      {

        rmfile( recognizedBundledataClassic ) ;

      }


      // if ( coverage_recobundles < coverageThreshold &&
      //                               adjacency_recobundles < adjacencyThreshold )
      if ( adjacency_recobundles < adjacencyThreshold )
      {

        if ( is_file( recognizedBundle ) )
        {

          rmfile( recognizedBundle ) ;

        }

        if ( is_file( recognizedBundledata ) )
        {

          rmfile( recognizedBundledata ) ;

        }

        if ( is_file( labelsRecognizedSBRPath ) )
        {

          rmfile( labelsRecognizedSBRPath ) ;

        }

        if ( is_file( labelsDictRecognizedSBRPath ) )
        {

          rmfile( labelsDictRecognizedSBRPath ) ;

        }

        return ;

      }

    }

  }
  else
  {

    if ( !is_file( recognizedBundle ) )
    {

      if ( is_file( recognizedBundledata ) )
      {

        rmfile( recognizedBundledata ) ;

      }

      if ( is_file( labelsRecognizedSBRPath ) )
      {

        rmfile( labelsRecognizedSBRPath ) ;

      }

      if ( is_file( labelsDictRecognizedSBRPath ) )
      {

        rmfile( labelsDictRecognizedSBRPath ) ;

      }

      return ;

    }


    float coverage_recobundles = getCoverageWithAtlas( recognizedBundle ) ;
    float adjacency_recobundles = getAdjacencyWithAtlas( recognizedBundle ) ;
    // if ( coverage_recobundles < coverageThreshold &&
    //                                 adjacency_recobundles < adjacencyThreshold )
    if ( adjacency_recobundles < adjacencyThreshold )
    {

      if ( is_file( recognizedBundle ) )
      {

        rmfile( recognizedBundle ) ;

      }

      if ( is_file( recognizedBundledata ) )
      {

        rmfile( recognizedBundledata ) ;

      }

      if ( is_file( labelsRecognizedSBRPath ) )
      {

        rmfile( labelsRecognizedSBRPath ) ;

      }

      if ( is_file( labelsDictRecognizedSBRPath ) )
      {

        rmfile( labelsDictRecognizedSBRPath ) ;

      }

      return ;

    }

  }

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
  index_sp = getFlagPosition( argc, argv, "-sp" ) ;
  index_force = getFlagPosition( argc, argv, "-force" ) ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_keep_tmp = getFlagPosition( argc, argv, "-ktf" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-i : Path to the input tractogram \n"
              << "-a : Path to the directory with the atlas \n"
              << "-ref : Path to the reference image where the tractogram is \n"
              << "-o : Path to the output directory \n"
              << "-nbPoints : Number of points per fiber (same number for all "
              << "fibers) \n"
              << "[-cc] : Path to the clientComputeCentroids.py file \n"
              << "[-rb] : Path to the clientRegisterBundles.py file \n"
              << "[-ods] : Path to the dipyServer.py file \n"
              << "[-cds] : Path to the clientCloseServer.py file \n"
              << "[-fa] : Path to full atlas tractogram (mandatory for global "
              << "SLR or if -anc is to given) \n"
              << "[-an] : Path to the atlas neighborhoods (must be given is -fa"
              << " is not given)\n"
              << "[-anc] : Path to the atlas neighborhoods centroids (must be "
              << "given is -fa is not given)\n"
              << "[-thrCov] : Threshold to keep bundles where coverage is "
              << "greater than thrCov (default : 0 -> keep all bundles ) \n"
              << "[-thrDBMP] : Threshold for maximum distance between medial "
              << "points (default = 50 mm) \n"
              << "[-tolP] : Tolerance for parameter p (for advanced users, "
              << "default = 0) \n"
              << "[-tolThr] : Tolerance for parameter thr (for advanced users, "
              << "default = 0) \n"
              << "[-tolMaxAng] : Tolerance for parameter max angle (for "
              << "advanced users,  default = 0) \n"
              << "[-tolMaxDirAng] : Tolerance for parameter max direction angle"
              << " (for advanced users,  default = 0) \n"
              << "[-tolMinShapeAng] : Tolerance for parameter min shape angle"
              << " (for advanced users,  default = 0) \n"
              << "[-tolMaxShapeAng] : Tolerance for parameter max shape angle"
              << " (for advanced users,  default = 0) \n"
              << "[-tolLenght] : Tolerance for parameter lenght (for advanced "
              << "users, default = 0) \n"
              << "[-tolThrCN] : tolerance for computeNeighborhood threshold "
              << "(default = 2.0) \n"
              << "[-tolDBMP] : Tolerance for distance between medial points "
              << "(for advanced users, default = 1.0) \n"
              << "[thrAdj] : keep bundle with adjacency greater than given value"
              << " (default : 0 -> keep all bundles ) \n"
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
              << "[-sp] : Save recognized bundles separetly (default : true)\n"
              << "[-force] : Force to overwrite files (default = false) \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads "
              << "(default : number of cores ) \n"
              << "[-ktf] : Keep temp files (default = false) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_atlas )
  {

    std::cout << "-a argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_fa && !index_anc )
  {

    std::cout << "At least -fa or -anc must be given ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_fa && !index_an )
  {

    std::cout << "At least -fa or -an must be given ..." << std::endl ;
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


  //////////////////////////////////////////////////////////////////////////////
  int nbCores = omp_get_num_procs() ;
  std::cout << "Number of cores in the system : " << nbCores << std::endl ;

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

    std::string outMessage = "ProjectAtlasRecoBundles : Only supported " \
                             "formats for input are .bundles/.bundlesdata, " \
                             ".trk, .tck \n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }


  //////////////////////////////////////////////////////////////////////////////
  char lastChar = inputBundlesFilename[ inputBundlesFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputBundlesFilename = inputBundlesFilename.substr( 0,
                                             inputBundlesFilename.size() - 1 ) ;

  }


  if ( !is_file( inputBundlesFilename ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ProjectAtlasRecoBundles : input file "
                  << inputBundlesFilename << "does not exists" << std::endl ;
    std::string outMessage = outMessageOss.str() ;
    throw( std::invalid_argument( outMessage ) ) ;


  }


  ////////////////////////////// Atlas directory ///////////////////////////////
  std::string atlasDirectory( argv[ index_atlas + 1 ] ) ;
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
  /////////////////////////////// Reference image //////////////////////////////
  std::string referenceFilename( argv[ index_reference + 1 ] ) ;
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

  ///////////////////////////////// Full atlas /////////////////////////////////
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

      std::string outMessage = "ConvertBundleFormat: Only supported formats " \
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

  ///////////////////////////// Atlas neighborhood /////////////////////////////
  std::string atlasNeighborhoodDirectory ;
  if ( index_an )
  {

    isAtlasNeighborhood = true ;

    atlasNeighborhoodDirectory = argv[ index_anc + 1 ] ;
    lastChar = atlasNeighborhoodDirectory[ atlasNeighborhoodDirectory.size()
                                                                         - 1 ] ;
    if ( lastChar != '/' )
    {

      atlasNeighborhoodDirectory = atlasNeighborhoodDirectory + "/" ;

    }
    if ( !is_dir( atlasNeighborhoodDirectory ) )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : Atlas neighborhood centroids directory "
                    << atlasNeighborhoodDirectory << " does not exists "
                    << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

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

    if ( toleranceP < -1 || toleranceP > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolP must be in [-1;1]" << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;


    }

  }

  if ( index_tolThr )
  {

    toleranceThr = std::stof( argv[ index_tolThr + 1 ] ) ;

    if ( toleranceThr < -1 || toleranceThr > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolThr must be in [-1;1]"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }

  if ( index_tolMaxAngle )
  {

    toleranceMaxAngle = std::stof( argv[ index_tolMaxAngle + 1 ] ) ;

    if ( toleranceMaxAngle < -1 || toleranceMaxAngle > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolMaxAng must be in [-1;1]"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }

  if ( index_tolMaxDirectionAngle )
  {

    toleranceMaxDirectionAngle = std::stof( argv[
                                            index_tolMaxDirectionAngle + 1 ] ) ;

    if ( toleranceMaxDirectionAngle < -1 || toleranceMaxDirectionAngle > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolMaxDirAng must be in [-1;1]"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }

  if ( index_tolMinShapeAngle )
  {

    toleranceMinShapeAngle = std::stof( argv[ index_tolMinShapeAngle + 1 ] ) ;

    if ( toleranceMinShapeAngle < -1 || toleranceMinShapeAngle > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolMinShapeAng must be in [-1;1]"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;



    }

  }

  if ( index_tolMaxShapeAngle )
  {

    toleranceMaxShapeAngle = std::stof( argv[ index_tolMaxShapeAngle + 1 ] ) ;

    if ( toleranceMaxShapeAngle < -1 || toleranceMaxShapeAngle > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolManShapeAng must be in [-1;1]"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }

  if ( index_tolLenght )
  {

    toleranceLenght = std::stof( argv[ index_tolLenght + 1 ] ) ;

    if ( toleranceLenght < -1 || toleranceLenght > 1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolLenght must be in [-1;1]"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }


  if ( index_tolThrCN )
  {

    toleranceThrComputeNeighborhood = std::stof( argv[ index_tolThrCN + 1 ] ) ;
    if ( toleranceThrComputeNeighborhood <= 0 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ERROR : argument -tolThrCN must be greater than 0"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    toleranceThrComputeNeighborhood = 2.0 ;

  }

  if ( index_tolDistBetMedPts )
  {

    toleranceDistanceBetweenMedialPoints = std::stof(
                                          argv[ index_tolDistBetMedPts + 1 ] ) ;

    if ( toleranceDistanceBetweenMedialPoints < -1 )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "Error argument : -tolDBMP must be greater than -1"
                                                                  << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }

  }
  else
  {

    toleranceDistanceBetweenMedialPoints = 1.0 ;

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

  /////////////////////////// Minimum number of fibers /////////////////////////
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

  /////////////////////////// Minimum number of fibers /////////////////////////
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

  ////////////////////////// computeNeighborhood file //////////////////////////
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

    // nbThreads = nbCores ;
    nbThreads = -1 ;

  }

  omp_set_num_threads( nbThreads ) ;

  #pragma omp parallel
  {

    nbThreadsUsed = omp_get_num_threads() ;

  }
  std::cout << "Number of threads : " << nbThreadsUsed << std::endl ;



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

  }


  /////////////////////////////// Getting format ///////////////////////////////
  std::string format ;
  if ( endswith( inputBundlesFilename, ".bundles" ) ||
                              endswith( inputBundlesFilename, ".bundlesdata" ) )
  {

    format = ".bundles" ;

  }
  else if ( endswith( inputBundlesFilename, ".trk" ) )
  {

    format = ".trk" ;

  }
  else if ( endswith( inputBundlesFilename, ".tck" ) )
  {

    format = ".tck" ;

  }
  // Already checked at the begginig of main if extension is one of those format


  //////////////// Preparing atlas bundles paths (tmpAtlasDir) /////////////////
  if ( verbose )
  {

    std::cout << "#########################################################\n" ;
    std::cout << "############# Preparing atlas bundles paths #############"
                                                                  << std::endl ;
    std::cout << "#########################################################\n" ;

  }

  checkAtlasDirectory( atlasDirectory, format ) ;
  std::vector<std::string> atlasBundleDirectories ;
  std::string tmpAtlasDir = getAtlasBunldesPaths( outputDirectory,
                                                  atlasDirectory,
                                                  format,
                                                  atlasBundleDirectories ) ;

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

  ///////////////// Getting atlas centroids if input is given //////////////////
  std::vector<std::string> atlasNeighborhoodCentroidsPaths ;
  if ( isAtlasNeighborhoodCentroids )
  {

    getAtlasNeighborhoodCentroids( atlasNeighborhoodCentroidsDirectory,
                                   atlasBundleDirectories,
                                   format,
                                   atlasNeighborhoodCentroidsPaths ) ;

  }


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
      outMessageOss << "ProjectAtlasRecoBundles : global SLR not compatible "
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
      std::cout << "Dipy SLR command :\n " << globalSLR << std::endl ;
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
                                << "-nbThreads " << nbCores << " "
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
                        << "-tolThr " << toleranceThrComputeNeighborhood << " "
                        << "-nbThreads " << nbCores << " "
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
  if ( verbose )
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

  if ( nbBundles == 0 )
  {

    std::cout << "ERROR : no valid bundles in atlas directory" << std::endl ;
    exit( 1 ) ;

  }

  #pragma omp parallel for schedule(dynamic)
  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    std::string atlasBundleDirectory = atlasBundleDirectories[ i ] ;
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
    atlasInfoPathOss << atlasBundleDirectory << bundleName ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      atlasInfoPathOss << ".bundles" ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      atlasInfoPathOss << ".minf" ;

    }
    std::string atlasInfoPath = atlasInfoPathOss.str() ;
    float thrDistanceBetweenMedialPointsBundle =
                                                thrDistanceBetweenMedialPoints ;

    if ( !index_thrDBMP )
    {

      thrDistanceBetweenMedialPointsBundle =
                  getAverageDistanceBetweenMedialPoints( atlasInfoPath ) / 2.0 ;
      thrDistanceBetweenMedialPointsBundle *= ( 1 +
                                        toleranceDistanceBetweenMedialPoints ) ;

    }


    std::ostringstream projectCommandOss ;
    projectCommandOss << projectAtlasFile << " "
                      << "-i " << movedTractogramNeighborhood << " "
                      << "-a " << atlasBundleDirectory << " "
                      << "-o " << outputDirectory << " "
                      << "-l " << bundleName << "_labels" << " "
                      << "-minNbFibers " << minimumNumberFibers << " "
                      << "-thrSim " << thrPercentageSimilarity << " "
                      << "-thrDBMP " << thrDistanceBetweenMedialPointsBundle
                                                                          << " "
		                  << "-tolP " << toleranceP << " "
                      << "-tolThr " << toleranceThr << " "
                      << "-tolMaxAng " << toleranceMaxAngle << " "
                      << "-tolMaxDirAng " << toleranceMaxDirectionAngle
                                                                      << " "
                      << "-tolMinShapeAng " << toleranceMinShapeAngle << " "
                      << "-tolMaxShapeAng " << toleranceMaxShapeAngle << " "
                      << "-tolLenght " << toleranceLenght << " "
                      << "-thrAdj " << adjacencyForCompareBundles << " "
                      << "-cb "
                      << "-nbThreads " << nbThreads << " "
                      << "-v 1" ;
    std::string projectCommand = projectCommandOss.str() ;

    std::ostringstream extractedBundlePathOss ;
    extractedBundlePathOss << outputDirectory << bundleName ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      extractedBundlePathOss << ".bundles" ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      extractedBundlePathOss << ".minf" ;

    }
    std::string extractedBundlePath = extractedBundlePathOss.str() ;

    std::ostringstream extractedBundleDataPathOss ;
    extractedBundleDataPathOss << outputDirectory << bundleName ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      extractedBundleDataPathOss << ".bundlesdata" ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      extractedBundleDataPathOss << format ;

    }
    std::string extractedBundleDataPath = extractedBundleDataPathOss.str() ;

    std::ostringstream extractedBundleLabelsPathOss ;
    extractedBundleLabelsPathOss << outputDirectory << bundleName << "_labels"
                                                                     << ".txt" ;
    std::string extractedBundleLabelsPath = extractedBundleLabelsPathOss.str() ;

    std::ostringstream extractedBundleLabelsDictPathOss ;
    extractedBundleLabelsDictPathOss << outputDirectory << bundleName
                                                        << "_labels" << ".txt" ;
    std::string extractedBundleLabelsDictPath =
                                        extractedBundleLabelsDictPathOss.str() ;

    int isProjectAtlasFail = 0 ;
    if ( is_file( extractedBundlePath ) && is_file( extractedBundleDataPath )
                              && is_file( extractedBundleLabelsPath )
                                   && is_file( extractedBundleLabelsDictPath ) )
    {

      if ( verbose > 1 )
      {

        std::cout << "WARNING : output files of extracted bundle " << bundleName
                  << " seem to exist and -force flag was not used, trying "
                  << "computations with existing files" << std::endl ;

      }

    }
    else
    {

      if ( doClassical )
      {

        isProjectAtlasFail = run_sh_process( projectCommand ) ;

      }
      else
      {

        isProjectAtlasFail = 0 ;

      }

    }
    if ( isProjectAtlasFail )
    {

      std::cout << "ERROR : could not project atlas without SBR, for bundle "
                << bundleName << "got exit code " << isProjectAtlasFail
                << " for command : \n" << projectCommand << std::endl ;

      // if ( is_dir( outputDirectory ) )
      // {
      //
      //   rmdir( outputDirectory ) ;
      //
      // }
      exit( 1 ) ;

    }


  }

  std::vector<std::string> labelsDictClassic ;
  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    std::string atlasBundleDirectory = atlasBundleDirectories[ i ] ;
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

    labelsDictClassic.push_back( bundleName ) ;

  }

  if ( labelsDictClassic.size() != nbBundles )
  {

    std::cout << "ERROR : the number of bundles in the atlas dir and the "
              << "classic labels dict is not the same" << std::endl ;
    exit( 1 ) ;

  }

  int nbFibersTractogram = getNbFibers( movedTractogram ) ;
  std::vector<std::vector<int16_t>> labelsClassic ;
  labelsClassic.resize( nbFibersTractogram ) ;
  for ( int i = 0 ; i < nbBundles ; i++ )
  {


    std::string movedTractogramNeighborhood = neighborhoodFilenames[ i ] ;

    std::string atlasNeighborhoodFile ;
    if ( isAtlasNeighborhoodCentroids )
    {

      atlasNeighborhoodFile = atlasNeighborhoodCentroidsPaths[ i ] ;


    }
    else
    {

      atlasNeighborhoodFile = neighborhoodAtlasFilenames[ i ] ;

    }

    std::string _tmpAtlasBundlePath = atlasNeighborhoodFile ;

    char lastChar = _tmpAtlasBundlePath[ _tmpAtlasBundlePath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      _tmpAtlasBundlePath = _tmpAtlasBundlePath.substr( 0,
                                              _tmpAtlasBundlePath.size() - 1 ) ;

    }

    std::experimental::filesystem::path tmpPath( _tmpAtlasBundlePath ) ;
    std::string bundleName = tmpPath.stem() ;


    std::ostringstream labelsRecognizedPathOss ;
    labelsRecognizedPathOss << outputDirectory << bundleName << "_labels"
                                                                     << ".txt" ;
    std::string labelsRecognizedPath = labelsRecognizedPathOss.str() ;

    std::ostringstream labelsDictRecognizedPathOss ;
    labelsDictRecognizedPathOss << outputDirectory << bundleName << "_labels"
                                                                    << ".dict" ;
    std::string labelsDictRecognizedPath = labelsDictRecognizedPathOss.str() ;


    if ( is_file( labelsRecognizedPath ) &&
                                           is_file( labelsDictRecognizedPath ) )
    {


      std::string neighborhoodTractogramBinPath = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "Index.bin" ) ;

      int nbFibersMovedTracotgramNeighborhood = getNbFibers(
                                                 movedTractogramNeighborhood ) ;

      std::vector<std::vector<int16_t>> labelsRecognized ;
      readingPredictedLabels( labelsRecognizedPath.c_str(),
                              labelsRecognized,
                              nbFibersMovedTracotgramNeighborhood ) ;

      std::vector<int64_t> indexInTractogram ;
      readIndexInTractogram( neighborhoodTractogramBinPath.c_str(),
                             indexInTractogram,
                             nbFibersMovedTracotgramNeighborhood ) ;

      int _tmpCounter = 0 ;
      for ( int _label = 0 ; _label < nbBundles ; _label++ )
      {

        if ( labelsDictClassic[ _label ] == bundleName )
        {

          for ( int _fiberIndex = 0 ;
                            _fiberIndex < nbFibersMovedTracotgramNeighborhood ;
                                                                 _fiberIndex++ )
          {

            int nbLabelsForFiberClassic =
                                        labelsRecognized[ _fiberIndex ].size() ;
            for ( int _labelIndex = 0 ; _labelIndex < nbLabelsForFiberClassic ;
                                                                 _labelIndex++ )
            {

              int labelRecognizedClassicTmp =
                                labelsRecognized[ _fiberIndex ][ _labelIndex ] ;

              if ( labelRecognizedClassicTmp == 0 )
              {

                int64_t _indexInTractogram = indexInTractogram[ _fiberIndex ] ;

                labelsClassic[ _indexInTractogram ].push_back( _label ) ;
                _tmpCounter++ ;

              }

            }

          }

          break ;

        }

      }

      if ( is_file( labelsRecognizedPath ) )
      {

        rmfile( labelsRecognizedPath ) ;

      }

      if ( is_file( labelsDictRecognizedPath ) )
      {

        rmfile( labelsDictRecognizedPath ) ;

      }

    }
    else
    {

      if ( is_file( labelsRecognizedPath ) )
      {

        rmfile( labelsRecognizedPath ) ;

      }

      if ( is_file( labelsDictRecognizedPath ) )
      {

        rmfile( labelsDictRecognizedPath ) ;

      }

    }

  }


  for ( int _fiberIndex = 0 ; _fiberIndex < nbFibersTractogram ; _fiberIndex++ )
  {

    if ( labelsClassic[ _fiberIndex ].empty() )
    {

      labelsClassic[ _fiberIndex ].push_back( -1 ) ;

    }

  }


  const std::chrono::duration< double > duration_no_sbr =
                          std::chrono::system_clock::now() - start_time_no_sbr ;

  if ( verbose )
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
  if ( index_ods )
  {

    launchDipyServiceOss << "python3 " << openDipyServerClientFile << " " ;

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

  std::vector<std::vector<int16_t>> labels ;
  labels.resize( nbFibersTractogram ) ;
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

      applyRecoBundles( movedTractogramNeighborhood,
                        atlasBundleDirectory,
                        atlasNeighborhoodFile,
                        atlasNeighborhoodCentroidsFile,
                        outputDirectory,
                        referenceFilename,
                        format,
                        nbPointsPerFiber,
                        portDipyServer,
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
    /*
    try
    {

      c.terminate() ;

    }
    catch(...){}
    */
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
  std::cout << "############# Saving labels in subject space ############"
                                                                  << std::endl ;
  std::cout << "#########################################################\n" ;
  for ( int i = 0 ; i < nbBundles ; i++ )
  {


    std::string movedTractogramNeighborhood = neighborhoodFilenames[ i ] ;

    std::string atlasNeighborhoodFile ;
    if ( isAtlasNeighborhoodCentroids )
    {

      atlasNeighborhoodFile = atlasNeighborhoodCentroidsPaths[ i ] ;


    }
    else
    {

      atlasNeighborhoodFile = neighborhoodAtlasFilenames[ i ] ;

    }

    std::string _tmpAtlasBundlePath = atlasNeighborhoodFile ;

    char lastChar = _tmpAtlasBundlePath[ _tmpAtlasBundlePath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      _tmpAtlasBundlePath = _tmpAtlasBundlePath.substr( 0,
                                              _tmpAtlasBundlePath.size() - 1 ) ;

    }

    std::experimental::filesystem::path tmpPath( _tmpAtlasBundlePath ) ;
    std::string bundleName = tmpPath.stem() ;


    std::ostringstream labelsRecognizedSBRPathOss ;
    labelsRecognizedSBRPathOss << outputDirectory << "labelsSBR_" << bundleName
                                                                     << ".txt" ;
    std::string labelsRecognizedSBRPath = labelsRecognizedSBRPathOss.str() ;

    std::ostringstream labelsDictRecognizedSBRPathOss ;
    labelsDictRecognizedSBRPathOss << outputDirectory << "labelsSBR_"
                                                      << bundleName << ".dict" ;
    std::string labelsDictRecognizedSBRPath =
                                          labelsDictRecognizedSBRPathOss.str() ;


    if ( is_file( labelsRecognizedSBRPath ) &&
                                        is_file( labelsDictRecognizedSBRPath ) )
    {


      std::string neighborhoodTractogramBinPath = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "Index.bin" ) ;

      int nbFibersMovedTracotgramNeighborhood = getNbFibers(
                                                 movedTractogramNeighborhood ) ;

      std::vector<std::vector<int16_t>> labelsRecognizedSBR ;
      readingPredictedLabels( labelsRecognizedSBRPath.c_str(),
                              labelsRecognizedSBR,
                              nbFibersMovedTracotgramNeighborhood ) ;

      std::vector<int64_t> indexInTractogram ;
      readIndexInTractogram( neighborhoodTractogramBinPath.c_str(),
                             indexInTractogram,
                             nbFibersMovedTracotgramNeighborhood ) ;

      int _tmpCounter = 0 ;
      for ( int _label = 0 ; _label < nbBundles ; _label++ )
      {

        if ( labelsDictClassic[ _label ] == bundleName )
        {

          for ( int _fiberIndex = 0 ;
                            _fiberIndex < nbFibersMovedTracotgramNeighborhood ;
                                                                 _fiberIndex++ )
          {

            int nbLabelsForFiberSBR = labelsRecognizedSBR[ _fiberIndex ].size() ;
            for ( int _labelIndex = 0 ; _labelIndex < nbLabelsForFiberSBR ;
                                                                 _labelIndex++ )
            {

              int labelRecognizedSBRtmp =
                             labelsRecognizedSBR[ _fiberIndex ][ _labelIndex ] ;

              if ( labelRecognizedSBRtmp == 0 )
              {

                int64_t _indexInTractogram = indexInTractogram[ _fiberIndex ] ;

                labels[ _indexInTractogram ].push_back( _label ) ;
                _tmpCounter++ ;

              }

            }

          }

          break ;

        }

      }

      if ( is_file( labelsRecognizedSBRPath ) )
      {

        rmfile( labelsRecognizedSBRPath ) ;

      }

      if ( is_file( labelsDictRecognizedSBRPath ) )
      {

        rmfile( labelsDictRecognizedSBRPath ) ;

      }

    }
    else
    {

      if ( is_file( labelsRecognizedSBRPath ) )
      {

        rmfile( labelsRecognizedSBRPath ) ;

      }

      if ( is_file( labelsDictRecognizedSBRPath ) )
      {

        rmfile( labelsDictRecognizedSBRPath ) ;

      }

      int _tmpCounter = 0 ;

      for ( int _label = 0 ; _label < nbBundles ; _label++ )
      {

        if ( labelsDictClassic[ _label ] == bundleName )
        {

          int nbLabeledFibersClassic = labelsClassic.size() ;

          for ( int _fiberIndex = 0; _fiberIndex < nbLabeledFibersClassic ;
                                                                 _fiberIndex++ )
          {

            int nbLabelsForFiberClassic = labelsClassic[ _fiberIndex ].size() ;
            for ( int _labelClassicIndex = 0 ;
                  _labelClassicIndex < nbLabelsForFiberClassic ;
                                                          _labelClassicIndex++ )
            {

              int _classicLabelTmp =
                            labelsClassic[ _fiberIndex ][ _labelClassicIndex ] ;

              if ( _classicLabelTmp == _label )
              {

                labels[ _fiberIndex ].push_back( _label ) ;
                _tmpCounter++ ;

              }

            }

          }

          break ;

        }

      }

    }

  }

  // Saving labels
  std::ostringstream labelsDictPathOss ;
  labelsDictPathOss << outputDirectory << "labels.dict" ;
  std::string labelsDictPath = labelsDictPathOss.str() ;
  saveLabelsDict( labelsDictPath.c_str(),
                  labelsDictClassic ) ;


  for ( int i = 0 ; i < nbFibersTractogram ; i++ )
  {

    if ( labels[ i ].empty() )
    {

      std::vector<int16_t> _tmpVectorNoLabel = { -1 } ;
      labels[ i ] = _tmpVectorNoLabel ;

    }

  }
  std::ostringstream labelsTxtPathOss ;
  labelsTxtPathOss << outputDirectory << "labels.txt" ;
  std::string labelsTxtPath = labelsTxtPathOss.str() ;
  saveLabels( labelsTxtPath.c_str(), labels ) ;



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
  std::ostringstream comparisonWithAtlasRecognizedSBRPathOss ;
  comparisonWithAtlasRecognizedSBRPathOss << outputDirectory
                                                  << "comparisonWithAtlas.tsv" ;
  std::string comparisonWithAtlasRecognizedSBRPath =
                                 comparisonWithAtlasRecognizedSBRPathOss.str() ;

  std::vector<std::vector<int16_t>> labelsSBR ;
  std::vector<std::string> labelsDictSBR ;
  std::vector<float> coveragesBundles ;
  std::vector<float> adjacencyBundles ;
  std::vector<float> overlapBundles ;
  std::vector<float> matrixTracksSBR ;
  std::vector<int32_t> pointsPerTrackSBR ;
  int curves_countSBR = 0 ;
  int _bundleLabel = 0 ;
  // #pragma omp parallel for
  for ( int _bundleIndex = 0 ; _bundleIndex < nbBundles ; _bundleIndex++ )
  {

    std::string movedTractogramNeighborhood = neighborhoodFilenames[
                                                                _bundleIndex ] ;

    std::string atlasNeighborhoodFile ;
    if ( isAtlasNeighborhoodCentroids )
    {

      atlasNeighborhoodFile = atlasNeighborhoodCentroidsPaths[ _bundleIndex ] ;


    }
    else
    {

      atlasNeighborhoodFile = neighborhoodAtlasFilenames[ _bundleIndex ] ;

    }

    std::string _tmpAtlasBundlePath = atlasNeighborhoodFile ;

    char lastChar = _tmpAtlasBundlePath[ _tmpAtlasBundlePath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      _tmpAtlasBundlePath = _tmpAtlasBundlePath.substr( 0,
                                              _tmpAtlasBundlePath.size() - 1 ) ;

    }

    std::experimental::filesystem::path tmpPath( _tmpAtlasBundlePath ) ;
    std::string bundleName = tmpPath.stem() ;


    std::ostringstream recognizedBundleOss ;
    recognizedBundleOss << outputDirectory << bundleName ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      recognizedBundleOss << ".bundles" ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      recognizedBundleOss << ".minf" ;

    }
    std::string recognizedBundlePath = recognizedBundleOss.str() ;

    std::ostringstream recognizedBundledataOss ;
    recognizedBundledataOss << outputDirectory << bundleName ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      recognizedBundledataOss << ".bundlesdata" ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      recognizedBundledataOss << format ;

    }
    std::string recognizedBundledataPath = recognizedBundledataOss.str() ;

    if ( !is_file( recognizedBundlePath )  ||
                                          !is_file( recognizedBundledataPath ) )
    {

      continue ;

    }

    BundlesData recognizedBundlesData( recognizedBundledataPath.c_str() ) ;

    std::vector<float>& recognizedMatrixTracks =
                                            recognizedBundlesData.matrixTracks ;
    std::vector<int32_t>& recognizedPointsPerTrack =
                                          recognizedBundlesData.pointsPerTrack ;
    int recognizedCurves_count = recognizedBundlesData.curves_count ;
    curves_countSBR += recognizedCurves_count ;
    int64_t offset = 0 ;
    for ( int fiber = 0 ; fiber < recognizedCurves_count ; fiber++ )
    {

      int32_t nbPointsFiber = recognizedPointsPerTrack[ fiber ] ;
      pointsPerTrackSBR.push_back( nbPointsFiber ) ;
      for ( int point = 0 ; point < nbPointsFiber ; point++ )
      {

        for ( int coord = 0 ; coord < 3 ; coord++ )
        {

          matrixTracksSBR.push_back( recognizedMatrixTracks[ 3 * point + coord +
                                                                    offset ] ) ;

        }

      }

      offset += 3 * nbPointsFiber ;

      std::vector<int16_t> _tmpVectorLabel = { (int16_t)_bundleLabel } ;
      labelsSBR.push_back( _tmpVectorLabel ) ;

    }

    labelsDictSBR.push_back( bundleName ) ;
    _bundleLabel++ ;

    float _tmpCoverage = getCoverageWithAtlas( recognizedBundlePath ) ;
    coveragesBundles.push_back( _tmpCoverage ) ;
    float _tmpAdjacency = getAdjacencyWithAtlas( recognizedBundlePath ) ;
    adjacencyBundles.push_back( _tmpAdjacency ) ;
    float _tmpOverlap = getOverlapWithAtlas( recognizedBundlePath ) ;
    overlapBundles.push_back( _tmpOverlap ) ;

    if ( !saveBundlesSeparetly )
    {

      if ( is_file( recognizedBundlePath ) )
      {

        rmfile( recognizedBundlePath ) ;

      }

      if ( is_file( recognizedBundledataPath ) )
      {

        rmfile( recognizedBundledataPath ) ;

      }

    }

  }

  //#########################################################################//
  //#########################################################################//
  //#########################################################################//
  // Do here selection of one label per fiber
  //#########################################################################//
  //#########################################################################//
  //#########################################################################//
  saveComparisonMeasuresWithAtlas(
                                coveragesBundles,
                                adjacencyBundles,
                                overlapBundles,
                                labelsDictSBR,
                                comparisonWithAtlasRecognizedSBRPath.c_str() ) ;

  saveLabelsDict( labelsDictRecognizedSBRPath.c_str(),
                  labelsDictSBR ) ;

  saveLabels( labelsRecognizedSBRPath.c_str(),
              labelsSBR ) ;


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


  if ( keepTmpFiles )
  {

    const std::chrono::duration< double > duration =
                                std::chrono::system_clock::now() - start_time ;

    if ( verbose )
    {

      std::cout << "Duration : " << duration.count() << std::endl ;

    }

  }
  std::cout << "Done" << std::endl ;


  ////////////////////////////////// Cleaning //////////////////////////////////
  if ( keepTmpFiles )
  {

    return 0 ;

  }
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

  if ( is_dir( tmpAtlasDir ) )
  {

    rmdir( tmpAtlasDir ) ;

  }

  const std::chrono::duration< double > duration =
                              std::chrono::system_clock::now() - start_time ;

  if ( verbose )
  {

    std::cout << "Duration : " << duration.count() << std::endl ;

  }

  std::cout << "Done" << std::endl ;

  return( 0 ) ;

}
