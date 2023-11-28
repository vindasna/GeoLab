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

#include "compareBundlesToAtlas.h"
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

  int index_input, index_atlas, index_reference, index_output, index_nbPoints,
                            index_thr, index_force, index_format, index_nbThreads, 
                                                        index_verbose, index_help ;

  index_input = getFlagPosition( argc, argv, "-i") ;
  index_atlas = getFlagPosition( argc, argv, "-a") ;
  index_reference = getFlagPosition( argc, argv, "-ref") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints") ;
  index_thr = getFlagPosition( argc, argv, "-thr" ) ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_format = getFlagPosition( argc, argv, "-f") ;
  index_force = getFlagPosition( argc, argv, "-force") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to compare bundles using to atlas : \n"
              << "-i : Path to the input directory with recognized bundles \n"
              << "-a : Path to the directory with the atlas \n"
              << "-f : bundles format among .bundles/.trk/.tck \n"
              << "[-o] : Path of the output file \n"
              << "[-nbPoints] : Number of points per fiber (same number for all"
              << " fibers) \n"
              << "[-thr] : Threshold for computing adjacency between fibers "
              << "(default = 10mm) \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads for parallel "
              << "regions (default : 1)\n"
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

  if ( !index_atlas )
  {

    std::cout << "-a argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_format )
  {

    std::cout << "-f argument required ..." << std::endl ;
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

  ///////////////// Get threshold for adjacency between fibers /////////////////
  if ( index_thr )
  {

    thresholdMetrics = std::stof( argv[ index_thr + 1 ] ) ;
    if ( thresholdMetrics <= 0 )
    {

      std::cout << "ERROR : -thr must be greater than 0" << std::endl ;
      exit( 1 ) ;

    }

  }

  ///////////////////////// If atlas is in .trk format /////////////////////////
  format = argv[ index_format + 1 ] ;
  if ( format == ".bundles" )
  {

    isBundlesFormat = true ;

  }
  else if ( format == ".trk" )
  {

    isTrkFormat = true ;

  }
  else if ( format == ".tck" )
  {

    isTckFormat = true ;

  }
  else
  {

    std::cout << "Argument -f must be either .bundles/.trk/.tck "
              << std::endl ;
    exit( 1 ) ;

  }



  ////////////////////////////// Input directoy //////////////////////////////
  std::string inputDirectoryPath( argv[ index_input + 1 ] ) ;

  char lastChar = inputDirectoryPath[ inputDirectoryPath.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    inputDirectoryPath = inputDirectoryPath + "/" ;

  }
  if ( !is_dir( inputDirectoryPath ) )
  {

    std::cout << "ERROR : Atlas directory " << inputDirectoryPath << " does not"
                                                    << " exists " << std::endl ;
    exit( 1 );

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

    std::cout << "ERROR : Atlas directory " << atlasDirectory << " does not "
                                                     << "exists " << std::endl ;
    exit( 1 );

  }

  ///////////////////////////////// Output file ////////////////////////////////
  std::string outputFilename ;
  if ( index_output )
  {

    outputFilename = argv[ index_output + 1 ] ;
    lastChar = outputFilename[ outputFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

    }

    if ( !endswith( outputFilename, ".tsv" ) )
    {

      std::ostringstream tmpOutputFilenameOss ;
      tmpOutputFilenameOss << outputFilename << ".tsv" ;
      outputFilename = tmpOutputFilenameOss.str() ;

    }

  }
  else
  {

    std::ostringstream tmpOutputFilenameOss ;
    tmpOutputFilenameOss << inputDirectoryPath << "comparisonWithAtlas.tsv" ;
    outputFilename = tmpOutputFilenameOss.str() ;

  }


  if ( is_file( outputFilename ) && !force )
  {

    std::cout << "ERROR : output file " << outputFilename << "already exists "
              << "and -force flag not used, exiting... "
              << std::endl ;
      exit( 1 ) ;

  }


  ///////////////////////// Number of points per fiber /////////////////////////
  if ( index_nbPoints )
  {

    nbPointsPerFiber = std::stoi( argv[ index_nbPoints + 1 ] ) ;

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

    nbThreads = 1 ;
    // nbThreads = -1 ;

  }

  omp_set_num_threads( nbThreads ) ;

  omp_set_nested( 1 ) ;
  
  
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

  ////////////////////// ///// Reading input directoy //////////////////////////
  AtlasBundles inputData( inputDirectoryPath.c_str(),
                          isBundlesFormat,
                          isTrkFormat,
                          isTckFormat,
                          0 ) ;

  //////////////////////////////// Reading atlas ///////////////////////////////
  AtlasBundles atlasData( atlasDirectory.c_str(),
                          isBundlesFormat,
                          isTrkFormat,
                          isTckFormat,
                          0 ) ;

  std::vector<std::string>& bundlesNames = inputData.bundlesNames ;
  int nbBundlesInput = bundlesNames.size() ;
  std::vector<float> coveragesBundles ;
  std::vector<float> adjacencyBundles ;
  std::vector<float> overlapBundles ;
  std::vector<float> disimilarityBundles ;
  std::vector<int> nbFibersBundles ;
  std::vector<std::string> outBundlesNames ;
  int tmpCounter = 1 ;
  #pragma omp parallel for num_threads(nbThreads) shared(nbFibersBundles) shared(disimilarityBundles) shared(overlapBundles) shared(adjacencyBundles) shared(coveragesBundles)
  for ( int _bundle = 0 ; _bundle < nbBundlesInput ; _bundle++ )
  {

    std::string _bundleName = bundlesNames[ _bundle ] ;

    if ( verbose == 1 )
    {

      #pragma omp critical
      {

        printf( "\rProcessing : [ %d  /  %d ]",
                                   tmpCounter, nbBundlesInput ) ;
        std::cout << "" << std::flush ;

      }

    }

    if ( verbose > 1 )
    {

      std::cout << "Comparing bundle " << _bundleName << std::endl ;

    }

    if ( _bundleName == "unlabeledFibers" )
    {

      continue ;

    }

    if ( _bundleName == "regroupedRecognized" )
    {

      continue ;

    }

    if ( verbose > 1 )
    {

      std::cout << "\tGetting bundle information and data from input folder and"
                << " atlas" << std::endl ;

    }

    int indexBundleInputData = inputData.findBundleIndexByName( _bundleName ) ;

    int indexBundleAtlasData = atlasData.findBundleIndexByName( _bundleName ) ;


    BundlesData& inputBundle = inputData.bundlesData[ indexBundleInputData ] ;


    BundlesData& atlasBundle = atlasData.bundlesData[ indexBundleAtlasData ] ;


    std::vector<float>& inputBundleTracks = inputBundle.matrixTracks ;
    std::vector<float>& atlasBundleTracks = atlasBundle.matrixTracks ;

    int nbCurvesInputBundle = inputBundle.curves_count ;
    int nbCurvesAtlasBundle = atlasBundle.curves_count ;


    std::vector<int32_t>& pointsPerTrackInputBundle =
                                                    inputBundle.pointsPerTrack ;
    int tmpNbPointsPerTrackInputBundle = pointsPerTrackInputBundle[ 0 ] ;
    for ( int _i = 1 ; _i < nbCurvesInputBundle ; _i++ )
    {

      if ( tmpNbPointsPerTrackInputBundle != pointsPerTrackInputBundle[ _i ] )
      {

        std::cout << "ERROR : in input directory, bundle " << _bundleName
                  << " has not the same number of points per track for all "
                  << "fibers"
                  << std::endl ;
        exit( 1 ) ;

      }

    }


    if ( nbPointsPerFiber == 0 )
    {

      nbPointsPerFiber = tmpNbPointsPerTrackInputBundle ;

    }


    if ( tmpNbPointsPerTrackInputBundle != nbPointsPerFiber )
    {

      std::cout << "ERROR : in input directory, bundle " << _bundleName
                << " has not the same number of points per track as the one "
                << "given by flag -nbPoints "
                << std::endl ;
      exit( 1 ) ;

    }

    std::vector<int32_t>& pointsPerTrackAtlasBundle =
                                                    atlasBundle.pointsPerTrack ;

    int tmpNbPointsPerTrackAtlasBundle = pointsPerTrackAtlasBundle[ 0 ] ;
    for ( int _i = 1 ; _i < nbCurvesAtlasBundle ; _i++ )
    {

      if ( tmpNbPointsPerTrackAtlasBundle != pointsPerTrackAtlasBundle[ _i ] )
      {

        std::cout << "ERROR : in atlas directory, bundle " << _bundleName
                  << " has not the same number of points per track for all "
                  << "fibers"
                  << std::endl ;
        exit( 1 ) ;

      }

    }


    if ( tmpNbPointsPerTrackAtlasBundle != nbPointsPerFiber )
    {

      std::cout << "ERROR : in atlas directory, bundle " << _bundleName
                << " has not the same number of points per track as the one "
                << "given by flag -nbPoints or as bundle in input directory "
                << std::endl ;
      exit( 1 ) ;

    }


    if ( verbose > 1 )
    {

      std::cout << "\tComputing metrics" << std::endl ;

    }
    float _tmpDisimilarity =  inputData.distanceBetweenBundles(
                                                         atlasBundle,
                                                         inputBundle,
                                                         nbThreads,
                                                         nbPointsPerFiber ) ;

    std::vector<int> nbAdjacentFibersInputToAtlas( nbCurvesInputBundle, 0 ) ;
    inputData.computeNumberAdjacentFibersBundle1ToBundle2(
                                                inputBundleTracks,
                                                atlasBundleTracks,
                                                nbCurvesInputBundle,
                                                nbCurvesAtlasBundle,
                                                nbPointsPerFiber,
                                                thresholdMetrics,
                                                nbThreads,
                                                nbAdjacentFibersInputToAtlas ) ;

    std::vector<int> nbAdjacentFibersAtlasToInput( nbCurvesAtlasBundle, 0 ) ;
    inputData.computeNumberAdjacentFibersBundle1ToBundle2(
                                                atlasBundleTracks,
                                                inputBundleTracks,
                                                nbCurvesAtlasBundle,
                                                nbCurvesInputBundle,
                                                nbPointsPerFiber,
                                                thresholdMetrics,
                                                nbThreads,
                                                nbAdjacentFibersAtlasToInput ) ;

    float coverageInputToAtlas = inputData.coverageRecognizedToAtlasBundles(
                                                nbAdjacentFibersInputToAtlas ) ;
    

    float coverageAtlasToInput = inputData.coverageRecognizedToAtlasBundles(
                                                nbAdjacentFibersAtlasToInput ) ;

    float _tmpOverlap = inputData.overlapRecognizedToAtlasBundles(
                                                nbAdjacentFibersInputToAtlas ) ;
    

    float _tmpAdjacency = inputData.bundlesAdjacency( coverageInputToAtlas,
                                                      coverageAtlasToInput ) ;
  

    if ( verbose > 1 )
    {

      std::cout << "\tAjacency : " << _tmpAdjacency << std::endl ;
      std::cout << "\nNumber of fibers : " << nbCurvesInputBundle << std::endl ;

    }

    #pragma omp critical
    {

      disimilarityBundles.push_back( _tmpDisimilarity ) ;
      coveragesBundles.push_back( coverageInputToAtlas ) ;
      overlapBundles.push_back( _tmpOverlap ) ;
      adjacencyBundles.push_back( _tmpAdjacency ) ;
      nbFibersBundles.push_back( nbCurvesInputBundle ) ;
      outBundlesNames.push_back( _bundleName ) ;

      tmpCounter += 1 ;

    }
    



  }


  saveComparisonMeasuresWithAtlas( coveragesBundles,
                                   adjacencyBundles,
                                   overlapBundles,
                                   disimilarityBundles,
                                   nbFibersBundles,
                                   outBundlesNames,
                                   outputFilename.c_str() ) ;

}
