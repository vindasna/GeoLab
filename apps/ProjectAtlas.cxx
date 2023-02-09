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

#include "ProjectAtlas.h"
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
///////////////////////////////////// Main /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  int index_input, index_atlas, index_output, index_labels, index_p, index_thr,
    index_minLen, index_maxLen, index_maxAngle, index_maxDirectionAngle,
    index_minShapeAngle, index_maxShapeAngle, index_thrDBMP, index_thrSim,
    index_minNbFibers, index_thrAdj, index_tolP, index_tolThr,
    index_tolMaxAngle, index_tolMaxDirectionAngle, index_tolMinShapeAngle,
    index_tolMaxShapeAngle, index_tolLenght, index_tolDistBetMedPts,
    index_useAvgThr, index_useMDF, index_useSimple, index_useMPAF,
    index_compare, index_saveExtractedBundles, index_save_unlabeled,
    index_nbThreads, index_verbose, index_help ;

  std::vector<std::string> possibleFlags{ "-i", "-a", "-o", "-l", "-p", "-thr",
                                          "-minLen", "-maxLen", "-maxAng",
                                          "-maxAngDir", "-minAngShape",
                                          "-maxAngShape", "-thrDBMP", "-thrSim",
                                          "-minNbFibers", "-thrAdj", "-tolP",
					                                "-tolThr", "-tolMaxAng",
					                                "-tolMaxDirAng", "-tolMinShapeAng",
					                                "-tolMaxShapeAng", "-tolLenght",
                                          "-tolDBMP", "-useAvgThr", "-useMPAF",
					                                "-useMDF", "-useSimple", "-cb",
					                                "-seb", "-su", "-nbThreads", "-v",
                                                                        "-h" } ;

  std::vector<bool> possibleFlagsNeedArgument{
                              true, true, true, true, true, true, true, true,
                              true, true, true, true, true, true, true, true,
                              true, true, true, true, true, true, true, true,
                              true, true, false, false, false, true,
			                        false, true, true, false } ;



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
  index_atlas = getFlagPosition( argc, argv, "-a" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_labels = getFlagPosition( argc, argv, "-l" ) ;
  index_p = getFlagPosition( argc, argv, "-p" ) ;
  index_thr = getFlagPosition( argc, argv, "-thr" ) ;
  index_minLen = getFlagPosition( argc, argv, "-minLen" ) ;
  index_maxLen = getFlagPosition( argc, argv, "-maxLen" ) ;
  index_maxAngle = getFlagPosition( argc, argv, "-maxAng" ) ;
  index_maxDirectionAngle = getFlagPosition( argc, argv, "-maxAngDir" ) ;
  index_minShapeAngle = getFlagPosition( argc, argv, "-minAngShape" ) ;
  index_maxShapeAngle = getFlagPosition( argc, argv, "-maxAngShape" ) ;
  index_thrDBMP = getFlagPosition( argc, argv, "-thrDBMP" ) ;
  index_thrSim = getFlagPosition( argc, argv, "-thrSim" ) ;
  index_minNbFibers = getFlagPosition( argc, argv, "-minNbFibers" ) ;
  index_thrAdj = getFlagPosition( argc, argv, "-thrAdj") ;
  index_tolP = getFlagPosition( argc, argv, "-tolP" ) ;
  index_tolThr = getFlagPosition( argc, argv, "-tolThr" ) ;
  index_tolMaxAngle = getFlagPosition( argc, argv, "-tolMaxAng" ) ;
  index_tolMaxDirectionAngle = getFlagPosition( argc, argv, "-tolMaxDirAng" ) ;
  index_tolMinShapeAngle = getFlagPosition( argc, argv, "-tolMinShapeAng" ) ;
  index_tolMaxShapeAngle = getFlagPosition( argc, argv, "-tolMaxShapeAng" ) ;
  index_tolLenght = getFlagPosition( argc, argv, "-tolLenght" ) ;
  index_tolDistBetMedPts = getFlagPosition( argc, argv, "-tolDBMP" ) ;
  index_useAvgThr = getFlagPosition( argc, argv, "-useAvgThr" ) ;
  index_useMPAF = getFlagPosition( argc, argv, "-useMPAF" ) ;
  index_useMDF = getFlagPosition( argc, argv, "-useMDF" ) ;
  index_useSimple = getFlagPosition( argc, argv, "-useSimple" ) ;
  index_compare = getFlagPosition( argc, argv, "-cb" ) ;
  index_saveExtractedBundles = getFlagPosition( argc, argv, "-seb" ) ;
  index_save_unlabeled = getFlagPosition( argc, argv, "-su" ) ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-i : Input tractogram from which to extract bundles \n"
              << "-a : Directory with the atlas (one file per bundle) \n"
              << "-o : Output directory to save the extracted bundles \n"
              << "-l : Name to give to the labels files ( name.txt and"
              << "name.dict) \n"
              << "[-p] : Threshold for the distance between the centers of the "
              << "track in the atlas and the track in the tractogram. It "
              << "impacts if the MMEA distance is computed for a fiber which "
              << "impacts the time per fiber (Default = 5 mm)\n"
              << "[-thr] : Threshold (in mm) of the MDA distance to decide if "
              << "the fiber belong to a bundle in the atlas (Default = 10 mm) \n"
              << "[-minLen] : Minimum length for fibers (Default = 20 mm) \n"
              << "[-maxLen] : Maximum length for fibers (Default = 200 mm) \n"
              << "[-maxAng] : Maximum angle between atlas fibers and tractogram"
              << " fibers in degrees (Default = 45 degrees) \n"
              << "[-maxAngDir] : Maximum angle between directions of atlas and "
              << "tractogram fibers in degrees (Default = 45 degrees) \n"
              << "[-minAngShape] : Minimum angle between vectors center fiber "
              << "- first point fiber and center fiber - last point fiber "
              << "(Default = 0 degrees) \n"
              << "[-maxAngShape] : Maximum angle between vectors center fiber "
              << "- first point fiber and center fiber - last point fiber "
              << "(Default = 180 degrees) \n"
              << "[-thrDBMP] : Threshold for maximum distance between medial "
              << "points (default = 5 mm) \n"
              << "[-thrSim] : Threshold for percentage of similarity in "
              << "projection i.e nbAtlasBundleFibersSimilar / "
              << "nbFibersAtlasBundle (Range [ 0 ; 1 ], default = 0.05) \n"
              << "[-minNbFibers] : Minimum numbers of fibers to validate an "
              << "extracted bundle (default = 20)\n"
              << "[-thrAdj] : Threshold for adjacency (Default = 10 mm) \n"
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
              << "[-tolDBMP] : Tolerance for distance between medial points "
              << "(for advanced users, default = 0) \n"
	            << "[-useAvgThr] : Use average instead of max for thr ( for "
	            << "advanced users, default : false ) \n"
              << "[-useMPAF] : Use medial point average fiber to determine "
              << "distance to center of bundle. Otherwise, use average of "
              << "medial points of fibers in bundle (deafault = true) \n"
              << "[-useMDF] : Use MDF distance for labeling\n"
              << "[-useSimple] : Use simple projection i.e do not use angle "
              << "between planes of bundles and angle between direction of "
              << "bundles for projection \n"
              << "[-cb] : Compute distance between recognized bundles and "
              << "atlas bundes \n"
              << "[-seb] : Save extracted .bundles (not only the labels file)\n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads "
              << "(default : number of cores ) \n"
              << "[-su] : Save bundle containing unlabeled fibers \n"
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

  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_labels )
  {

    std::cout << "-l argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_p )
  {

    p = std::stof( argv[ index_p + 1 ] ) ;
    if ( p <= 0 )
    {

      std::cout << "Error argument: p must be positive ( p > 0 ), got p = "
                << p << std::endl ;
      exit( 1 ) ;

    }
    useDefautlP = true ;

  }

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

  if ( index_minLen )
  {

    minimumLenghtFiber = std::stof( argv[ index_minLen + 1 ] ) ;
    if ( minimumLenghtFiber <= 0 )
    {

      std::cout << "Error argument: minimumLenghtFiber must be positive ( "
                << "minimumLenghtFiber > 0 ), got minimumLenghtFiber = "
                 << minimumLenghtFiber << std::endl ;
      exit( 1 ) ;

    }

    useDefaultMinLen = true ;

  }

  if ( index_maxLen )
  {

    maximumLenghtFiber = std::stof( argv[ index_maxLen + 1 ] ) ;
    if ( maximumLenghtFiber <= 0 )
    {

      std::cout << "Error argument: maximumLenghtFiber must be positive ( "
                << "maximumLenghtFiber > 0 ), got maximumLenghtFiber = "
                << maximumLenghtFiber << std::endl ;
      exit( 1 ) ;

    }
    else if ( maximumLenghtFiber <= minimumLenghtFiber )
    {

      std::cout << "Error argument: maximumLenghtFiber must be superior to "
                << "manimumLenghtFiber, got maximumLenghtFiber = "
                << maximumLenghtFiber << " and minimumLenghtFiber = "
                << minimumLenghtFiber << std::endl ;
      exit( 1 ) ;

    }

    useDefaultMaxLen = true ;

  }

  if ( index_maxAngle )
  {

    maxAngle = std::stof( argv[ index_maxAngle + 1 ] ) ;

    if ( maxAngle < 0 || maxAngle > 90 )
    {

      std::cout << "Error argument : maxAng has to be between 0 and 90 degrees"
                << std::endl ;
      exit( 1 ) ;

    }

    useDefaultMaxAngle = true ;

  }

  if ( index_maxDirectionAngle )
  {

    maxDirectionAngle = std::stof( argv[ index_maxDirectionAngle + 1 ] ) ;

    if ( maxDirectionAngle < 0 || maxDirectionAngle > 90 )
    {

      std::cout << "Error argument : maxAngDir must be between 0 and 90 degrees"
                << std::endl ;
      exit( 1 ) ;

    }

    useDefautlMaxDirectionAngle = true ;

  }

  if ( index_minShapeAngle )
  {

    minShapeAngle = std::stof( argv[ index_minShapeAngle + 1 ] ) ;

    if ( minShapeAngle < 0 || minShapeAngle > 180 )
    {

      std::cout << "Error argument : maxAngShape must be between 0 and 90 "
                << "degrees" << std::endl ;
      exit( 1 ) ;

    }

    useDefautlMinShapeAngle = true ;

  }

  if ( index_maxShapeAngle )
  {

    maxShapeAngle = std::stof( argv[ index_maxShapeAngle + 1 ] ) ;

    if ( maxShapeAngle < 0 || maxShapeAngle > 180 )
    {

      std::cout << "Error argument : maxAngShape must be between 0 and 90 "
                << "degrees" << std::endl ;
      exit( 1 ) ;

    }

    useDefautlMaxShapeAngle = true ;

  }

  if ( index_thrDBMP )
  {

    thrDistanceBetweenMedialPoints = std::stof( argv[ index_thrDBMP + 1 ] ) ;

    if ( thrDistanceBetweenMedialPoints < 0 )
    {

      std::cout << "Error argument : thrDBMP must be positive " << std::endl ;
      exit( 1 ) ;

    }

    useDefaultThrDistanceBetweenMedialPoints = true ;


  }

  if ( index_thrSim )
  {

    thrPercentageSimilarity = std::stof( argv[ index_thrSim + 1 ] ) ;

    if ( thrPercentageSimilarity < 0 || thrPercentageSimilarity > 1 )
    {

      std::cout << "Error argument : thrSim must be between 0 and 1 "
                << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_minNbFibers )
  {

    minimumNumberFibers = std::stoi( argv[ index_minNbFibers + 1 ] ) ;

    if ( minimumNumberFibers <= 0 )
    {

      std::cout << "Error argument : minNbFibers must be greater than 0"
                << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_thrAdj )
  {

    thresholdAdjacency = std::stof( argv[ index_thrAdj + 1 ] ) ;

    if ( thresholdAdjacency <= 0 )
    {

      std::cout << "Error argument : thrAdj must be greater than 0"
                << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_tolP )
  {

    toleranceP = std::stof( argv[ index_tolP + 1 ] ) ;

    // if ( toleranceP < -1 || toleranceP > 1 )
    // {
    //
    //   std::cout << "Error argument : -tolP must be in [-1;1]" << std::endl ;
    //   exit( 1 ) ;
    //
    // }

  }

  if ( index_tolThr )
  {

    toleranceThr = std::stof( argv[ index_tolThr + 1 ] ) ;

    // if ( toleranceThr < -1 || toleranceThr > 1 )
    // {
    //
    //   std::cout << "Error argument : -tolThr must be in [-1;1]" << std::endl ;
    //   exit( 1 ) ;
    //
    // }

  }

  if ( index_tolMaxAngle )
  {

    toleranceMaxAngle = std::stof( argv[ index_tolMaxAngle + 1 ] ) ;

    // if ( toleranceMaxAngle < -1 || toleranceMaxAngle > 1 )
    // {
    //
    //   std::cout << "Error argument : -tolMaxAng must be in [-1;1]" << std::endl ;
    //   exit( 1 ) ;
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
    //   std::cout << "Error argument : -tolMaxDirAng must be in [-1;1]"
    //                                                               << std::endl ;
    //   exit( 1 ) ;
    //
    // }

  }

  if ( index_tolMinShapeAngle )
  {

    toleranceMinShapeAngle = std::stof( argv[ index_tolMinShapeAngle + 1 ] ) ;

    // if ( toleranceMinShapeAngle < -1 || toleranceMinShapeAngle > 1 )
    // {
    //
    //   std::cout << "Error argument : -tolMinShapeAng must be in [-1;1]"
    //                                                               << std::endl ;
    //   exit( 1 ) ;
    //
    // }

  }

  if ( index_tolMaxShapeAngle )
  {

    toleranceMaxShapeAngle = std::stof( argv[ index_tolMaxShapeAngle + 1 ] ) ;

    // if ( toleranceMaxShapeAngle < -1 || toleranceMaxShapeAngle > 1 )
    // {
    //
    //   std::cout << "Error argument : -tolMaxShapeAng must be in [-1;1]"
    //                                                               << std::endl ;
    //   exit( 1 ) ;
    //
    // }

  }

  if ( index_tolLenght )
  {

    toleranceLenght = std::stof( argv[ index_tolLenght + 1 ] ) ;

    // if ( toleranceLenght < -1 || toleranceLenght > 1 )
    // {
    //
    //   std::cout << "Error argument : -tolLenght must be in [-1;1]"
    //                                                               << std::endl ;
    //   exit( 1 ) ;
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
    //   std::cout << "Error argument : -tolDBMP must be greater than -1"
    //                                                               << std::endl ;
    //   exit( 1 ) ;
    //
    // }

  }

  if ( index_useMPAF )
  {

    std::string tmpParsed( argv[ index_useMPAF + 1 ] ) ;

    if ( tmpParsed == "True" || tmpParsed == "TRUE" || tmpParsed == "true" )
    {

      useMedialPointAverageFiber = true ;

    }
    else if ( tmpParsed == "False" || tmpParsed == "FALSE" ||
                                                         tmpParsed == "false" )
    {

      useMedialPointAverageFiber = false ;

    }
    else
    {

      std::cout << "ERROR : Invalid argument for '-useMPAF', only posible "
                << "options are {TRUE, True, true, FALSE, False, false } \n" ;
      exit( 1 ) ;

    }


  }
  else
  {

    useMedialPointAverageFiber = true ;

  }

  if ( index_useAvgThr )
  {

    std::string tmpParsed( argv[ index_useAvgThr + 1 ] ) ;

    if ( tmpParsed == "True" || tmpParsed == "TRUE" || tmpParsed == "true" )
    {

      useAvgThr = true ;

    }
    else if ( tmpParsed == "False" || tmpParsed == "FALSE" ||
                                                         tmpParsed == "false" )
    {

      useAvgThr = false ;

    }
    else
    {

      std::cout << "ERROR : Invalid argument for '-useAvgThr', only posible "
                << "options are {TRUE, True, true, FALSE, False, false } \n" ;
      exit( 1 ) ;

    }


  }
  else
  {

    useAvgThr = false ;

  }


  if ( index_useMDF )
  {

    useMDFDistance = true ;

  }
  else
  {

    useMDFDistance = false ;

  }

  if ( index_useSimple )
  {

    useSimpleProjection = true ;

  }
  else
  {

    useSimpleProjection = false ;

  }

  if ( index_compare )
  {

    compareRecognizedWithAtlas = true ;

  }
  else
  {

    compareRecognizedWithAtlas = false ;

  }


  if ( index_saveExtractedBundles )
  {

    std::string tmpParsed( argv[ index_saveExtractedBundles + 1 ] ) ;

    if ( tmpParsed == "True" || tmpParsed == "TRUE" || tmpParsed == "true" )
    {

      saveExtractedBundles = true ;

    }
    else if ( tmpParsed == "False" || tmpParsed == "FALSE" ||
                                                         tmpParsed == "false" )
    {

      saveExtractedBundles = false ;

    }
    else
    {

      std::cout << "ERROR : Invalid argument for '-seb', only posible "
                << "options are {TRUE, True, true, FALSE, False, false } \n" ;
      exit( 1 ) ;

    }


  }
  else
  {

    saveExtractedBundles = true ;

  }


  if ( index_save_unlabeled )
  {

    saveUnlabeled = true ;

  }
  else
  {

    saveUnlabeled = false ;

  }

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

    // nbThreads = omp_get_num_procs() ;
    nbThreads = -1 ;

  }

  omp_set_num_threads( nbThreads ) ;

  #pragma omp parallel
  {

    nbThreadsUsed = omp_get_num_threads() ;

  }
  std::cout << "Number of threads : " << nbThreadsUsed << std::endl ;


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

  std::string inputBundlesFilename( argv[ index_input + 1 ] ) ;
  char lastChar = inputBundlesFilename[ inputBundlesFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputBundlesFilename = inputBundlesFilename.substr( 0,
                                             inputBundlesFilename.size() - 1 ) ;

  }

  std::string atlasDirectory( argv[ index_atlas + 1 ] ) ;
  lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    atlasDirectory = atlasDirectory + "/" ;

  }

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

  std::string labelsName( argv[ index_labels + 1 ] ) ;
  lastChar = labelsName[ labelsName.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    labelsName = labelsName.substr( 0, labelsName.size() - 1 ) ; ;

  }


  //   xxxxxxxxxxxxxxxxxxxxx Reading Input tractogram xxxxxxxxxxxxxxxxxxxxx   //
  if ( !is_file( inputBundlesFilename ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ProjectAtlas : input tractogram " << inputBundlesFilename
                  << " does not exists" << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }



  BundlesData inputTractogram( inputBundlesFilename.c_str() ) ;
  BundlesMinf inputTractogramInfo( inputBundlesFilename.c_str() ) ;

  //   xxxxxxxxxxxxxxxxxxxxxxxxxx Reading Atlas xxxxxxxxxxxxxxxxxxxxxxxxxx   //

  AtlasBundles atlasData( atlasDirectory.c_str(),
                          inputTractogramInfo.isBundles,
                          inputTractogramInfo.isTrk,
                          inputTractogramInfo.isTck,
                          verbose ) ;

  //////////////////////////////// Sanity cheks ////////////////////////////////
  int nbFibersTractogram = inputTractogram.curves_count ;
  int nbPoints = inputTractogram.pointsPerTrack[ 0 ] ;

  for ( int fiber = 1 ; fiber < nbFibersTractogram ; fiber++ )
  {

    if ( inputTractogram.pointsPerTrack[ fiber]  != nbPoints )
    {

      std::cout << "Error tractogram : the number of points in each fiber "
                << "of the tractogram has to be the same, got "
                << inputTractogram.pointsPerTrack[ fiber ]
                << " for fiber " << fiber << " and " << nbPoints
                << " for fiber 0 " << std::endl ;
      exit( 1 ) ;

    }

  }

  int nbBundlesAtlas = atlasData.bundlesData.size() ;

  for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
  {

    int nbFibersInBundle = atlasData.bundlesData[ bundle ].curves_count ;
    for ( int fiber = 0 ; fiber < nbFibersInBundle ; fiber++ )
    {

      if ( atlasData.bundlesData[ bundle ].pointsPerTrack[ fiber ] !=
                                                                     nbPoints )
      {

        std::cout << "Error atlas : the number of points in each fiber of the "
                  << "atlas has to be the sames as the number of points in "
                  << "the tractogram, got nbPointsAtlas = "
                  << atlasData.bundlesData[ bundle ].pointsPerTrack[ fiber ]
                  << " and nbPointsTractogram = " << nbPoints << std::endl ;
        exit( 1 ) ;

      }

    }

  }



  ///////////////////////////////// Projection /////////////////////////////////
  RecognizedBundles recognized( verbose,
                                p,
                                thrDistance ,
                                minimumLenghtFiber,
                                maximumLenghtFiber,
                                maxAngle,
                                maxDirectionAngle,
                                minShapeAngle,
                                maxShapeAngle,
                                compareRecognizedWithAtlas,
                                useMDFDistance,
                                useMedialPointAverageFiber,
                                useSimpleProjection,
                                toleranceP,
                                toleranceThr,
                                toleranceMaxAngle,
                                toleranceMaxDirectionAngle,
                                toleranceMinShapeAngle,
                                toleranceMaxShapeAngle,
                                toleranceLenght,
                                toleranceDistanceBetweenMedialPoints,
                                thrPercentageSimilarity,
                                thrDistanceBetweenMedialPoints,
                                minimumNumberFibers,
                                thresholdAdjacency,
                                useDefautlP,
                                useDefaultThr,
                                useDefaultMaxAngle,
                                useDefautlMaxDirectionAngle,
                                useDefautlMinShapeAngle,
                                useDefautlMaxShapeAngle,
                                useDefaultThrDistanceBetweenMedialPoints,
                                useDefaultMinLen,
                                useDefaultMaxLen,
                                saveExtractedBundles,
                                saveUnlabeled,
			                          useAvgThr,
                                nbThreads ) ;

    recognized.projectAtlas( atlasData,
                             inputTractogram,
                             thresholdAdjacency,
                             outputDirectory,
                             labelsName ) ;


  return 0 ;

}
