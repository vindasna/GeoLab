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

#include "statisticalFiltering.h"
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
void filterBundle( const BundlesData& atlasBundleData,
                   const BundlesMinf& atlasBundleInfo,
                   float tolerance,
                   std::vector<bool>& fibersToKeep )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;
  int nbPoints = atlasBundleData.pointsPerTrack[ 0 ] ;

  fibersToKeep.resize( nbFibersAtlasBundle ) ;


  //////////////////////////////////////////////////////////////////////////////

  if ( verbose )
  {

    std::cout << "Computing features for bundle " << atlasBundleInfo.bundleName
                                                                  << std::endl ;

    if ( verbose > 1 )
    {

      std::cout << "Number of points : " << nbPoints << "\nNumber of fibers in "
                << "input tractogram : " << nbFibersAtlasBundle << std::endl ;

    }

  }


  // Getting average features bundle
  float maxLength = atlasBundleInfo.maxLength * ( 1 + tolerance );
  float minLength = atlasBundleInfo.minLength * ( 1 - tolerance );
  float maxDistanceMedialPointFiberToAtlas = atlasBundleInfo.maxRadius * ( 1 + tolerance );
  float maxShapeAngle = atlasBundleInfo.maxShapeAngle * ( 1 + tolerance );
  float minShapeAngle = atlasBundleInfo.minShapeAngle * ( 1 - tolerance );
  float maxAngleBetweenPlanes = atlasBundleInfo.maxAngle * ( 1 + tolerance );
  float maxDirectionAngle = atlasBundleInfo.maxDirectionAngle * ( 1 + tolerance );
  float maxDMDA = atlasBundleInfo.maxDisimilarity * ( 1 + tolerance );

  //////////////////////////////////////////////////////////////////////////////
  std::vector<float> centerBundleAtlas = atlasBundleInfo.centerBundle ;

  #pragma omp paraelle for schedule(dynamic)
  for ( int fiberIndex1 = 0; fiberIndex1 < nbFibersAtlasBundle ; fiberIndex1++ )
  {

    // length fiber
    float lengthFiber = atlasBundleData.computeLengthFiber(
                                                   atlasBundleData.matrixTracks,
                                                   fiberIndex1,
                                                   nbPoints ) ;
    
    if ( lengthFiber > maxLength || lengthFiber < minLength )
    {

        fibersToKeep[ fiberIndex1 ] = false ;
        continue ;
    
    }

    

    // Distance to center of gravity and medial point fiber 1
    std::array<float, 3> medialPointFiber1{0, 0, 0} ;
    atlasBundleData.computeMedialPointFiberWithDistance( fiberIndex1,
                                                         medialPointFiber1 ) ;
    float tmpDistanceToCenterBundle = 0 ;
    for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
    {

        tmpDistanceToCenterBundle += pow( medialPointFiber1[ tmpCoord ] - 
                                          centerBundleAtlas[ tmpCoord ], 2 ) ;

    }
    tmpDistanceToCenterBundle = sqrt( tmpDistanceToCenterBundle ) ;

    if ( tmpDistanceToCenterBundle > maxDistanceMedialPointFiberToAtlas )
    {

        fibersToKeep[ fiberIndex1 ] = false ;
        continue ;

    }
    


    // Fiber points
    std::array<float, 3> endPoint1{0, 0, 0} ;
    std::array<float, 3> endPoint2{0, 0, 0} ;
    std::vector<float> fiber1Points( 3 * nbPoints, 0 ) ;
    int64_t offsetTractogram = 3 * nbPoints * fiberIndex1 ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      endPoint1[ i ] = atlasBundleData[ 3 * 0 + i + offsetTractogram ] ;
      endPoint2[ i ] = atlasBundleData[ 3 * ( nbPoints - 1 ) + i
                                                          + offsetTractogram ] ;

    }

    float distanceToOriginEndPoint1 = 0 ;
    float distanceToOriginEndPoint2 = 0 ;
    for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
    {

      distanceToOriginEndPoint1 += pow( endPoint1[ tmpCoord ] - 
                                                   endPoint1[ tmpCoord ], 2 ) ;
      distanceToOriginEndPoint2 += pow( endPoint2[ tmpCoord ] - 
                                                   endPoint2[ tmpCoord ], 2 ) ;

    }

    distanceToOriginEndPoint1 = sqrt( distanceToOriginEndPoint1 ) ;
    distanceToOriginEndPoint2 = sqrt( distanceToOriginEndPoint2 ) ;

    if ( distanceToOriginEndPoint1 < distanceToOriginEndPoint2 )
    {

      for ( int tmpPoint = 0 ; tmpPoint < nbPoints ; tmpPoint++ )
      {

        for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
        {

          fiber1Points[ 3 * tmpPoint + tmpCoord ] = 
                              atlasBundleData[ 3 * tmpPoint + tmpCoord +
                                                           offsetTractogram ] ;
          

        }

      }
    }
    else
    {

      for ( int tmpPoint = 0 ; tmpPoint < nbPoints ; tmpPoint++ )
      {

        for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
        {

          fiber1Points[ 3 * tmpPoint + tmpCoord ] = 
                              atlasBundleData[ 3 * ( nbPoints - tmpPoint - 1 ) 
                                             + tmpCoord + offsetTractogram  ] ;
          

        }

      }

    }
    


    // Shape angle
    std::array<float, 3> point1{0, 0, 0} ;
    std::array<float, 3> point2{0, 0, 0} ;

    std::array<float, 3> vector1{0, 0, 0} ;
    std::array<float, 3> vector2{0, 0, 0} ;
    for ( int i = 0 ; i < 3 ; i++ )
    {

      point1[ i ] = atlasBundleData[ 3 * 0 + i + offsetTractogram ] ;
      point2[ i ] = atlasBundleData[ 3 * ( nbPoints - 1 ) + i
                                                          + offsetTractogram ] ;

      vector1[ i ] = point1[ i ] - medialPointFiber1[ i ] ;
      vector2[ i ] = point2[ i ] - medialPointFiber1[ i ] ;

    }

    float shapeAngle = atlasBundleData.computeAngleBetweenVectors(
                                                            vector1, vector2 ) ;
    
    if ( shapeAngle > maxShapeAngle || shapeAngle < minShapeAngle )
    {

        fibersToKeep[ fiberIndex1 ] = false ;
        continue ;


    }

    // Normal vector fiber 1
    std::array<float, 3> normalVectorFiber1{0, 0, 0} ;
    atlasBundleData.computeNormalVectorFiberTractogram(
                                                atlasBundleData.matrixTracks,
                                                medialPointFiber1,
                                                fiberIndex1,
                                                nbPoints,
                                                normalVectorFiber1 ) ;
    
    // Direction vector fiber 1
    std::array<float, 3> directionVectorFiber1{0, 0, 0} ;
    atlasBundleData.computeDirectionVectorFiberTractogram(
                                                atlasBundleData.matrixTracks,
                                                medialPointFiber1,
                                                normalVectorFiber1,
                                                fiberIndex1,
                                                nbPoints,
                                                directionVectorFiber1 ) ;
    
    bool isBreak = false ;
    for ( int fiberIndex2 = fiberIndex1 + 1 ; fiberIndex2 < nbFibersAtlasBundle ; fiberIndex2++ )
    {

        // No need to compute length, shape angle or distance to center of gravity for fiber 2 
        // as througth the loop fiber 1 will go through all the fibers in the bundle

        // Medial point fiber 2
        std::array<float, 3> medialPointFiber2{0, 0, 0} ;
        atlasBundleData.computeMedialPointFiberWithDistance( fiberIndex1,
                                                             medialPointFiber2 ) ;

        // Normal vector fiber 2
        std::array<float, 3> normalVectorFiber2{0, 0, 0} ;
        atlasBundleData.computeNormalVectorFiberTractogram(
                                                    atlasBundleData.matrixTracks,
                                                    medialPointFiber2,
                                                    fiberIndex2,
                                                    nbPoints,
                                                    normalVectorFiber2 ) ;
        
        float angleBetweenPlanes = atlasBundleData.computeAngleBetweenPlanes( 
                                                             normalVectorFiber1, normalVectorFiber2 ) ;
        if ( angleBetweenPlanes > maxAngleBetweenPlanes )
        {

            fibersToKeep[ fiberIndex1 ] = false ;
            isBreak = true ;
            break ;

        }
        

        

        // Direction vector fiber 2
        std::array<float, 3> directionVectorFiber2{0, 0, 0} ;
        atlasBundleData.computeDirectionVectorFiberTractogram(
                                                    atlasBundleData.matrixTracks,
                                                    medialPointFiber2,
                                                    normalVectorFiber2,
                                                    fiberIndex2,
                                                    nbPoints,
                                                    directionVectorFiber2 ) ;
        float angleBetweenDirections = atlasBundleData.computeAngleBetweenDirections( 
                                                       directionVectorFiber1, directionVectorFiber2 ) ;
        if ( angleBetweenDirections > maxDirectionAngle )
        {

            fibersToKeep[ fiberIndex1 ] = false ;
            isBreak = true ;
            break ;

        }
        


        std::vector<float> fiber2ToFiber1( 3 * nbPoints, 0 ) ;
        std::array<float, 3> normalVectorFiber2ToFiber1{0, 0, 0} ;
        atlasBundleData.registerFiber(
                                 atlasBundleData.matrixTracks,
                                 atlasBundleData.matrixTracks,
                                 fiberIndex1,
                                 fiberIndex2,
                                 nbPoints,
                                 fiber2ToFiber1,
                                 normalVectorFiber2ToFiber1 ) ;
        
        
        float tmpDistanceMDAD = 
                 atlasBundleData.computeMDADBetweenTwoFibersAfterAlignement(
                                            fiber2ToFiber1,
                                            atlasBundleData.matrixTracks,
                                            0,
                                            fiberIndex1,
                                            true,
                                            nbPoints ) ;
        if ( tmpDistanceMDAD > maxDMDA )
        {

            fibersToKeep[ fiberIndex1 ] = false ;
            isBreak = true ;
            break ;

        }
        

    }

    if ( !isBreak )
    {

        fibersToKeep[ fiberIndex1 ] = true ;

    }

  }

}



////////////////////////////////////////////////////////////////////////////////
void saveFeaturesFibersBundle(
                     const std::vector<std::vector<float>>& fibersFeturesValues,
                     int nbPoints,
                     const char* savePath )
{

  std::string savePathStr = savePath ;

  std::ofstream file( savePathStr.c_str() ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path "
              << std::endl ;
    exit( 1 );

  }

  for ( int tmpPoint = 0 ; tmpPoint < nbPoints ; tmpPoint++ )
  {

    for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
    {

      file << "point_" << tmpPoint << "_" << tmpCoord << "\t" ;

    }

  }

  file << "lengthFiber\t"
       << "shapeAngle\t" ;
  
  for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
  {

    file << "normal_vector_" << tmpCoord << "\t" ;

  }

  for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
  {

    file << "direction_vector_" << tmpCoord << "\t" ;
    
  }

  file << "angleToDirectionX\t"
       << "angleToPlaneYZ" ;

  file << "\n" ;

  for ( std::vector<float>tmpFeatures : fibersFeturesValues )
  {

    for ( int featureIndex = 0 ; featureIndex < tmpFeatures.size() ; 
                                                               featureIndex++ )
    {

      float tmpValue = tmpFeatures[ featureIndex ] ;

      file << tmpValue << "\t" ;

    }

    file << "\n" ;

  }

  file.close() ;

}
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  int index_input, index_tractogram, index_labels, index_dict, index_output, 
                  index_format, index_tolerance, index_nbThreads, index_verbose, 
                                                                    index_help ;

  std::vector<std::string> possibleFlags{
                                        "-i", "-t", "-l", "-d", "-o", "-format", 
                                              "-tol", "-nbThreads", "-v", "-h" } ;

  std::vector<bool> possibleFlagsNeedArgument{ true, true, true, true, true, true,
                                                       true, true, true, false } ;

  if ( possibleFlagsNeedArgument.size() != possibleFlags.size() )
  {

    std::cout << "getFeaturesFIbersAtlas CODING ERROR : possibleFlags and "
              << "possibleFlagsNeedArgument must have the same size"
              << std::endl ;
    exit( 1 ) ;

  }



  for ( int k = 1 ; k < argc ; k++ )
  {

    std::string arg = argv[ k ] ;
    auto it = std::find( possibleFlags.begin(), possibleFlags.end(), arg ) ;
    if( it == possibleFlags.end() )
    {

      if ( k == 1 )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for getFeaturesFIbersAtlas command" << std::endl ;
        exit( 1 ) ;

      }

      std::string arg2 = argv[ k - 1 ] ;
      it = std::find( possibleFlags.begin(), possibleFlags.end(), arg2 ) ;
      if( it == possibleFlags.end() )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for getFeaturesFIbersAtlas command" << std::endl ;
        exit( 1 ) ;

      }
      int indexFound = it - possibleFlags.begin() ;
      if ( !possibleFlagsNeedArgument[ indexFound ] )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for getFeaturesFIbersAtlas command" << std::endl ;
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
  index_tractogram = getFlagPosition( argc, argv, "-t" ) ;
  index_labels = getFlagPosition( argc, argv, "-l" ) ;
  index_dict = getFlagPosition( argc, argv, "-d" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_format = getFlagPosition( argc, argv, "-format") ;
  index_tolerance = getFlagPosition( argc, argv, "-tol" ) ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to get features of fibers in atlas : \n"
              << "-i : Input directory with the bundles (two file per bundle)\n"
              << "-t : Input tractogram used for labeling \n"
              << "-l : Labels file (.txt)\n"
              << "-d : dictionary file (.dict)\n"
              << "-o : Output directory to save output files\n"
              << "-format : format of bundles in atlas among .bundles/.trk/"
              << ".tck \n"
              << "[-tol] : Tolerance (default = 0.5)\n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads "
              << "(default : number of cores ) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_tractogram )
  {

    std::cout << "-t argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_labels )
  {

    std::cout << "-l argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_dict )
  {

    std::cout << "-d argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_format )
  {

    std::cout << "-format argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_tolerance )
  {

    tolerance = std::stof( argv[ index_tolerance + 1 ] ) ;

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

  // Atlas directory

  std::string atlasDirectory = argv[ index_input + 1 ] ;
  char lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    atlasDirectory = atlasDirectory + "/" ;

  }
  if ( !is_dir( atlasDirectory ) )
  {

    std::cout << "Input atlas is a file" << std::endl ;

    lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      atlasDirectory = atlasDirectory.substr( 0,
                                              atlasDirectory.size() - 1 ) ;

    }

    if ( !is_file( atlasDirectory ) )
    {

      std::stringstream outMessageOss ;

      outMessageOss << "ERROR : Atlas directory/file " << atlasDirectory
                                        << " does not exists " << std::endl ;
      std::string outMessage = outMessageOss.str() ;
      throw( std::invalid_argument( outMessage ) ) ;

    }
    else
    {

      if ( !endswith( atlasDirectory, format ) )
      {

        std::cout << "ERROR : atlas bundle file must be in " << format
                  << " format" << std::endl ;
        exit( 1 ) ;

      }

    }

  }
  else
  {

    std::cout << "Input atlas is a directory" << std::endl ;

  }

  // Input tractogram

  std::string inputFilename( argv[ index_tractogram + 1 ] ) ;
  lastChar = inputFilename[ inputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputFilename = inputFilename.substr( 0, inputFilename.size() - 1 ) ;

  }

  std::string inputBundlesFilename ;
  std::string inputBundlesDataFilename ;
  std::string inputTRKFilename ;
  bool isBundlesFormat = false ;
  bool isTrkFormat = false ;
  bool isTckFormat = false ;

  if ( endswith( inputFilename, ".bundlesdata" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".bundles" ) ;

    isBundlesFormat = true ;

  }
  else if ( endswith( inputFilename, ".bundles" ) )
  {

    inputBundlesFilename = inputFilename ;

    inputBundlesDataFilename = replaceExtension( inputFilename,
                                                              ".bundlesdata" ) ;

    isBundlesFormat = true ;

  }
  else if ( endswith( inputFilename, ".trk" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".minf" ) ;

    isTrkFormat = true ;


  }
  else if ( endswith( inputFilename, ".tck" ) )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = replaceExtension( inputFilename, ".minf" ) ;

    isTckFormat = true ;


  }
  else
  {

    std::cout << "The only tractogram format supported are .bundles/.trk/.tck"
              << std::endl ;
    exit( 1 ) ;

  }


  std::string format ;
  if ( isBundlesFormat )
  {

    format = ".bundles" ;

  }
  else if ( isTrkFormat )
  {

    format = ".trk" ;

  }
  else if ( isTckFormat )
  {

    format = ".tck" ;

  }


  BundlesData inputTractogram( inputBundlesDataFilename.c_str() ) ;
  BundlesMinf inputTractogramInfo( inputBundlesFilename.c_str() ) ;


  // Reading labels
  std::string predictedLabelsFilename( argv[ index_labels + 1 ] ) ;
  lastChar = predictedLabelsFilename[ predictedLabelsFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    predictedLabelsFilename = predictedLabelsFilename.substr( 0,
                                       predictedLabelsFilename.size() - 1 ) ;

  }
  std::vector<std::vector<int16_t>> streamlinesLabels ;
  readingPredictedLabels( predictedLabelsFilename.c_str(), streamlinesLabels, 
                                              inputTractogram.curves_count ) ;


  std::vector<std::string> dictionaryLabels ;
  if ( index_dict )
  {

    std::string dictionaryLabelsFilename( argv[ index_dict + 1 ] ) ;
    lastChar = dictionaryLabelsFilename[ dictionaryLabelsFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      dictionaryLabelsFilename = dictionaryLabelsFilename.substr( 0,
                                         dictionaryLabelsFilename.size() - 1 ) ;

    }

    readLabelsDict( dictionaryLabelsFilename.c_str(),
                    dictionaryLabels ) ;


  }


  // Getting output directory

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


  // Getting format 
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

    std::cout << "Argument -if must be either .bundles/.trk/.tck "
              << std::endl ;
    exit( 1 ) ;

  }



  //   xxxxxxxxxxxxxxxxxxxxxxxxxx Reading Atlas xxxxxxxxxxxxxxxxxxxxxxxxxx   //

  AtlasBundles atlasData ;
  if ( is_dir( atlasDirectory ) )
  {

    atlasData.read( atlasDirectory.c_str(),
                    isBundlesFormat,
                    isTrkFormat,
                    isTckFormat,
                    verbose ) ;

  }
  else if ( is_file( atlasDirectory ) )
  {

    BundlesData tmpBundlesData( atlasDirectory.c_str() ) ;
    BundlesMinf tmpBundlesMinf( atlasDirectory.c_str() ) ;

    std::string tmpBundleName = basenameNoExtension( atlasDirectory ) ;

    atlasData.bundlesMinf.push_back( tmpBundlesMinf ) ;
    atlasData.bundlesData.push_back( tmpBundlesData ) ;
    atlasData.bundlesNames.push_back( tmpBundleName ) ;
    atlasData.isBundles = tmpBundlesData.isBundles ;
    atlasData.isTrk = tmpBundlesData.isTrk ;
    atlasData.isTck = tmpBundlesData.isTck ;

  }



  //////////////////////////////// Sanity cheks ////////////////////////////////
  int nbBundlesAtlas = atlasData.bundlesData.size() ;

  int nbPoints = atlasData.bundlesData[ 0 ].pointsPerTrack[ 0 ] ;
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


  //////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<int16_t>> newStreamlinesLabels( streamlinesLabels.size(), 
                                                  std::vector<int16_t>( 1, -1 ) ) ;
  std::vector<std::string> newDictionaryLabels = dictionaryLabels ;
  for ( int bundleIndex = 0 ; bundleIndex < nbBundlesAtlas ; bundleIndex++ )
  {

    std::cout << "Processing : [" << bundleIndex + 1 << "/" << nbBundlesAtlas
                                                    << "]            " << "\r" ;

    BundlesMinf& atlasBundleMinf = atlasData.bundlesMinf[ bundleIndex ] ;
    BundlesData& atlasBundleData = atlasData.bundlesData[ bundleIndex ] ;
    std::string atlasBundleName = atlasData.bundlesNames[ bundleIndex ] ;

    std::vector<bool> fibersToKeep ;
    filterBundle( atlasBundleData, atlasBundleMinf, tolerance, fibersToKeep ) ;

    BundlesMinf& newAtlasBundleMinf = atlasData.bundlesMinf[ bundleIndex ] ;
    BundlesData& newAtlasBundleData = atlasData.bundlesData[ bundleIndex ] ;

    int tmpLabel = getLabelFromName( newDictionaryLabels, atlasBundleName ) ;
    // Get fibers with labels
    std::vector<int64_t> tmpFibersWithLabel ;
    for ( int64_t tmpI = 0 ; tmpI < streamlinesLabels.size() ; tmpI++ )
    {

        std::vector<int16_t>& tmpLabelsFiber = streamlinesLabels[ tmpI ] ;
        if ( std::find( tmpLabelsFiber.begin(), tmpLabelsFiber.end(), 
                                                            tmpLabel ) != tmpLabelsFiber.end() )
        {

            tmpFibersWithLabel.push_back( tmpI ) ;

        }

    }

    if ( newAtlasBundleData.curves_count != tmpFibersWithLabel.size() )
    {

        std::cout << "ERROR : the number of streamlines in the atlas bundle has to be the same "
                  << "as the number of streamlines with the bundle label in the labels.txt " 
                  << "file. Got " << newAtlasBundleData.curves_count << " for the bunlde and "
                  << tmpFibersWithLabel.size() << " for the labels" << std::endl ;
        exit( 1 ) ;
        
    }


    std::vector<float> newStreamlines ;
    std::vector<int32_t> newPointsPerTrack ;
    int64_t newCurveCount = 0 ;
    std::vector<int64_t> newTmpFibersWithLabel ;
    
    for ( int fiberIndex = 0 ; fiberIndex < atlasBundleData.curves_count ; fiberIndex++ )
    {

        int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

        if ( fibersToKeep[ fiberIndex ] )
        {


            newCurveCount++ ;
            newPointsPerTrack.push_back( nbPoints ) ;


            for ( int point = 0 ; point < nbPoints ; point++ )
            {

                for ( int coord = 0 ; coord < 3 ; coord++ )
                {

                    newStreamlines.push_back( atlasBundleData.matrixTracks[ 3 * point + coord + offsetTractogram ] ) ;

                }

            }


            newTmpFibersWithLabel.push_back( tmpFibersWithLabel[ fiberIndex ] ) ;

        }

    }
    
    for ( int64_t tmpIndex : newTmpFibersWithLabel )
    {

        if ( newStreamlinesLabels[ tmpIndex ].size() == 1 && newStreamlinesLabels[ tmpIndex ][ 0 ] == -1 )
        {

            newStreamlinesLabels[ tmpIndex ] = std::vector<int16_t>( 1, tmpLabel ) ;

        }
        else
        {

            newStreamlinesLabels[ tmpIndex ].push_back( tmpLabel ) ;

        }


    }


    newAtlasBundleData.matrixTracks = newStreamlines ;
    newAtlasBundleData.curves_count = newCurveCount ;
    newAtlasBundleData.pointsPerTrack = newPointsPerTrack ;

    newAtlasBundleMinf.curves_count = newCurveCount ;


    std::stringstream outPathOss ;

    outPathOss << outputDirectory << atlasBundleName << format ;
    std::string outPath = outPathOss.str() ;

    newAtlasBundleData.write( outPath.c_str(), newAtlasBundleMinf ) ;


  }

  std::stringstream outPathLabelsOss ;

  outPathLabelsOss << outputDirectory << "labels.txt" ;
  std::string outPathLabels = outPathLabelsOss.str() ;

  saveLabels( outPathLabels.c_str(), newStreamlinesLabels ) ;

  std::stringstream outPathLabelsDictOss ;

  outPathLabelsDictOss << outputDirectory << "labels.dict" ;
  std::string outPathLabelsDict = outPathLabelsDictOss.str() ;

  saveLabelsDict( outPathLabelsDict.c_str(), newDictionaryLabels ) ;


  std::cout << "\nDone" << std::endl ;

}
