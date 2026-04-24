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

#include "getFeaturesFibersAtlas.h"
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
void computeFeaturesFibersBundle(
                          const BundlesData& atlasBundleData,
                          const BundlesMinf& atlasBundleInfo,
                          std::vector<std::vector<float>>& fibersFeturesValues )
{

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;
  int nbPoints = atlasBundleData.pointsPerTrack[ 0 ] ;

  fibersFeturesValues.resize( nbFibersAtlasBundle ) ;


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
  float averageLength = atlasBundleInfo.averageLength ;
  float averageDistanceMedialPointFiberToAtlas = atlasBundleInfo.averageRadius ;
  float averageDistanceMedialPoints = 
                            atlasBundleInfo.averageDistanceBetweenMedialPoints ;
  float averageDistanceEndpoints1 = atlasBundleInfo.averageDisimilarity / 2 ;
  float averageDistanceEndpoints2 = atlasBundleInfo.averageDisimilarity / 2 ;
  float averageShapeAngle = atlasBundleInfo.averageShapeAngle ;
  float averageAngleBetweenPlanes = atlasBundleInfo.averageAngle ;
  float averageDirectionAngle = atlasBundleInfo.averageDirectionAngle ;
  float averageDMDA = atlasBundleInfo.averageDisimilarity ;

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

      distanceToOriginEndPoint1 += endPoint1[ tmpCoord ] * 
                                                        endPoint1[ tmpCoord ] ;
      distanceToOriginEndPoint2 += endPoint2[ tmpCoord ] * 
                                                        endPoint2[ tmpCoord ] ;

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
    std::array<float, 3> medialPointFiber1{0, 0, 0} ;
    atlasBundleData.computeMedialPointFiberWithDistance( fiberIndex1,
                                                         medialPointFiber1 ) ;
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


    // Normal vector
    std::array<float, 3> normalVectorFiber1{0, 0, 0} ;
    atlasBundleData.computeNormalVectorFiberTractogram(
                                                   atlasBundleData.matrixTracks,
                                                   medialPointFiber1,
                                                   fiberIndex1,
                                                   nbPoints,
                                                   normalVectorFiber1 ) ;
    std::array<float, 3> xVector = {1, 0, 0 } ;
    float angleToPlaneYZ = 
       atlasBundleData.computeAngleBetweenPlanes( xVector, normalVectorFiber1 ) ;

    

    // Direction vector
    std::array<float, 3> directionVectorFiber1{0, 0, 0} ;
    atlasBundleData.computeDirectionVectorFiberTractogram(
                                                   atlasBundleData.matrixTracks,
                                                   medialPointFiber1,
                                                   normalVectorFiber1,
                                                   fiberIndex1,
                                                   nbPoints,
                                                   directionVectorFiber1 ) ;
    float angleToDirectionX = 
       atlasBundleData.computeAngleBetweenDirections( xVector, directionVectorFiber1 ) ;

    
    // Number of feauters : 3 * nbPoints (coord fiber) + 3 (normalVector)
    // + 3 (direction vector) + 1 (shape angle) + 1 (length) 
    // + 1 (angle to direction X) + 1 angle to plane ZY
    std::vector<float> tmpFeaturesFiberAtlas( 3 * nbPoints + 10, -1.0 ) ;
    for ( int tmpPoint = 0 ; tmpPoint < nbPoints ; tmpPoint++ )
    {

      for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
      {

        tmpFeaturesFiberAtlas[ 3 * tmpPoint + tmpCoord ] = 
                                      fiber1Points[ 3 * tmpPoint + tmpCoord ] ;

      }

    }
    tmpFeaturesFiberAtlas[ 3 * nbPoints + 0 ] = lengthFiber ;
    tmpFeaturesFiberAtlas[ 3 * nbPoints + 1 ] = shapeAngle ;
    for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
    {

      tmpFeaturesFiberAtlas[ 3 * nbPoints + 2 + tmpCoord ] = 
                                               normalVectorFiber1[ tmpCoord ] ;

    }
    for ( int tmpCoord = 0 ; tmpCoord < 3 ; tmpCoord++ )
    {

      tmpFeaturesFiberAtlas[ 3 * nbPoints + 5 + tmpCoord ] = 
                                            directionVectorFiber1[ tmpCoord ] ;

    }

    tmpFeaturesFiberAtlas[ 3 * nbPoints + 8 + 0 ] = angleToDirectionX ;

    tmpFeaturesFiberAtlas[ 3 * nbPoints + 8 + 1 ] = angleToPlaneYZ ;

    
    fibersFeturesValues[ fiberIndex1 ] = tmpFeaturesFiberAtlas ;

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

  int index_atlas, index_output, index_format, index_nbThreads, 
                                                     index_verbose, index_help ;

  std::vector<std::string> possibleFlags{
                                        "-a", "-o", "-format", "-nbThreads", 
                                                                  "-v", "-h" } ;

  std::vector<bool> possibleFlagsNeedArgument{ true, true, true, true, true, 
                                                                       false } ;

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

  index_atlas = getFlagPosition( argc, argv, "-a" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_format = getFlagPosition( argc, argv, "-format") ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to get features of fibers in atlas : \n"
              << "-a : Directory with the atlas (one file per bundle), must "
              << "be given if other than ESBA atlas. You can also give a .tck "
              << "file of the atlas bundle (just one) you want to extract. \n"
              << "-o : Output directory to save the features .tsv files\n"
              << "-format : format of bundles in atlas among .bundles/.trk/"
              << ".tck \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads "
              << "(default : number of cores ) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
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

  if ( !index_format )
  {

    std::cout << "-fomat argument required ..." << std::endl ;
    exit( 1 ) ;

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

  std::string atlasDirectory = argv[ index_atlas + 1 ] ;
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
  for ( int bundleIndex = 0 ; bundleIndex < nbBundlesAtlas ; bundleIndex++ )
  {

    std::cout << "Processing : [" << bundleIndex + 1 << "/" << nbBundlesAtlas
                                                    << "]            " << "\r" ;

    BundlesMinf& atlasBundleMinf = atlasData.bundlesMinf[ bundleIndex ] ;
    BundlesData& atlasBundleData = atlasData.bundlesData[ bundleIndex ] ;
    std::string atlasBundleName = atlasData.bundlesNames[ bundleIndex ] ;

    std::vector<std::vector<float>> fibersFeturesValues ;
    computeFeaturesFibersBundle( atlasBundleData, atlasBundleMinf,
                                                         fibersFeturesValues ) ;


    std::stringstream outPathOss ;

    outPathOss << outputDirectory << atlasBundleName << ".tsv" ;
    std::string outPath = outPathOss.str() ;

    saveFeaturesFibersBundle( fibersFeturesValues, nbPoints, 
                                                            outPath.c_str() ) ;

  }


  std::cout << "\nDone" << std::endl ;

}
