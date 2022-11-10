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
#include <random>
#include <experimental/filesystem>

#include "computeNeighborhood.h"




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
//////////////////// Save index in tractogram of neighbors /////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void saveIndexInTractogram( const char* labelsBinaryFilename,
                                            const std::vector<int64_t>& labels )
{

  std::ofstream file ;
  file.open( labelsBinaryFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << labelsBinaryFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  int n_curves = labels.size() ;

  for ( int track = 0 ; track < n_curves ; track++ )
  {

    int64_t label = labels[ track ] ;

    file.write( reinterpret_cast<char*>( &label ), sizeof( int64_t ) ) ;

  }

  file.close() ;


}


////////////////////////////////////////////////////////////////////////////////
//////////////////// Read index in tractogram neighbors ////////////////////////
////////////////////////////////////////////////////////////////////////////////
void readIndexInTractogram( const char* predictedLabelsFilename,
                             std::vector<int64_t>& predictedLabels,
                             int nbFibers,
                             int verbose )
{

  predictedLabels.resize( nbFibers ) ;

  std::ifstream file ;
  file.open( predictedLabelsFilename, std::ios::binary | std::ios::in ) ;
  if ( file.fail() )
  {

    std::cout << "Problem reading file : " << predictedLabelsFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  for ( int fiber = 0 ; fiber < nbFibers ; fiber++ )
  {

    file.read( reinterpret_cast<char*>( &( predictedLabels[ fiber ] ) ),
                                                           sizeof( int64_t ) ) ;

  }

  file.close() ;

}


////////////////////////////////////////////////////////////////////////////////
/////////////////// Compute distance of a fiber to a bundle ////////////////////
////////////////////////////////////////////////////////////////////////////////
float computeDistanceFiberToBundle( BundlesDataFormat& atlasBundleData,
                                    const std::vector<float>& fiber )
{

  int nbPointsAtlasFibers = atlasBundleData.pointsPerTrack[ 0 ] ;
  int nbPointsFiber = ( int )( fiber.size() / 3 ) ;
  if ( nbPointsFiber != nbPointsAtlasFibers )
  {

    std::cout << "ERROR in computeNeighborhood::computeDistanceFiberToBundle : "
              << "the number of points in the fiber and in tha bundle fibers "
              << "must be the same "
              << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPointsFiber ;

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  // float distance = 0 ;
  std::vector<float> distances( nbFibersAtlasBundle, 0 ) ;

  #pragma omp parallel for
  for ( int atlasFiberIndex = 0 ; atlasFiberIndex < nbFibersAtlasBundle ;
                                                             atlasFiberIndex++ )
  {

    // distance += atlasBundleData.computeMDADBetweenTwoFibers(
    //                                                atlasBundleData.matrixTracks,
    //                                                fiber,
    //                                                atlasFiberIndex,
    //                                                0,
    //                                                nbPoints ) ;
    distances[ atlasFiberIndex ] =
                             atlasBundleData.computeMDFBetweenTwoFibers(
                                                   atlasBundleData.matrixTracks,
                                                   fiber,
                                                   atlasFiberIndex,
                                                   0,
                                                   nbPoints ) ;

  }

  int minIndex = std::min_element( distances.begin(), distances.end() ) -
                                                             distances.begin() ;

  // distance /= nbFibersAtlasBundle ;
  float distance = distances[ minIndex ] ;
  return distance ;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  const auto start_time = std::chrono::system_clock::now() ;

  int index_input, index_atlas, index_output, index_thr, index_minLen,
                          index_maxLen, index_tolThr, index_verbose, index_help ;
  index_input = getFlagPosition( argc, argv, "-i" ) ;
  index_atlas = getFlagPosition( argc, argv, "-a" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_thr = getFlagPosition( argc, argv, "-thr" ) ;
  index_minLen = getFlagPosition( argc, argv, "-minLen" ) ;
  index_maxLen = getFlagPosition( argc, argv, "-maxLen" ) ;
  index_tolThr = getFlagPosition( argc, argv, "-tolThr" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-i : Input tractogram from which to extract neighborhood to "
              << "atlas bundles \n"
              << "-a : Atlas directory \n"
              << "-o : Output directory \n"
              << "[-minLen] : Only use fibers longer than minLen for the "
              << "neighborhood (default = 10 ) \n"
              << "[-maxLen] : Only use fibers shorter than maxLen for the "
              << "neighborhood (default = atlasBundleMaxLength * 1.2) \n"
              << "[-tolThr] : multipicative value for distance threshold "
              << "(default = 1.2) \n"
              << "[-thr] : Threshold distance (Default = maxRadius) \n"
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

  if ( index_thr )
  {

    thresholdDistance = std::stof( argv[ index_thr + 1 ] ) ;

    if ( thresholdDistance <= 0 )
    {

      std::cout << "Error argument : thr must be greater than 0"
                << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_minLen )
  {

    minLength = std::stof( argv[ index_minLen + 1 ] ) ;

    if ( minLength <= 0 )
    {

      std::cout << "Error argument : minLen must be greater than 0"
                << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_maxLen )
  {

    maxLength = std::stof( argv[ index_maxLen + 1 ] ) ;

    if ( maxLength <= 0 )
    {

      std::cout << "Error argument : minLen must be greater than 0"
                << std::endl ;
      exit( 1 ) ;

    }

    useDefaultMaxLength = false ;

  }
  else
  {

    useDefaultMaxLength = true ;

  }

  if ( index_tolThr )
  {

    toleranceThr = std::stof( argv[ index_tolThr + 1 ] ) ;

    if ( toleranceThr <= 0 )
    {

      std::cout << "Error argument : -tolThr must be greater than 0"
                << std::endl ;
      exit( 1 ) ;

    }

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

  std::string inputTractogramPath( argv[ index_input + 1 ] ) ;
  char lastChar = inputTractogramPath[ inputTractogramPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputTractogramPath = inputTractogramPath.substr( 0,
                                              inputTractogramPath.size() - 1 ) ;

  }
  std::string inputBundlesDataFilename ;
  std::string inputBundlesFilename ;
  std::string inputTRKFilename ;

  if ( inputTractogramPath.find( ".bundlesdata" ) != std::string::npos )
  {

    inputBundlesDataFilename = inputTractogramPath ;

    inputBundlesFilename = inputTractogramPath ;
    size_t index = inputTractogramPath.find( ".bundlesdata" ) ;
    inputBundlesFilename.replace( index, 12, ".bundles") ;

    isBundlesFormat = true ;

  }
  else if ( inputTractogramPath.find( ".bundles" ) != std::string::npos )
  {

    inputBundlesFilename = inputTractogramPath ;

    inputBundlesDataFilename = inputTractogramPath ;
    size_t index = inputTractogramPath.find( ".bundles" ) ;
    inputBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

    isBundlesFormat = true ;

  }
  else if ( inputTractogramPath.find( ".trk" ) != std::string::npos )
  {

    inputTRKFilename = inputTractogramPath ;

    isTRKFormat = true ;

  }
  else
  {

    std::cout << "ERROR : Only output format supported is "
              << ".bundles/.bundlesdata " << std::endl ;
    exit( 1 ) ;

  }

  /// ********************************************************************** ///

  std::string atlasDirectory( argv[ index_atlas + 1 ] ) ;
  lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    atlasDirectory = atlasDirectory + '/' ;

  }

  /// ********************************************************************** ///

  std::string outputDirectory( argv[ index_output + 1 ] ) ;
  lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectory = outputDirectory + '/' ;

  }


  //////////////////////////////////////////////////////////////////////////////
  // isBundlesFormat = true ;
  // isTRKFormat = false ;

  //   xxxxxxxxxxxxxxxxxxxxx Reading inputTractogram xxxxxxxxxxxxxxxxxxxxxx   //

  BundlesFormat inputTractogramInfo ;
  BundlesDataFormat inputTractogram ;
  TrkFormat trkData ;


  if ( isBundlesFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << inputBundlesFilename << std::endl ;

    }

    inputTractogramInfo.bundlesReading(
                                       inputBundlesFilename.c_str(), verbose ) ;


    if ( verbose )
    {

      std::cout << "Reading : " << inputBundlesDataFilename << std::endl ;

    }


    inputTractogram.bundlesdataReading( inputBundlesDataFilename.c_str(),
                                       inputBundlesFilename.c_str(),
                                       verbose ) ;

  }
  else if ( isTRKFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << inputTRKFilename << std::endl ;

    }

    trkData.trkReading( inputTRKFilename.c_str(), verbose ) ;
    inputTractogram.matrixTracks = trkData.matrixTracks ;
    inputTractogram.pointsPerTrack = trkData.pointsPerTrack ;
    inputTractogram.curves_count = trkData.curves_count ;

  }
  else
  {

    std::cout << "ERROR : The only supported format are .bundles/.bundlesdata "
              << "and .trk" << std::endl ;
    exit( 1 ) ;

  }

  //   xxxxxxxxxxxxxxxxxxxxxx Reading atlas Bundle xxxxxxxxxxxxxxxxxxxxxxxx   //
  bool isBundlesFormatAtlas = false ;
  bool isTRKFormatAtlas = false ;
  if ( atlasDirectory.find( "-trk" ) == std::string::npos )
  {

    isBundlesFormatAtlas = true ;
    isTRKFormatAtlas = false ;

  }
  else
  {

    isBundlesFormatAtlas = false ;
    isTRKFormatAtlas = true ;

  }
  AtlasBundles atlasBundles( atlasDirectory.c_str(),
                             isBundlesFormatAtlas,
                             isTRKFormatAtlas,
                             verbose ) ;

  // BundlesFormat atlasBundleInfo( atlasBundle.c_str(), verbose ) ;


  //////////////////////////////////////////////////////////////////////////////
  if ( verbose )
  {

    if ( !isBundlesFormatAtlas || index_thr )
    {

      std::cout << "Threshold distance : " << thresholdDistance << " mm \n" ;

    }
    else
    {

      std::cout << "Threshold distance : maxRadius" << std::endl ;

    }

    std::cout << "Min length accepted for fibers : " << minLength << "\n" ;

    if ( useDefaultMaxLength )
    {

      std::cout << "Max length accepted for fibers : atlasBundleMaxLength * "
                                                        "tolThr"  << std::endl ;

    }
    else
    {

      std::cout << "Max length accepted for fibers : " << maxLength << "\n" ;

    }


  }

  int nbBundlesAtlas = atlasBundles.bundles.size() ;

  if ( verbose )
  {

    std::cout << "Processing bundles..." << std::endl ;

  }

  #pragma omp parallel for
  for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
  {

    /*
    if ( verbose < 2 )
    {

      printf( "Processing bundle [ %d / %d ] \r", bundle + 1, nbBundlesAtlas ) ;
      fflush( stdout ) ;

    }
    */



    BundlesFormat& atlasBundleInfo = atlasBundles.bundles[ bundle ] ;
    BundlesDataFormat& atlasBundleData = atlasBundles.bundlesData[ bundle ] ;

    std::vector<float> medialPointAtlasBundle = atlasBundleInfo.centerBundle ;


    float thresholdDistanceBundle = thresholdDistance * toleranceThr ;
    if ( !index_thr )
    {


      thresholdDistanceBundle = atlasBundleInfo.maxRadius * toleranceThr ;
      // float thresholdDistanceBundle = atlasBundleInfo.averageRadius ;
      if ( thresholdDistanceBundle == 0 )
      {

        std::cout << "WARNING : maxRadio equals to 0 for bundle "
                  << atlasBundles.bundlesNames[ bundle ] << "... Using default "
                  << "threshold of 10 mm " << std::endl ;
        thresholdDistanceBundle = 10 ;

      }

    }

    float maxLengthBundle = maxLength ;
    if ( useDefaultMaxLength )
    {

      maxLengthBundle = atlasBundleInfo.maxLength * 1.2 ;

    }


    if ( verbose > 1 )
    {

      #pragma omp critical
      {

        std::cout << atlasBundles.bundlesNames[ bundle ] << " threshold "
                  << "distance : " << thresholdDistanceBundle << std::endl ;

        std::cout << atlasBundles.bundlesNames[ bundle ] << " max length "
                  << "allowed : " << maxLengthBundle << std::endl ;

      }

    }

    int64_t nbFibersTractogram = inputTractogram.curves_count ;

    int nbPoints = inputTractogram.pointsPerTrack[ 0 ] ;
    // std::vector to create large arrays
    std::vector< int > indexNeighborFibers( nbFibersTractogram, 0 ) ;
    int curveCountNeighborhood = 0 ;
    int64_t nbElementsExtractedNeighborhood = 0 ;

    bool isOK = false ;

    while ( !isOK )
    {

      curveCountNeighborhood = 0 ;
      nbElementsExtractedNeighborhood = 0 ;


      for ( int fiberIndex = 0 ; fiberIndex < nbFibersTractogram ; fiberIndex++ )
      {

        indexNeighborFibers[ fiberIndex ] = 0 ;

        int nbPointsFiber = inputTractogram.pointsPerTrack[ fiberIndex ] ;

        int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;

        if ( nbPointsFiber != nbPoints )
        {

          std::cout << "ERROR : The number of points of each fiber in the input "
                    << "tractogram has to be the same, got " << nbPoints
                    << " for fiber 0 and " << nbPointsFiber << " for fiber "
                    << fiberIndex << std::endl ;
          exit( 1 ) ;

        }

        std::vector<float> medialPointFiberTractogram( 3, 0 ) ;

        // inputTractogram.computeMedialPointFiberTractogram(
        //                                           fiberIndex,
        //                                           medialPointFiberTractogram ) ;
        inputTractogram.computeMedialPointFiberWithDistance(
                                                  fiberIndex,
                                                  medialPointFiberTractogram ) ;

        float distance = 0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {

          distance += pow( medialPointFiberTractogram[ i ] -
                                              medialPointAtlasBundle[ i ], 2 ) ;

        }
        distance = sqrt( distance ) ;

        float lengthFiber = inputTractogram.computeLengthFiber( fiberIndex ) ;

        if ( distance < thresholdDistanceBundle && lengthFiber > minLength &&
                                                lengthFiber < maxLengthBundle  )
        {

          indexNeighborFibers[ fiberIndex ] = 1 ;

          nbElementsExtractedNeighborhood += 3 * nbPointsFiber ;
          curveCountNeighborhood += 1 ;

        }

      }


      if ( curveCountNeighborhood > minNbCurvesNeighborhood )
      {

        isOK = true ;

      }
      else
      {

        thresholdDistanceBundle += 5 ;

      }

    }


    if ( verbose > 1 )
    {

      #pragma omp critical
      {

        std::cout << "Number of fibers in neighborhood for bundle "
                  << atlasBundles.bundlesNames[ bundle ] << " : "
                  << curveCountNeighborhood << std::endl ;

      }

    }


    std::vector<std::vector<float>> extractedNeighborhood ;
    std::vector<float> extractedNeighborhoodTotal(
                                          nbElementsExtractedNeighborhood, 0 ) ;
    std::vector<int32_t> pointsPerTrackNeighborhood ;
    int32_t fiberIndexNeighborhood = 0 ;

    std::vector<int64_t> indexInTractogramOfNeighbors ;

    for ( int fiberIndex = 0 ; fiberIndex < nbFibersTractogram ; fiberIndex++ )
    {

      if ( indexNeighborFibers[ fiberIndex ] == 1 )
      {

        int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;
        int offsetExtractedNeighborhood = 3 * nbPoints * fiberIndexNeighborhood ;

        pointsPerTrackNeighborhood.push_back( inputTractogram.pointsPerTrack[
                                                                fiberIndex ] ) ;

        std::copy( inputTractogram.matrixTracks.begin() + offsetTractogram,
                   inputTractogram.matrixTracks.begin() + offsetTractogram +
                   3 * nbPoints, extractedNeighborhoodTotal.begin() +
                                                 offsetExtractedNeighborhood ) ;

        indexInTractogramOfNeighbors.push_back( ( int64_t )fiberIndex ) ;

        fiberIndexNeighborhood += 1 ;

      }

    }




    if ( isBundlesFormat )
    {




      //------------------------- Saving neighborhood ------------------------//

      BundlesFormat extractedNeighborhoodInfo( inputBundlesFilename.c_str(),
                                               verbose ) ;

      std::string outputBundlesFilename = outputDirectory +
                                          atlasBundles.bundlesNames[ bundle ] +
                                          ".bundles" ;

      if ( verbose > 1 )
      {

        std::cout << "Saving " << outputBundlesFilename << std::endl ;

      }

      extractedNeighborhoodInfo.curves_count = curveCountNeighborhood ;
      // extractedNeighborhoodInfo.curves_count = atlasBundleData.curves_count ;
      // extractedNeighborhoodInfo.curves_count = 10 ;

      extractedNeighborhoodInfo.bundlesWriting( outputBundlesFilename.c_str(),
                                                                     verbose ) ;



      BundlesDataFormat extractedNeighborhoodData( extractedNeighborhoodTotal,
                                                   pointsPerTrackNeighborhood,
                                                   curveCountNeighborhood ) ;


      std::string outputBundlesDataFilename =
                                           outputDirectory +
                                           atlasBundles.bundlesNames[ bundle ] +
                                           ".bundlesdata" ;

      // Saving index in tractogram of neighbors fibers
      std::string indexInTractogramOfNeighborsFilename = outputDirectory +
                                          atlasBundles.bundlesNames[ bundle ] +
                                          "Index.bin" ;

      saveIndexInTractogram( indexInTractogramOfNeighborsFilename.c_str(),
                                                indexInTractogramOfNeighbors ) ;


      //
      std::vector<int64_t> indexInTractogramRead ;
      int32_t testNbFibers = indexInTractogramOfNeighbors.size() ;
      readIndexInTractogram( indexInTractogramOfNeighborsFilename.c_str(),
                             indexInTractogramRead,
                             testNbFibers,
                             verbose ) ;


      // std::cout << "Test read index : " << std::endl ;
      // for ( int i = 0 ; i < testNbFibers ; i++ )
      // {
      //
      //   std::cout << indexInTractogramRead[ i ] << "\t|\t" <<
      //                           indexInTractogramOfNeighbors[ i ] << std::endl ;
      //
      // }
      // std::cout << "Done" << std::endl ;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
      // // Selecting maxNbFibers fibers in moving bundle that resemble the most to
      // // the reference bundle
      // // int maxNbFibers = atlasBundleData.curves_count  ;
      // int maxNbFibers = 10  ;
      //
      // std::vector<float> selectedMovingFibers( 3 * nbPoints * maxNbFibers, 0 ) ;
      // std::vector<int> pointsPerTrackSelectedFibers( maxNbFibers, 0 ) ;
      //
      //
      // std::vector<float> distancesMovingFibersToReference( curveCountNeighborhood, 0 ) ;
      // #pragma omp parallel for
      // for ( int fiberIndex = 0 ; fiberIndex < curveCountNeighborhood ; fiberIndex++ )
      // {
      //
      //   std::vector<float> fiber( 3 * nbPoints, 0 ) ;
      //   extractedNeighborhoodData.getFiberFromTractogram(
      //                                           extractedNeighborhoodData.matrixTracks,
      //                                           fiberIndex,
      //                                           nbPoints,
      //                                           fiber ) ;
      //   distancesMovingFibersToReference[ fiberIndex ] =
      //                computeDistanceFiberToBundle( atlasBundleData, fiber ) ;
      //
      // }
      //
      // std::vector<int> indices( curveCountNeighborhood, 0 ) ;
      // std::iota( indices.begin(), indices.end(), 0 ) ;
      // std::partial_sort( indices.begin(), indices.begin() + maxNbFibers,
      //                    indices.end(),
      //                    [&distancesMovingFibersToReference]( int i, int j )
      //                    { return distancesMovingFibersToReference[ i ] >
      //                             distancesMovingFibersToReference[ j ] ; } ) ;
      //
      //
      // for ( int fiberSelectedIndex = 0 ; fiberSelectedIndex < maxNbFibers ;
      //                                                       fiberSelectedIndex++ )
      // {
      //
      //   int offsetExtractedNeighborhood = 3 * nbPoints * fiberSelectedIndex ;
      //
      //   // Find minimum element index
      //   int minIndex = indices[ fiberSelectedIndex ] ;
      //
      //   std::vector<float> fiber( 3 * nbPoints, 0 ) ;
      //   // extractedNeighborhoodData.getFiberFromTractogram(
      //   //                                         extractedNeighborhoodData.matrixTracks,
      //   //                                         fiberSelectedIndex,
      //   //                                         nbPoints,
      //   //                                         fiber ) ;
      //   extractedNeighborhoodData.getFiberFromTractogram(
      //                                           extractedNeighborhoodData.matrixTracks,
      //                                           minIndex,
      //                                           nbPoints,
      //                                           fiber ) ;
      //
      //   std::copy( fiber.begin(),
      //              fiber.begin() + 3 * nbPoints,
      //              selectedMovingFibers.begin() + offsetExtractedNeighborhood ) ;
      //
      //   pointsPerTrackSelectedFibers[ fiberSelectedIndex ] =
      //                               extractedNeighborhoodData.pointsPerTrack[ minIndex ] ;
      //
      // }
      //
      // extractedNeighborhoodData.matrixTracks = std::vector<float>() ; // deallocation
      // extractedNeighborhoodData.matrixTracks = selectedMovingFibers ;
      //
      // extractedNeighborhoodData.pointsPerTrack = std::vector<int>() ; // deallocation
      // extractedNeighborhoodData.pointsPerTrack = pointsPerTrackSelectedFibers ;
      //
      // extractedNeighborhoodData.curves_count = maxNbFibers ;
      //
      //
      // if ( verbose )
      // {
      //
      //   std::cout << " Done " << std::endl ;
      //
      // }
      //
      // float lengthTest = extractedNeighborhoodData.computeLengthFiber( 0 ) ;
      // std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" ;
      // std::cout <<  lengthTest << std::endl ;
      // std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" ;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

      if ( verbose > 1 )
      {

        std::cout << "Saving " << outputBundlesDataFilename << std::endl ;

      }

      extractedNeighborhoodData.bundlesdataWriting(
                                             outputBundlesDataFilename .c_str(),
                                             verbose ) ;




    }
    else
    {

      std::cout << "ERROR : Only output files supported are .bundles/"
                << ".bundlesdata " << std::endl ;
      exit( 1 ) ;

    }

  }

  std::cout << "\n" ;

  //
  const std::chrono::duration< double > duration =
                              std::chrono::system_clock::now() - start_time ;

  if ( verbose )
  {

    std::cout << "Duration : " << duration.count() << std::endl ;

  }

  return( 0 ) ;

}
