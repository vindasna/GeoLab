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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include "applyPredictedLabels.h"
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
/////////////// Function to read the predicted labels file ////////////////////
///////////////////////////////////////////////////////////////////////////////

void readingPredictedLabels( const char* predictedLabelsFilename,
                             std::vector<std::vector<int16_t>>& predictedLabels,
                             int nbFibers,
                             int verbose )
{

  if ( verbose )
  {

    std::cout << "Reading labels : " << predictedLabelsFilename << std::endl ;

  }

  predictedLabels.resize( nbFibers ) ;

  const char delim = ':' ;
  std::string line ;
  std::ifstream dictFile ;
  dictFile.open( predictedLabelsFilename ) ;
  if ( dictFile.fail() )
  {

    std::cout << "Problem reading file : " << predictedLabelsFilename
                                           << std::endl ;
    exit( 1 ) ;

  }
  while ( std::getline( dictFile, line ) )
  {

    std::vector< std::string > out ;
    std::stringstream ss( line ) ;
    std::string s ;
    while ( std::getline( ss, s, delim ) )
    {

      s.erase( std::remove( s.begin(), s.end(), ' ' ), s.end() ) ;
      out.push_back( s ) ;

    }

    predictedLabels[ stoi( out[ 0 ] ) ].push_back( stoi( out[ 1 ] ) ) ;

  }

  dictFile.close() ;

}
////////////////////////////////////////////////////////////////////////////////
////// Function to read the predicted labels binary for SupWMA file ////////////
////////////////////////////////////////////////////////////////////////////////

void readingPredictedLabelsSupWMA(
                             const char* predictedLabelsFilename,
                             std::vector<std::vector<int16_t>>& predictedLabels,
                             int nbFibers,
                             int verbose )
{

  if ( verbose )
  {

    std::cout << "Reading SupWMA labels : " << predictedLabelsFilename
                                                                  << std::endl ;

  }

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

    int16_t _label ;

    file.read( reinterpret_cast<char*>( &( _label ) ), sizeof( int16_t ) ) ;
    predictedLabels[ fiber ].push_back( _label ) ;

  }

  file.close() ;

}


///////////////////////////////////////////////////////////////////////////////
/////////////// Function to read the dictionary labels file ///////////////////
///////////////////////////////////////////////////////////////////////////////

void readingDictionaryLabels( const char* dictionaryLabelsFilename,
                              std::vector<std::string>& dictionaryLabels,
                              int verbose )
{

  const char delim = ':' ;
  std::string line ;
  std::ifstream dictFile ;
  dictFile.open( dictionaryLabelsFilename ) ;
  if ( dictFile.fail() )
  {

    std::cout << "Problem reading file : " << dictionaryLabelsFilename
                                           << std::endl ;
    exit( 1 ) ;

  }
  while ( std::getline( dictFile, line ) )
  {

    std::vector< std::string > out ;
    std::stringstream ss( line ) ;
    std::string s ;
    while ( std::getline( ss, s, delim ) )
    {

      s.erase( std::remove( s.begin(), s.end(), ' ' ), s.end() ) ;
      out.push_back( s ) ;

    }

    dictionaryLabels.push_back( out[ 0 ] ) ;

  }

  dictFile.close() ;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

  int index_input, index_labels, index_output, index_dict, index_save_unlabeled,
                                     index_isSupWMA, index_verbose, index_help ;

  std::vector<std::string> possibleFlags{ "-i", "-l", "-o", "-d", "-su",
                                                       "-supWMA", "-v", "-h" } ;
  std::vector<bool> possibleFlagsNeedArgument{ true, true, true, true, false,
                                                          false, true, false } ;


  for ( int k = 1 ; k < argc ; k++ )
  {

    std::string arg = argv[ k ] ;
    auto it = std::find( possibleFlags.begin(), possibleFlags.end(), arg ) ;
    if( it == possibleFlags.end() )
    {

      if ( k == 1 )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for applyPredictedLabels command" << std::endl ;
        exit( 1 ) ;

      }

      std::string arg2 = argv[ k - 1 ] ;
      it = std::find( possibleFlags.begin(), possibleFlags.end(), arg2 ) ;
      if( it == possibleFlags.end() )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for applyPredictedLabels command" << std::endl ;
        exit( 1 ) ;

      }
      int indexFound = it - possibleFlags.begin() ;
      if ( !possibleFlagsNeedArgument[ indexFound ] )
      {

        std::cout << "ERROR in arguments : " << arg << " is not an argument "
                  << "for applyPredictedLabels command" << std::endl ;
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
  index_labels = getFlagPosition( argc, argv, "-l" ) ;
  index_output = getFlagPosition( argc, argv, "-o" ) ;
  index_dict = getFlagPosition( argc, argv, "-d" ) ;
  index_save_unlabeled = getFlagPosition( argc, argv, "-su" ) ;
  index_isSupWMA = getFlagPosition( argc, argv, "-supWMA" ) ;
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : \n"
              << "-i : Input tractogram from which to extract bundles \n"
              << "-l : Predicted labels (.txt) \n"
              << "-o : Output directory to save the extracted bundles \n"
              << "[-d] : Dictionary containing the names associated to each "
              << "label (.dict) \n"
              << "[-su] : Save bundle containing unlabeled fibers \n"
              << "[-supWMA] : Use this flag is SupWMA .bin predicted labels \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_input )
  {

    std::cout << "-i argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_labels )
  {

    std::cout << "-l argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  if ( index_save_unlabeled )
  {

    saveUnlabeled = true ;

  }
  else
  {

    saveUnlabeled = false ;

  }

  if ( index_isSupWMA )
  {

    isSupWMA = true ;

  }
  else
  {

    isSupWMA = false ;

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

  std::string predictedLabelsFilename( argv[ index_labels + 1 ] ) ;
  lastChar = predictedLabelsFilename[ predictedLabelsFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    predictedLabelsFilename = predictedLabelsFilename.substr( 0,
                                       predictedLabelsFilename.size() - 1 ) ;

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

    readingDictionaryLabels( dictionaryLabelsFilename.c_str(),
                             dictionaryLabels,
                             verbose ) ;


  }


  //   xxxxxxxxxxxxxxxxxxxxx Reading Input tractogram xxxxxxxxxxxxxxxxxxxxx   //

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


  //   xxxxxxxxxxxxxxxxxxxxxxxxxx Reading label xxxxxxxxxxxxxxxxxxxxxxxxxx   //
  if ( verbose )
  {

    std::cout << "Reading predicted labels : " << predictedLabelsFilename
                                               << std::endl ;
  }

  if ( isSupWMA )
  {

    readingPredictedLabelsSupWMA( predictedLabelsFilename.c_str(),
                                  predictedLabels,
                                  inputTractogram.curves_count,
                                  verbose ) ;



  }
  else
  {

    readingPredictedLabels( predictedLabelsFilename.c_str(),
                            predictedLabels,
                            inputTractogram.curves_count,
                            verbose ) ;

  }


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

  int nbLabels = predictedLabels.size() ;

  if ( nbLabels != nbFibersTractogram )
  {

    std::cout << "ERROR labels : the number of predicted labels has to be the "
              << "same as the number of fibers in the tractogram, got "
              << nbLabels << " labels and " << nbFibersTractogram
              << " fibers" << std::endl ;
    exit( 1 ) ;

  }

  if ( verbose )
  {

    std::cout << "Minimum number of fibers per bundle : " <<
                                              minimumNumberFibers << std::endl ;

  }




  //////////////////////////////// Applt labels ////////////////////////////////
  if ( verbose )
  {

    std::cout << "Labeling fibers in tractogram... " << std::endl ;

  }

  // Getting the number of bundles
  int16_t nbBundles = 0 ;

  if ( index_dict )
  {

    nbBundles = dictionaryLabels.size() ;

    if ( verbose )
    {

      std::cout << "Number of clusters according to dictionary with labels : "
                << nbBundles << std::endl ;

    }

    // The following commented lines where used when predictedLabels was a
    // std::vector<int64_t> where each fiber had just one label
    // int16_t nbRecognizedBundles = *std::max_element( predictedLabels.begin(),
    //                                                predictedLabels.end() ) + 1 ;
    int nbRecognizedBundles = -1 ;
    for ( int fiberIndex = 0 ; fiberIndex < nbFibersTractogram ; fiberIndex++ )
    {

      int nbLabelsFiber = predictedLabels[ fiberIndex ].size() ;
      for ( int _label = 0 ; _label < nbLabelsFiber ; _label++ )
      {

        if ( nbRecognizedBundles < predictedLabels[ fiberIndex ][ _label ] )
        {

          nbRecognizedBundles = predictedLabels[ fiberIndex ][ _label ] ;

        }

      }

    }

    if ( nbRecognizedBundles > nbBundles )
    {

      std::cout << "ERROR : the number of recognize bundles according to the "
                << "binary is superior to the number of bundles in the "
                << "dictionary file " << std::endl ;
      exit( 1 ) ;

    }

  }
  else
  {

    nbBundles = -1 ;
    for ( int fiberIndex = 0 ; fiberIndex < nbFibersTractogram ; fiberIndex++ )
    {

      int nbLabelsFiber = predictedLabels[ fiberIndex ].size() ;
      for ( int _label = 0 ; _label < nbLabelsFiber ; _label++ )
      {

        if ( nbBundles < predictedLabels[ fiberIndex ][ _label ] )
        {

          nbBundles = predictedLabels[ fiberIndex ][ _label ] ;

        }

      }

    }
    if ( nbBundles == -1 )
    {

      std::cout << "ERROR : could not get the number of bundles from the "
                << "labels file " << std::endl ;
      exit( 1 ) ;

    }
    else
    {

      nbBundles += 1 ; // Because in labels indexation starts at 0

    }

    if ( verbose )
    {

      std::cout << "Number of clusters according to labels files : "
                << nbBundles << std::endl ;

    }


  }


  if ( predictedLabels.size() != nbFibersTractogram )
  {

    std::cout << "ERROR : the number of fibers in predictedLabels is different "
              << "than the number of fibers in the input tractogram, got "
              << predictedLabels.size() << " and " << nbFibersTractogram
              << " respectively" << std::endl ;
    exit( 1 ) ;

  }


  std::vector< BundlesData > labeledFibers ;
  labeledFibers.resize( nbBundles ) ;
  int64_t nbUnlabelledFibers = 0 ;
  for ( int i = 0 ; i < nbFibersTractogram ; i++ )
  {

    int nbLabelsFiber = predictedLabels[ i ].size() ;
    for ( int _label = 0 ; _label < nbLabelsFiber ; _label++ )
    {

      int labelIndex = predictedLabels[ i ][ _label ] ;
      if ( labelIndex > -1 )
      {

        labeledFibers[ labelIndex ].curves_count += 1 ;

      }
      else
      {

        nbUnlabelledFibers += 1 ;

      }

    }

  }

  int nbBundlesRecognized = 0 ;
  for ( int bundle = 0 ; bundle < nbBundles ; bundle++ )
  {

    if ( labeledFibers[ bundle ].curves_count > minimumNumberFibers )
    {

      nbBundlesRecognized += 1 ;

    }

  }

  if ( verbose )
  {

    printf( "Number of bundles recognized : %d \n", nbBundlesRecognized ) ;


  }

  BundlesData unlabeledFibersData ;
  if ( nbUnlabelledFibers > 0 )
  {

    unlabeledFibersData.curves_count = nbUnlabelledFibers ;
    unlabeledFibersData.matrixTracks.resize( 3 * nbPoints * nbUnlabelledFibers,
                                                                           0 ) ;
    unlabeledFibersData.pointsPerTrack.resize( nbUnlabelledFibers, 0 ) ;

  }

  if ( nbBundlesRecognized )
  {

    std::vector<int64_t> offsetLabeledFibers( nbBundles, 0 ) ;
    std::vector<int64_t> offsetpointsPerTrackLabeledFiber( nbBundles, 0 ) ;
    int64_t offsetUnlablledFibers = 0 ;
    int64_t offsetPointsPerTrackUnlablledFibers = 0 ;

    for ( int i = 0 ; i < nbFibersTractogram ; i++ )
    {

      int64_t offsetTractogram = 3 * nbPoints * i  ;

      int nbLabelsFiber = predictedLabels[ i ].size() ;

      for ( int _label = 0 ; _label < nbLabelsFiber ; _label++ )
      {

        int labelIndex = predictedLabels[ i ][ _label ] ;
        if ( labelIndex < 0 )
        {

          std::copy( inputTractogram.matrixTracks.begin() + offsetTractogram,
                     inputTractogram.matrixTracks.begin() + offsetTractogram +
                                                                   3 * nbPoints,
                     unlabeledFibersData.matrixTracks.begin() +
                                                       offsetUnlablledFibers ) ;


          unlabeledFibersData.pointsPerTrack[
                                         offsetPointsPerTrackUnlablledFibers ] =
                                           inputTractogram.pointsPerTrack[ i ] ;

          offsetUnlablledFibers += 3 * nbPoints ;
          offsetPointsPerTrackUnlablledFibers += 1 ;

          continue ;

        }

        if ( labeledFibers[ labelIndex ].curves_count > 0 )
        {

          if ( offsetLabeledFibers[ labelIndex ] == 0 )
          {

            labeledFibers[ labelIndex ].matrixTracks.resize(
                  labeledFibers[ labelIndex ].curves_count * 3 * nbPoints, 0 ) ;
            labeledFibers[ labelIndex ].pointsPerTrack.resize(
                                 labeledFibers[ labelIndex ].curves_count, 0 ) ;

          }


          std::copy( inputTractogram.matrixTracks.begin() + offsetTractogram,
                     inputTractogram.matrixTracks.begin() + offsetTractogram +
                                                                   3 * nbPoints,
                     labeledFibers[ labelIndex ].matrixTracks.begin() +
                                           offsetLabeledFibers[ labelIndex ] ) ;

          labeledFibers[ labelIndex ].pointsPerTrack[
                              offsetpointsPerTrackLabeledFiber[ labelIndex ] ] =
                                           inputTractogram.pointsPerTrack[ i ] ;

          offsetLabeledFibers[ labelIndex ] += nbPoints * 3 ;

          offsetpointsPerTrackLabeledFiber[ labelIndex ] += 1 ;

        }

      }

    }

  }
  else
  {

    std::cout << "No bundles of the atlas found in the input tractogram \n" ;
    return 0 ;

  }


  lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectory = outputDirectory + "/" ;

  }

  if ( verbose )
  {

    std::cout << "Saving extracted bundles in : " << outputDirectory
                                                  << std::endl ;

  }

  int countProcessedRecognizedBundles = 0 ;
  // #pragma omp parallel for reduction( + : countProcessedRecognizedBundles )
  for ( int bundle = 0 ; bundle < nbBundles ; bundle++ )
  {

    if ( verbose && ( countProcessedRecognizedBundles % 1 == 0 ||
          ( countProcessedRecognizedBundles + 1 ) == nbBundlesRecognized ) &&
          labeledFibers[ bundle ].curves_count > minimumNumberFibers )
    {

      printf("\rProcessing fiber [ %10d / %10d ]",
                                        countProcessedRecognizedBundles + 1 ,
                                                        nbBundlesRecognized ) ;
      std::cout << "" << std::flush ;

    }

    if ( labeledFibers[ bundle ].curves_count > minimumNumberFibers )
    {

      std::string labelName = "" ;
      if ( index_dict )
      {

        labelName += dictionaryLabels[ bundle ] ;

      }
      else
      {

        labelName = "cluster_" + std::to_string( bundle ) ;

      }

      BundlesMinf outBundlesInfo = inputTractogramInfo ;

      outBundlesInfo.curves_count = labeledFibers[ bundle ].curves_count ;

      std::stringstream outBundlesFilenameOss ;
      outBundlesFilenameOss << outputDirectory << labelName << format ;
      std::string outBundlesFilename = outBundlesFilenameOss.str() ;

      if ( isBundlesFormat )
      {

        labeledFibers[ bundle ].isBundles = true ;
        labeledFibers[ bundle ].isTrk = false ;
        labeledFibers[ bundle ].isTck = false ;

      }
      else if ( isTrkFormat )
      {

        labeledFibers[ bundle ].isBundles = false ;
        labeledFibers[ bundle ].isTrk = true ;
        labeledFibers[ bundle ].isTck = false ;

      }
      else if ( isTckFormat )
      {

        labeledFibers[ bundle ].isBundles = false ;
        labeledFibers[ bundle ].isTrk = false ;
        labeledFibers[ bundle ].isTck = true ;

      }

      labeledFibers[ bundle ].write( outBundlesFilename.c_str(),
                                                              outBundlesInfo ) ;


      countProcessedRecognizedBundles += 1 ;

    }

  }

  if ( saveUnlabeled && nbUnlabelledFibers > 0 )
  {

    std::cout << "\nSaving unlabeled fibers..." << std::endl ;


    std::string labelName = "unlabeled" ;

    BundlesMinf outBundlesInfo = inputTractogramInfo ;

    outBundlesInfo.curves_count = unlabeledFibersData.curves_count ;


    std::stringstream outBundlesFilenameOss ;
    outBundlesFilenameOss << outputDirectory << labelName << format ;
    std::string outBundlesFilename = outBundlesFilenameOss.str() ;

    unlabeledFibersData.write( outBundlesFilename.c_str(), outBundlesInfo ) ;

    std::cout << "Done" << std::endl ;


  }


  if ( verbose )
  {

    std::cout << "\nDone" << std::endl ;

  }
  return 0 ;

}
