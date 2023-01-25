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

#include <random>

#include <boost/process.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "computeScoresFromLabels.h"
#include "ioWrapper.h"



////////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ////////////////
////////////////////////////////////////////////////////////////////////////////

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
////////////////////// Function to save confusion matrix ///////////////////////
////////////////////////////////////////////////////////////////////////////////
void saveConfusionMatrix(
                      const char* confusionMatrixPath,
                      const std::vector<std::vector<int32_t>>& confusionMatrix )
{

  std::cout << "Saving confusion matrix in : " << confusionMatrixPath
                                               << std::endl ;

  std::ofstream file ;
  file.open( confusionMatrixPath ) ;
  if ( file.fail() )
  {

    std::cout << "   Problem opening file to write : " << confusionMatrixPath <<
                                                                     std::endl ;

    exit( 1 ) ;

  }

  int sizeConfusionMatrix = confusionMatrix.size() ;

  std::cout << "   Size confusion matrix : " << sizeConfusionMatrix
                                                                  << std::endl ;


  for ( int _line_i = 0 ; _line_i < sizeConfusionMatrix ; _line_i++ )
  {

    for ( int _column_i = 0 ; _column_i < sizeConfusionMatrix ; _column_i++ )
    {

      file << confusionMatrix[ _line_i ][ _column_i ] << "\t" ;

    }

    file << "\n" ;

  }

  file.close() ;

  std::cout << "Done" << std::endl ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////// Function to save scores per bundle ///////////////////////
////////////////////////////////////////////////////////////////////////////////
void saveScorePerBundle( const char* scoresPerBundlePath,
                         const std::vector<float>& precisionPerBundle,
                         const std::vector<float>& sensitivityPerBundle,
                         const std::vector<float>& accuracyPerBundle,
                         const std::vector<float>& jaccardPerBundle,
                         const std::vector<std::string>& bundleNames,
                         float averagePrecision,
                         float averageSensitivity,
                         float averageAccuracy,
                         float averageJaccard )
{

  std::cout << "Saving scores per bundle in : " << scoresPerBundlePath
                                                    << std::endl ;
  std::ofstream file ;
  file.open( scoresPerBundlePath ) ;
  if ( file.fail() )
  {

    std::cout << "   Problem opening file to write : " << scoresPerBundlePath <<
                                                                   std::endl ;

    exit( 1 ) ;

  }

  file << "Bundle\tPrecision\tSensitivity\tAccuracy\tJaccard\n" ;

  int nbBundles = precisionPerBundle.size() - 1 ;
  for ( int _bundleIndex = 0 ; _bundleIndex < nbBundles + 1 ; _bundleIndex++ )
  {

    if ( _bundleIndex == nbBundles )
    {

      file << "Unlabeled\t" ;

    }
    else
    {

      file << bundleNames[ _bundleIndex ] << "\t" ;

    }
    file << precisionPerBundle[ _bundleIndex ] << "\t" ;
    file << sensitivityPerBundle[ _bundleIndex ] << "\t" ;
    file << accuracyPerBundle[ _bundleIndex ] << "\t" ;
    file << jaccardPerBundle[ _bundleIndex ] << "\t" ;

    file << "\n" ;

  }

  file << "Average\t" << averagePrecision << "\t" << averageSensitivity << "\t"
                      << averageAccuracy << "\t" << averageJaccard << "\n" ;

  file.close() ;

  std::cout << "Done" << std::endl ;

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  const auto start_time = std::chrono::system_clock::now() ;

  int index_pl, index_pd, index_tl, index_td, index_st, index_output,
                       index_force, index_nbThreads, index_verbose, index_help ;

  index_pl = getFlagPosition( argc, argv, "-pl") ;
  index_pd = getFlagPosition( argc, argv, "-pd") ;
  index_tl = getFlagPosition( argc, argv, "-tl") ;
  index_td = getFlagPosition( argc, argv, "-td") ;
  index_st = getFlagPosition( argc, argv, "-st") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_force = getFlagPosition( argc, argv, "-force") ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-pl : Path to the predicted labels .txt \n"
              << "-pd : Path to the predicted labels dictionary .dict \n"
              << "-tl : Path to the true labels .txt \n"
              << "-td : Path to the true labels dictionary .dict \n"
              << "-st : Path to the input subject's tractogram of the "
              << "ProjectAtlasGeoLab command \n"
              << "-o : Path to the output directory \n"
              << "[-force] : Force to overwrite files (default = false) \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads before "
              << "applyRecoBundles (default : let openMP decide ) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_pl )
  {

    std::cout << "-pl argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_pd )
  {

    std::cout << "-pd argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_tl )
  {

    std::cout << "-tl argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_td )
  {

    std::cout << "-td argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_st )
  {

    std::cout << "-st argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }



  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Checking arguments /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  /////////////////////////// Predicted labels path ////////////////////////////
  predictedLabelsPath = argv[ index_pl + 1 ] ;
  char lastChar = predictedLabelsPath[ predictedLabelsPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    predictedLabelsPath = predictedLabelsPath.substr( 0,
                                              predictedLabelsPath.size() - 1 ) ;

  }
  if ( !is_file( predictedLabelsPath ) )
  {

    std::cout << "ERROR : Predicted labels " << predictedLabelsPath << " does "
                                                 << "not exists " << std::endl ;
    exit( 1 );

  }

  //////////////////////////// Predicted dict path /////////////////////////////
  predictedDictPath = argv[ index_pd + 1 ] ;
  lastChar = predictedDictPath[ predictedDictPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    predictedDictPath = predictedDictPath.substr( 0,
                                              predictedDictPath.size() - 1 ) ;

  }
  if ( !is_file( predictedDictPath ) )
  {

    std::cout << "ERROR : Predicted dict " << predictedDictPath << " does not "
                                                     << "exists " << std::endl ;
    exit( 1 );

  }

  ////////////////////////////// True labels path //////////////////////////////
  trueLabelsPath = argv[ index_tl + 1 ] ;
  lastChar = trueLabelsPath[ trueLabelsPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    trueLabelsPath = trueLabelsPath.substr( 0, trueLabelsPath.size() - 1 ) ;

  }
  if ( !is_file( trueLabelsPath ) )
  {

    std::cout << "ERROR : True labels " << trueLabelsPath << " does not exists"
                                                                  << std::endl ;
    exit( 1 );

  }

  /////////////////////////////// True dict path ///////////////////////////////
  trueDictPath = argv[ index_td + 1 ] ;
  lastChar = trueDictPath[ trueDictPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    trueDictPath = trueDictPath.substr( 0, trueDictPath.size() - 1 ) ;

  }
  if ( !is_file( trueDictPath ) )
  {

    std::cout << "ERROR : True dict " << trueDictPath << " does not exists"
                                                                  << std::endl ;
    exit( 1 );

  }


  ////////////////////////// Subject's tracogram path //////////////////////////
  subjectTractogramPath = argv[ index_st + 1 ] ;
  lastChar = subjectTractogramPath[ subjectTractogramPath.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    subjectTractogramPath = subjectTractogramPath.substr( 0,
                                            subjectTractogramPath.size() - 1 ) ;

  }
  if ( !is_file( subjectTractogramPath ) )
  {

    std::cout << "ERROR : Subject's tractogram " << subjectTractogramPath
                                            << " does not exists" << std::endl ;
    exit( 1 );

  }

  if ( !endswith( subjectTractogramPath, ".bundles" ) &&
                            !endswith( subjectTractogramPath, ".bundesdata" ) &&
                                   !endswith( subjectTractogramPath, ".trk" ) &&
                                    !endswith( subjectTractogramPath, ".tck" ) )
  {

    std::cout << "ERROR : subject's tractogram must be .bundles/.trk/.tck"
              << std::endl ;
    exit( 1 ) ;

  }


  /////////////////////////////////// Force ////////////////////////////////////
  // Must be done before outputPath because if !force and
  // countFilesDirectory( outputPath ) > 0 then exit command
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


  ////////////////////////////// Output directory //////////////////////////////
  outputDirectory = argv[ index_output + 1 ]  ;
  lastChar = outputDirectory[ outputDirectory.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectory = outputDirectory + "/" ;

  }



  if ( !is_dir( outputDirectory ) )
  {

    mkdir( outputDirectory ) ;

  }



  ///////////////////////////////// nbThreads //////////////////////////////////
  if ( index_nbThreads )
  {

    nbThreads = std::stoi( argv[ index_nbThreads + 1 ] ) ;
    if ( nbThreads < 0 )
    {

      std::cout << "Invalid argument for -nbThreads : you must give a postive "
                << "integer " << std::endl ;
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
  BundlesMinf subjectTractogramInfo( subjectTractogramPath.c_str() ) ;

  int nbFibers = subjectTractogramInfo.curves_count ;

  std::vector<std::vector<std::string>> predictedLabelsByName ;
  std::cout << "Reading predicted labels... " ;
  readLabelsWithDict( predictedDictPath.c_str(), predictedLabelsPath.c_str(),
                                             predictedLabelsByName, nbFibers ) ;
  std::cout << "Done" << std::endl ;

  std::cout << "Reading predicted dictionary... " ;
  std::vector<std::string> predictedDict ;
  readLabelsDict( predictedDictPath.c_str(), predictedDict ) ;
  int nbPredictedLabels = predictedDict.size() ;
  std::cout << "Done" << std::endl ;

  std::cout << "Reading true labels... " ;
  std::vector<std::vector<std::string>> trueLabelsByName ;
  readLabelsWithDict( trueDictPath.c_str(), trueLabelsPath.c_str(),
                                                  trueLabelsByName, nbFibers ) ;
  std::cout << "Done" << std::endl ;

  std::cout << "Reading true dict... " ;
  std::vector<std::string> trueDict ;
  readLabelsDict( trueDictPath.c_str(), trueDict ) ;
  int nbTrueLabels = trueDict.size() ;
  std::cout << "Done" << std::endl ;

  std::vector<std::string>& bundlesNames = trueDict ;
  int nbLabels = nbTrueLabels ;
  // int nbLabels = -1 ;
  // if ( nbTrueLabels > nbPredictedLabels )
  // {
  //
  //   nbLabels = nbTrueLabels ;
  //   bundlesNames = trueDict ;
  //
  // }
  // else
  // {
  //
  //   nbLabels = nbPredictedLabels ;
  //   bundlesNames = predictedDict ;
  //
  // }


  std::cout << "Computing confusion matrix... " ;
  std::vector<std::vector<int32_t>> confusionMatrix( nbLabels + 1,
                                     std::vector<int32_t>( nbLabels + 1, 0 ) ) ;
  // Race condition
  // #pragma omp parallel for
  for ( int fiberIndex = 0 ; fiberIndex < nbFibers ; fiberIndex++ )
  {

    std::vector<std::string>& predictedLabels =
                                           predictedLabelsByName[ fiberIndex ] ;
    std::vector<std::string>& trueLabels = trueLabelsByName[ fiberIndex ] ;

    for ( int predictedLabelIndex = 0 ;
                                predictedLabelIndex < predictedLabels.size() ;
                                                       predictedLabelIndex++ )
    {

      for ( int trueLabelIndex = 0 ; trueLabelIndex < trueLabels.size() ;
                                                            trueLabelIndex++ )
      {

        int line = -1 ;
        int column = -1 ;

        if ( predictedLabels[ predictedLabelIndex ] == "Unlabeled" )
        {

          line =  nbLabels ;

        }
        else
        {

          line = getLabelFromName( bundlesNames,
                                  predictedLabels[ predictedLabelIndex ] ) ;
          if ( line == -1 )
          {

            if ( verbose > 2 )
            {

              std::cout << "WARNING : predicted bundle "
                        << predictedLabels[ predictedLabelIndex ]
                        << " not found in true labels dictionary"
                        << std::endl ;

            }

            continue ;

            // std::stringstream outMessageOss ;
            // outMessageOss << "ERROR : while computing confusion matrix could "
            //               << "not get the line number " << std::endl ;
            // std::string outMessage = outMessageOss.str() ;
            //
            // throw( std::invalid_argument( outMessage ) ) ;

          }

        }

        if ( trueLabels[ trueLabelIndex ] == "Unlabeled" )
        {

          column =  nbLabels ;

        }
        else
        {

          column = getLabelFromName( bundlesNames,
                                            trueLabels[ trueLabelIndex ] ) ;
          if ( column == -1 )
          {

            std::stringstream outMessageOss ;
            outMessageOss << "ERROR : while computing confusion matrix could "
                          << "not get the column number " << std::endl ;
            std::string outMessage = outMessageOss.str() ;

            throw( std::invalid_argument( outMessage ) ) ;

          }

        }

        confusionMatrix[ line ][ column ] += 1 ;

        // if ( trueLabels[ trueLabelIndex ] ==
        //                               predictedLabels[ predictedLabelIndex ] )
        // {
        //
        //   break ;
        //
        // }


      }

    }

  }
  std::cout << "Done" << std::endl ;

  // Computing TP, TN, FP and FN
  std::cout << "Computing TP, TN, FP, FN... " ;
  std::vector<int> truePositives( nbLabels + 1, 0 ) ;
  std::vector<int> trueNegatives( nbLabels + 1, 0 ) ;
  std::vector<int> falsePositives( nbLabels + 1, 0 ) ;
  std::vector<int> falseNegatives( nbLabels + 1, 0 ) ;
  for ( int labelIndex = 0 ; labelIndex < nbLabels + 1 ; labelIndex++ )
  {

    truePositives[ labelIndex ] = confusionMatrix[ labelIndex ][ labelIndex ] ;

    for ( int i = 0 ; i < nbLabels + 1 ; i++ )
    {

      if ( i != labelIndex )
      {

        trueNegatives[ labelIndex ] += confusionMatrix[ i ][ i ] ;

      }

      if ( i != labelIndex )
      {

        falsePositives[ labelIndex ] += confusionMatrix[ i ][ labelIndex ] ;

      }

      if ( i != labelIndex )
      {

        falseNegatives[ labelIndex ] += confusionMatrix[ labelIndex ][ i ] ;

      }


    }


  }
  std::cout << "Done" << std::endl ;


  std::vector<float> precisionPerBundle( nbLabels + 1, 0 ) ;
  std::vector<float> sensitivityPerBundle( nbLabels + 1, 0 ) ;
  std::vector<float> accuracyPerBundle( nbLabels + 1, 0 ) ;
  std::vector<float> jaccardPerBundle( nbLabels + 1, 0 ) ;
  float averagePrecision = 0 ;
  float averageSensitivity = 0 ;
  float averageAccuracy = 0 ;
  float averageJaccard = 0 ;
  for ( int labelIndex = 0 ; labelIndex < nbLabels + 1 ; labelIndex++ )
  {

    float _tp = truePositives[ labelIndex ] ;
    float _tn = trueNegatives[ labelIndex ] ;
    float _fp = falsePositives[ labelIndex ] ;
    float _fn = falseNegatives[ labelIndex ] ;

    float precisionTmp = 0 ;
    float sensitivityTmp = 0 ;
    float jaccardTmp = 0 ;
    if ( _tp != 0 )
    {

      precisionTmp = _tp / ( _tp + _fp ) ;

      sensitivityTmp = _tp / ( _tp + _fn ) ;

      jaccardTmp = ( _tp ) / ( _tp + _fn + _fp ) ;

    }

    float accuracyTmp = 0 ;
    if ( ( _tp + _tn ) != 0 )
    {

      accuracyTmp = ( _tp + _tn ) / ( _tp + _tn + _fp + _fn ) ;

    }

    precisionPerBundle[ labelIndex ] = precisionTmp ;
    sensitivityPerBundle[ labelIndex ] = sensitivityTmp ;
    accuracyPerBundle[ labelIndex ] = accuracyTmp ;
    jaccardPerBundle[ labelIndex ] = jaccardTmp ;

    averagePrecision += precisionTmp ;
    averageSensitivity += sensitivityTmp ;
    averageAccuracy += accuracyTmp ;
    averageJaccard += jaccardTmp ;

  }

  averagePrecision /= nbLabels ;
  averageSensitivity /= nbLabels ;
  averageAccuracy /= nbLabels ;
  averageJaccard /= nbLabels ;


  // Saving results
  std::string tractogramName = basenameNoExtension( subjectTractogramPath ) ;

  std::stringstream confusionMatrixPathOss ;
  confusionMatrixPathOss << outputDirectory << tractogramName
                                            << "_confusionMatrix.tsv" ;
  confusionMatrixPath = confusionMatrixPathOss.str() ;
  if ( !force && is_file( confusionMatrixPath ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "Confusion matrix output file " << confusionMatrixPath
                  << " already exists and the -force flag was not used"
                  << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }
  saveConfusionMatrix( confusionMatrixPath.c_str(), confusionMatrix ) ;


  std::stringstream scoresPathOss ;
  scoresPathOss << outputDirectory << tractogramName << "_scores.tsv" ;
  scoresPath = scoresPathOss.str() ;
  if ( !force && is_file( scoresPath ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "Scores per bundle output file " << scoresPath
                  << " already exists and the -force flag was not used"
                  << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }
  saveScorePerBundle( scoresPath.c_str(),
                      precisionPerBundle,
                      sensitivityPerBundle,
                      accuracyPerBundle,
                      jaccardPerBundle,
                      bundlesNames,
                      averagePrecision,
                      averageSensitivity,
                      averageAccuracy,
                      averageJaccard ) ;

}
