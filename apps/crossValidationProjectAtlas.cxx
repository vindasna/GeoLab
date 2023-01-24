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

#include "crossValidationProjectAtlas.h"
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

  if ( verbose > 1 )
  {

    std::cout << "\tSaving confusion matrix in : " << confusionMatrixPath
                                                 << std::endl ;
  }
  std::ofstream file ;
  file.open( confusionMatrixPath ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << confusionMatrixPath <<
                                                                     std::endl ;

    exit( 1 ) ;

  }

  int sizeConfusionMatrix = confusionMatrix.size() ;
  std::cout << "Size confusion matrix : " << sizeConfusionMatrix << std::endl ;

  for ( int _line_i = 0 ; _line_i < sizeConfusionMatrix ; _line_i++ )
  {

    for ( int _column_i = 0 ; _column_i < sizeConfusionMatrix ; _column_i++ )
    {

      file << confusionMatrix[ _line_i ][ _column_i ] << "\t" ;

    }

    file << "\n" ;

  }

  file.close() ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////// Function to save scores per bundle ///////////////////////
////////////////////////////////////////////////////////////////////////////////
void saveScorePerBundle( const char* scoresPerBundlePath,
                         const std::vector<float>& precisionPerBundle,
                         const std::vector<float>& recallPerBundle,
                         const std::vector<float>& accuracyPerBundle,
                         const std::vector<float>& weightsBundles,
                         const std::vector<std::string>& bundleNames,
                         float averagePrecision,
                         float averageRecall,
                         float averageAccuracy )
{

  if ( verbose > 1 )
  {

    std::cout << "\tSaving scores per bundle in : " << scoresPerBundlePath
                                                    << std::endl ;
  }
  std::ofstream file ;
  file.open( scoresPerBundlePath ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << scoresPerBundlePath <<
                                                                   std::endl ;

    exit( 1 ) ;

  }

  file << "Bundle\tPrecision\tRecall\tAccuracy\tWeights\n" ;

  int nbBundles = precisionPerBundle.size() ;
  for ( int _bundleIndex = 0 ; _bundleIndex < nbBundles ; _bundleIndex++ )
  {

    file << bundleNames[ _bundleIndex ] << "\t" ;
    file << precisionPerBundle[ _bundleIndex ] << "\t" ;
    file << recallPerBundle[ _bundleIndex ] << "\t" ;
    file << accuracyPerBundle[ _bundleIndex ] << "\t" ;
    file << weightsBundles[ _bundleIndex ] << "\t" ;

    file << "\n" ;

  }

  file << "Average\t" << averagePrecision << "\t" << averageRecall << "\t"
                                                  << averageAccuracy << "\n" ;

  file.close() ;

}
////////////////////////////////////////////////////////////////////////////////
/////////////////////// Generates true with X probability //////////////////////
////////////////////////////////////////////////////////////////////////////////
bool generateTrueWithXprobability( float probability )
{

  // boost::random::mt19937 gen ;
  //
  // // Generator of integers between [1, 10]
  // boost::random::uniform_int_distribution<> dist( 1, 10 ) ;
  // int randomInteger = dist( gen ) ;

  std::random_device rd ;
  std::mt19937 gen( rd() ) ;
  std::uniform_real_distribution<> dist( 0, 100 ) ;

  int randomInteger = dist( gen ) ;

  if ( randomInteger <= probability * 100 )
  {

    return( true ) ;

  }

  return( false ) ;

}

////////////////////////////////////////////////////////////////////////////////
//////////////// Function to get number of fibers atlas bundle /////////////////
////////////////////////////////////////////////////////////////////////////////
int getNbFibers( const std::string& bundleFilename )
{

  BundlesMinf bundle( bundleFilename.c_str() ) ;
  if ( bundle.curves_count <= 0 )
  {

    std::cout << "ERROR : got invalid fiber count of " << bundle.curves_count
                                                       << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.curves_count ) ;

}




////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  const auto start_time = std::chrono::system_clock::now() ;

  int index_atlas, index_reference, index_output, index_fa, index_lfa,
      index_mlfa, index_anc, index_thrCov, index_thrAdj, index_pa, index_cv,
      index_cn, index_cc, index_rb, index_ods, index_cds, index_nbPoints,
      index_k, index_format, index_ana, index_slr, index_sp, index_force,
                                    index_verbose, index_nbThreads, index_help ;

  index_atlas = getFlagPosition( argc, argv, "-a") ;
  index_reference = getFlagPosition( argc, argv, "-ref") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_cc = getFlagPosition( argc, argv, "-cc") ;
  index_rb = getFlagPosition( argc, argv, "-rb") ;
  index_ods = getFlagPosition( argc, argv, "-ods") ;
  index_cds = getFlagPosition( argc, argv, "-cds") ;
  index_ana = getFlagPosition( argc, argv, "-ana") ;
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints") ;
  index_k = getFlagPosition( argc, argv, "-k") ;
  index_format = getFlagPosition( argc, argv, "-format") ;
  index_fa = getFlagPosition( argc, argv, "-fa") ;
  index_lfa = getFlagPosition( argc, argv, "-lfa") ;
  index_mlfa = getFlagPosition( argc, argv, "-mlfa") ;
  index_anc = getFlagPosition( argc, argv, "-anc") ;
  index_thrCov = getFlagPosition( argc, argv, "-thrCov") ;
  index_thrAdj = getFlagPosition( argc, argv, "-thrAdj") ;
  index_pa = getFlagPosition( argc, argv, "-pa") ;
  index_cv = getFlagPosition( argc, argv, "-cv") ;
  index_cn = getFlagPosition( argc, argv, "-cn") ;
  index_slr = getFlagPosition( argc, argv, "-slr") ;
  index_sp = getFlagPosition( argc, argv, "-sp") ;
  index_force = getFlagPosition( argc, argv, "-force") ;
  index_nbThreads = getFlagPosition( argc, argv, "-nbThreads") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-a : path to atlas directoy \n"
              << "-fa : Path to full atlas tractogram \n"
              << "-lfa : labels per fiber full atlas \n"
              << "-mlfa : multi-labels per fiber full atlas \n"
              << "-ref : Path to the reference image where the tractogram is \n"
              << "-o : Path to the output directory \n"
              << "-ana : Path to the analyseAtlasBundle.py file \n"
              << "-nbPoints : Number of points per fiber (same number for all "
              << "fibers) \n"
              << "-k : Value of k for k-fold cross-validation \n"
              << "-format : format used among .tck/.trk/.bundles \n :"
              << "[-cc] : Path to the clientComputeCentroids.py file \n"
              << "[-rb] : Path to the clientRegisterBundles.py file \n"
              << "[-ods] : Path to the dipyServer.py file \n"
              << "[-cds] : Path to the clientCloseServer.py file \n"
              << "[-anc] : Path to the atlas neighborhoods centroids (must be "
              << "given is -fa is not given)\n"
              << "[-thrCov] : Threshold to keep bundles where coverage is "
              << "greater than thrCov (default : 0 -> keep all bundles ) \n"
              << "[-thrAdJ] : Threshold to keep bundles where adjacency is "
              << "greater than thrAdj (default : 0 -> keep all bundles ) \n"
              << "[-pa] : Path to the ProjectAtlas file \n"
              << "[-cv] : Path to the ConvertBundleFormat file \n"
              << "[-cn] : Path to the computeNeighborhood file \n"
              << "[-slr] : Do global SLR step (default : false)\n"
              << "[-sp] : Save recognized bundles separetly (default : true)\n"
              << "[-force] : Force to overwrite files (default = false) \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads before "
              << "applyRecoBundles (default : let openMP decide ) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_atlas )
  {

    std::cout << "-a argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_fa )
  {

    std::cout << "-fa argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_lfa )
  {

    std::cout << "-lfa argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_mlfa )
  {

    std::cout << "-mlfa argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_fa && !index_anc )
  {

    std::cout << "At least -fa or -anc must be given ..." << std::endl ;
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


  if ( !index_ana )
  {

    std::cout << "-ana argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_nbPoints )
  {

    std::cout << "-nbPoints argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_k )
  {

    std::cout << "-k argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_nbPoints )
  {

    std::cout << "-nb argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_format )
  {

    std::cout << "-format argument required ..." << std::endl ;
    exit( 1 ) ;

  }


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Checking arguments /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  ////////////////////////////// Atlas directory ///////////////////////////////
  atlasDirectory = argv[ index_atlas + 1 ] ;
  char lastChar = atlasDirectory[ atlasDirectory.size() - 1 ] ;
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

  /////////////////////////////// Reference image //////////////////////////////
  referenceFilename = argv[ index_reference + 1 ]  ;
  lastChar = referenceFilename[ referenceFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    referenceFilename = referenceFilename.substr( 0,
                                                referenceFilename.size() - 1 ) ;

  }
  if ( referenceFilename.find( ".nii.gz" ) != std::string::npos )
  {

    std::cout << "ERROR : .nii.gz format not supported for the reference image "
              << std::endl ;
    exit( 1 ) ;

  }
  else if ( referenceFilename.find( ".nii" ) != std::string::npos )
  {

    std::cout << "Reference image : OK " << std::endl ;

  }
  else
  {

    std::cout << "The only reference image format supported is .nii"
              << std::endl ;
    exit( 1 ) ;

  }
  if ( !is_file( referenceFilename ) )
  {

    std::cout << "ERROR : Reference image file " << referenceFilename
                                            << " does not exists" << std::endl ;

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

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }

      throw( std::invalid_argument( outMessage ) ) ;


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

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }

      throw( std::invalid_argument( outMessage ) ) ;

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

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }

      throw( std::invalid_argument( outMessage ) ) ;

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

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }

      throw( std::invalid_argument( outMessage ) ) ;

    }


  }
  else
  {

    closeDipyServerClientFile = "clientCloseServer.py" ;

  }

  /////////////////////// Analyse atlas bundles command ////////////////////////
  analyseAtlasBundleFile = argv[ index_ana + 1 ] ;
  lastChar = analyseAtlasBundleFile[ analyseAtlasBundleFile.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    analyseAtlasBundleFile = analyseAtlasBundleFile.substr( 0,
                                              analyseAtlasBundleFile.size() - 1 ) ;

  }
  if ( !is_file( analyseAtlasBundleFile ) )
  {

    std::cout << "ERROR : analyseAtlasBundle.py file " << analyseAtlasBundleFile
              << " does not exists " << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

  ///////////////////////// Number of points per fiber /////////////////////////
  nbPointsPerFiber = std::stoi( argv[ index_nbPoints + 1 ] ) ;

  //////////////////////////// Format of tractograms ///////////////////////////
  format = argv[ index_format + 1 ] ;
  if ( format != ".bundles" && format != ".bundlesdata" && format != ".trk"
                                                           && format != ".tck" )
  {

    std::cout << "ERROR : input -format must be .bundles/.trk/.tck" << "\n" ;
    exit( 1 ) ;

  }

  ///////////////////////////////// Full atlas /////////////////////////////////
  if ( index_fa )
  {

    std::string tmpFullAtlasPath( argv[ index_fa + 1 ] ) ;
    isFullAtlas = true ;

    char lastChar = tmpFullAtlasPath[ tmpFullAtlasPath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      tmpFullAtlasPath = tmpFullAtlasPath.substr( 0,
                                                tmpFullAtlasPath.size() - 1 ) ;

    }

    if ( !endswith( tmpFullAtlasPath, format ) )
    {

      std::cout << "ERROR : fullAtlaas given to -fa and format given to -format"
                << "must be the same" << std::endl ;
      exit( 1 ) ;

    }

    if ( format == ".bundles" )
    {

      fullAtlasBundlesFilename = tmpFullAtlasPath ;
      fullAtlasBundlesDataFilename = replaceExtension( tmpFullAtlasPath,
                                                              ".bundlesdata" ) ;

    }
    else if ( format == ".bundlesdata" )
    {

      fullAtlasBundlesDataFilename = tmpFullAtlasPath ;
      fullAtlasBundlesFilename = replaceExtension( tmpFullAtlasPath,
                                                                  ".bundles" ) ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      fullAtlasBundlesDataFilename = tmpFullAtlasPath ;
      fullAtlasBundlesFilename = replaceExtension( tmpFullAtlasPath, ".minf" ) ;

    }
    else
    {

      std::cout << "The only full atlas formate supported are .bundles/.trk/"
                << ".tck" << std::endl ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }

      exit( 1 ) ;

    }
    if ( !is_file( fullAtlasBundlesFilename ) )
    {

      std::cout << "ERROR : Full atlas file " << fullAtlasBundlesFilename
                                            << " does not exists" << std::endl ;

    }
    else
    {

      std::cout << "Full atlas tractogram : OK " << std::endl ;

    }

  }

  ////////////////////////////// Full atlas labels /////////////////////////////
  if ( index_lfa )
  {

    fullAtlasLabelsFilename = argv[ index_lfa + 1 ] ;
    isFullAtlas = true ;

    char lastChar = fullAtlasLabelsFilename[ fullAtlasLabelsFilename.size() - 1
                                                                             ] ;
    if ( lastChar == '/' )
    {

      fullAtlasLabelsFilename = fullAtlasLabelsFilename.substr( 0,
                                          fullAtlasLabelsFilename.size() - 1 ) ;

    }

    if ( !is_file( fullAtlasLabelsFilename ) )
    {

      std::cout << "ERROR : Full atlas labels file " << fullAtlasLabelsFilename
                                            << " does not exists" << std::endl ;

    }
    else
    {

      std::cout << "Full atlas labels : OK " << std::endl ;

    }

    fullAtlasLabelsDictFilename = replaceExtension( fullAtlasLabelsFilename,
                                                                     ".dict" ) ;

  }

  /////////////////////////// Full atlas multi labels //////////////////////////
  if ( index_mlfa )
  {

    fullAtlasMultiLabelsFilename = argv[ index_mlfa + 1 ] ;
    isFullAtlas = true ;

    char lastChar = fullAtlasMultiLabelsFilename[
                                     fullAtlasMultiLabelsFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      fullAtlasMultiLabelsFilename = fullAtlasMultiLabelsFilename.substr( 0,
                                     fullAtlasMultiLabelsFilename.size() - 1 ) ;

    }

    if ( !is_file( fullAtlasMultiLabelsFilename ) )
    {

      std::cout << "ERROR : Full atlas multi labels file "
                << fullAtlasMultiLabelsFilename << " does not exists"
                << std::endl ;

    }
    else
    {

      std::cout << "Full atlas multi labels : OK " << std::endl ;

    }

    fullAtlasMultiLabelsDictFilename = replaceExtension(
                                                   fullAtlasMultiLabelsFilename,
                                                    ".dict" ) ;

  }

  //////////////////////// Atlas neighborhood centroids ////////////////////////
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

      std::cout << "ERROR : Atlas neighborhood centroids directory "
                << atlasNeighborhoodCentroidsDirectory << " does not exists "
                << std::endl ;
      exit( 1 );

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

  ///////////////////////////// Project atlas file /////////////////////////////
  if ( index_pa )
  {

    std::string tmpFile( argv[ index_pa + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::cout << "ERROR : projecAtlas file " << tmpFile
                << " does not exists " << std::endl ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }

      exit( 1 ) ;

    }

    projectAtlasFile = tmpFile ;

  }

  ////////////////////////// ConvertBundleFormat file //////////////////////////
  if ( index_cv )
  {

    std::string tmpFile( argv[ index_cv + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::cout << "ERROR : ConvertBundleFormat file " << tmpFile
                << " does not exists " << std::endl ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }

    convertBundleFormatsFile = tmpFile ;

  }

  ////////////////////////// computeNeighborhood file //////////////////////////
  if ( index_cn )
  {

    std::string tmpFile( argv[ index_cn + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::cout << "ERROR : computeNeighborhood file " << tmpFile
                << " does not exists " << std::endl ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
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

      std::cout << "Argument of -slr must be either \"true\" or \"false\" "
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
    if ( nbThreads < 0 )
    {

      std::cout << "Invalid argument for -nbThreads : you must give a postive "
                << "integer " << std::endl ;
      exit( 1 ) ;

    }

  }

  ////////////////////////////////// k value ///////////////////////////////////
  if ( index_k )
  {

    k_value = std::stoi( argv[ index_k + 1 ] ) ;
    if ( k_value < 0 )
    {

      std::cout << "Invalid argument for -k : you must give a postive "
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

    std::cout << "ERROR : the only formats supported for -format are "
              << ".bundles/.trk./tck" << std::endl ;
    exit( 1 ) ;

  }
  percentageTest = 1.0 / k_value ;

  std::vector<float> accuraciesKfold( k_value ) ;
  std::vector<float> f1ScoreKfold( k_value ) ;

  AtlasBundles atlasData( atlasDirectory.c_str(),
                          isBundles,
                          isTrk,
                          isTck,
                          0 ) ;

  BundlesData fullAtlasData( fullAtlasBundlesDataFilename.c_str() ) ;

  int nbFibers = fullAtlasData.curves_count ;

  std::vector< std::vector<int16_t>> fullAtlasLabels ;
  readingPredictedLabels( fullAtlasLabelsFilename.c_str(),
                          fullAtlasLabels,
                          nbFibers ) ;
  std::vector<std::string> fullAtlasLabelsDict ;
  readLabelsDict( fullAtlasLabelsDictFilename.c_str(),
                  fullAtlasLabelsDict ) ;

  int nbBundles = fullAtlasLabelsDict.size() ;

  std::vector< std::vector<int16_t>> fullAtlasMultiLabels ;
  readingPredictedLabels( fullAtlasMultiLabelsFilename.c_str(),
                          fullAtlasMultiLabels,
                          nbFibers ) ;
  std::vector<std::string> fullAtlasMultiLabelsDict ;
  readLabelsDict( fullAtlasMultiLabelsDictFilename.c_str(),
                  fullAtlasMultiLabelsDict ) ;

  if ( nbBundles != fullAtlasMultiLabelsDict.size() )
  {

    std::cout << "ERROR : the number of bundles in labels dict has to be the "
              << "same as in multi-label dict, got " << nbBundles
              << " for full atlas labels and "
              << fullAtlasMultiLabelsDict.size()
              << " for full atlas multi-labels" << std::endl ;
    exit( 1 ) ;

  }

  if ( verbose )
  {

    std::cout << "Separating set in " << k_value << " subsets... "
                                                                  << std::endl ;

  }


  bool isDivisionDone = false ;
  for ( int crossValidationIndex = 0 ; crossValidationIndex < k_value ;
                                                        crossValidationIndex++ )
  {

    std::ostringstream dirSubset_iOss ;
    dirSubset_iOss << outputDirectory << "subset_"
                                            << crossValidationIndex + 1 << "/" ;
    std::string dirSubset_i = dirSubset_iOss.str() ;
    if ( !is_dir( dirSubset_i ) )
    {

      mkdir( dirSubset_i ) ;

    }

    if ( countFilesDirectory( dirSubset_i ) > 0 && !force )
    {

      std::cout << "\tWARNING : it seems the set has already been divided "
                << "and -force flag was not used, skipping separation step "
                << std::endl ;
      isDivisionDone = true ;
      break ;

    }

  }

  for ( int bundleIndex = 0 ; bundleIndex < nbBundles ; bundleIndex++ )
  {

    printf( "\rProcessing bundle : [%d/%d ]", bundleIndex + 1, nbBundles ) ;
    std::cout << "" << std::flush ;

    if ( isDivisionDone )
    {

      break ;

    }

    BundlesMinf& tmpBundles = atlasData.bundlesMinf[ bundleIndex ] ;
    BundlesData& tmpBundlesData = atlasData.bundlesData[ bundleIndex ] ;
    std::string bundleName = atlasData.bundlesNames[ bundleIndex ] ;

    int labelInFullAtlas = -2 ;
    for ( int i = 0 ; i < fullAtlasLabelsDict.size() ; i++ )
    {

      std::string _tmpBundleName = fullAtlasLabelsDict[ i ] ;
      if ( _tmpBundleName == bundleName )
      {

        labelInFullAtlas = i ;
        break ;

      }

    }
    if ( labelInFullAtlas == -2 )
    {

      std::cout << "ERROR : could not find label for " << bundleName
                << " in full atlas labels dictionary " << std::endl ;
      exit( 1 ) ;

    }

    std::vector<float>& tmpMatrixTracks = tmpBundlesData.matrixTracks ;
    std::vector<int32_t>& tmpPointsPerTrack = tmpBundlesData.pointsPerTrack ;
    int nbFibersBundle = tmpPointsPerTrack.size() ;
    int nbFibersTest = nbFibersBundle / k_value ;

    std::vector<int> indexCurrentBundleInFullAtlas ;

    for ( int fiberIndex = 0 ; fiberIndex < nbFibers ; fiberIndex++ )
    {

      for ( int _tmpIndex = 0 ;
                      _tmpIndex < fullAtlasLabels[ fiberIndex ].size() ;
                                                                   _tmpIndex++ )
      {

        if ( fullAtlasLabels[ fiberIndex ][ _tmpIndex ] == labelInFullAtlas )
        {

          indexCurrentBundleInFullAtlas.push_back( fiberIndex ) ;

        }

      }

      if ( indexCurrentBundleInFullAtlas.size() == nbFibersBundle )
      {

        break ;

      }

    }

    std::vector<std::vector<int>> fiberIndexPerSet( k_value ) ;
    for ( int crossValidationIndex = 0 ; crossValidationIndex < k_value ;
                                                        crossValidationIndex++ )
    {

      std::ostringstream dirSubset_iOss ;
      dirSubset_iOss << outputDirectory << "subset_"
                                            << crossValidationIndex + 1 << "/" ;
      std::string dirSubset_i = dirSubset_iOss.str() ;

      if ( crossValidationIndex == k_value - 1 )
      {

        nbFibersTest += nbFibersBundle % k_value ;

      }

      std::vector<float> testMatrixTracks ;
      std::vector<int32_t> testPointsPerTracks ;
      std::vector<float> trainMatrixTracks ;
      std::vector<int32_t> trainPointsPerTracks ;

      std::vector<int>& fiberIndexTest = fiberIndexPerSet[
                                                        crossValidationIndex ] ;

      if ( indexCurrentBundleInFullAtlas.size() != nbFibersBundle )
      {

        std::cout << "ERROR : the number of fibers with label "
                  << bundleName << " in the full atlas is different than "
                  << " in bundle file in atlas directory, got "
                  << indexCurrentBundleInFullAtlas.size()
                  << " in full atlas and " << nbFibersBundle
                  << " in atlas directoy "
                  << std::endl ;
        exit( 1 ) ;

      }

      bool _stop = false ;
      while( !_stop )
      {

        for ( int fiberIndex = 0 ; fiberIndex < nbFibersBundle ;
                                                              fiberIndex++ )
        {

          int fiberIndexBundleInTractogram = indexCurrentBundleInFullAtlas[
                                                                fiberIndex ] ;

          bool isIndexAlreadySelected = false ;
          for ( int _tmpFoldIndex = 0 ;
                                       _tmpFoldIndex <= crossValidationIndex ;
                                                             _tmpFoldIndex++ )
          {

            if ( !fiberIndexPerSet[ _tmpFoldIndex ].empty() &&
                         std::find( fiberIndexPerSet[ _tmpFoldIndex ].begin(),
                                      fiberIndexPerSet[ _tmpFoldIndex ].end(),
                                             fiberIndexBundleInTractogram ) !=
                                     fiberIndexPerSet[ _tmpFoldIndex ].end() )
            {

              isIndexAlreadySelected = true ;
              break ;

            }

          }
          if ( isIndexAlreadySelected )
          {

            continue ;

          }


          for ( int _tmpIndex = 0 ; _tmpIndex <
                    fullAtlasLabels[ fiberIndexBundleInTractogram ].size() ;
                                                                 _tmpIndex++ )
          {

            if ( fullAtlasLabels[ fiberIndexBundleInTractogram ][
                                             _tmpIndex ] == labelInFullAtlas )
            {

              if ( generateTrueWithXprobability( percentageTest ) )
              {

                fiberIndexTest.push_back( fiberIndexBundleInTractogram ) ;
                if ( fiberIndexTest.size() == nbFibersTest )
                {

                  _stop = true ;
                  break ;

                }

              }

            }

          }

          if ( _stop )
          {

            break ;

          }

        }

      }

      if ( crossValidationIndex == k_value - 1 )
      {

        for ( int _tmpInd = 0 ; _tmpInd < k_value ; _tmpInd++ )
        {

          for ( int _tmpInd2 = 0 ; _tmpInd2 < k_value ; _tmpInd2++ )
          {

            if ( _tmpInd == _tmpInd2 )
            {

              continue ;

            }

            for ( int _iTmp = 0 ; _iTmp < fiberIndexPerSet[ _tmpInd ].size() ;
                                                                       _iTmp++ )
            {

              if ( std::find( fiberIndexPerSet[ _tmpInd2 ].begin(),
                              fiberIndexPerSet[ _tmpInd2 ].end(),
                              fiberIndexPerSet[ _tmpInd ][ _iTmp ] ) !=
                                            fiberIndexPerSet[ _tmpInd2 ].end() )
              {

                std::cout << "Repeated index in separation of bundle "
                          << bundleName << std::endl ;
                std::cout << _tmpInd << "-fold" << std::endl ;
                for ( int _iTmp2 ; _iTmp2 < fiberIndexPerSet[ _tmpInd ].size() ;
                                                                      _iTmp2++ )
                {

                  std::cout << fiberIndexPerSet[ _tmpInd ][ _iTmp2 ]
                                                                  << std::endl ;

                }
                std::cout << _tmpInd2 << "-fold" << std::endl ;
                for ( int _iTmp2 ; _iTmp2 < fiberIndexPerSet[ _tmpInd ].size() ;
                                                                      _iTmp2++ )
                {

                  std::cout << fiberIndexPerSet[ _tmpInd2 ][ _iTmp2 ]
                                                                  << std::endl ;

                }
                exit( 1 ) ;

              }

            }

          }

        }

      }




      std::vector<std::vector<int16_t>> multiLabelsBundle ;
      std::vector<std::string>& multiLabelsDictBundle =
                                                      fullAtlasMultiLabelsDict ;
      for ( int i = 0 ; i < nbFibersTest ; i++ )
      {

        // Only works if all the fibers have the same number of points
        int fiberIndex = fiberIndexTest[ i ] ;
        int64_t offset = 3 * nbPointsPerFiber * fiberIndex ;
        testPointsPerTracks.push_back(
                                fullAtlasData.pointsPerTrack[ fiberIndex ] ) ;
        for ( int _point = 0 ; _point < nbPointsPerFiber ; _point++ )
        {

          for ( int coord = 0 ; coord < 3 ; coord++ )
          {

            testMatrixTracks.push_back( fullAtlasData.matrixTracks[
                                               3 * _point + coord + offset ] ) ;

          }

        }

        std::vector<int16_t> _tmpVectorLabel ;
        for ( int _tmpIndex = 0 ;
                       _tmpIndex < fullAtlasMultiLabels[ fiberIndex ].size() ;
                                                                 _tmpIndex++ )
        {

          int16_t _tmpLabelFiber =
                               fullAtlasMultiLabels[ fiberIndex ][ _tmpIndex ] ;

          if ( _tmpLabelFiber == -1 )
          {

            _tmpVectorLabel.push_back( fullAtlasLabels[ fiberIndex ][ 0 ] ) ;

          }
          else
          {

            _tmpVectorLabel.push_back( _tmpLabelFiber ) ;

          }

        }

        multiLabelsBundle.push_back( _tmpVectorLabel ) ;

      }


      for ( int i = 0 ; i < indexCurrentBundleInFullAtlas.size() ; i++ )
      {

        int _tmpFiberIndex = indexCurrentBundleInFullAtlas[ i ] ;

        if ( std::find( fiberIndexTest.begin(), fiberIndexTest.end(),
                                    _tmpFiberIndex ) == fiberIndexTest.end() )
        {

          // Only works if all the fibers have the same number of points
          int fiberIndex = _tmpFiberIndex ;
          int64_t offset = 3 * nbPointsPerFiber * fiberIndex ;
          trainPointsPerTracks.push_back(
                                fullAtlasData.pointsPerTrack[ fiberIndex ] ) ;
          for ( int _point = 0 ; _point < nbPointsPerFiber ; _point++ )
          {

            for ( int coord = 0 ; coord < 3 ; coord++ )
            {

              trainMatrixTracks.push_back( fullAtlasData.matrixTracks[
                                               3 * _point + coord + offset ] ) ;

            }

          }

        }

      }


      // Saving test
      std::ostringstream bundlesDataTestnOss ;
      bundlesDataTestnOss << dirSubset_i << bundleName << format ;
      std::string bundlesDataTest = bundlesDataTestnOss.str() ;
      if ( format == ".bundles" )
      {

        bundlesDataTest = replaceExtension( bundlesDataTest, ".bundlesdata" ) ;

      }

      std::ostringstream bundlesTestOss ;
      bundlesTestOss << dirSubset_i << bundleName << format ;
      std::string bundlesTest = bundlesTestOss.str() ;
      if ( format == ".trk" || format == ".tck" )
      {

        bundlesTest = replaceExtension( bundlesTest , ".minf" ) ;

      }
      else if ( format == ".bundlesdata" )
      {

        bundlesTest = replaceExtension( bundlesTest, ".bundles" ) ;

      }

      int curvesCountTest = testPointsPerTracks.size() ;

      BundlesMinf testBundleBundles( tmpBundles ) ;
      testBundleBundles.curves_count = curvesCountTest ;
      testBundleBundles.write( bundlesTest.c_str() ) ;

      BundlesData testBundleBundlesData( tmpBundlesData ) ;
      testBundleBundlesData.matrixTracks = testMatrixTracks,
      testBundleBundlesData.pointsPerTrack = testPointsPerTracks,
      testBundleBundlesData.curves_count = curvesCountTest ;
      testBundleBundlesData.write( bundlesDataTest.c_str(),
                                                           testBundleBundles ) ;




      std::ostringstream bundlesTestMultiLabelsOss ;
      bundlesTestMultiLabelsOss << dirSubset_i << bundleName << "_labels.txt" ;
      std::string bundlesTestMultiLabels = bundlesTestMultiLabelsOss.str() ;

      saveLabels( bundlesTestMultiLabels.c_str(),
                  multiLabelsBundle ) ;

      std::ostringstream bundlesTestMultiLabelsDictOss ;
      bundlesTestMultiLabelsDictOss << dirSubset_i << bundleName
                                                           << "_labels.dict" ;
      std::string bundlesTestMultiLabelsDict =
                                         bundlesTestMultiLabelsDictOss.str() ;

      saveLabelsDict( bundlesTestMultiLabelsDict.c_str(),
                      multiLabelsDictBundle ) ;


    }

  }

  std::cout << "\nDone" << std::endl ;




  for ( int crossValidationIndex = 0 ; crossValidationIndex < k_value ;
                                                       crossValidationIndex ++ )
  {

    std::ostringstream foldDirectoryOss ;
    foldDirectoryOss << outputDirectory << "fold_" << crossValidationIndex
                                                                        << "/" ;
    std::string foldDirectory = foldDirectoryOss.str() ;
    if ( !is_dir( foldDirectory ) )
    {

      mkdir( foldDirectory ) ;

    }


    if ( verbose )
    {

      std::cout << "Creating train set for fold " << crossValidationIndex
                                                                  << std::endl ;

    }

    std::ostringstream trainDirectoryOss ;
    trainDirectoryOss << foldDirectory << "train" << "/" ;
    std::string trainDirectory = trainDirectoryOss.str() ;
    if ( !is_dir( trainDirectory ) )
    {

      mkdir( trainDirectory ) ;

    }

    bool isTrainDone = false ;
    if ( countFilesDirectory( trainDirectory ) > 5 )
    {

      std::cout << "WARNING : for fold " << crossValidationIndex << " train set"
                << " seems to already exists, skipping computations..."
                << std::endl ;
      isTrainDone = true ;

    }

    std::vector<std::string>& bundlesNames = fullAtlasLabelsDict ;
    std::vector<BundlesMinf> trainBundles( nbBundles ) ;
    std::vector<BundlesData> trainBundlesData( nbBundles ) ;
    std::vector<std::vector<std::vector<int16_t>>> labelsBundles( nbBundles ) ;
    if ( !isTrainDone )
    {

      std::vector<int> foldsForTrain ;
      for ( int _fold = 0 ; _fold < k_value ; _fold++ )
      {

        if ( _fold == crossValidationIndex )
        {

          continue ;

        }

        std::ostringstream dirSubset_iOss ;
        dirSubset_iOss << outputDirectory << "subset_" << _fold + 1 << "/" ;
        std::string dirSubset_i = dirSubset_iOss.str() ;


        std::vector<std::string> filesInSubsetDir ;
        listDir( dirSubset_i, filesInSubsetDir ) ;


        for ( int _fileIndex = 0 ; _fileIndex < filesInSubsetDir.size() ;
                                                                  _fileIndex++ )
        {

          std::string _filePath = filesInSubsetDir[ _fileIndex ] ;
          if ( !endswith( _filePath , format ) )
          {

            continue ;

          }

          std::string bundlesPath = _filePath ;
          if ( format == ".bundlesdata" )
          {

            bundlesPath = replaceExtension( bundlesPath, ".bundles" ) ;

          }
          else if ( format == ".trk" || format == ".tck" )
          {

            bundlesPath = replaceExtension( bundlesPath, ".minf" ) ;

          }


          std::string bundlesDataPath = replaceExtension( bundlesPath,
                                                                      format ) ;
          if ( format == ".bundles" )
          {

            bundlesDataPath = replaceExtension( bundlesPath, ".bundlesdata" ) ;

          }


          std::string bundleName = basenameNoExtension( bundlesPath ) ;
          int bundleIndexInFullAtlasDict = -1 ;
          for ( int i = 0 ; i < nbBundles ; i++ )
          {

            if ( bundlesNames[ i ] == bundleName )
            {

              bundleIndexInFullAtlasDict = i ;
              break ;

            }

          }


          if ( !is_file( bundlesPath ) || !is_file( bundlesDataPath ) )
          {

            std::cout << "ERROR : could not find .bundles/.bundlesdata or "
                      << ".minf/.trk or .minf/.tck couple for  bundle "
                      << bundleName << " in subset directory "
                      << dirSubset_i << std::endl ;
            exit( 1 ) ;

          }


          BundlesMinf _tmpBundles( bundlesPath.c_str() ) ;
          BundlesData _tmpBundlesData( bundlesDataPath.c_str() ) ;


          if ( trainBundles[ bundleIndexInFullAtlasDict ].curves_count > 0 )
          {

            trainBundles[ bundleIndexInFullAtlasDict ].curves_count +=
                                                      _tmpBundles.curves_count ;

          }
          else
          {

            trainBundles[ bundleIndexInFullAtlasDict ] = _tmpBundles ;


          }

          // Works for different number of points per fiber
          if ( trainBundlesData[ bundleIndexInFullAtlasDict ].curves_count > 0 )
          {

            std::vector<int32_t>& _tmpPointsPerTrack = trainBundlesData[
                                   bundleIndexInFullAtlasDict ].pointsPerTrack ;
            std::vector<float>& _tmpMatrixTracks = trainBundlesData[
                                     bundleIndexInFullAtlasDict ].matrixTracks ;

            int64_t offset = 0 ;
            for ( int _fiber = 0 ; _fiber < _tmpBundlesData.curves_count ;
                                                                      _fiber++ )
            {

              int _tmpNbPointsFiber = _tmpBundlesData.pointsPerTrack[ _fiber ] ;
              _tmpPointsPerTrack.push_back( _tmpNbPointsFiber ) ;
              for ( int _point = 0 ; _point < _tmpNbPointsFiber ; _point++ )
              {

                for ( int coord = 0 ; coord < 3 ; coord++ )
                {

                  _tmpMatrixTracks.push_back( _tmpBundlesData[ offset +
                                                        3 * _point + coord ] ) ;

                }

              }

              offset += 3 * _tmpNbPointsFiber ;

            }

            trainBundlesData[ bundleIndexInFullAtlasDict ].curves_count +=
                                                  _tmpBundlesData.curves_count ;

          }
          else
          {

            trainBundlesData[ bundleIndexInFullAtlasDict ] = _tmpBundlesData ;

          }

          std::ostringstream _tmpDictBundlePathOss ;
          _tmpDictBundlePathOss << dirSubset_i << bundleName
                                                             << "_labels.dict" ;
          std::string _tmpDictBundlePath =_tmpDictBundlePathOss.str() ;
          std::vector<std::string> _tmpDictBundle ;
          readLabelsDict( _tmpDictBundlePath.c_str(),
                          _tmpDictBundle ) ;

          std::ostringstream _tmpLabelsBundlePathOss ;
          _tmpLabelsBundlePathOss << dirSubset_i << bundleName
                                                              << "_labels.txt" ;
          std::string _tmpLabelsBundlePath =_tmpLabelsBundlePathOss.str() ;
          std::vector<std::vector<int16_t>> _tmpLabelsBundle ;
          readingPredictedLabels( _tmpLabelsBundlePath.c_str(),
                                  _tmpLabelsBundle,
                                  _tmpBundlesData.curves_count ) ;


          for ( int i = 0 ; i < _tmpBundlesData.curves_count ; i++ )
          {

            std::vector<int16_t> _tmpLabels ;
            for ( int _k = 0 ; _k < _tmpLabelsBundle[ i ].size() ; _k++ )
            {

              int16_t _tmpLabelInFiberBundle = _tmpLabelsBundle[ i ][ _k ] ;
              std::string _tmpBundleNameInBundleDict = _tmpDictBundle[
                                                    _tmpLabelInFiberBundle ] ;

              int16_t _tmpLabelInFullAtlasDict = -2 ;
              for ( int _tmpDictIndex = 0 ;  _tmpDictIndex < nbBundles ;
                                                             _tmpDictIndex++ )
              {

                if ( _tmpBundleNameInBundleDict == bundlesNames[
                                                             _tmpDictIndex ] )
                {

                  _tmpLabelInFullAtlasDict  = _tmpDictIndex ;
                  break ;

                }

              }
              if ( _tmpLabelInFullAtlasDict == -2 )
              {

                std::cout << "ERROR : did not find label for bundle "
                          << _tmpBundleNameInBundleDict << " in fullAtlasDict"
                          << std::endl ;
                exit( 1 ) ;

              }

              _tmpLabels.push_back( _tmpLabelInFullAtlasDict ) ;

            }

            labelsBundles[ bundleIndexInFullAtlasDict ].push_back(
                                                                _tmpLabels ) ;

          }

        }

      }


      AtlasBundles trainAtlas( trainBundles, trainBundlesData, bundlesNames,
                                                     isBundles, isTrk, isTck ) ;
      trainAtlas.write( trainDirectory.c_str(), isBundles, isTrk, isTck, 0 ) ;
      for ( int i = 0 ; i < bundlesNames.size() ; i++ )
      {

        std::string bundleName = bundlesNames[ i ] ;
        std::vector<std::vector<int16_t>> bundleLabels = labelsBundles[ i ] ;

        std::ostringstream _tmpBundleLabelsPathOss ;
        _tmpBundleLabelsPathOss << trainDirectory << bundleName
                                                              << "_labels.txt" ;
        std::string _tmpBundleLabelsPath = _tmpBundleLabelsPathOss.str() ;

        std::string _tmpBundleDictPath = replaceExtension( _tmpBundleLabelsPath,
                                                                     ".dict" ) ;

        saveLabels( _tmpBundleLabelsPath.c_str(), bundleLabels ) ;
        saveLabelsDict( _tmpBundleDictPath.c_str(), bundlesNames ) ;

      }

    }


    if ( verbose )
    {

      std::cout << "Creating validation set for fold " << crossValidationIndex
                                                                  << std::endl ;

    }

    std::ostringstream valDirectoryOss ;
    valDirectoryOss << foldDirectory << "val" << "/" ;
    std::string valDirectory = valDirectoryOss.str() ;
    if ( !is_dir( valDirectory ) )
    {

      mkdir( valDirectory ) ;

    }

    bool isValDone = false ;
    if ( countFilesDirectory( valDirectory ) > 5 )
    {

      std::cout << "WARNING : for fold " << crossValidationIndex
                << " validation ser seems to already exists, skipping "
                << "computations..." << std::endl ;
      isValDone = true ;

    }

    if ( !isValDone )
    {

      std::ostringstream dirSubset_iOss ;
      dirSubset_iOss << outputDirectory << "subset_" << crossValidationIndex + 1
                                                                        << "/" ;
      std::string dirSubset_i = dirSubset_iOss.str() ;

      copytree( dirSubset_i, valDirectory, true ) ;

    }



    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////// Fuse test bundles /////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Fusing test set " << std::endl ;

    }

    std::string& testDir = valDirectory ;
    std::string& trainDir = trainDirectory ;

    std::ostringstream fullTestdBundlesDataPathOss ;
    fullTestdBundlesDataPathOss << testDir << "fullTest" << format ;
    std::string fullTestdBundlesDataPath = fullTestdBundlesDataPathOss.str() ;
    if ( format == ".bundles" )
    {

      fullTestdBundlesDataPath = replaceExtension( fullTestdBundlesDataPath,
                                                              ".bundlesdata" ) ;

    }

    std::ostringstream fullTestBundlesPathOss ;
    fullTestBundlesPathOss << testDir << "fullTest" << format ;
    std::string fullTestBundlesPath = fullTestBundlesPathOss.str() ;
    if ( format == ".trk" || format == ".tck" )
    {

      fullTestBundlesPath = replaceExtension( fullTestBundlesPath, ".minf" ) ;

    }

    std::ostringstream fullTestLabelsPathOss ;
    fullTestLabelsPathOss << testDir << "fullTestLabels.txt" ;
    std::string fullTestLabelsPath = fullTestLabelsPathOss.str() ;

    std::ostringstream fullTestLabelsDictPathOss ;
    fullTestLabelsDictPathOss << testDir << "fullTestLabels.dict" ;
    std::string fullTestLabelsDictPath = fullTestLabelsDictPathOss.str() ;

    int curves_countFullTest = -1 ;

    if ( !force && is_file( fullTestdBundlesDataPath ) &&
                                                is_file( fullTestBundlesPath ) )
    {

      std::cout << "\tWARNING : it seems the test set has already been "
                << " fuse into one file and -force flag was not used, skipping "
                << "fusing step " << std::endl ;

      curves_countFullTest = getNbFibers( fullTestBundlesPath ) ;

    }
    else
    {

      AtlasBundles testAtlasData( testDir.c_str(),
                                  isBundles,
                                  isTrk,
                                  isTck,
                                  0 ) ;

      std::vector<std::vector<int16_t>> labelsTest ;
      std::vector<std::string>& dictTest = fullAtlasLabelsDict ;
      std::vector<float> fullTestMatrixTracks ;
      std::vector<int32_t> fullTestPointsPerTrack ;

      for ( int16_t bundleIndex = 0 ; bundleIndex < nbBundles ; bundleIndex++ )
      {

        BundlesMinf& tmpBundles = testAtlasData.bundlesMinf[ bundleIndex ] ;
        BundlesData& tmpBundlesData = testAtlasData.bundlesData[ bundleIndex ] ;
        std::string bundleName = testAtlasData.bundlesNames[ bundleIndex ] ;
        std::vector<float>& tmpMatrixTracks = tmpBundlesData.matrixTracks ;
        std::vector<int32_t>& tmpPointsPerTrack =
                                                 tmpBundlesData.pointsPerTrack ;
        int nbFibersBundleTest = tmpPointsPerTrack.size() ;


        std::ostringstream bundlesTestMultiLabelsOss ;
        bundlesTestMultiLabelsOss << testDir << bundleName << "_labels.txt" ;
        std::string bundlesTestMultiLabels = bundlesTestMultiLabelsOss.str() ;

        std::vector<std::vector<int16_t>> multiLabelsBundle ;
        readingPredictedLabels( bundlesTestMultiLabels.c_str(),
                                multiLabelsBundle,
                                nbFibersBundleTest ) ;

        for ( int fiberIndex = 0 ; fiberIndex < nbFibersBundleTest ;
                                                                  fiberIndex++ )
        {

          labelsTest.push_back( multiLabelsBundle[ fiberIndex ] ) ;

          fullTestPointsPerTrack.push_back( tmpPointsPerTrack[ fiberIndex ] ) ;

          for ( int _point = 0 ; _point < nbPointsPerFiber ; _point++ )
          {

            int offset = 3 * nbPointsPerFiber * fiberIndex ;
            for ( int coord = 0 ; coord < 3 ; coord++ )
            {

              fullTestMatrixTracks.push_back( tmpMatrixTracks[ 3 * _point +
                                                            coord + offset ] ) ;

            }

          }

        }

      }

      curves_countFullTest = fullTestPointsPerTrack.size() ;

      BundlesMinf fullTestBundles( tmpBundles ) ;
      fullTestBundles.curves_count = curves_countFullTest ;
      // fullTestBundles.bundles = "[ '255', 0 ]" ;
      fullTestBundles.write( fullTestBundlesPath.c_str() ) ;


      BundlesData fullTestdBundlesData( tmpBundlesData ) ;
      fullTestdBundlesData.matrixTracks = fullTestMatrixTracks ;
      fullTestdBundlesData.pointsPerTrack = fullTestPointsPerTrack ;
      fullTestdBundlesData.curves_count = curves_countFullTest ;
      fullTestdBundlesData.write( fullTestdBundlesDataPath.c_str(),
                                                             fullTestBundles ) ;


      saveLabels( fullTestLabelsPath.c_str(), labelsTest ) ;
      saveLabelsDict( fullTestLabelsDictPath.c_str(), dictTest ) ;

    }

    if ( curves_countFullTest == -1 )
    {

      std::cout << "ERROR : could not get the number of fiber in the full "
                << "atlas test " << std::endl ;
      exit( 1 ) ;

    }

    // continue ;


    ////////////////////////////////////////////////////////////////////////////
    /////////////////////////// Analyse train atlas ////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Analysing train set " << std::endl ;

    }
    std::ostringstream analyseAtlasCommandOss ;
    analyseAtlasCommandOss << "python3 " << analyseAtlasBundleFile << " "
                           << "-i " << trainDir << " "
                           << "-r " << referenceFilename << " "
                           << "-force "
                           << "-saveFig false "
                           << "-v 0" ;
    std::string analyseAtlasCommand = analyseAtlasCommandOss.str() ;
    std::ostringstream analysisDirectoryPathOss ;
    analysisDirectoryPathOss << trainDir << "Analysis/" ;
    std::string analysisDirectoryPath = analysisDirectoryPathOss.str() ;
    if ( !force && is_dir( analysisDirectoryPath ) )
    {

      std::cout << "\tWARNING : analysis directory : "
                << analysisDirectoryPath << " exists, it seems atlas analysis "
                << "has been done and -force flag was not use, skipping "
                << "atlas analysis step " << std::endl ;

    }
    else
    {

      run_sh_process( analyseAtlasCommand ) ;

    }



    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////////// Projection ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Applying method to test set " << std::endl ;

    }
    std::ostringstream resultsTestDir_Oss ;
    resultsTestDir_Oss << foldDirectory << "resultsTest" << "/" ;
    std::string resultsTestDir = resultsTestDir_Oss.str() ;

    std::ostringstream projectionCommandOss ;
    projectionCommandOss << "ProjectAtlasGeoLab "
                         << "-i " << fullTestBundlesPath << " "
                         << "-a " << trainDir << " "
                         << "-fa " << fullAtlasBundlesFilename << " "
                         << "-ref " << referenceFilename << " "
                         << "-o " << resultsTestDir << " " ;
    if ( index_cc )
    {

      projectionCommandOss << "-cc " << computeCentroidsClientFilename << " " ;

    }

    if ( index_rb )
    {

      projectionCommandOss << "-rb " << registerBundlesClientFile << " " ;

    }

    if ( index_ods )
    {

      projectionCommandOss << "-ods " << openDipyServerClientFile << " " ;

    }

    if ( index_cds )
    {

      projectionCommandOss << "-cds " << closeDipyServerClientFile << " " ;

    }
    projectionCommandOss << "-nbPoints " << nbPointsPerFiber << " "
                         << "-minNbFibers 1 "
                         // << "-thrSim 0.00000001 "
                         << "-thrSim 0.3 "
                         << "-sp true "
                         << "-v 1 " ;
    std::string projectionCommand = projectionCommandOss.str() ;

    std::ostringstream regroupedRecognizedBundleOss ;
    regroupedRecognizedBundleOss << resultsTestDir
                                              << "regroupedRecognized.bundles" ;
    std::string regroupedRecognizedBundle = regroupedRecognizedBundleOss.str() ;

    if ( is_file( regroupedRecognizedBundle ) )
    {

      std::cout << "\tWARNING : file : " << regroupedRecognizedBundle
                << " already exists, skipping applying method to test set "
                << std::endl ;

    }
    else
    {

      run_sh_process( projectionCommand ) ;

    }



    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////// Building confusion matrix ////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Building confusion matrix " << std::endl ;

    }

    if ( verbose > 1 )
    {

      std::cout << "\tReading recognized labels" <<  std::endl ;

    }

    std::ostringstream labelsRecognizedTestPathOss ;
    labelsRecognizedTestPathOss << resultsTestDir << "labels.txt" ;
    std::string labelsRecognizedTestPath = labelsRecognizedTestPathOss.str() ;
    std::vector<std::vector<int16_t>> labelsRecognizedTest ;
    readingPredictedLabels( labelsRecognizedTestPath.c_str(),
                            labelsRecognizedTest,
                            curves_countFullTest ) ;


    if ( verbose > 1 )
    {

      std::cout << "\tReading recognized labels dictionary" <<  std::endl ;

    }

    std::ostringstream dictRecognizedTestPathOss ;
    dictRecognizedTestPathOss << resultsTestDir << "labels.dict" ;
    std::string dictRecognizedTestPath = dictRecognizedTestPathOss.str() ;
    std::vector<std::string> dictRecognizedTest ;
    readLabelsDict( dictRecognizedTestPath.c_str(),
                    dictRecognizedTest ) ;


    if ( verbose > 1 )
    {

      std::cout << "\tReading actual labels" <<  std::endl ;

    }
    std::vector<std::vector<int16_t>> labelsFullTest ;
    readingPredictedLabels( fullTestLabelsPath.c_str(),
                            labelsFullTest,
                            curves_countFullTest ) ;


    if ( verbose > 1 )
    {

      std::cout << "\tReading actual dictionary labels" << std::endl ;

    }
    std::vector<std::string> dictFullTest ;
    readLabelsDict( fullTestLabelsDictPath.c_str(),
                    dictFullTest ) ;


    int nbBundlesTest = dictFullTest.size() ;
    std::vector<std::vector<int32_t>> confusionMatrix ;
    int sizeConfusionMatrix = nbBundlesTest + 1 ; // Because theres a class for
                                                  // Unlabeled (last)
    // int sizeConfusionMatrix = nbBundlesTest ;
    confusionMatrix.resize( sizeConfusionMatrix ) ;
    for ( int lineIndex = 0 ; lineIndex < sizeConfusionMatrix ; lineIndex++ )
    {

      confusionMatrix[ lineIndex ].resize( sizeConfusionMatrix, 0 ) ;

    }

    // We are building the confusion matrix using the predicted dictionary
    // bundles names order
    if ( verbose > 1 )
    {

      std::cout << "\tBuilding confusion matrix" << std::endl ;

    }

    for ( int fiberIndex = 0 ; fiberIndex < curves_countFullTest ;
                                                                  fiberIndex++ )
    {

      for ( int _tmpIndex = 0 ;
                         _tmpIndex < labelsRecognizedTest[ fiberIndex ].size() ;
                                                                   _tmpIndex++ )
      {

        int _predictedLabel = labelsRecognizedTest[ fiberIndex ][ _tmpIndex ] ;

        int columnIndex ;
        if ( _predictedLabel == -1 )
        {

          columnIndex = nbBundlesTest ;

        }
        else
        {

          columnIndex = _predictedLabel ;

        }


        int lineIndex ;
        for ( int _tmpIndex2 = 0 ;
                              _tmpIndex2 < labelsFullTest[ fiberIndex ].size() ;
                                                                  _tmpIndex2++ )
        {

          int _trueLabel = labelsFullTest[ fiberIndex ][ _tmpIndex2 ] ;
          if ( _trueLabel == -1 )
          {

            std::cout << "ERROR : in labelsFullTest, found label with value -1"
                      << std::endl ;
            exit( 1 ) ;
          }

          std::string _trueLabelBundleName = dictFullTest[ _trueLabel ] ;

          lineIndex = -1 ;
          for ( int _predictedDictIndex = 0 ;
                                           _predictedDictIndex < nbBundlesTest ;
                                                         _predictedDictIndex++ )
          {

            std::string _tmpBundleName = dictRecognizedTest[
                                                         _predictedDictIndex ] ;

            if ( _tmpBundleName == _trueLabelBundleName )
            {

              lineIndex = _predictedDictIndex  ;
              break ;

            }

          }
          if ( lineIndex == -1 )
          {

            std::cout << "ERROR : while building confusion matrix, could not "
                      << "get line index " << std::endl ;
            exit( 1 ) ;

          }

          if ( lineIndex == columnIndex )
          {

            break ;

          }

        }

        if ( lineIndex != columnIndex )
        {

          int _trueLabel = labelsFullTest[ fiberIndex ][ 0 ] ;
          if ( _trueLabel == -1 )
          {

            std::cout << "ERROR : in labelsFullTest, found label with value -1"
                      << std::endl ;
            exit( 1 ) ;
          }

          std::string _trueLabelBundleName = dictFullTest[ _trueLabel ] ;

          lineIndex = -1 ;
          for ( int _predictedDictIndex = 0 ;
                                           _predictedDictIndex < nbBundlesTest ;
                                                         _predictedDictIndex++ )
          {

            std::string _tmpBundleName = dictRecognizedTest[
                                                           _predictedDictIndex ] ;

            if ( _tmpBundleName == _trueLabelBundleName )
            {

              lineIndex = _predictedDictIndex  ;
              break ;

            }

          }
          if ( lineIndex == -1 )
          {

            std::cout << "ERROR : while building confusion matrix, could not "
                      << "get line index " << std::endl ;
            exit( 1 ) ;

          }

        }

        confusionMatrix[ lineIndex ][ columnIndex ] += 1 ;

      }

    }



    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////// Saving confusion matrix //////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    std::ostringstream confusionMatrixPathOss ;
    confusionMatrixPathOss << resultsTestDir << "confusionMatrix.tsv" ;
    std::string confusionMatrixPath = confusionMatrixPathOss.str() ;
    saveConfusionMatrix( confusionMatrixPath.c_str(),
                         confusionMatrix ) ;


    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////// Computing accuracy ////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Computing accuracy " << std::endl ;

    }
    float  sumDiagonal = 0 ;
    float  sumAllElements = 0 ;
    for ( int lineIndex = 0 ; lineIndex < sizeConfusionMatrix ; lineIndex++ )
    {

      for ( int columnIndex = 0 ; columnIndex < sizeConfusionMatrix ;
                                                                 columnIndex++ )
      {

        if ( columnIndex == lineIndex )
        {

          sumDiagonal += ( float )confusionMatrix[ lineIndex ][ columnIndex ] ;

        }

        sumAllElements += ( float )confusionMatrix[ lineIndex ][ columnIndex ] ;

      }

    }

    float accuracyTest = sumDiagonal / sumAllElements ;
    // float accuracyTest = sumDiagonal / curves_countFullTest ;

    std::cout << "Accuracy = " << accuracyTest << std::endl ;

    ////////////////////////////////////////////////////////////////////////////
    /////////////////////// Computing weights per bundle ///////////////////////
    ////////////////////////////////////////////////////////////////////////////
    AtlasBundles trainAtlasData( trainDir.c_str(),
                                true,
                                false,
                                0 ) ;

    int nbBundlesTrain = trainAtlasData.bundles.size() ;
    int64_t totalNbFibersTrain = 0 ;
    std::vector<float> weightsBundles( nbBundlesTrain, 0 ) ;
    // Same order as recognized bundles
    for ( int _bundleTrainIndex = 0 ; _bundleTrainIndex < nbBundlesTrain ;
                                                           _bundleTrainIndex++ )
    {

      std::string _bundleTrainName = trainAtlasData.bundlesNames[
                                                           _bundleTrainIndex ] ;
      int _indexInRecognized = -1 ;
      int _nbBundlesInRecognized = dictRecognizedTest.size() ;
      for ( int _recognizedIndex = 0 ;
                                     _recognizedIndex < _nbBundlesInRecognized ;
                                                            _recognizedIndex++ )
      {

        if ( dictRecognizedTest[ _recognizedIndex ] == _bundleTrainName )
        {

          _indexInRecognized = _recognizedIndex ;
          break ;

        }

      }

      if ( _indexInRecognized == -1 )
      {

        std::cout << "ERROR : while computing weights could not find bundle "
                  << _bundleTrainName << " in dictionary with names of "
                  << "recognized bundles " << std::endl ;
        exit( 1 ) ;

      }

      BundlesMinf& _bundleTrain = trainAtlasData.bundles[ _bundleTrainIndex ] ;
      int tmpNbFibers = _bundleTrain.curves_count ;
      totalNbFibersTrain += tmpNbFibers ;
      weightsBundles[ _indexInRecognized ] = ( float )tmpNbFibers ;

    }

    float sumWeights = 0 ;
    for ( int _bundleTrainIndex = 0 ; _bundleTrainIndex < nbBundlesTrain ;
                                                           _bundleTrainIndex++ )
    {

      weightsBundles[ _bundleTrainIndex ] /= ( float )totalNbFibersTrain ;
      sumWeights += weightsBundles[ _bundleTrainIndex ] ;

    }

    if ( verbose > 1 )
    {

      std::cout << "\tSum weights : " << sumWeights << std::endl ;
    }




    ////////////////////////////////////////////////////////////////////////////
    ///////////////////// Computing scores for each bundle /////////////////////
    ////////////////////////////////////////////////////////////////////////////
    std::vector<int> tpPerBundle( sizeConfusionMatrix, 0 ) ;
    std::vector<int> tnPerBundle( sizeConfusionMatrix, 0 ) ;
    std::vector<int> fpPerBundle( sizeConfusionMatrix, 0 ) ;
    std::vector<int> fnPerBundle( sizeConfusionMatrix, 0 ) ;
    for ( int _bundleIndex = 0 ; _bundleIndex < nbBundlesTest ; _bundleIndex++ )
    {

      for ( int lineIndex = 0 ; lineIndex < sizeConfusionMatrix ; lineIndex++ )
      {

        for ( int columnIndex = 0 ; columnIndex < sizeConfusionMatrix ;
                                                                 columnIndex++ )
        {

          if ( columnIndex == lineIndex )
          {

            if ( _bundleIndex == lineIndex )
            {

              tpPerBundle[ _bundleIndex ] += confusionMatrix[ lineIndex ][
                                                                 columnIndex ] ;

            }
            else
            {

              tnPerBundle[ _bundleIndex ] += confusionMatrix[ lineIndex ][
                                                                 columnIndex ] ;


            }

          }
          else
          {

            if ( lineIndex == _bundleIndex )
            {


              fnPerBundle[ _bundleIndex ] += confusionMatrix[ lineIndex ][
                                                                 columnIndex ] ;

            }

            if ( columnIndex == _bundleIndex )
            {

              fpPerBundle[ _bundleIndex ] += confusionMatrix[ lineIndex ][
                                                                 columnIndex ] ;


            }

          }

        }

      }

    }

    std::vector<float> precisionPerBundle( sizeConfusionMatrix, 0 ) ;
    std::vector<float> recallPerBundle( sizeConfusionMatrix, 0 ) ;
    std::vector<float> accuracyPerBundle( sizeConfusionMatrix, 0 ) ;
    float averagePrecision = 0 ;
    float averageRecall = 0 ;
    float averageAccuracy = 0 ;
    for ( int _bundleIndex = 0 ; _bundleIndex < nbBundlesTest ; _bundleIndex++ )
    {

      if ( tpPerBundle[ _bundleIndex ] > 0 )
      {

        recallPerBundle[ _bundleIndex ] = ( float )tpPerBundle[ _bundleIndex ] /
                 ( tpPerBundle[ _bundleIndex ] + fpPerBundle[ _bundleIndex ] ) ;

        precisionPerBundle[ _bundleIndex ] =
                                          ( float )tpPerBundle[ _bundleIndex ] /
                 ( tpPerBundle[ _bundleIndex ] + fnPerBundle[ _bundleIndex ] ) ;



      }
      else
      {

        precisionPerBundle[ _bundleIndex ] = 0 ;

        recallPerBundle[ _bundleIndex ] = 0 ;

      }

      if ( tpPerBundle[ _bundleIndex ] > 0 && tnPerBundle[ _bundleIndex ] > 0 )
      {

        accuracyPerBundle[ _bundleIndex ] =
                 ( float )
                 ( tpPerBundle[ _bundleIndex ] + tnPerBundle[ _bundleIndex ] ) /
                 ( tpPerBundle[ _bundleIndex ] + tnPerBundle[ _bundleIndex ] +
                   fpPerBundle[ _bundleIndex ] + fnPerBundle[ _bundleIndex ] ) ;

      }
      else
      {

        accuracyPerBundle[ _bundleIndex ] = 0 ;

      }

      averagePrecision += precisionPerBundle[ _bundleIndex ] ;
      averageRecall += recallPerBundle[ _bundleIndex ] ;
      averageAccuracy += accuracyPerBundle[ _bundleIndex ] ;

    }

    averagePrecision /= nbBundlesTest ;
    averageRecall /= nbBundlesTest ;
    averageAccuracy /= nbBundlesTest ;

    std::ostringstream scoresPerBundlePathOss ;
    scoresPerBundlePathOss << resultsTestDir << "scoresPerBundle.tsv" ;
    std::string scoresPerBundlePath = scoresPerBundlePathOss.str() ;
    saveScorePerBundle( scoresPerBundlePath.c_str(),
                        precisionPerBundle,
                        recallPerBundle,
                        accuracyPerBundle,
                        weightsBundles,
                        dictRecognizedTest,
                        averagePrecision,
                        averageRecall,
                        averageAccuracy ) ;



    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////// Computing macro F1-score /////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Computing macro F1-score " << std::endl ;

    }

    // Getting sum of True Positives for each class
    if ( verbose )
    {

      std::cout << "\tGetting true positives " << std::endl ;

    }

    std::vector<int> truePositives( sizeConfusionMatrix, 0 ) ;
    for ( int _bundleTestIndex = 0 ; _bundleTestIndex < sizeConfusionMatrix ;
                                                            _bundleTestIndex++ )
    {

      truePositives[ _bundleTestIndex ] =
                       confusionMatrix[ _bundleTestIndex ][ _bundleTestIndex ] ;

    }

    // Getting sum of False positives for each class
    if ( verbose )
    {

      std::cout << "\tGetting false positives " << std::endl ;

    }

    std::vector<int> falsePositives( sizeConfusionMatrix, 0 ) ;
    for ( int lineIndex = 0 ; lineIndex < sizeConfusionMatrix ; lineIndex++ )
    {

      for ( int columnIndex = 0 ; columnIndex < sizeConfusionMatrix ;
                                                                 columnIndex++ )
      {

        if ( columnIndex != lineIndex )
        {

          falsePositives[ columnIndex ] +=
                                    confusionMatrix[ lineIndex ][columnIndex ] ;

        }

      }

    }

    // Getting sum of False Negatives for each class
    if ( verbose )
    {

      std::cout << "\tGetting false negatives " << std::endl ;

    }

    std::vector<int> falseNegatives( sizeConfusionMatrix, 0 ) ;
    for ( int lineIndex = 0 ; lineIndex < sizeConfusionMatrix ; lineIndex++ )
    {

      for ( int columnIndex = 0 ; columnIndex < sizeConfusionMatrix ;
                                                                 columnIndex++ )
      {

        if ( columnIndex != lineIndex )
        {

          falseNegatives[ lineIndex ] +=
                                    confusionMatrix[ lineIndex ][columnIndex ] ;

        }

      }

    }

    // Getting sum of True Negatives for each class
    if ( verbose )
    {

      std::cout << "\tGetting true negatives " << std::endl ;

    }

    std::vector<int> trueNegatives( sizeConfusionMatrix, 0 ) ;
    for ( int _bundleTestIndex = 0 ; _bundleTestIndex < sizeConfusionMatrix ;
                                                            _bundleTestIndex++ )
    {

      trueNegatives[ _bundleTestIndex ] = sumAllElements -
                                          truePositives[ _bundleTestIndex ] -
                                          falsePositives[ _bundleTestIndex ] -
                                          falseNegatives[ _bundleTestIndex ] ;

    }


    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //
    ////////////////////////////////////////////////////////////////////////////
    /////////////////////// Computing balanced accuracy ////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Computing balanced accuracy " << std::endl ;

    }
    // Getting sum of rows for each class
    std::vector<int> sumRows( sizeConfusionMatrix, 0 ) ;
    if ( verbose )
    {

      std::cout << "\tGetting sum of rows " << std::endl ;

    }

    for ( int lineIndex = 0 ; lineIndex < sizeConfusionMatrix ; lineIndex++ )
    {

      for ( int columnIndex = 0 ; columnIndex < sizeConfusionMatrix ;
                                                                 columnIndex++ )
      {

        sumRows[ lineIndex ] +=  confusionMatrix[ lineIndex ][ columnIndex ] ;

      }

    }

    if ( verbose )
    {

      std::cout << "\tComputing balanced accuracy " << std::endl ;

    }

    float balancedAccuracy = 0 ;
    int nbExtractedBundles = 0 ;
    for ( int _bundleTestIndex = 0 ; _bundleTestIndex < sizeConfusionMatrix ;
                                                            _bundleTestIndex++ )
    {

      if ( sumRows[ _bundleTestIndex ] == 0 )
      {

        continue ;

      }
      nbExtractedBundles++ ;

      balancedAccuracy += ( ( float )truePositives[ _bundleTestIndex ] /
                                       ( float )sumRows[  _bundleTestIndex ] ) ;

    }
    // balancedAccuracy /= ( float )nbExtractedBundles ;
    balancedAccuracy /= ( float )sizeConfusionMatrix ;
    std::cout << "Balanced Accuracy : " << balancedAccuracy << std::endl ;


    ////////////////////////////////////////////////////////////////////////////
    ////////////////// Computing balanced weighted accuracy ////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if ( verbose )
    {

      std::cout << "Computing weighted balanced accuracy " << std::endl ;

    }
    float balancedWeightedAccuracy = 0 ;
    for ( int _bundleTestIndex = 0 ; _bundleTestIndex < sizeConfusionMatrix ;
                                                            _bundleTestIndex++ )
    {

      if ( sumRows[ _bundleTestIndex ] == 0 )
      {

        continue ;

      }

      balancedWeightedAccuracy += weightsBundles[ _bundleTestIndex ] * (
                                    ( float )truePositives[ _bundleTestIndex ] /
                                       ( float )sumRows[  _bundleTestIndex ] ) ;

    }
    // balancedWeightedAccuracy /= ( ( float )sizeConfusionMatrix * sumWeights )  ;
    balancedWeightedAccuracy /= ( ( float )sumWeights )  ;
    std::cout << "Weihgted Balanced Accuracy : " << balancedWeightedAccuracy
                                                                  << std::endl ;
    exit( 1 ) ;
    continue ;
    // exit( 1 ) ;


    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------------------- //

    // // Checking values
    // for ( int _bundleTestIndex = 0 ; _bundleTestIndex < sizeConfusionMatrix ;
    //                                                         _bundleTestIndex++ )
    // {
    //
    //   std::string _bundleName = dictFullTest[ _bundleTestIndex ] ;
    //
    //   float _tmpValue1 = truePositives[ _bundleTestIndex ] +
    //                                         falsePositives[ _bundleTestIndex ] ;
    //   float _tmpValue2 = truePositives[ _bundleTestIndex ] +
    //                                         falseNegatives[ _bundleTestIndex ] ;
    //
    //   if ( _tmpValue1 == 0 )
    //   {
    //
    //     std::cout << "ERROR : for bundle " << _bundleName << " the sum of "
    //               << "TP and FP is zero " << std::endl ;
    //     exit( 1 ) ;
    //
    //   }
    //
    //
    //   if ( _tmpValue2 == 0 )
    //   {
    //
    //     std::cout << "ERROR : for bundle " << _bundleName << " the sum of "
    //               << "TP and FN is zero " << std::endl ;
    //     exit( 1 ) ;
    //
    //   }
    //
    //
    // }



    // Computing macro precision and macro recall
    if ( verbose )
    {

      std::cout << "\tComputing macro precision and recall " << std::endl ;

    }

    float macroPrecisionTest = 0 ;
    float macroRecallTest = 0 ;
    for ( int _bundleTestIndex = 0 ; _bundleTestIndex < sizeConfusionMatrix ;
                                                            _bundleTestIndex++ )
    {

      float precisionBundle = 0 ;
      float recallBundle =0 ;
      if ( truePositives[ _bundleTestIndex ] == 0 )
      {


        precisionBundle = 0 ;

        recallBundle = 0 ;

      }
      else
      {

        precisionBundle = truePositives[ _bundleTestIndex ] /
                                ( truePositives[ _bundleTestIndex ] +
                                          falsePositives[ _bundleTestIndex ] ) ;

        recallBundle = truePositives[ _bundleTestIndex ] /
                                ( truePositives[ _bundleTestIndex ] +
                                          falseNegatives[ _bundleTestIndex ] ) ;

      }


      macroPrecisionTest += precisionBundle ;
      macroRecallTest += recallBundle ;

    }
    macroPrecisionTest /= sizeConfusionMatrix ;
    macroRecallTest /= sizeConfusionMatrix ;

    if ( macroRecallTest == 0 )
    {

      std::cout << "ERROR : the value of macro recall is 0 " << std::endl ;
      exit( 1 ) ;

    }

    if ( macroPrecisionTest == 0 )
    {

      std::cout << "ERROR : the value of macro precision is 0 " << std::endl ;
      exit( 1 ) ;

    }

    // Computing macro F1-score
    if ( verbose )
    {

      std::cout << "\tComputing macro F1-score " << std::endl ;

    }
    float macroF1ScoreTest = 2 * ( macroPrecisionTest * macroRecallTest ) /
                              ( 1 / macroPrecisionTest + 1 / macroRecallTest ) ;


    // Saving in vector
    accuraciesKfold[ crossValidationIndex ] = accuracyTest ;
    f1ScoreKfold[ crossValidationIndex ] = macroF1ScoreTest ;

    if ( verbose )
    {

      std::cout << "Accuracy : " << accuracyTest << std::endl ;
      std::cout << "F1-score : " << macroF1ScoreTest << std::endl ;

    }

    exit( 1 ) ;


  }

  if ( verbose )
  {

    std::cout << "\nDone " << std::endl ;

  }

  std::cout << "Accuracies : " << std::endl ;
  for ( int _index = 0 ; _index < k_value ; _index++ )
  {

    std::cout << accuraciesKfold[ _index ] << std::endl ;

  }

  std::cout << "Macro F1-scores : " << std::endl ;
  for ( int _index = 0 ; _index < k_value ; _index++ )
  {

    std::cout << f1ScoreKfold[ _index ] << std::endl ;

  }




}
