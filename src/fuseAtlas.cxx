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

#include "fuseAtlas.h"
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

  int index_input, index_output, index_format, index_nbPoints, index_force,
                                                     index_verbose, index_help ;

  index_input = getFlagPosition( argc, argv, "-i") ;
  index_output = getFlagPosition( argc, argv, "-o") ;
  index_format = getFlagPosition( argc, argv, "-f") ;
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints") ;
  index_force = getFlagPosition( argc, argv, "-force") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-i : Path to the input directory containing the bundles \n"
              << "-o : Path to the output directory \n"
              << "-f : bundles format among .bundles/.trk/.tck \n"
              << "[-nbPoints] : Number of points per fiber (same number for all "
              << "fibers) \n"
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

  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
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

  ////////////////////////////// Input Directory ///////////////////////////////
  inputDirectoryPath = argv[ index_input + 1 ] ;
  char lastChar = inputDirectoryPath[ inputDirectoryPath.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    inputDirectoryPath = inputDirectoryPath + "/" ;

  }

  if ( !is_dir( inputDirectoryPath ) )
  {

    std::cout << "ERROR : input directory path " << inputDirectoryPath
              << " does not exists" << std::endl ;
    exit( 1 ) ;

  }

  ////////////////////////////// Output directory //////////////////////////////
  outputDirectoryPath = argv[ index_output + 1 ] ;
  lastChar = outputDirectoryPath[ outputDirectoryPath.size() - 1 ] ;
  if ( lastChar != '/' )
  {

    outputDirectoryPath = outputDirectoryPath + "/" ;


  }

  if ( !is_dir( outputDirectoryPath ) )
  {

    mkdir( outputDirectoryPath ) ;

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



  ///////////////////////// Number of points per fiber /////////////////////////
  if ( index_nbPoints )
  {

    nbPointsPerFiber = std::stoi( argv[ index_nbPoints + 1 ] ) ;

    if ( nbPointsPerFiber <= 0 )
    {

      std::cout << "ERROR : -nbPoints must be greater than 0 " << std::endl ;
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
  // Reading Atlas
  AtlasBundles atlasData( inputDirectoryPath.c_str(),
                          isBundlesFormat,
                          isTrkFormat,
                          isTckFormat,
                          0 ) ;



  // Building fused atlas labels dictionary
  if ( verbose )
  {

    std::cout << "\n--------------------------------------------" << std::endl ;
    std::cout << "Building fused atlas labels dictionary... " ;
    std::cout << "--------------------------------------------" << std::endl ;

  }
  std::vector<std::string> labelsDictFusedAtlas ;
  int64_t nbFibersFusedAtlas = 0 ;
  int nbBundles = atlasData.bundlesNames.size() ;
  for ( int bundleIndex = 0 ; bundleIndex < nbBundles ; bundleIndex++ )
  {

    std::string& bundleName = atlasData.bundlesNames[ bundleIndex ] ;

    nbFibersFusedAtlas += atlasData.bundlesData[ bundleIndex ].curves_count ;

    // Reading dictionary
    std::ostringstream labelsDictBundlePathOss ;
    labelsDictBundlePathOss << inputDirectoryPath << bundleName
                                                             << "_labels.dict" ;
    std::string labelsDictBundlePath = labelsDictBundlePathOss.str() ;
    std::vector<std::string> labelsDict ;
    if ( !is_file( labelsDictBundlePath ) )
    {

      labelsDict.push_back( bundleName ) ;
      std::cout << "\nWARNING : labels dictionary file for bundle "
                << bundleName << " not found in " << labelsDictBundlePath
                << ", creating new label : " << bundleName << std::endl ;

    }
    else
    {

      readLabelsDict( labelsDictBundlePath.c_str(), labelsDict ) ;

    }

    for ( int labelIndex = 0 ; labelIndex < labelsDict.size() ; labelIndex++ )
    {

      if ( std::find( labelsDictFusedAtlas.begin(), labelsDictFusedAtlas.end(),
                                                  labelsDict[ labelIndex ] ) ==
                                                    labelsDictFusedAtlas.end() )
      {

        labelsDictFusedAtlas.push_back( labelsDict[ labelIndex ] ) ;

      }

    }

  }
  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }

  if ( verbose )
  {

    std::cout << "\n--------------------------------------------" << std::endl ;
    std::cout << "Fusing atlas" << std::endl ;
    std::cout << "--------------------------------------------" << std::endl ;

  }
  std::vector<float> fusedAtlasMatrixTracks ;
  std::vector<int32_t> fusedAtlasPointsPerTrack ;
  std::vector<std::vector<int16_t>> labelsFusedAtlas( nbFibersFusedAtlas ) ;
  int64_t fiberIndexInFusedAtlas = 0 ;
  for ( int bundleIndex = 0 ; bundleIndex < nbBundles ; bundleIndex++ )
  {

    std::string& bundleName = atlasData.bundlesNames[ bundleIndex ] ;

    BundlesData& atlasBundleData = atlasData.bundlesData[ bundleIndex ] ;

    int nbCurvesBundle = atlasBundleData.curves_count ;

    // Reading labels
    std::ostringstream labelsBundlePathOss ;
    labelsBundlePathOss << inputDirectoryPath << bundleName << "_labels.txt" ;
    std::string labelsBundlePath = labelsBundlePathOss.str() ;
    std::vector<std::vector<int16_t>> labelsBundle ;
    if ( !is_file( labelsBundlePath ) )
    {

      std::cout << "\nWARNING : labels file for bundle " << bundleName
                << " not found in " << labelsBundlePath
                << ", creting new label" << std::endl ;


      int tmpLabel = -2 ;
      for ( int i = 0 ; i < labelsDictFusedAtlas.size() ; i++ )
      {

        if ( labelsDictFusedAtlas[ i ] == bundleName )
        {

          tmpLabel = i ;
          break ;

        }

      }
      if ( tmpLabel == -2 )
      {

        std::cout << "ERROR : could not find bundle " << bundleName
                  << " in labelsDictFusedAtlas" << std::endl ;
        exit( 1 ) ;

      }

      labelsBundle.resize( nbCurvesBundle,
                                         std::vector<int16_t>( 1, tmpLabel ) ) ;

    }
    else
    {

      readingPredictedLabels( labelsBundlePath.c_str(),
                              labelsBundle,
                              nbCurvesBundle ) ;

    }



    // Reading dictionary
    std::ostringstream labelsDictBundlePathOss ;
    labelsDictBundlePathOss << inputDirectoryPath << bundleName
                                                             << "_labels.dict" ;
    std::string labelsDictBundlePath = labelsDictBundlePathOss.str() ;

    std::vector<std::string> labelsDict ;
    bool isBundleLabelsDict = false ;
    if ( !is_file( labelsDictBundlePath ) )
    {

      labelsDict.push_back( bundleName ) ;
      std::cout << "\nWARNING : labels dictionary file for bundle "
                << bundleName << " not found in " << labelsDictBundlePath
                << ", creating new label : " << bundleName << std::endl ;

    }
    else
    {

      readLabelsDict( labelsDictBundlePath.c_str(), labelsDict ) ;
      isBundleLabelsDict = true ;

    }




    for ( int fiberIndex = 0 ; fiberIndex < nbCurvesBundle ; fiberIndex++ )
    {

      if ( fiberIndex%1000 == 0 || fiberIndex == nbCurvesBundle - 1 )
      {

        printf( "\rBundle : [%d/%d]\t|\tFiber : [%d/%d]\t|\t%s",
                                     bundleIndex + 1, nbBundles, fiberIndex + 1,
                                     nbCurvesBundle, bundleName.c_str() ) ;
        std::cout << "" << std::flush ;

      }


      if ( labelsDictFusedAtlas.size() > 0 )
      {

        std::vector<int16_t>& labelsFiberInBundle = labelsBundle[ fiberIndex ] ;
        for ( int _tmpIndex = 0 ; _tmpIndex < labelsFiberInBundle.size() ;
                                                                   _tmpIndex++ )
        {

          std::string labelName ;
          if ( labelsFiberInBundle[ _tmpIndex ] == -1 )
          {

            // labelName = "unlabeledFibers" ;
            labelsFusedAtlas[ fiberIndexInFusedAtlas ].push_back( -1 ) ;
            continue ;

          }
          else
          {

            if ( isBundleLabelsDict )
            {

              labelName = labelsDict[ labelsFiberInBundle[ _tmpIndex ] ] ;

            }
            else
            {

              labelName = labelsDict[ 0 ] ;


            }


          }

          int16_t labelInFusedAtlas = -2 ;
          for ( int labelIndex = 0 ; labelIndex < labelsDictFusedAtlas.size() ;
                                                                  labelIndex++ )
          {

            if ( labelsDictFusedAtlas[ labelIndex ] == labelName )
            {

              labelInFusedAtlas = labelIndex ;
              break ;

            }

          }

          if ( labelInFusedAtlas == -2 )
          {

            std::cout << "\nWARNING : could not find bundle " << labelName
                      << " in labelsDictFusedAtlas  " << std::endl ; ;

          }

          labelsFusedAtlas[ fiberIndexInFusedAtlas ].push_back(
                                                           labelInFusedAtlas ) ;

        }

      }

      fusedAtlasPointsPerTrack.push_back( atlasBundleData.pointsPerTrack[
                                                                fiberIndex ] ) ;

      int tmpOffset = 3 * nbPointsPerFiber * fiberIndex  ;
      for ( int point = 0 ; point < nbPointsPerFiber ; point++ )
      {

        for ( int coord = 0 ; coord < 3 ; coord++ )
        {

          fusedAtlasMatrixTracks.push_back( atlasBundleData.matrixTracks[
                                             tmpOffset + 3 * point + coord ] ) ;

        }

      }

      fiberIndexInFusedAtlas++ ;

    }

  }


  // Saving labels
  bool isLabelsFusedAtlas = true ;
  for ( std::vector<int16_t> _labels : labelsFusedAtlas )
  {

    if ( _labels.size() == 0 )
    {

      isLabelsFusedAtlas = false ;
      break ;

    }

  }
  if ( verbose && isLabelsFusedAtlas )
  {

    std::cout << "\nSaving labels fused atlas... " ;

  }

  std::ostringstream outLabelsBundlePathOss ;
  outLabelsBundlePathOss << outputDirectoryPath << "labels.txt" ;
  std::string outLabelsBundlePath = outLabelsBundlePathOss.str() ;

  if ( isLabelsFusedAtlas )
  {

    saveLabels( outLabelsBundlePath.c_str(),
                labelsFusedAtlas ) ;

  }



  if ( verbose && isLabelsFusedAtlas )
  {

    std::cout << "Done" << std::endl ;

  }


  // Saving dictionary
  if ( verbose && labelsDictFusedAtlas.size() > 0 )
  {

    std::cout << "Saving labels dictionary fused atlas... " ;

  }

  std::ostringstream outLabelsDictBundlePathOss ;
  outLabelsDictBundlePathOss << outputDirectoryPath << "labels.dict" ;
  std::string outLabelsDictBundlePath = outLabelsDictBundlePathOss.str() ;

  if ( labelsDictFusedAtlas.size() > 0 )
  {

    saveLabelsDict( outLabelsDictBundlePath.c_str(),
                    labelsDictFusedAtlas ) ;

  }


  if ( verbose && labelsDictFusedAtlas.size() > 0 )
  {

    std::cout << "Done" << std::endl ;

  }

  // Saving .bundles/.trk/.tck files
  if ( verbose )
  {

    std::cout << "Saving fused atlas... " ;

  }

  BundlesMinf fusedAtlasInfo = atlasData.bundlesMinf[ 0 ] ;
  fusedAtlasInfo.curves_count = fusedAtlasPointsPerTrack.size() ;


  std::ostringstream outBundlesDataPathOss ;
  outBundlesDataPathOss << outputDirectoryPath << "fusedAtlas" << format ;
  std::string outBundlesDataPath = outBundlesDataPathOss.str() ;
  BundlesData fusedAtlasData = atlasData.bundlesData[ 0 ] ;
  fusedAtlasData.matrixTracks = fusedAtlasMatrixTracks ;
  fusedAtlasData.pointsPerTrack = fusedAtlasPointsPerTrack ;
  fusedAtlasData.curves_count = fusedAtlasPointsPerTrack.size() ;
  fusedAtlasData.write( outBundlesDataPath.c_str(), fusedAtlasInfo ) ;

  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }


}
