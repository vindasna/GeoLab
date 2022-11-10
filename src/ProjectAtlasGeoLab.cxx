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
//////////////////// Function to check if file exists //////////////////////////
////////////////////////////////////////////////////////////////////////////////

inline bool is_file( const std::string& path )
{

  struct stat buffer ;
  return( stat ( path.c_str(), &buffer ) == 0 ) ;

}

////////////////////////////////////////////////////////////////////////////////
////////////////// Function to check if directory exists ///////////////////////
////////////////////////////////////////////////////////////////////////////////

inline bool is_dir( const std::string& path )
{

  return( std::experimental::filesystem::is_directory( path ) ) ;

}


////////////////////////////////////////////////////////////////////////////////
/////////////////////// Function to create directory ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

inline bool mkdir( const std::string& path )
{

  return( std::experimental::filesystem::create_directory( path ) ) ;

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////// Function to delete file /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

inline bool rmfile( const std::string& path )
{

  if ( is_file( path ) )
  {

    return( std::experimental::filesystem::remove( path ) ) ;

  }

  return( false ) ;

}

////////////////////////////////////////////////////////////////////////////////
/////////////////////// Function to delete directory ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline bool rmdir( const std::string& path )
{

  if ( is_dir( path ) )
  {

    if ( std::experimental::filesystem::remove_all( path ) )
    {

      return( true ) ;

    }
    else
    {

      return( false ) ;
    }

  }
  else
  {

    return( false ) ;

  }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////// Function to copy file ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline bool copy( const std::string& source,
                  const std::string& destination )
{

  return( std::experimental::filesystem::copy_file( source, destination ) ) ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////// Function to rename file //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline bool rename( const std::string& source, const std::string& destination )
{

  try
  {

    std::experimental::filesystem::rename( source, destination ) ;
    if ( !is_file( source ) && is_file( destination ) )
    {

      return( true ) ;

    }
    return( false ) ;

  }
  catch ( ... )
  {

    return( false ) ;

  }

}


////////////////////////////////////////////////////////////////////////////////
/////////////// Function to check if string ends with sub-string ///////////////
////////////////////////////////////////////////////////////////////////////////
inline bool endswith( const std::string& input,
                      const std::string& substring )
{

  std::string tmpString ;
  char lastChar = input[ input.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    tmpString = input.substr( 0, input.size() - 1 ) ;

  }
  else
  {

    tmpString = input ;

  }

  std::string endSubstring = input.substr( input.size() - substring.size(),
                                                                input.size() ) ;

  if ( endSubstring == substring )
  {

    return( true ) ;

  }
  else
  {

    return( false ) ;

  }

}

////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to get parent directory //////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline std::string dirname( const std::string& path )
{

  std::string tmpString = path ;
  char lastChar = tmpString[ tmpString.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    tmpString = tmpString.substr( 0, tmpString.size() - 1 ) ;

  }

  std::experimental::filesystem::path p( tmpString ) ;

  std::string dirname = p.parent_path() ;
  lastChar = dirname[ dirname.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    dirname = dirname.substr( 0, dirname.size() - 1 ) ;

  }
  dirname += "/" ;
  return( dirname ) ;

}

////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to replace extension file ////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline std::string replaceExtension( const std::string& path,
                                     const std::string& newExtension )
{

  std::string tmpString = path ;
  char lastChar = tmpString[ tmpString.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    tmpString = tmpString.substr( 0, tmpString.size() - 1 ) ;

  }

  std::experimental::filesystem::path p( tmpString ) ;
  std::string filenameNoExtension = p.stem() ;
  std::ostringstream _tmpOss ;
  _tmpOss << dirname( path ) << filenameNoExtension ;
  std::string newPath = _tmpOss.str() ;
  newPath += newExtension ;

  return( newPath ) ;

}


////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to check atlas directory /////////////////////////
////////////////////////////////////////////////////////////////////////////////
int countFilesDirectory( const std::string& path )
{

  if ( !is_dir( path ) )
  {

    return( -1 ) ;

  }

  int nbFilesInDirectory = 0 ;

  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                                        path ) )
  {

    nbFilesInDirectory += 1 ;

  }

  return( nbFilesInDirectory ) ;

}



////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to count files directory /////////////////////////
////////////////////////////////////////////////////////////////////////////////
void checkAtlasDirectory( const std::string& path,
                          std::string& outputDirectory )
 // Not const for outputDirectory
{

  int nbCorrectFiles = 0 ;

  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                                        path ) )
  {

    std::string tmpBundlesDataFilename ;
    std::string tmpBundlesFilename ;

    tmpBundlesDataFilename = file.path() ;
    tmpBundlesFilename = tmpBundlesDataFilename ;
    std::string key (".bundlesdata") ;

    if ( tmpBundlesDataFilename.find( ".bundlesdata" ) != std::string::npos )
    {

      std::size_t found = tmpBundlesFilename.rfind( key ) ;
      tmpBundlesFilename.replace( found, key.length(), ".bundles" ) ;

      if ( is_file( tmpBundlesDataFilename ) && is_file( tmpBundlesFilename ) )
      {

        nbCorrectFiles += 1 ;

      }
      else
      {

        if (  is_file( tmpBundlesDataFilename ) &&
                                                !is_file( tmpBundlesFilename ) )
        {

          std::cout << "ERROR : File " << tmpBundlesDataFilename << " found "
                    << "but complementary file " << tmpBundlesFilename
                    << " not found " << std::endl ;

          if ( is_dir( outputDirectory ) )
          {

            rmdir( outputDirectory ) ;

          }
          exit( 1 ) ;

        }

        if (  !is_file( tmpBundlesDataFilename ) &&
                                                is_file( tmpBundlesFilename ) )
        {

          std::cout << "ERROR : File " << tmpBundlesFilename << " found "
                    << "but complementary file " << tmpBundlesDataFilename
                    << " not found " << std::endl ;

          if ( is_dir( outputDirectory ) )
          {

            rmdir( outputDirectory ) ;

          }
          exit( 1 ) ;

        }

      }

    }

  }

  if ( nbCorrectFiles > 0 )
  {

    std::cout << "Atlas directory : OK " << std::endl ;

  }
  else
  {

    std::cout << "ERROR : the atlas directory " << path << " does not contain "
              << "any bundles in .bundles format " << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

}


////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to launch process with timeout ///////////////////
////////////////////////////////////////////////////////////////////////////////
std::string run_sh_process_timeout( const std::string& command, int timeout )
{

  try
  {

    boost::process::ipstream out ; // To not pipe output in main process
    boost::process::ipstream err ; // To not pipe error in main process

    boost::process::group g ;
    boost::process::child c( command.c_str(), boost::process::std_out > out,
                                              boost::process::std_err > err, g ) ;
    if ( !g.wait_for( std::chrono::seconds( timeout ) ) )
    {

      g.terminate() ;

    }

    if ( c.running() )
    {

      c.wait() ;

    }


    std::string line;
    std::ostringstream outStringOss ;
    while ( std::getline( out, line ) )
    {

      outStringOss << line << "\n" ;

    }
    std::string outString = outStringOss.str() ;
    return( outString ) ;

  }
  catch ( ... )
  {

    std::ostringstream outStringOss ;
    outStringOss << "ERROR boost " ;
    std::string outString = outStringOss.str() ;
    return( outString ) ;

  }

}

////////////////////////////////////////////////////////////////////////////////
////////////////// Function to launch process without timeout //////////////////
////////////////////////////////////////////////////////////////////////////////
int run_sh_process( const std::string& command )
{

  boost::process::ipstream out ; // To not pipe output in main process
  boost::process::ipstream err ; // To not pipe error in main process

  boost::process::group g ;
  boost::process::child c( command.c_str(), boost::process::std_out > out,
                                            boost::process::std_err > err, g ) ;


  if ( verbose > 1 )
  {

    std::string line;
    while ( std::getline( out, line ) )
    {

      std::cout << line << std::endl ;

    }

  }
  else
  {

    // Using this because c.wait() bugs with long process (why?)

    std::string line;
    while ( std::getline( out, line ) )
    {

      if ( verbose > 1 )
      {

        std::ostringstream _tmpOss ;
        _tmpOss << "\rRunning command : " << command ;
        std::string _tmp = _tmpOss.str() ;
        std::cout << _tmp ;

      }

    }

  }

  try
  {

    c.wait() ;

  }
  catch ( ... ){}

  return( c.exit_code() ) ;

}



///////////////////////////////////////////////////////////////////////////////
/////////////// Function to read the predicted labels file ////////////////////
///////////////////////////////////////////////////////////////////////////////
void readingPredictedLabels( const char* predictedLabelsFilename,
                             std::vector<std::vector<int16_t>>& predictedLabels,
                             const std::string& outputDirectory,
                             int nbFibers,
                             int verbose )
{

  predictedLabels.resize( nbFibers ) ;

  const char delim = ':' ;
  std::string line ;
  std::ifstream dictFile ;
  dictFile.open( predictedLabelsFilename ) ;
  if ( dictFile.fail() )
  {

    std::cout << "Problem reading file : " << predictedLabelsFilename
                                           << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
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


///////////////////////////////////////////////////////////////////////////////
/////////////// Function to save the predicted labels file ////////////////////
///////////////////////////////////////////////////////////////////////////////
void saveLabels( const char* labelsBinaryFilename,
                 const std::vector<std::vector<int16_t>>& labels,
                 const std::string& outputDirectory )
{

  std::ofstream file ;
  file.open( labelsBinaryFilename ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << labelsBinaryFilename <<
                                                                     std::endl ;
    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

  int n_curves = labels.size() ;

  for ( int track = 0 ; track < n_curves ; track++ )
  {

    int nbLabels = labels[ track ].size() ;
    for ( int _labelIndex = 0 ; _labelIndex < nbLabels ; _labelIndex++ )
    {

      int label = labels[ track ][ _labelIndex ] ;

      file << track << " : " << label << std::endl ;


    }

  }

  file.close() ;


}



///////////////////////////////////////////////////////////////////////////////
/////////////// Function to read the dictionary labels file ///////////////////
///////////////////////////////////////////////////////////////////////////////
void readLabelsDict( const char* labelsDictFilename,
                     std::vector<std::string>& bundlesNames,
                     const std::string& outputDirectory )
{

  const char delim = ':' ;
  std::string line ;
  std::ifstream dictFile ;
  dictFile.open( labelsDictFilename ) ;
  if ( dictFile.fail() )
  {

    std::cout << "Problem reading file : " << labelsDictFilename
                                           << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
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

    bundlesNames.push_back( out[ 0 ] ) ;

  }

  dictFile.close() ;


}

///////////////////////////////////////////////////////////////////////////////
/////////////// Function to save the dictionary labels file ///////////////////
///////////////////////////////////////////////////////////////////////////////
void saveLabelsDict( const char* labelsDictFilename,
                     const std::vector<std::string>& bundlesNames,
                     const std::string& outputDirectory )
{

  std::ofstream file( labelsDictFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path " ;
    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 );

  }

  int nbBundles = bundlesNames.size() ;

  for ( int bundle = 0 ; bundle < nbBundles ; bundle++ )
  {

    file << bundlesNames[ bundle ] << " : " << bundle << std::endl ;

  }

  file.close() ;


}



////////////////////////////////////////////////////////////////////////////////
//////////////////// Read index in tractogram neighbors ////////////////////////
////////////////////////////////////////////////////////////////////////////////
void readIndexInTractogram( const char* predictedLabelsFilename,
                            std::vector<int64_t>& predictedLabels,
                            const std::string& outputDirectory,
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

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
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
/////////////// Function to get vector with altas bundles paths ////////////////
////////////////////////////////////////////////////////////////////////////////
std::string getAtlasBunldesPaths( const std::string& outputDirectory,
                                  const std::string& atlasDirectory,
                                  std::vector<std::string>& atlasBundlesPaths )
{

  std::ostringstream oss ;
  oss << outputDirectory << "tmpAtlasDir/" ;
  std::string tmpAtlasDir = oss.str() ;

  if ( !is_dir( tmpAtlasDir ) )
  {

    mkdir( tmpAtlasDir ) ;

  }


  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                              atlasDirectory ) )
  {

    std::string tmpBundlesDataFilename ;
    std::string tmpBundlesFilename ;

    tmpBundlesDataFilename = file.path() ;
    tmpBundlesFilename = tmpBundlesDataFilename ;
    std::string key(".bundlesdata") ;

    if ( tmpBundlesDataFilename.find( ".bundlesdata" ) != std::string::npos )
    {

      std::size_t found = tmpBundlesFilename.rfind( key ) ;
      tmpBundlesFilename.replace( found, key.length(), ".bundles" ) ;

      if ( is_file( tmpBundlesDataFilename ) && is_file( tmpBundlesFilename ) )
      {

        std::experimental::filesystem::path tmpPath( tmpBundlesFilename ) ;
        std::string bundleName = tmpPath.stem() ;

        std::ostringstream tmpBundleDirOss ;
        tmpBundleDirOss << tmpAtlasDir << bundleName << "/" ;
        std::string tmpBundleDir = tmpBundleDirOss.str( ) ;
        if ( !is_dir( tmpBundleDir ) )
        {

          mkdir( tmpBundleDir ) ;

        }

        std::ostringstream _bundlesPathOss ;
        _bundlesPathOss << tmpBundleDir << bundleName << ".bundles" ;
        std::string _bundlesPath = _bundlesPathOss.str() ;

        std::ostringstream _bundlesdataPathOss ;
        _bundlesdataPathOss << tmpBundleDir << bundleName << ".bundlesdata" ;
        std::string _bundlesdataPath = _bundlesdataPathOss.str() ;

        if ( is_file( _bundlesPath ) )
        {

          rmfile( _bundlesPath ) ;

        }
        bool isCopy1 = copy( tmpBundlesFilename, _bundlesPath ) ;
        if ( !isCopy1 )
        {

          std::cout << "ERROR : unable to copy " << tmpBundlesFilename << "to "
                    << _bundlesPath << std::endl ;

          if ( is_dir( outputDirectory ) )
          {

            rmdir( outputDirectory ) ;

          }
          exit( 1 ) ;

        }

        if ( is_file( _bundlesdataPath ) )
        {

          rmfile( _bundlesdataPath ) ;

        }
        bool isCopy2 = copy( tmpBundlesDataFilename, _bundlesdataPath ) ;
        if ( !isCopy2 )
        {

          std::cout << "ERROR : unable to copy " << tmpBundlesDataFilename
                    << " to " << _bundlesdataPath << std::endl ;

          if ( is_dir( outputDirectory ) )
          {

            rmdir( outputDirectory ) ;

          }
          exit( 1 ) ;

        }

        atlasBundlesPaths.push_back( tmpBundleDir ) ;

      }

    }

  }

  return( tmpAtlasDir ) ;


}


////////////////////////////////////////////////////////////////////////////////
/////////////// Function to get vector with altas bundles paths ////////////////
////////////////////////////////////////////////////////////////////////////////
void getAtlasNeighborhoodCentroids(
                              const std::string& outputDirectory,
                              const std::string& inputDirectory,
                              const std::vector<std::string>& atlasBundlesPaths,
                              std::vector<std::string>& bundlesPaths )
{


  int nbBundles = atlasBundlesPaths.size() ;

  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    std::string _tmpAtlasBundlePath = atlasBundlesPaths[ i ] ;

    char lastChar = _tmpAtlasBundlePath[ _tmpAtlasBundlePath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      _tmpAtlasBundlePath = _tmpAtlasBundlePath.substr( 0,
                                              _tmpAtlasBundlePath.size() - 1 ) ;

    }

    std::experimental::filesystem::path tmpPath( _tmpAtlasBundlePath ) ;
    std::string bundleName = tmpPath.stem() ;

    std::ostringstream tmpBundlesPathOss ;
    tmpBundlesPathOss << inputDirectory << bundleName << ".bundles" ;
    std::string tmpBundlesPath = tmpBundlesPathOss.str() ;

    std::ostringstream tmpBundlesDataPathOss ;
    tmpBundlesDataPathOss << inputDirectory << bundleName << ".bundlesdata" ;
    std::string tmpBundlesDataPath = tmpBundlesDataPathOss.str() ;

    if ( is_file( tmpBundlesPath ) && is_file( tmpBundlesDataPath ) )
    {

      bundlesPaths.push_back( tmpBundlesPath ) ;

    }
    else
    {

      std::cout << "ERROR : atlas neighborhood centroids file "
                << tmpBundlesPath << " or file " << tmpBundlesDataPath
                << " does not exists " << std::endl ;
      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }

  }

  if ( bundlesPaths.size() != nbBundles )
  {

    std::cout << "ERROR : could not find centroids in " << inputDirectory
              << " for all bundles in atlas " << std::endl ;
    exit( 1 ) ;

  }

}

////////////////////////////////////////////////////////////////////////////////
/////////////// Function to get vector with nieghborhood paths /////////////////
////////////////////////////////////////////////////////////////////////////////
void getNeighborhoodFilenames(
                             const std::string& tmpNeighborhoodDir,
                             const std::vector<std::string>& atlasBundlesPaths,
                             std::vector<std::string>& neighborhoodFilenames,
                             const std::string& outputDirectory )
{

  int nbBundles = atlasBundlesPaths.size() ;

  if ( !is_dir( tmpNeighborhoodDir ) )
  {

    mkdir( tmpNeighborhoodDir ) ;

  }

  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    std::string _tmpAtlasBundlePath = atlasBundlesPaths[ i ] ;

    char lastChar = _tmpAtlasBundlePath[ _tmpAtlasBundlePath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      _tmpAtlasBundlePath = _tmpAtlasBundlePath.substr( 0,
                                              _tmpAtlasBundlePath.size() - 1 ) ;

    }

    std::experimental::filesystem::path tmpPath( _tmpAtlasBundlePath ) ;
    std::string bundleName = tmpPath.stem() ;

    std::ostringstream tmpNeighborhoodPathOss ;
    tmpNeighborhoodPathOss << tmpNeighborhoodDir << bundleName << ".bundles" ;
    std::string tmpNeighborhoodPath = tmpNeighborhoodPathOss.str() ;

    if ( is_file( tmpNeighborhoodPath ) )
    {

      neighborhoodFilenames.push_back( tmpNeighborhoodPath ) ;

    }
    else
    {

      std::cout << "ERROR : neighborhood " << tmpNeighborhoodPath
                << " does not exists " << std::endl ;
      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
///////////////////// Function to convert .bundles -> trk //////////////////////
////////////////////////////////////////////////////////////////////////////////
void convertBundlesFormat( const std::string& inputBundles,
                           const std::string& outputTrk,
                           const std::string& referenceImage,
                           const std::string& outputDirectory )
{

  std::ostringstream commandOss ;
  commandOss << convertBundleFormatsFile << "  "
             << "-i " << inputBundles << " "
             << "-o " << outputTrk << " "
             << "-r " << referenceImage << " "
             << "-v " ;
  std::string command = commandOss.str() ;
  if ( is_file( outputTrk ) && !force )
  {

    if ( verbose > 1 )
    {

      std::cout << "WARNING : output file of convertBundlesFormat : "
                << outputTrk << "already exists and -force flag was not used, "
                << "trying computations with existing file" << std::endl ;
      return ;

    }

  }
  int isFail = run_sh_process( command ) ;
  if ( isFail )
  {

    std::cout << "ERROR : could not convert " << inputBundles << " to "
              << outputTrk << ", got exit code " << isFail << std::endl ;

    if ( is_dir( outputDirectory ) )
    {


      rmdir( outputDirectory ) ;
    }
    exit( 1 ) ;

  }

}


////////////////////////////////////////////////////////////////////////////////
//////// Function to save adjacencies, coverages and overlaps in a .tsv ////////
////////////////////////////////////////////////////////////////////////////////
void saveComparisonMeasuresWithAtlas(
                                    const std::vector<float>& coveragesBundles,
                                    const std::vector<float>& adjacencyBundles,
                                    const std::vector<float>& overlapBundles,
                                    const std::vector<std::string>& labelsDict,
                                    const char* fileName,
                                    const std::string& outputDirectory )
{

  std::ofstream file ;
  file.open( fileName ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << fileName << std::endl ;
    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

  int nbBundles = coveragesBundles.size() ;
  if ( adjacencyBundles.size() != nbBundles )
  {

    std::cout << "ERROR in saveComparisonMeasuresWithAtlas : the number of "
              << "bundles with coverage is different than the number of bundles"
              << " with adjacency " << std::endl ;
    exit( 1 ) ;

  }
  if ( overlapBundles.size() != nbBundles )
  {

    std::cout << "ERROR in saveComparisonMeasuresWithAtlas : the number of "
              << "bundles with coverage is different than the number of bundles"
              << " with overlap " << std::endl ;
    exit( 1 ) ;

  }
  if ( labelsDict.size() != nbBundles )
  {

    std::cout << "ERROR in saveComparisonMeasuresWithAtlas : the number of "
              << "bundles in labelsDict is different than the number of bundles"
              << " with overlap " << std::endl ;
    exit( 1 ) ;

  }

  file << "Bundle_Name\tCoverage\tAdjacency\tOverlap" << std::endl ;

  for ( int _bundle = 0 ; _bundle < nbBundles ; _bundle++ )
  {

    std::string _bundleName = labelsDict[ _bundle ] ;
    float _coverage = coveragesBundles[ _bundle ] ;
    float _adjacency = adjacencyBundles[ _bundle ] ;
    float _overlap = overlapBundles[ _bundle ] ;

    file << _bundleName << "\t" << _coverage << "\t" << _adjacency << "\t"
                                                      << _overlap << std::endl ;

  }

  file.close() ;


}

////////////////////////////////////////////////////////////////////////////////
/////////// Function to get coverage with atlas of extracted bundle ////////////
////////////////////////////////////////////////////////////////////////////////
float getCoverageWithAtlas( const std::string& bundleFilename )
{


  BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.coverageWithAtlas < 0 )
  {

    std::cout << "ERROR : got invalid coverage of " << bundle.coverageWithAtlas
                                                    << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.coverageWithAtlas ) ;

}

////////////////////////////////////////////////////////////////////////////////
/////////// Function to get adjacency with atlas of extracted bundle ///////////
////////////////////////////////////////////////////////////////////////////////
float getAdjacencyWithAtlas( const std::string& bundleFilename )
{


  BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.adjacencyWithAtlas < 0 )
  {

    std::cout << "ERROR : got invalid coverage of " << bundle.adjacencyWithAtlas
                                                    << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.adjacencyWithAtlas ) ;

}

////////////////////////////////////////////////////////////////////////////////
/////////// Function to get adjacency with atlas of extracted bundle ///////////
////////////////////////////////////////////////////////////////////////////////
float getOverlapWithAtlas( const std::string& bundleFilename )
{


  BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.overlapWithAtlas < 0 )
  {

    std::cout << "ERROR : got invalid coverage of " << bundle.overlapWithAtlas
                                                    << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.overlapWithAtlas ) ;

}



////////////////////////////////////////////////////////////////////////////////
/////////////// Function to get average radius of atlas  bundle ////////////////
////////////////////////////////////////////////////////////////////////////////
float getAverageRadiusAtlasBundle( const std::string& bundleFilename )
{

  BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.averageRadius <= 0 )
  {

    std::cout << "ERROR : got invalid radius of " << bundle.averageRadius
                                                  << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.averageRadius ) ;

}

////////////////////////////////////////////////////////////////////////////////
/// Function to get average distance between medial points of atlas  bundle ////
////////////////////////////////////////////////////////////////////////////////
float getAverageDistanceBetweenMedialPoints( const std::string& bundleFilename )
{

  BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.averageDistanceBetweenMedialPoints <= 0 )
  {

    std::cout << "ERROR : got invalid radius of "
              << bundle.averageDistanceBetweenMedialPoints << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.averageDistanceBetweenMedialPoints ) ;

}



////////////////////////////////////////////////////////////////////////////////
//////////////// Function to get number of fibers atlas bundle /////////////////
////////////////////////////////////////////////////////////////////////////////
int getNbFibers( const std::string& bundleFilename )
{

  BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.curves_count <= 0 )
  {

    std::cout << "ERROR : got invalid fiber count of " << bundle.curves_count
                                                       << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.curves_count ) ;

}



////////////////////////////////////////////////////////////////////////////////
//////////////////////// Function to apply Recobundles /////////////////////////
////////////////////////////////////////////////////////////////////////////////
void applyRecoBundles( const std::string& movedTractogramNeighborhood,
                       const std::string& atlasBundleDirectory,
                       const std::string& atlasNeighborhoodFile,
                       const std::string& outputDirectory,
                       const std::string& referenceImage,
                       int nbPointsPerFiber,
                       int verbose )
{

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
  recognizedBundleOss << outputDirectory << bundleName << ".bundles" ;
  std::string recognizedBundle = recognizedBundleOss.str() ;

  std::ostringstream recognizedBundledataOss ;
  recognizedBundledataOss << outputDirectory << bundleName << ".bundlesdata" ;
  std::string recognizedBundledata = recognizedBundledataOss.str() ;

  std::ostringstream recognizedBundleClassicOss ;
  recognizedBundleClassicOss << outputDirectory << bundleName
                                         << "_classic.bundles" ;
  std::string recognizedBundleClassic = recognizedBundleClassicOss.str() ;

  std::ostringstream recognizedBundledataClassicOss ;
  recognizedBundledataClassicOss << outputDirectory << bundleName
                                             << "_classic.bundlesdata" ;
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

  float coverage_classic = -1 ;
  float adjacency_classic = -1 ;
  // if ( false )
  if ( is_file( recognizedBundle ) )
  {

    coverage_classic = getCoverageWithAtlas( recognizedBundle ) ;
    adjacency_classic = getAdjacencyWithAtlas( recognizedBundle ) ;
    // if ( coverage_classic > 0.80 && adjacency_classic > 0.80 )
    if ( adjacency_classic > 0.90 )
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

        if ( is_dir( outputDirectory ) )
        {

          rmdir( outputDirectory ) ;

        }
        exit( 1 ) ;

      }

      isRename = rename( recognizedBundledata, recognizedBundledataClassic ) ;
      if ( !isRename )
      {

        std::cout << "ERROR : could not copy " << recognizedBundledata << " to "
                  << recognizedBundledataClassic << std::endl ;

        if ( is_dir( outputDirectory ) )
        {

          rmdir( outputDirectory ) ;

        }
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
  if ( !( endswith( atlasNeighborhoodFile, ".bundles" ) ) )
  {

    std::ostringstream _tmpOss ;
    _tmpOss << atlasNeighborhoodFile << bundleName << ".bundles" ;
    atlasNeighborhood = _tmpOss.str() ;

  }
  else
  {

    atlasNeighborhood = atlasNeighborhoodFile ;

  }

  if ( !is_file( atlasNeighborhood ) )
  {

    std::cout << "ERROR : in applyRecoBundles, atlas neighborhood file "
              << atlasNeighborhood << " does not exists " << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }


  // Computing centroids neighborhood atlas
  std::ostringstream atlasBundleFileOss ;
  atlasBundleFileOss << atlasBundleDirectory << bundleName << ".bundles" ;
  std::string atlasBundleFile = atlasBundleFileOss.str() ;
  float averageRadius = getAverageRadiusAtlasBundle( atlasBundleFile ) ;
  std::string atlasNeighborhoodCentroids ;
  if ( !isAtlasNeighborhoodCentroids )
  {

    atlasNeighborhoodCentroids = replaceExtension( atlasNeighborhood,
                                                   "_centroids.bundles" ) ;
    std::ostringstream computeCentroidsCommandOss ;
    computeCentroidsCommandOss << "python3 "
                               << computeCentroidsFilename << " "
                               << "-i " << atlasNeighborhood << " "
                               << "-o " << atlasNeighborhoodCentroids << " "
                               << "-r " << referenceImage << " "
                               << "-thr " << averageRadius << " "
                               << "-nbPoints " << nbPointsPerFiber << " "
                               << "-v 1 " ;
    std::string computeCentroidsCommand = computeCentroidsCommandOss.str() ;
    int isFailCentroids = 0 ;
    if ( countFilesDirectory( atlasNeighborhoodCentroids ) > 5 && !force )
    {

      if ( verbose > 1 )
      {

        std::cout << "WARNING : output directory of "
                  << computeCentroidsFilename << " : "
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

        isFailCentroids = run_sh_process( computeCentroidsCommand ) ;

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

      std::cout << "ERROR WITH COMMAND : \n " << computeCentroidsCommand
                                              << std::endl ;
      exit( 1 ) ;
      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }

  }
  else
  {

    atlasNeighborhoodCentroids = atlasNeighborhoodFile ;

  }

  // Compute centroids neighborhood moved tractogram
  std::string movedTractogramNeighborhoodCentroids = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "_centroids.bundles" ) ;
  std::ostringstream computeCentroidsCommand2Oss ;
  computeCentroidsCommand2Oss << "python3 "
                         << computeCentroidsFilename << " "
                         << "-i " << movedTractogramNeighborhood << " "
                         << "-o " << movedTractogramNeighborhoodCentroids << " "
                         << "-r " << referenceImage << " "
                         << "-thr " << averageRadius << " "
                         << "-nbPoints " << nbPointsPerFiber << " "
                         << "-v 1 " ;
  std::string computeCentroidsCommand2 = computeCentroidsCommand2Oss.str() ;
  bool isFailCentroids = 0 ;
  if ( countFilesDirectory( movedTractogramNeighborhoodCentroids ) > 5 &&
                                                                        !force )
  {

    if ( verbose > 1 )
    {

      std::cout << "WARNING : output directory of " << computeCentroidsFilename
                << " : " << movedTractogramNeighborhoodCentroids << " already "
                << "exists with more than 5 files and -force flag was not used,"
                << " trying computations  with existing directory"
                << std::endl ;

    }

  }
  else
  {

    int _tmpNbFibers = getNbFibers( movedTractogramNeighborhood ) ;
    if ( _tmpNbFibers > 500 )
    {

      isFailCentroids = run_sh_process( computeCentroidsCommand2 ) ;

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
    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }


  // Registering centroids
  std::string neighborhoodRegistered = replaceExtension(
                                                    movedTractogramNeighborhood,
                                                    "_moved.bundles" ) ;

  std::ostringstream registerBundlesCommadOss ;
  registerBundlesCommadOss << "python3 "
                         << registerBundlesFile << " "
                         << "-s " << atlasNeighborhoodCentroids << " "
                         << "-m " << movedTractogramNeighborhoodCentroids << " "
                         << "-ra " << referenceImage << " "
                         << "-b " << movedTractogramNeighborhood << " "
                         << "-o " << neighborhoodRegistered << " "
                         << "-n " << nbPointsPerFiber << " "
                         << "-xfm " << "rigid" << " "
                         << "-v 1" ;
  std::string registerBundlesCommad = registerBundlesCommadOss.str() ;
  // int timeout = 100 ; // In s
  int timeout = 500 ; // In s
  std::string _tmpError1 ;
  if ( is_file( neighborhoodRegistered ) && !force )
  {

    if ( verbose > 1 )
    {

      std::cout << "WARNING : output file of " << registerBundlesFile
                << " : " << neighborhoodRegistered << " already exists and "
                << "the -force flag was not used, trying computations with "
                << "existing file" << std::endl ;

    }

  }
  else
  {

    _tmpError1 = run_sh_process_timeout( registerBundlesCommad, timeout ) ;

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
                << registerBundlesCommad << "\nFail to register "
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
  atlasInfoPathOss << atlasBundleDirectory << bundleName << ".bundles" ;
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
                         << "-tolDBMP " << toleranceDistanceBetweenMedialPoints
                                                                         << " "
                         << "-thrAdj " << adjacencyForCompareBundles << " "
                         << "-cb "
                         << "-v 1" ;
  // projectAtlasCommandOss << projectAtlasFile << " "
  //                        << "-i " << neighborhoodRegistered << " "
  //                        << "-a " << atlasBundleDirectory << " "
  //                        << "-o " << outputDirectory << " "
  //                        << "-l " << "labelsSBR_" << bundleName << " "
  //                        << "-minNbFibers " << minimumNumberFibers << " "
  //                        << "-thrSim " << thrPercentageSimilarity << " "
  //                        << "-useSimple "
  //                        << "-cb "
  //                        << "-v 1" ;
  std::string projectAtlasCommand = projectAtlasCommandOss.str() ;
  // Here we do not use force because we want to always do the RecoBundles
  // projection
  int isFailProjection = run_sh_process( projectAtlasCommand ) ;
  if ( isFailProjection )
  {

    std::cout << "ERROR : could not project " << neighborhoodRegistered
              << " to " << atlasBundleDirectory << ", got exit code "
              << isFailProjection << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

  // Getting labels
  if ( is_file( recognizedBundleClassic ) &&
                                        is_file( recognizedBundledataClassic ) )
  {

    if ( coverage_classic == -1 )
    {

      std::cout << "ERROR : invalid coverage_classic : -1 " << std::endl ;
      exit( 1 ) ;
      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
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
        exit( 1 ) ;

      }

      bool _tmpIsRename2 = rename( recognizedBundledataClassic,
                                                        recognizedBundledata ) ;
      if ( !_tmpIsRename2 )
      {

        std::cout << "ERROR : could not move file "
                  << recognizedBundledataClassic << " to "
                  << recognizedBundledata << std::endl ;
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
  index_nbPoints = getFlagPosition( argc, argv, "-nbPoints" ) ;
  index_fa = getFlagPosition( argc, argv, "-fa" ) ;
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
  index_verbose = getFlagPosition( argc, argv, "-v" ) ;
  index_help = getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-i : Path to the input tractogram \n"
              << "-a : Path to the directory with the atlas \n"
              << "-ref : Path to the reference image where the tractogram is \n"
              << "-o : Path to the output directory \n"
              << "-cc : Path to the computeCentroids.py file \n"
              << "-rb : Path to the registerBundles.py file \n"
              << "-nbPoints : Number of points per fiber (same number for all "
              << "fibers) \n"
              << "[-fa] : Path to full atlas tractogram (mandatory for global "
              << "SLR or if -anc is to given) \n"
              << "[-anc] : Path to the atlas neighborhoods centroids (must be "
              << "given is -fa is not given)\n"
              << "[-thrCov] : Threshold to keep bundles where coverage is "
              << "greater than thrCov (default : 0 -> keep all bundles ) \n"
              << "[-thrAdJ] : Threshold for adjacency (Default = 10 mm) \n"
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
              << "(default = 1.2) \n"
              << "[-tolDBMP] : Tolerance for distance between medial points "
              << "(for advanced users, default = 0) \n"
              << "thrAdj : keep bundle with adjacency greater than given value"
              << " (default : 0 -> keep all bundles ) \n"
              << "[-minNbFibers] : Minimum number of fiber to consider a bundle"
              << " recognized ( default : 20 )\n"
              << "[-thrSim] : Threshold for percentage of similarity in "
              << "projection i.e nbAtlasBundleFibersSimilar / "
              << "nbFibersAtlasBundle (Range [ 0 ; 1 ], default = 0.05) \n"
              << "[-adjCB] : adjacency for computing comparison with atlas "
              << "(default = 5 mm) \n"
              << "[-pa] : Path to the ProjectAtlas file \n"
              << "[-cv] : Path to the ConvertBundleFormat file \n"
              << "[-cn] : Path to the computeNeighborhood file \n"
              << "[-slr] : Do global SLR step (default : false)\n"
              << "[-cp] : Do first a classical projection without SBR (default "
	            << ": false)\n"
              << "[-sp] : Save recognized bundles separetly (default : true)\n"
              << "[-force] : Force to overwrite files (default = false) \n"
              << "[-nbThreads] : Sets the value of omp_set_num_threads before "
              << "applyRecoBundles (default : let openMP decide ) \n"
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

  if ( !index_cc )
  {

    std::cout << "-cc argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_rb )
  {

    std::cout << "-rb argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_nbPoints )
  {

    std::cout << "-nb argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Checking arguments /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  ////////////////////////////// Input tractogram //////////////////////////////
  std::string inputFilename( argv[ index_input + 1 ] ) ;
  std::string inputBundlesFilename ;
  std::string inputBundlesDataFilename ;
  std::string inputTRKFilename ;
  char lastChar = inputFilename[ inputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    inputFilename = inputFilename.substr( 0, inputFilename.size() - 1 ) ;

  }
  if ( inputFilename.find( ".bundlesdata" ) != std::string::npos )
  {

    inputBundlesDataFilename = inputFilename ;

    inputBundlesFilename = inputFilename ;
    size_t index = inputFilename.find( ".bundlesdata" ) ;
    inputBundlesFilename.replace( index, 12, ".bundles" ) ;

    inputTRKFilename = inputFilename ;
    index = inputFilename.find( ".bundlesdata" ) ;
    inputTRKFilename.replace( index, 12, ".trk" ) ;

  }
  else if ( inputFilename.find( ".bundles" ) != std::string::npos )
  {

    inputBundlesFilename = inputFilename ;

    inputBundlesDataFilename = inputFilename ;
    size_t index = inputFilename.find( ".bundles" ) ;
    inputBundlesDataFilename.replace( index, 8, ".bundlesdata" ) ;

    inputTRKFilename = inputFilename ;
    inputTRKFilename.replace( index, 8, ".trk" ) ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }
  if ( !is_file( inputFilename ) )
  {

    std::cout << "ERROR : Input tratogram file " << inputFilename
                                            << " does not exists" << std::endl ;

  }
  else
  {

    std::cout << "Input tractogram : OK " << std::endl ;

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
  /////////////////////////////// Reference image //////////////////////////////
  std::string referenceFilename( argv[ index_reference + 1 ] ) ;
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
  computeCentroidsFilename = argv[ index_cc + 1 ] ;
  lastChar = computeCentroidsFilename[ computeCentroidsFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    computeCentroidsFilename = computeCentroidsFilename.substr( 0,
                                         computeCentroidsFilename.size() - 1 ) ;

  }
  if ( !is_file( computeCentroidsFilename ) )
  {

    std::cout << "ERROR : computeCentroids.py file " << computeCentroidsFilename
              << " does not exists " << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

  ////////////////////////// Register bundles command //////////////////////////
  registerBundlesFile = argv[ index_rb + 1 ] ;
  lastChar = registerBundlesFile[ registerBundlesFile.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    registerBundlesFile = registerBundlesFile.substr( 0,
                                              registerBundlesFile.size() - 1 ) ;

  }
  if ( !is_file( registerBundlesFile ) )
  {

    std::cout << "ERROR : computeCentroids.py file " << registerBundlesFile
              << " does not exists " << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

  ///////////////////////// Number of points per fiber /////////////////////////
  nbPointsPerFiber = std::stoi( argv[ index_nbPoints + 1 ] ) ;

  ///////////////////////////////// Full atlas /////////////////////////////////
  std::string fullAtlasFilename ;
  std::string fullAtlasTrkFilename ;
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
    if ( tmpFullAtlasPath.find( ".bundlesdata" ) != std::string::npos )
    {

      fullAtlasFilename = tmpFullAtlasPath ;
      size_t index = tmpFullAtlasPath.find( ".bundlesdata" ) ;
      fullAtlasFilename.replace( index, 12, ".bundles") ;
      fullAtlasTrkFilename = tmpFullAtlasPath ;
      fullAtlasTrkFilename.replace( index, 12, ".trk") ;

    }
    else if ( tmpFullAtlasPath.find( ".bundles" ) != std::string::npos )
    {

      fullAtlasFilename = tmpFullAtlasPath ;

      fullAtlasTrkFilename = tmpFullAtlasPath ;
      size_t index = tmpFullAtlasPath.find( ".bundles" ) ;
      fullAtlasTrkFilename.replace( index, 8, ".trk") ;

    }
    else
    {

      std::cout << "The only full atlas format supported is .bundles"
                << std::endl ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }

      exit( 1 ) ;

    }
    if ( !is_file( fullAtlasFilename ) )
    {

      std::cout << "ERROR : Full atlas file " << fullAtlasFilename
                                            << " does not exists" << std::endl ;

    }
    else
    {

      std::cout << "Full atlas tractogram : OK " << std::endl ;

    }

  }

  //////////////////////// Atlas neighborhood centroids ////////////////////////
  std::string atlasNeighborhoodCentroidsDirectory ;
  if ( index_anc )
  {

    isAtlasNeighborhoodCentroids = true ;

    atlasNeighborhoodCentroidsDirectory = argv[ index_anc + 1 ] ;
    lastChar = atlasNeighborhoodCentroidsDirectory[ atlasDirectory.size() - 1 ] ;
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

  ////////////////// Threshold distance between medial points //////////////////
  if ( index_thrDBMP )
  {

    thrDistanceBetweenMedialPoints = std::stof( argv[ index_thrDBMP + 1 ] ) ;

    if ( thrDistanceBetweenMedialPoints < 0 )
    {

      std::cout << "Error argument : thrDBMP must be positive " << std::endl ;
      exit( 1 ) ;

    }


  }

  /////////////////////////////// Tolerances ///////////////////////////////////

  if ( index_tolP )
  {

    toleranceP = std::stof( argv[ index_tolP + 1 ] ) ;

    if ( toleranceP < -1 || toleranceP > 1 )
    {

      std::cout << "Error argument : -tolP must be in [-1;1]" << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_tolThr )
  {

    toleranceThr = std::stof( argv[ index_tolThr + 1 ] ) ;

    if ( toleranceThr < -1 || toleranceThr > 1 )
    {

      std::cout << "Error argument : -tolThr must be in [-1;1]" << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_tolMaxAngle )
  {

    toleranceMaxAngle = std::stof( argv[ index_tolMaxAngle + 1 ] ) ;

    if ( toleranceMaxAngle < -1 || toleranceMaxAngle > 1 )
    {

      std::cout << "Error argument : -tolMaxAng must be in [-1;1]" << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_tolMaxDirectionAngle )
  {

    toleranceMaxDirectionAngle = std::stof( argv[
                                            index_tolMaxDirectionAngle + 1 ] ) ;

    if ( toleranceMaxDirectionAngle < -1 || toleranceMaxDirectionAngle > 1 )
    {

      std::cout << "Error argument : -tolMaxDirAng must be in [-1;1]"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_tolMinShapeAngle )
  {

    toleranceMinShapeAngle = std::stof( argv[ index_tolMinShapeAngle + 1 ] ) ;

    if ( toleranceMinShapeAngle < -1 || toleranceMinShapeAngle > 1 )
    {

      std::cout << "Error argument : -tolMinShapeAng must be in [-1;1]"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_tolMaxShapeAngle )
  {

    toleranceMaxShapeAngle = std::stof( argv[ index_tolMaxShapeAngle + 1 ] ) ;

    if ( toleranceMaxShapeAngle < -1 || toleranceMaxShapeAngle > 1 )
    {

      std::cout << "Error argument : -tolMaxShapeAng must be in [-1;1]"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_tolLenght )
  {

    toleranceLenght = std::stof( argv[ index_tolLenght + 1 ] ) ;

    if ( toleranceLenght < -1 || toleranceLenght > 1 )
    {

      std::cout << "Error argument : -tolLenght must be in [-1;1]"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

  }


  if ( index_tolThrCN )
  {

    toleranceThrComputeNeighborhood = std::stof( argv[ index_tolThrCN + 1 ] ) ;
    if ( toleranceThrComputeNeighborhood <= 0 )
    {

      std::cout << "ERROR : argument -tolThrCN must be greater than 0"
                                                                  << std::endl ;
      exit( 1 ) ;
    }

  }

  if ( index_tolDistBetMedPts )
  {

    toleranceDistanceBetweenMedialPoints = std::stof(
                                          argv[ index_tolDistBetMedPts + 1 ] ) ;

    if ( toleranceDistanceBetweenMedialPoints < -1 )
    {

      std::cout << "Error argument : -tolDBMP must be greater than -1"
                                                                  << std::endl ;
      exit( 1 ) ;

    }

  }


  /////////////////////////// Minimum number of fibers /////////////////////////
  if ( index_minNbFibers )
  {

    minimumNumberFibers = std::stoi( argv[ index_minNbFibers + 1 ] ) ;
    if ( minimumNumberFibers < 1 )
    {

      std::cout << "ERROR : argument -minNbFibers must be > 0 " << std::endl ;
      exit( 1 ) ;

    }

  }

  /////////////////////////// Minimum number of fibers /////////////////////////
  if ( index_thrSim )
  {

    thrPercentageSimilarity = std::stof( argv[ index_thrSim + 1 ] ) ;
    if ( thrPercentageSimilarity <= 0 || thrPercentageSimilarity > 1 )
    {

      std::cout << "ERROR : argument -thrSim must be greater than 0 and "
                << "lower or equal to 1 " << std::endl ;
      exit( 1 ) ;
    }

  }

  /////////////////////////// Minimum number of fibers /////////////////////////
  if ( index_adjCB )
  {

    adjacencyForCompareBundles = std::stof( argv[ index_adjCB + 1 ] ) ;
    if ( adjacencyForCompareBundles <= 0 )
    {

      std::cout << "ERROR : argument -adjCB must be greater than 0 "
                                                                  << std::endl ;
      exit( 1 ) ;
    }

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

  ////////////////////////// computeNeighborhood file //////////////////////////
  if ( index_rb )
  {

    std::string tmpFile( argv[ index_rb + 1 ] ) ;
    if ( !is_file( tmpFile ) )
    {

      std::cout << "ERROR : registerBunldes file " << tmpFile
                << " does not exists " << std::endl ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }

    registerBundlesFile = tmpFile ;

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
  //////////////// Preparing atlas bundles paths (tmpAtlasDir) /////////////////
  if ( verbose )
  {

    std::cout << "#########################################################\n" ;
    std::cout << "############# Preparing atlas bundles paths #############"
                                                                  << std::endl ;
    std::cout << "#########################################################\n" ;

  }
  checkAtlasDirectory( atlasDirectory, outputDirectory ) ;
  std::vector<std::string> atlasBundleDirectories ;
  std::string tmpAtlasDir = getAtlasBunldesPaths( outputDirectory,
                                                  atlasDirectory,
                                                  atlasBundleDirectories ) ;

  if ( verbose )
  {

    std::cout << "Done " << std::endl ;

  }


  ///////////////// Getting atlas centroids if input is given //////////////////
  std::vector<std::string> atlasNeighborhoodCentroidsPaths ;
  if ( isAtlasNeighborhoodCentroids )
  {

    getAtlasNeighborhoodCentroids( outputDirectory,
                                   atlasNeighborhoodCentroidsDirectory,
                                   atlasBundleDirectories,
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

    std::ostringstream tmpSLRdirOss ;
    tmpSLRdirOss << outputDirectory << "tmpSLR/" ;
    std::string tmpSLRdir = tmpSLRdirOss.str() ;
    if ( !is_dir( tmpSLRdir ) )
    {

      mkdir( tmpSLRdir ) ;

    }

    convertBundlesFormat( fullAtlasFilename,
                          fullAtlasTrkFilename,
                          referenceFilename,
                          outputDirectory ) ;


    convertBundlesFormat( inputBundlesFilename,
                          inputTRKFilename,
                          referenceFilename,
                          outputDirectory ) ;

    std::ostringstream movedTractogramTrkOss ;
    movedTractogramTrkOss << tmpSLRdir << "moved.trk" ;
    std::string movedTractogramTrk = movedTractogramTrkOss.str() ;

    std::ostringstream globalSLROss ;
    // globalSLROss << "dipy_slr "
    //              << fullAtlasTrkFilename << " "
    //              << inputTRKFilename << " "
    //              << "--greater_than 10 "
    //              << "--less_than 200 "
    //              << "--out_dir " << tmpSLRdir << " "
    //              << "--out_moved " << movedTractogramTrk << " "
    //              << "--force " ;
    globalSLROss << "dipy_slr "
                 << fullAtlasTrkFilename << " "
                 << inputTRKFilename << " "
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
      movedTractogramOss << tmpSLRdir << "moved.bundles" ;
      std::string movedTractogram = movedTractogramOss.str() ;

      convertBundlesFormat( movedTractogramTrk,
                            movedTractogram,
                            referenceFilename,
                            outputDirectory ) ;

      // Cleaning
      rmfile( fullAtlasTrkFilename ) ;
      rmfile( inputTRKFilename ) ;
      rmfile( movedTractogramTrk ) ;

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
				                        << "-maxLen " << 100 << " "
                                << "-tolThr " << toleranceThrComputeNeighborhood
                                                                          << " "
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

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }


  std::vector<std::string> neighborhoodFilenames ;
  getNeighborhoodFilenames( tmpNeighborhoodDir,
                            atlasBundleDirectories,
                            neighborhoodFilenames,
                            outputDirectory ) ;


  if ( verbose )
  {

    std::cout << "Done" << std::endl ;

  }

  ///////////////////////// Computing neighborhood atlas ///////////////////////
  std::vector<std::string> neighborhoodAtlasFilenames ;
  std::string tmpNeighborhoodAtlasDir ;
  if ( !isAtlasNeighborhoodCentroids && isFullAtlas )
  {

    if ( verbose )
    {

      std::cout << "#########################################################"
                << std::endl ;
      std::cout << "Computing atlas neighborhood... " << std::endl ;
      std::cout << "#########################################################"
                << std::endl ;

    }

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

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }

    getNeighborhoodFilenames( tmpNeighborhoodAtlasDir,
                              atlasBundleDirectories,
                              neighborhoodAtlasFilenames,
                              outputDirectory ) ;


    if ( verbose )
    {

      std::cout << "Done" << std::endl ;

    }

  }

  ////////////////////////// Projecting atlas without SBR ////////////////////////
  /*
  // Not necessary becaus done later on neighborhoods (maybe remove this part)
  if ( verbose && doClassical )
  {

    std::cout << "#########################################################\n" ;
    std::cout << "############## Projecting atlas without SBR #############"
                                                                  << std::endl ;
    std::cout << "#########################################################\n" ;

  }

  std::ostringstream projectCommandOss ;
  projectCommandOss << projectAtlasFile << " "
                    << "-i " << movedTractogram << " "
                    << "-a " << atlasDirectory << " "
                    << "-o " << outputDirectory << " "
                    << "-l " << "labels" << " "
                    << "-tolP " << toleranceP << " "
                    << "-tolThr " << toleranceThr << " "
                    << "-tolMaxAng " << toleranceMaxAngle << " "
                    << "-tolMaxDirAng " << toleranceMaxDirectionAngle
                                                                    << " "
                    << "-tolMinShapeAng " << toleranceMinShapeAngle << " "
                    << "-tolMaxShapeAng " << toleranceMaxShapeAngle << " "
                    << "-tolLenght " << toleranceLenght << " "
                    << "-cb "
                    << "-v 1" ;
  std::string projectCommand = projectCommandOss.str() ;
  int isProjectAtlasFail = 0 ;
  if ( countFilesDirectory( outputDirectory ) > 5 && !force )
  {
    if ( verbose > 1 )
    {
      std::cout << "WARNING : output directory of " << projectAtlasFile
                << " : " <<  outputDirectory << " already exists with more than"
                << " 5 files and -force flag was not used, trying computations "
                << "with existing directory" << std::endl ;

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

    std::cout << "ERROR : could not project atlas without SBR, got exit code "
              << isProjectAtlasFail << std::endl ;

    if ( is_dir( outputDirectory ) )
    {

      rmdir( outputDirectory ) ;

    }
    exit( 1 ) ;

  }

  if ( verbose && doClassical )
  {

    std::cout << "Done" << std::endl ;

  }

  std::ostringstream labelsDictPathOss ;
  labelsDictPathOss << outputDirectory << "labels.dict" ;
  std::string labelsDictPath = labelsDictPathOss.str() ;

  std::ostringstream labelsTxtPathOss ;
  labelsTxtPathOss << outputDirectory << "labels.txt" ;
  std::string labelsTxtPath = labelsTxtPathOss.str() ;

  std::ostringstream labelsDictClassicPathOss ;
  labelsDictClassicPathOss << outputDirectory << "labelsClassic.dict" ;
  std::string labelsDictClassicPath = labelsDictClassicPathOss.str() ;

  std::ostringstream labelsTxtClassicPathOss ;
  labelsTxtClassicPathOss << outputDirectory << "labelsClassic.txt" ;
  std::string labelsTxtClassicPath = labelsTxtClassicPathOss.str() ;

  if ( doClassical )
  {

    bool isRename = rename( labelsDictPath, labelsDictClassicPath ) ;
    if ( !isRename )
    {

      std::cout << "ERROR : could not copy " << labelsDictPath << " to "
                << labelsDictClassicPath << std::endl ;

      exit( 1 ) ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }

    isRename = rename( labelsTxtPath, labelsTxtClassicPath ) ;
    if ( !isRename )
    {

      std::cout << "ERROR : could not copy " << labelsTxtPath << " to "
                << labelsTxtClassicPath << std::endl ;

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
      exit( 1 ) ;

    }
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
  */


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

  if ( nbThreads )
  {

    omp_set_num_threads( nbThreads ) ;

  }
  #pragma omp parallel for
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
    atlasInfoPathOss << atlasBundleDirectory << bundleName << ".bundles" ;
    std::string atlasInfoPath = atlasInfoPathOss.str() ;
    float thrDistanceBetweenMedialPointsBundle =
                                                thrDistanceBetweenMedialPoints ;
    if ( !index_thrDBMP )
    {

      thrDistanceBetweenMedialPointsBundle =
                        getAverageDistanceBetweenMedialPoints( atlasInfoPath ) ;
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
                      << "-tolDBMP " << toleranceDistanceBetweenMedialPoints
                                                                          << " "
                      << "-thrAdj " << adjacencyForCompareBundles << " "
                      << "-cb "
                      << "-v 1" ;
    std::string projectCommand = projectCommandOss.str() ;

    std::ostringstream extractedBundlePathOss ;
    extractedBundlePathOss << outputDirectory << bundleName << ".bundles" ;
    std::string extractedBundlePath = extractedBundlePathOss.str() ;

    std::ostringstream extractedBundleDataPathOss ;
    extractedBundleDataPathOss << outputDirectory << bundleName
                                                             << ".bundlesdata" ;
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

      if ( is_dir( outputDirectory ) )
      {

        rmdir( outputDirectory ) ;

      }
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
                              outputDirectory,
                              nbFibersMovedTracotgramNeighborhood,
                              1 ) ;

      std::vector<int64_t> indexInTractogram ;
      readIndexInTractogram( neighborhoodTractogramBinPath.c_str(),
                             indexInTractogram,
                             outputDirectory,
                             nbFibersMovedTracotgramNeighborhood,
                             1 ) ;

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


  //////////////////
  // std::ostringstream labelsDictClassicPathOss ;
  // labelsDictClassicPathOss << outputDirectory << "labelsClassic.dict" ;
  // std::string labelsDictClassicPath = labelsDictClassicPathOss.str() ;
  // saveLabelsDict( labelsDictClassicPath.c_str(),
  //                 labelsDictClassic,
  //                 outputDirectory ) ;
  //
  // std::ostringstream labelsTxtClassicPathOss ;
  // labelsTxtClassicPathOss << outputDirectory << "labelsClassic.txt" ;
  // std::string labelsTxtClassicPath = labelsTxtClassicPathOss.str() ;
  // saveLabels( labelsTxtClassicPath.c_str(),
  //             labelsClassic,
  //                  const std::string& outputDirectory ) ;
  //////////////////

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

  std::vector<std::vector<int16_t>> labels ;
  labels.resize( nbFibersTractogram ) ;
  int countProcessSubjects = 1 ;
  if ( nbThreads )
  {

    omp_set_num_threads( nbThreads ) ;

  }
  #pragma omp parallel for
  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    std::string atlasBundleDirectory = atlasBundleDirectories[ i ] ;
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

    applyRecoBundles( movedTractogramNeighborhood,
                      atlasBundleDirectory,
                      atlasNeighborhoodFile,
                      outputDirectory,
                      referenceFilename,
                      nbPointsPerFiber,
                      verbose ) ;

    // Ensure that the section {} will be only access by a thread at the time
    // #pragma omp critical
    // {
    //
    //   if ( verbose > 1 )
    //   {
    //
    //     std::cout << "-----------------------------------------------------\n" ;
    //
    //   }
    //
    //   if ( verbose == 1 )
    //   {
    //
    //     printf( "\rProjecting atlas bundles with SBR : [ %d  /  %d ]",
    //                                          countProcessSubjects, nbBundles ) ;
    //     std::cout << "" << std::flush ;
    //
    //   }
    //   else if ( verbose > 1 )
    //   {
    //
    //     printf( "Projecting atlas bundles with SBR : [ %d / %d ]\n",
    //                                          countProcessSubjects, nbBundles ) ;
    //
    //   }
    //
    //   if ( verbose > 1 )
    //   {
    //
    //     std::cout << "-----------------------------------------------------\n" ;
    //
    //   }
    //
    //   countProcessSubjects++ ;
    //
    // }



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
                              outputDirectory,
                              nbFibersMovedTracotgramNeighborhood,
                              1 ) ;

      std::vector<int64_t> indexInTractogram ;
      readIndexInTractogram( neighborhoodTractogramBinPath.c_str(),
                             indexInTractogram,
                             outputDirectory,
                             nbFibersMovedTracotgramNeighborhood,
                             1 ) ;

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
                  labelsDictClassic,
                  outputDirectory ) ;


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
  saveLabels( labelsTxtPath.c_str(), labels, outputDirectory ) ;



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
    recognizedBundleOss << outputDirectory << bundleName << ".bundles" ;
    std::string recognizedBundlePath = recognizedBundleOss.str() ;

    std::ostringstream recognizedBundledataOss ;
    recognizedBundledataOss << outputDirectory << bundleName << ".bundlesdata" ;
    std::string recognizedBundledataPath = recognizedBundledataOss.str() ;

    if ( !is_file( recognizedBundlePath )  ||
                                          !is_file( recognizedBundledataPath ) )
    {

      continue ;

    }

    BundlesDataFormat recognizedBundlesData( recognizedBundledataPath.c_str(),
                                             recognizedBundlePath.c_str(),
                                             0  ) ;
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
  saveComparisonMeasuresWithAtlas( coveragesBundles,
                                   adjacencyBundles,
                                   overlapBundles,
                                   labelsDictSBR,
                                   comparisonWithAtlasRecognizedSBRPath.c_str(),
                                   outputDirectory ) ;

  saveLabelsDict( labelsDictRecognizedSBRPath.c_str(),
                  labelsDictSBR,
                  outputDirectory ) ;

  saveLabels( labelsRecognizedSBRPath.c_str(),
              labelsSBR,
              outputDirectory ) ;


  BundlesDataFormat regroupedRecognizedBundlesData( matrixTracksSBR,
                                                    pointsPerTrackSBR,
                                                    curves_countSBR ) ;
  std::ostringstream regroupedRecognizedBundledataOss ;
  regroupedRecognizedBundledataOss << outputDirectory
                                   << "regroupedRecognized.bundlesdata" ;
  std::string regroupedRecognizedBundledataPath =
                                      regroupedRecognizedBundledataOss.str() ;
  regroupedRecognizedBundlesData.bundlesdataWriting(
                                    regroupedRecognizedBundledataPath.c_str(),
                                    0 ) ;


  BundlesFormat regroupedRecognizedBundles( movedTractogram.c_str(), 0 ) ;
  std::ostringstream regroupedRecognizedBundleOss ;
  regroupedRecognizedBundleOss << outputDirectory
                                   << "regroupedRecognized.bundles" ;
  std::string regroupedRecognizedBundlePath =
                                      regroupedRecognizedBundleOss.str() ;
  regroupedRecognizedBundles.curves_count = curves_countSBR ;
  regroupedRecognizedBundles.bundlesWriting(
                                        regroupedRecognizedBundlePath.c_str(),
                                        0 ) ;


  std::cout << "Done" << std::endl ;


  ////////////////////////////////// Cleaning //////////////////////////////////
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
