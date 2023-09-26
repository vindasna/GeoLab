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

#include "ioWrapper.h"
#include "bundlesMinf.h"
// #include "BundlesFormat.h"


////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to check if file exists //////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool is_file( const std::string& path )
{

  struct stat buffer ;
  return( stat ( path.c_str(), &buffer ) == 0 ) ;

}

////////////////////////////////////////////////////////////////////////////////
////////////////// Function to check if directory exists ///////////////////////
////////////////////////////////////////////////////////////////////////////////

bool is_dir( const std::string& path )
{

  return( std::experimental::filesystem::is_directory( path ) ) ;

}


////////////////////////////////////////////////////////////////////////////////
/////////////////////// Function to create directory ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool mkdir( const std::string& path )
{

  return( std::experimental::filesystem::create_directory( path ) ) ;

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////// Function to delete file /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool rmfile( const std::string& path )
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
bool rmdir( const std::string& path )
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
bool copy( const std::string& source,
                  const std::string& destination )
{

  return( std::experimental::filesystem::copy_file( source, destination ) ) ;

}

////////////////////////////////////////////////////////////////////////////////
/////////////////////// Function to copy directory /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool copytree( const std::string& source,
               const std::string& destination,
               bool forceOverride )
{

  if ( !is_dir( source ) )
  {

    std::cout << "ERROR : in copytree(), source " << source
              << " is not a directoy " << std::endl ;
    exit( 1 ) ;

  }

  if ( is_dir( destination ) && !forceOverride )
  {

    std::cout << "ERROR : in copytree(), destination " << destination
              << " already exists and forceOverride is set to false"
              << std::endl ;
    exit( 1 ) ;

  }


  std::experimental::filesystem::copy( source, destination ) ;

  if ( !is_dir( destination ) )
  {

    return( false ) ;

  }
  else
  {

    if ( countFilesDirectory( source ) != countFilesDirectory( destination ) )
    {

      return( false ) ;

    }
    else
    {

      return( true ) ;

    }

  }

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////// Function to rename file //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool rename( const std::string& source, const std::string& destination )
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
//////////////////// Function to get parent directory //////////////////////////
////////////////////////////////////////////////////////////////////////////////
std::string dirname( const std::string& path )
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
/////////////////////// Function to get basename of path ///////////////////////
////////////////////////////////////////////////////////////////////////////////
std::string basename( const std::string& path )
{

  std::string tmpString = path ;
  char lastChar = tmpString[ tmpString.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    tmpString = tmpString.substr( 0, tmpString.size() - 1 ) ;

  }

  std::experimental::filesystem::path p( tmpString ) ;

  std::ostringstream _tmpOss ;
  _tmpOss << p.filename() ;
  std::string _tmp = _tmpOss.str() ;
  _tmp = _tmp.substr( 1, _tmp.size() ) ;
  _tmp = _tmp.substr( 0, _tmp.size() - 1 ) ;

  return( _tmp ) ;

}


////////////////////////////////////////////////////////////////////////////////
////////////// Function to get basename of path without extension //////////////
////////////////////////////////////////////////////////////////////////////////
std::string basenameNoExtension( const std::string& path )
{

  std::string tmpString = path ;
  char lastChar = tmpString[ tmpString.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    tmpString = tmpString.substr( 0, tmpString.size() - 1 ) ;

  }

  std::experimental::filesystem::path p( tmpString ) ;

  std::ostringstream _tmpOss ;
  _tmpOss << p.stem() ;
  std::string _tmp = _tmpOss.str() ;
  _tmp = _tmp.substr( 1, _tmp.size() ) ;
  _tmp = _tmp.substr( 0, _tmp.size() - 1 ) ;

  return( _tmp ) ;

}


////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to replace extension file ////////////////////////
////////////////////////////////////////////////////////////////////////////////
std::string replaceExtension( const std::string& path,
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
//////////////////// Function to get path without extension ////////////////////
////////////////////////////////////////////////////////////////////////////////
std::string getFilenameNoExtension( const std::string& path )
{

  std::string tmpString = path ;
  char lastChar = tmpString[ tmpString.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    tmpString = tmpString.substr( 0, tmpString.size() - 1 ) ;

  }

  std::experimental::filesystem::path p( tmpString ) ;
  std::string filenameNoExtension = p.stem() ;

  return( filenameNoExtension ) ;

}


////////////////////////////////////////////////////////////////////////////////
/////////////// Function to check if string ends with sub-string ///////////////
////////////////////////////////////////////////////////////////////////////////
bool endswith( const std::string& input,
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
//////// Function to get files in directory with specific extension ////////////
////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> getFilesInDirectoryWithExtension(
                                                  const std::string& path,
                                                  const std::string& extension )
{

  if ( !is_dir( path ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ioWrapper->getFilesInDirectoryWithExtension : directory "
                  << path << " does not exists" ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  std::vector<std::string> outFiles ;


  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                                        path ) )
  {

    std::string tmpFilename = file.path() ;

    if ( endswith( tmpFilename, extension ) )
    {

      outFiles.push_back( tmpFilename ) ;

    }

  }

  if ( outFiles.size() == 0 )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ioWrapper->getFilesInDirectoryWithExtension : the "
                  << "directory " << path << " does not contain any files in "
                  << extension << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    // throw( std::invalid_argument( outMessage ) ) ;
    std::cout << "WARNING : " << outMessage ;

  }

  return( outFiles ) ;

}


////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to get files in directoy /////////////////////////
////////////////////////////////////////////////////////////////////////////////
void listDir( const std::string& path,
              std::vector<std::string>& dirList )
{

  if ( !is_dir( path ) )
  {

    std::cout << "ERROR : in listDir, directoy " << path << "does not exists"
                                                         << std::endl ;

    exit( 1 ) ;

  }

  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                                        path ) )
  {

    std::ostringstream fileNameOss ;
    fileNameOss << file ;
    std::string fileName = fileNameOss.str() ;
    fileName = fileName.substr( 1, fileName.size() ) ;
    fileName = fileName.substr( 0, fileName.size() - 1 ) ;

    dirList.push_back( fileName ) ;

  }


}

////////////////////////////////////////////////////////////////////////////////
//////////////////// Function to count files directory /////////////////////////
////////////////////////////////////////////////////////////////////////////////
void checkAtlasDirectory( const std::string& path, const std::string& format )
{

  int nbCorrectFiles = 0 ;

  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                                        path ) )
  {

    std::string tmpBundlesDataFilename ;
    std::string tmpBundlesMinfFilename ;

    tmpBundlesDataFilename = file.path() ;
    tmpBundlesMinfFilename = tmpBundlesDataFilename ;

    bool isFileWithExtension = false ;

    if ( format == ".bundles" &&
                                endswith( tmpBundlesDataFilename, ".bundles" ) )
    {

      tmpBundlesDataFilename = replaceExtension( tmpBundlesDataFilename,
                                                            ".bundlesdata" ) ;

      if ( !is_file( tmpBundlesDataFilename ) ||
                                            !is_file( tmpBundlesMinfFilename ) )
      {

        std::stringstream outMessageOss ;
        outMessageOss << "ioWrapper->checkAtlasDirectory : for .bundles/ "
                      << ".bundlesdata format each .bundles file must have "
                      << "its .bundlesdata file which is not the case for "
                      << tmpBundlesMinfFilename << " and "
                      << tmpBundlesDataFilename << std::endl ;
        std::string outMessage = outMessageOss.str() ;

        throw( std::invalid_argument( outMessage ) ) ;

      }

      nbCorrectFiles += 1 ;

    }
    else if ( format == ".bundlesdata" &&
                            endswith( tmpBundlesDataFilename, ".bundlesdata" ) )
    {

      tmpBundlesMinfFilename = replaceExtension( tmpBundlesDataFilename,
                                                                  ".bundles" ) ;

      if ( !is_file( tmpBundlesDataFilename ) ||
                                            !is_file( tmpBundlesMinfFilename ) )
      {

        std::stringstream outMessageOss ;
        outMessageOss << "ioWrapper->checkAtlasDirectory : for .bundles/ "
                      << ".bundlesdata format each .bundles file must have "
                      << "its .bundlesdata file which is not the case for "
                      << tmpBundlesMinfFilename << " and "
                      << tmpBundlesDataFilename << std::endl ;
        std::string outMessage = outMessageOss.str() ;

        throw( std::invalid_argument( outMessage ) ) ;

      }

      nbCorrectFiles += 1 ;


    }
    else if ( format == ".trk" && endswith( tmpBundlesDataFilename, ".trk" ) )
    {

      tmpBundlesMinfFilename = replaceExtension( tmpBundlesDataFilename,
                                                                     ".minf" ) ;

      if ( !is_file( tmpBundlesDataFilename ) ||
                                            !is_file( tmpBundlesMinfFilename ) )
      {

        std::stringstream outMessageOss ;
        outMessageOss << "ioWrapper->checkAtlasDirectory : for .trk format each"
                      << " .trk file must have its .minf file which is not the "
                      << "case for " << tmpBundlesMinfFilename << " and "
                      << tmpBundlesDataFilename << std::endl ;
        std::string outMessage = outMessageOss.str() ;

        throw( std::invalid_argument( outMessage ) ) ;

      }

      nbCorrectFiles += 1 ;


    }
    else if ( format == ".tck" && endswith( tmpBundlesDataFilename, ".tck" ) )
    {

      tmpBundlesMinfFilename = replaceExtension( tmpBundlesDataFilename,
                                                                     ".minf" ) ;

      if ( !is_file( tmpBundlesDataFilename ) ||
                                            !is_file( tmpBundlesMinfFilename ) )
      {

        std::stringstream outMessageOss ;
        outMessageOss << "ioWrapper->checkAtlasDirectory : for .tck format each"
                      << " .trk file must have its .minf file which is not the "
                      << "case for " << tmpBundlesMinfFilename << " and "
                      << tmpBundlesDataFilename << std::endl ;
        std::string outMessage = outMessageOss.str() ;

        throw( std::invalid_argument( outMessage ) ) ;

      }

      nbCorrectFiles += 1 ;


    }
    else if ( format != ".bundles" && format != ".bundlesdata" &&
              format != ".trk" && format != ".tck" )
    {

      std::string outMessage = "The only atlas supported format ara .bundles/" \
                               ".bundlesdata, .trk and .tck \n" ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

  }

  if ( nbCorrectFiles > 0 )
  {

    std::cout << "Atlas directory : " << path << std::endl
              << "Number of bundles with " << format << " format : "
              << nbCorrectFiles << std::endl ;

  }
  else
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ioWrapper->checkAtlasDirectory : the atlas directory "
                  << path << " does not contain any bundles in " << format
                  << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

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
int run_sh_process( const std::string& command, int verbose )
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
    while ( std::getline( out, line ) ){}

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
                             int nbFibers )
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
                 const std::vector<std::vector<int16_t>>& labels )
{

  std::ofstream file ;
  file.open( labelsBinaryFilename ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << labelsBinaryFilename <<
                                                                     std::endl ;
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
                     std::vector<std::string>& bundlesNames )
{

  const char delim = ':' ;
  std::string line ;
  std::ifstream dictFile ;
  dictFile.open( labelsDictFilename ) ;
  if ( dictFile.fail() )
  {

    std::cout << "Problem reading file : " << labelsDictFilename
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

    bundlesNames.push_back( out[ 0 ] ) ;

  }

  dictFile.close() ;


}

///////////////////////////////////////////////////////////////////////////////
/////////////// Function to save the dictionary labels file ///////////////////
///////////////////////////////////////////////////////////////////////////////
void saveLabelsDict( const char* labelsDictFilename,
                     const std::vector<std::string>& bundlesNames )
{

  std::ofstream file( labelsDictFilename ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path " ;
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
////// Function to read labels into a list of vecotr of bundles names //////////
////////////////////////////////////////////////////////////////////////////////
void readLabelsWithDict( const char* labelsDictFilename,
                         const char* labelsBinaryFilename,
                         std::vector<std::vector<std::string>>& labelsByName,
                         int nbFibers )
{

  // Read dictionary
  std::vector<std::string> bundlesNames ;

  const char delim = ':' ;
  std::string line ;
  std::ifstream dictFile ;
  dictFile.open( labelsDictFilename ) ;
  if ( dictFile.fail() )
  {

    std::cout << "Problem reading file : " << labelsDictFilename
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

    bundlesNames.push_back( out[ 0 ] ) ;

  }

  dictFile.close() ;

  int nbLabels = bundlesNames.size() ;

  // Read predicted labels
  labelsByName.resize( nbFibers, std::vector<std::string>() ) ;

  std::ifstream labelsFile ;
  labelsFile.open( labelsBinaryFilename ) ;
  if ( labelsFile.fail() )
  {

    std::cout << "Problem reading file : " << labelsBinaryFilename
                                           << std::endl ;

    exit( 1 ) ;

  }
  while ( std::getline( labelsFile, line ) )
  {

    std::vector< std::string > out ;
    std::stringstream ss( line ) ;
    std::string s ;
    while ( std::getline( ss, s, delim ) )
    {

      s.erase( std::remove( s.begin(), s.end(), ' ' ), s.end() ) ;
      out.push_back( s ) ;

    }

    int labelFiber = stoi( out[ 1 ] ) ;


    if ( labelFiber > nbLabels )
    {

      std::stringstream outMessageOss ;
      outMessageOss << "ioWrapper->readLabelsWithDict : labelFiber = "
                    << labelFiber << " is greater than the labels in the "
                    << "dictionary = " << bundlesNames.size() << std::endl ;
      std::string outMessage = outMessageOss.str() ;

      throw( std::invalid_argument( outMessage ) ) ;

    }

    std::string labelName ;
    if ( labelFiber == -1 )
    {

      labelName = "Unlabeled" ;


    }
    else
    {

      labelName = bundlesNames[ labelFiber ] ;

    }


    labelsByName[ stoi( out[ 0 ] ) ].push_back( labelName ) ;

  }

  labelsFile.close() ;

}


////////////////////////////////////////////////////////////////////////////////
////// Function to read labels into a list of vecotr of bundles names //////////
////////////////////////////////////////////////////////////////////////////////
void readLabelsWithDictSupWMA(
                            const char* labelsDictFilename,
                            const char* labelsBinaryFilename,
                            std::vector<std::vector<std::string>>& labelsByName,
                            int nbFibers )
{

  // Read dictionary
  std::vector<std::string> bundlesNames ;

  const char delim = ':' ;
  std::string line ;
  std::ifstream dictFile ;
  dictFile.open( labelsDictFilename ) ;
  if ( dictFile.fail() )
  {

    std::cout << "Problem reading file : " << labelsDictFilename
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

    bundlesNames.push_back( out[ 0 ] ) ;

  }

  dictFile.close() ;

  int nbLabels = bundlesNames.size() ;

  // Read predicted labels
  labelsByName.resize( nbFibers, std::vector<std::string>() ) ;

  std::ifstream file ;
  file.open( labelsBinaryFilename, std::ios::binary | std::ios::in ) ;
  if ( file.fail() )
  {

    std::cout << "Problem reading file : " << labelsBinaryFilename <<
                                                                     std::endl ;
    exit( 1 ) ;

  }

  for ( int fiber = 0 ; fiber < nbFibers ; fiber++ )
  {

    int16_t _label ;

    file.read( reinterpret_cast<char*>( &( _label ) ), sizeof( int16_t ) ) ;
    labelsByName[ fiber ].push_back( bundlesNames[ _label ] ) ;

  }

  file.close() ;

}

////////////////////////////////////////////////////////////////////////////////
/////////////////// Get label number from name and dict ////////////////////////
////////////////////////////////////////////////////////////////////////////////
int getLabelFromName( const std::vector<std::string>& bundlesDict,
                      const std::string& bundleName )
{


  int nbBundles = bundlesDict.size() ;

  if ( nbBundles <= 0 )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ioWrapper->getLabelFromName : the size of bundlesDict "
                  << " must be strictly positive " << nbBundles << std::endl ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    if ( bundlesDict[ i ] == bundleName )
    {

      return i ;

    }

  }

  return -1 ;

}

////////////////////////////////////////////////////////////////////////////////
//////////////////// Read index in tractogram neighbors ////////////////////////
////////////////////////////////////////////////////////////////////////////////
void readIndexInTractogram( const char* predictedLabelsFilename,
                            std::vector<int64_t>& predictedLabels,
                            int nbFibers )
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
/////////////// Function to get vector with altas bundles paths ////////////////
////////////////////////////////////////////////////////////////////////////////
std::string getAtlasBunldesPaths( const std::string& outputDirectory,
                                  const std::string& atlasDirectory,
                                  const std::string& format,
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

    if ( endswith( tmpBundlesDataFilename, format ) )
    {

      if ( format == ".bundles" )
      {

        tmpBundlesDataFilename = replaceExtension( tmpBundlesFilename,
                                                              ".bundlesdata" ) ;


      }
      else if ( format == ".bundlesdata" )
      {

        tmpBundlesFilename = replaceExtension( tmpBundlesFilename,
                                                                  ".bundles" ) ;


      }
      else if ( format == ".trk" || format == ".tck" )
      {

        tmpBundlesFilename = replaceExtension( tmpBundlesFilename, ".minf" ) ;

      }

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
        _bundlesPathOss << tmpBundleDir << bundleName ;
        if ( format == ".bundles" || format == ".bundlesdata" )
        {

          _bundlesPathOss << ".bundles" ;


        }
        else if ( format == ".trk"  || format == ".tck" )
        {

          _bundlesPathOss << ".minf" ;


        }

        std::string _bundlesPath = _bundlesPathOss.str() ;

        std::ostringstream _bundlesdataPathOss ;
        _bundlesdataPathOss << tmpBundleDir << bundleName ;
        if ( format == ".bundles" || format == ".bundlesdata" )
        {

          _bundlesdataPathOss << ".bundlesdata" ;


        }
        else if ( format == ".trk" || format == ".tck" )
        {

          _bundlesdataPathOss << format ;


        }

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
                              const std::string& inputDirectory,
                              const std::vector<std::string>& atlasBundlesPaths,
                              const std::string& format,
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
    tmpBundlesPathOss << inputDirectory << bundleName ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      tmpBundlesPathOss << ".bundles" ;

    }
    else if ( format == ".trk" || format == ".tck" )
    {

      tmpBundlesPathOss << ".minf" ;

    }


    std::string tmpBundlesPath = tmpBundlesPathOss.str() ;

    std::ostringstream tmpBundlesDataPathOss ;
    tmpBundlesDataPathOss << inputDirectory << bundleName ;
    if ( format == ".bundles" || format == ".bundlesdata" )
    {

      tmpBundlesDataPathOss << ".bundlesdata" ;
    }
    else if ( format == ".trk" || format == ".tck" )
    {

      tmpBundlesDataPathOss << format ;
    }

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
                             const std::string& format,
                             std::vector<std::string>& neighborhoodFilenames )
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
    tmpNeighborhoodPathOss << tmpNeighborhoodDir << bundleName << format ;
    std::string tmpNeighborhoodPath = tmpNeighborhoodPathOss.str() ;

    if ( is_file( tmpNeighborhoodPath ) )
    {

      neighborhoodFilenames.push_back( tmpNeighborhoodPath ) ;

    }
    else
    {

      std::cout << "ERROR : neighborhood " << tmpNeighborhoodPath
                << " does not exists " << std::endl ;
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
                           bool flip_x,
                           bool flip_y,
                           bool flip_z,
                           bool force,
                           int verbose )
{

  std::string convertBundleFormatsFile = "ConvertBundleFormat" ;

  std::ostringstream commandOss ;
  commandOss << convertBundleFormatsFile << "  "
             << "-i " << inputBundles << " "
             << "-o " << outputTrk << " "
             << "-r " << referenceImage << " " ;
  
  if ( flip_x )
  {

    commandOss << "-x " ;

  }

  if ( flip_y )
  {

    commandOss << "-y " ;

  }

  if ( flip_x )
  {

    commandOss << "-z " ;

  }

  commandOss << "-v " ;
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
                                   const std::vector<float>& disimilarities,
                                   const std::vector<int>& nbFibersBundles,
                                   const std::vector<std::string>& bundlesNames,
                                   const char* fileName )
{

  std::ofstream file ;
  file.open( fileName ) ;
  if ( file.fail() )
  {

    std::cout << "Problem opening file to write : " << fileName << std::endl ;
    exit( 1 ) ;

  }

  int bundlesNamesSize = bundlesNames.size() ;
  bool isUnlabeled = false ;
  bool isRegrouped = false ;
  for ( int _bundle = 0 ; _bundle < bundlesNames.size() ; _bundle++ )
  {

    if ( bundlesNames[ _bundle ] == "unlabeledFibers"  )
    {

      bundlesNamesSize -= 1 ;
      isUnlabeled = true ;

    }

    if ( bundlesNames[ _bundle ]  == "regroupedRecognized" )
    {

       bundlesNamesSize -= 1 ;
       isRegrouped = true ;

    }

    if ( isUnlabeled && isRegrouped )
    {

      break ;

    }

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
  if ( bundlesNamesSize != nbBundles )
  {

    std::cout << "ERROR in saveComparisonMeasuresWithAtlas : the number of "
              << "bundles in bundlesNames is different than the number of bundles"
              << " with coverage " << std::endl ;
    exit( 1 ) ;

  }

  file << "Bundle_Name\tCoverage\tAdjacency\tOverlap\tDisimilarity\tNbFibers"
                                                                  << std::endl ;

  int _bundle = 0 ;
  for ( int i = 0 ; i < bundlesNames.size() ; i++ )
  {

    std::string _bundleName = bundlesNames[ i ] ;
    if ( _bundleName != "unlabeledFibers" &&
		                                     _bundleName != "regroupedRecognized"  )
    {

      float _coverage = coveragesBundles[ _bundle ] ;
      float _adjacency = adjacencyBundles[ _bundle ] ;
      float _overlap = overlapBundles[ _bundle ] ;
      float _disimilarity = disimilarities[ _bundle ] ;
      float _nbFibers = nbFibersBundles[ _bundle ] ;

      file << _bundleName << "\t" << _coverage << "\t" << _adjacency << "\t"
           << _overlap << "\t" << _disimilarity << "\t" << _nbFibers
                                                                  << std::endl ;

      _bundle += 1 ;

    }

  }

  file.close() ;


}



////////////////////////////////////////////////////////////////////////////////
///////////////// Function to read comparison with atlas .tsv //////////////////
////////////////////////////////////////////////////////////////////////////////
void readComparisonMeasuresWithAtlas(
                                       const char* fileName,
                                        std::vector<float>& coveragesBundles,
                                        std::vector<float>& adjacencyBundles,
                                        std::vector<float>& overlapBundles,
                                        std::vector<float>& disimilarities,
                                        std::vector<int>& nbFibersBundles,
                                        std::vector<std::string>& bundlesNames )
{

  std::string fileNameStr = fileName ;

  if ( !( endswith( fileName, ".tsv" ) ) )
  {

    std::stringstream outMessageOss ;
    outMessageOss << "ioWrapper->readComparisonMeasuresWithAtlas : fileName "
                  << fileName << " is not .tsv" ;
    std::string outMessage = outMessageOss.str() ;

    throw( std::invalid_argument( outMessage ) ) ;

  }

  const char delim = '\t' ;
  std::string line ;
  std::ifstream tsvFile ;
  tsvFile.open( fileName ) ;
  if ( tsvFile.fail() )
  {

    std::cout << "Problem reading file : " << fileName
                                           << std::endl ;

    exit( 1 ) ;

  }
  while ( std::getline( tsvFile, line ) )
  {

    std::vector< std::string > out ;
    std::stringstream ss( line ) ;
    std::string s ;
    while ( std::getline( ss, s, delim ) )
    {

      s.erase( std::remove( s.begin(), s.end(), ' ' ), s.end() ) ;
      s.erase( std::remove( s.begin(), s.end(), '\n' ), s.end() ) ;
      out.push_back( s ) ;

    }

    bundlesNames.push_back( out[ 0 ] ) ;
    coveragesBundles.push_back( stof( out[ 1 ] ) ) ;
    adjacencyBundles.push_back( stof( out[ 2 ] ) ) ;
    overlapBundles.push_back( stof( out[ 3 ] ) ) ;
    disimilarities.push_back( stof( out[ 4 ] ) ) ;
    nbFibersBundles.push_back( stof( out[ 5 ] ) ) ;


  }

  tsvFile.close() ;


}

////////////////////////////////////////////////////////////////////////////////
/////////// Function to get coverage with atlas of extracted bundle ////////////
////////////////////////////////////////////////////////////////////////////////
float getCoverageWithAtlas( const std::string& bundleFilename )
{


  // BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  BundlesMinf bundle( bundleFilename.c_str() ) ;
  if ( bundle.coverageWithAtlas < 0 )
  {

    std::cout << "ERROR : got invalid coverage of " << bundle.coverageWithAtlas
              << " for file " << bundleFilename << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.coverageWithAtlas ) ;

}

////////////////////////////////////////////////////////////////////////////////
/////////// Function to get adjacency with atlas of extracted bundle ///////////
////////////////////////////////////////////////////////////////////////////////
float getAdjacencyWithAtlas( const std::string& bundleFilename )
{


  BundlesMinf bundle( bundleFilename.c_str() ) ;
  // BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.adjacencyWithAtlas < 0 )
  {

    std::cout << "ERROR : got invalid adjacency of "
              << bundle.adjacencyWithAtlas << " for file " << bundleFilename
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


  BundlesMinf bundle( bundleFilename.c_str() ) ;
  // BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.overlapWithAtlas < 0 )
  {

    std::cout << "ERROR : got invalid overlap of " << bundle.overlapWithAtlas
              << " for file " << bundleFilename << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.overlapWithAtlas ) ;

}



////////////////////////////////////////////////////////////////////////////////
/////////////// Function to get average radius of atlas  bundle ////////////////
////////////////////////////////////////////////////////////////////////////////
float getAverageRadiusAtlasBundle( const std::string& bundleFilename )
{

  BundlesMinf bundle( bundleFilename.c_str() ) ;
  // BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.averageRadius <= 0 )
  {

    std::cout << "ERROR : got invalid radius of " << bundle.averageRadius
              << " for file " << bundleFilename << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.averageRadius ) ;

}

////////////////////////////////////////////////////////////////////////////////
/// Function to get average distance between medial points of atlas  bundle ////
////////////////////////////////////////////////////////////////////////////////
float getAverageDistanceBetweenMedialPoints( const std::string& bundleFilename )
{

  BundlesMinf bundle( bundleFilename.c_str() ) ;
  // BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.averageDistanceBetweenMedialPoints <= 0 )
  {

    std::cout << "ERROR : got invalid radius of "
              << bundle.averageDistanceBetweenMedialPoints
              << " for file " << bundleFilename << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.averageDistanceBetweenMedialPoints ) ;

}

////////////////////////////////////////////////////////////////////////////////
/// Function to get maximum distance between medial points of atlas  bundle ////
////////////////////////////////////////////////////////////////////////////////

float getMaximumDistanceBetweenMedialPoints( const std::string& bundleFilename )
{

  BundlesMinf bundle( bundleFilename.c_str() ) ;
  // BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.maxDistanceBetweenMedialPoints <= 0 )
  {

    std::cout << "ERROR : got invalid radius of "
              << bundle.maxDistanceBetweenMedialPoints
              << " for file " << bundleFilename << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.maxDistanceBetweenMedialPoints ) ;

}



////////////////////////////////////////////////////////////////////////////////
//////////////// Function to get number of fibers atlas bundle /////////////////
////////////////////////////////////////////////////////////////////////////////
int getNbFibers( const std::string& bundleFilename )
{

  BundlesMinf bundle( bundleFilename.c_str() ) ;
  // BundlesFormat bundle( bundleFilename.c_str(), 0 ) ;
  if ( bundle.curves_count <= 0 )
  {

    std::cout << "ERROR : got invalid fiber count of " << bundle.curves_count
              << " for file " << bundleFilename << std::endl ;
    exit( 1 ) ;

  }
  return( bundle.curves_count ) ;

}


////////////////////////////////////////////////////////////////////////////////
///////////////////////// Function to get memory usage /////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Based on https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0 * 1e-6 ;
   resident_set = rss * page_size_kb;
}
