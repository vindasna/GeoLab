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

#include <map>

#include <boost/algorithm/string.hpp>


#include "getFileForPlotsRegression.h"




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
////////////////////////// Fonction to read .tsv info //////////////////////////
////////////////////////////////////////////////////////////////////////////////
void readTsvInfo( const char* tsvFilePath,
                  std::map< std::string, std::vector< std::string > >& tsvInfo,
                  int verbose )
{

  if ( verbose )
  {

    std::cout << "Reading " << tsvFilePath << std::endl ;

  }

  std::ifstream inFile( tsvFilePath ) ;
  std::string line;
  while ( getline( inFile, line ) )
  {
    // Split line into tab-separated parts
    std::vector< std::string > lineSplit;
    split( lineSplit, line, boost::is_any_of( "\t" ) ) ;

    // tsvInfo.push_back( lineSplit ) ;
    std::vector< std::string >::iterator newEnd;
    newEnd = remove_if( lineSplit.begin(), lineSplit.end(),
                                  []( std::string s ){ return( s == "" ) ; } ) ;
    lineSplit.erase( newEnd, lineSplit.end( ) );
    std::string participantID = lineSplit[ 0 ] ;
    std::vector< std::string > participantInfo =
                                    { lineSplit.begin() + 1, lineSplit.end() } ;
    tsvInfo[ participantID ] = participantInfo ;

  }

  inFile.close() ;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////// Fonction to read .tsv info //////////////////////////
////////////////////////////////////////////////////////////////////////////////
void saveOutMapDict(
       const char* savePath,
       const std::map< std::string, std::vector< std::vector< std::string > > >&
                                                                         outMap,
       const std::vector< std::string >& bundlesNames,
       int verbose )
{

  std::ofstream file( savePath ) ;
  if ( !file )
  {

    std::cout << "Cannot save file, there's a problem with the saving path : "
              << savePath << std::endl ;
    exit( 1 );

  }

  int nbBundles = outMap.size() ;
  if ( nbBundles != bundlesNames.size() )
  {

    std::cout << "ERROR saveOutMapDict : the number of bundles in outMap ("
              << nbBundles << ") is different than the number of bundles in "
              << "bundlesNames (" << bundlesNames.size() << ")" << std::endl ;
    exit( 1 ) ;

  }

  for ( int i = 0 ; i < nbBundles ; i++ )
  {

    std::string bundleName = bundlesNames[ i ] ;
    if ( ( outMap.find( bundleName ) == outMap.end() ) )
    {
      std::cout << "ERROR saveOutMapDict : the bundle " << bundleName
                << " is not in the outMap " << std::endl ;
      exit( 1 ) ;

    }

    file << bundleName << ":\n" ;

    std::vector< std::vector< std::string > > bundleInfo =
                                                       outMap.at( bundleName ) ;

    int nbSubjectsWithBundle = std::stoi( bundleInfo[ 3 ][ 0 ] ) ;
    for ( int j = 0 ; j < nbSubjectsWithBundle ; j++ )
    {

      std::string subjectAge = bundleInfo[ 0 ][ j ] ;
      std::string meanMeasureValue = bundleInfo[ 1 ][ j ] ;
      std::string subject = bundleInfo[ 2 ][ j ] ;

      file << subjectAge << "\t" << meanMeasureValue << "\t" << subject
                                                                       << "\n" ;

    }

    file << "\n\n" ;


  }

  file.close() ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////// Fonction to get age subject //////////////////////////
////////////////////////////////////////////////////////////////////////////////
std::string getAgeSubject(
             const std::string subjectId,
             const std::map< std::string, std::vector< std::string > >& tsvInfo,
             int verbose )
{
  std::string ageSubject ;

  // Cannot use the [] operator because map passed as const
  if ( !( tsvInfo.find( subjectId ) == tsvInfo.end() ) )
  {

    ageSubject = tsvInfo.at( subjectId )[ 2 ] ;

  }
  else
  {

    std::cout << "ERROR in getAgeSubject : subjectID -> " << subjectId
              << " not in .tsv info file " << std::endl ;
    exit( 1 ) ;

  }

  return( ageSubject ) ;

}

///////////////////////////////////// Main /////////////////////////////////////
int main( int argc, char* argv[] )
{

  int index_ukb_dti, index_output, index_measure, index_tsv, index_verbose,
                                                                    index_help ;
  index_ukb_dti =   getFlagPosition( argc, argv, "-ukb-dti" ) ;
  index_measure =   getFlagPosition( argc, argv, "-m" ) ;
  index_tsv =   getFlagPosition( argc, argv, "-tsv" ) ;
  index_output =   getFlagPosition( argc, argv, "-o" ) ;
  index_verbose =   getFlagPosition( argc, argv, "-v" ) ;
  index_help =   getFlagPosition( argc, argv, "-h" ) ;

  if ( index_help > 0 )
  {

    std::cout << "Function to convert bundles format : " << std::endl
              << "-ukb-dti : UKBiobank directory with dti measures "
              << std::endl
              << "-m : Name of the measure to plot (options : FA, MD, OD, "
              << "ICVF, ISOVF) " << std::endl
              << "-tsv : tsv file with info about subjects " << std::endl
              << "-o : Out file name " << std::endl
              << "[-v] : set verbosity level at 1 " << std::endl
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }


  std::string ukbDtiDirectory ;
  if ( !index_ukb_dti )
  {

    std::cout << "-ukb-dti argument required ..." << std::endl ;
    exit( 1 ) ;

  }
  else
  {

    ukbDtiDirectory = argv[ index_ukb_dti + 1 ] ;
    char lastChar = ukbDtiDirectory[ ukbDtiDirectory.size() - 1 ] ;
    if ( lastChar != '/' )
    {

      ukbDtiDirectory = ukbDtiDirectory + "/" ;

    }

  }

  std::string measureName ;
  if ( !index_measure )
  {

    std::cout << "-m argument required ..." << std::endl ;
    exit( 1 ) ;

  }
  else
  {

    measureName = argv[ index_measure + 1 ] ;
    if ( measureName != "FA" && measureName != "MD" && measureName != "OD"
         && measureName != "ICVF" && measureName != "ISOVF" )
    {

      std::cout << "ERROR : -m can only take the following values -> FA, OD, "
                << "MD, ICVF, ISOVF" << std::endl ;
      exit( 1 ) ;
    }

  }

  std::string tsvFilePath ;
  if ( !index_tsv )
  {

    std::cout << "-tsv argument required ..." << std::endl ;
    exit( 1 ) ;

  }
  else
  {

    tsvFilePath = argv[ index_tsv + 1 ] ;
    char lastChar = tsvFilePath[ tsvFilePath.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      tsvFilePath = tsvFilePath.substr( 0, tsvFilePath.size() - 1 ) ;

    }

  }

  std::string outputFilename ;
  if ( !index_output )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }
  else
  {

    outputFilename = argv[ index_output + 1 ] ;
    char lastChar = outputFilename[ outputFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

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

  /////////////////////////// Reading .tsv info file ///////////////////////////
  // std::vector< std::vector< std::string > > tsvInfo ;
  std::map< std::string, std::vector< std::string > > tsvInfo ;
  readTsvInfo( tsvFilePath.c_str(), tsvInfo, verbose ) ;


  ////////////////////// Getting subjects in ukb-dti dir ///////////////////////
  std::vector< std::string > subjectsIdUkbDti ;
  std::vector< std::string > subjectsPathUkbDti ;
  for ( const auto & file : std::experimental::filesystem::directory_iterator(
                                                    ukbDtiDirectory.c_str() ) )
  {

    std::string tmpSubjectPath = file.path() ;
    std::string key( "sub-" ) ;

    if ( tmpSubjectPath.find( key ) != std::string::npos )
    {

      subjectsPathUkbDti.push_back( tmpSubjectPath ) ;

      std::size_t found = tmpSubjectPath.rfind( key ) ;
      std::string tmpSubjectID = { tmpSubjectPath.begin() + found +key.length(),
                                                        tmpSubjectPath.end() } ;

      subjectsIdUkbDti.push_back( tmpSubjectID ) ;

    }

  }

  int nbSubjects = subjectsIdUkbDti.size() ;


  // Map to store output
  std::map< std::string, std::vector< std::vector< std::string > > > outMap ;

  int count_subjects_with_data = 0 ;
  std::vector< std::string > bundlesNames ;
  std::cout << "Getting data..." << std::endl ;
  for( int i = 0; i < nbSubjects ; i++ )
  {
    std::string subjectId = subjectsIdUkbDti[ i ] ;
    std::string measurePath = ukbDtiDirectory + "/sub-" + subjectId + "/" +
                                                          measureName + ".txt" ;
    printf( "\rProcessing : [ %d / %d ]", i + 1, nbSubjects ) ;
    fflush( stdout ) ;

    std::string subjectAge = getAgeSubject( subjectId, tsvInfo, verbose ) ;

    const char delim = '\t' ;
    std::string line ;
    std::ifstream measureFile ;
    measureFile.open( measurePath ) ;
    if ( measureFile.fail() )
    {

      std::cout << "Problem reading file : " << measurePath << std::endl ;
      exit( 1 ) ;

    }
    while ( std::getline( measureFile, line ) )
    {

      std::vector< std::string > out ;
      std::stringstream ss( line ) ;
      std::string s ;
      while ( std::getline( ss, s, delim ) )
      {

        out.push_back( s ) ;


      }

      std::string key( ".txt" ) ;
      std::size_t found  = out[ 0 ].rfind( key ) ;
      if ( found != std::string::npos )
      {

        out[ 0 ].replace( found, key.length(), "" ) ;

      }


      if ( std::find( bundlesNames.begin(), bundlesNames.end(), out[ 0 ] ) ==
                                                            bundlesNames.end() )
      {

        bundlesNames.push_back( out[ 0 ] ) ;

      }

      if ( outMap.find( out[ 0 ] ) == outMap.end() )
      {

        outMap[ out[ 0 ] ] = { { subjectAge }, { out[ 1 ] },
                                             { "sub-" + subjectId }, { "1" } } ;

      }
      else
      {

        outMap[ out[ 0 ] ][ 0 ].push_back( subjectAge ) ;
        outMap[ out[ 0 ] ][ 1 ].push_back( out[ 1 ] ) ;
        outMap[ out[ 0 ] ][ 2 ].push_back( "sub-" + subjectId ) ;
        std::string _count = outMap[ out[ 0 ] ][ 3 ][ 0 ] ;
        outMap[ out[ 0 ] ][ 3 ][ 0 ] = std::to_string(
                                                     std::stoi( _count ) + 1 ) ;

      }

    }

  }


  saveOutMapDict( outputFilename.c_str(), outMap, bundlesNames, verbose ) ;

  std::cout << "\nDone" << std::endl ;


}
