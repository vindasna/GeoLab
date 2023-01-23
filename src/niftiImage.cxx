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

#include <stdexcept>

#include "niftiImage.h"
#include "ioWrapper.h"


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Constructors /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
NiftiImage::NiftiImage(){}


NiftiImage::NiftiImage( const char* path )
{

  read( path ) ;

}


NiftiImage::NiftiImage( const NiftiImage& niftiImageInfo )
{

  this->vox_to_ras = niftiImageInfo.vox_to_ras ;
  this->size = niftiImageInfo.size ;
  this->resolution = niftiImageInfo.resolution ;

}


NiftiImage::NiftiImage( std::vector<std::vector<float>> vox_to_ras,
                        std::vector<short> size,
                        std::vector<float> resolution )
{

  this->vox_to_ras = vox_to_ras ;
  this->size = size ;
  this->resolution = resolution ;

}


////////////////////////////////// Destructor //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
NiftiImage::~NiftiImage() {}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Methods ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void NiftiImage::read( const char* path )
{

  short dim[ 8 ] ;
  float pixdim[ 8 ] ;
  std::string pathStr = path ;
  if ( endswith( pathStr, ".nii" ) )
  {

    std::ifstream file ;
    file.open( path, std::ios::binary | std::ios::in ) ;
    if ( file.fail() )
    {

      std::ostringstream outMessageOss ;
      outMessageOss << "Problem reading file : " << path << std::endl ;

      throw ( std::invalid_argument( outMessageOss.str() ) ) ;

    }

    file.seekg( 40, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &dim ), sizeof( dim ) ) ;

    file.seekg( 76, std::ios_base::beg ) ;
    file.read( reinterpret_cast<char*>( &pixdim ), sizeof( pixdim ) ) ;

    file.seekg( 280, std::ios_base::beg ) ;
    for ( int i = 0 ; i < 3 ; i ++ )
    {

      this->vox_to_ras[ i ].resize( 4, 0 ) ;

      float buffer[ 4 ] ;
      file.read( reinterpret_cast<char*>( &buffer ), sizeof( buffer ) ) ;
      for ( int j = 0 ; j < 4 ; j++ )
      {

        this->vox_to_ras[ i ][ j ] = buffer[ j ] ;

      }

    }
    this->vox_to_ras[ 3 ].resize( 4, 0.0f ) ;
    this->vox_to_ras[ 3 ][ 3 ] = 1.0f ;

    file.close() ;

  }
  else
  {

    const char* outMessage = "NiftiImage::read : The only format supported is "
                                                                      ".nii\n" ;
    throw( std::invalid_argument( outMessage ) ) ;

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    this->resolution[ i ] = pixdim[ i + 1 ] ;

  }

  for ( int i = 0 ; i < 3 ; i++ )
  {

    this->size[ i ] = dim[ i + 1 ] ;

  }


}
