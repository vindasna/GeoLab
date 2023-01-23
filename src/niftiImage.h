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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class NiftiImage
{

  public :
  /////////////////////////////// Public Fields ////////////////////////////////
  std::vector<std::vector<float>> vox_to_ras = std::vector<std::vector<float>>(
                                               4, std::vector<float>( 4, 0 ) ) ;
  std::vector<short int> size = std::vector<short int >( 3, 0 ) ;
  std::vector<float> resolution = std::vector<float>( 3, 0 ) ;

  //////////////////////////////// Constructors ////////////////////////////////
  NiftiImage() ;

  NiftiImage( const char* path ) ;

  NiftiImage( const NiftiImage& niftiImageInfo ) ;

  NiftiImage( std::vector<std::vector<float>> vox_to_ras,
              std::vector<short> size,
              std::vector<float> resolution ) ;

  ///////////////////////////////// Destructor /////////////////////////////////
  virtual ~NiftiImage() ;


  ////////////////////////////////// Methods ///////////////////////////////////
  void read( const char* path ) ;

} ;
