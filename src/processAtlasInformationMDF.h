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

#include "RecognizedBundles.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct niiFormat
{

  int sizeof_hdr = 348 ; // Must be 348 (bytes)
  char data_type[ 10 ] = "None" ;
  char db_name[ 18 ] = "None" ;
  int extents = 0 ;
  short session_error = 0 ;
  char regular ;
  char dim_info ;
  short dim[ 8 ] = { 0 } ;
  float intent_p1 = 0 ;
  float intent_p2 = 0 ;
  float intent_p3 = 0 ;
  short intent_code = 0 ;
  short datatype = 0 ;
  short bitpix = 0 ;
  short slice_start = 0 ;
  float pixdim[ 8 ] = { 0 } ;
  float vox_offset = 0 ;
  float scl_slope = 0 ;
  float scl_inter = 0 ;
  short slice_end = 0 ;
  char slice_code ;
  char xyzt_units ;
  float cal_max = 0 ;
  float cal_min = 0 ;
  float slice_duration = 0 ;
  float toffset = 0 ;
  int glmax = 0 ;
  int gmin = 0 ;
  char descrip[ 80 ] = "None" ;
  char aux_file[ 24 ] = "None" ;
  short qform_code = 0 ;
  short sform_code = 0 ;
  float quatern_b = 0 ;
  float quatern_c = 0 ;
  float quatern_d = 0 ;
  float qoffset_x = 0 ;
  float qoffset_y = 0 ;
  float qoffset_z = 0 ;
  float srow_x[ 4 ] = { 0 } ;
  float srow_y[ 4 ] = { 0 } ;
  float srow_z[ 4 ] = { 0 } ;
  char intent_name[ 16 ] = "None" ;
  char magic[ 4 ] = "N" ;

};

int verbose = 0 ;
bool useMDFDistance = false ;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;


void computeCenterAtlasBundleFibers(
                           BundlesDataFormat& atlasBundleData,
                           std::vector<float>& medialPointsAtlasBundleFibers ) ;

void computeAverageFiberBundle(
                        BundlesDataFormat& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::vector<float>& medialPointAtlasBundle ) ;

void computeDistancesToCenterBundle(
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        const std::vector<float>& medialPointAtlasBundle,
                        int nbFibersAtlasBundle,
                        std::vector<float>& distancesToCenterAtlasBundle ) ;

void computeNormalVectorFibersAtlasBundle(
                                     BundlesDataFormat& atlasBundleData,
                                     std::vector<float>& normalVectorsBundle ) ;

void computeAnglesBundle( BundlesDataFormat& atlasBundleData,
                          const std::vector<float>& normalVectorsBundle,
                          std::vector<float>& anglesAtlasBundle ) ;


void computeDirectionAnglesBundle(
                               BundlesDataFormat& atlasBundleData,
                               const std::vector<float>& normalVectorsBundle,
                               int nbPoints,
                               std::vector<float>& directionAnglesAtlasBundle ) ;

void computeShapeAnglesBundle(
                        BundlesDataFormat& atlasBundleData,
                        int nbPoints,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        std::vector<float>& shapeAnglesAtlasBundle ) ;

void computeAverageDisimilarity(
                        BundlesDataFormat& atlasBundleData,
                        int nbPoints,
                        std::vector<float>& disimilaritiesAtlasBundle ) ;
