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

#include "RecognizedBundles.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int verbose = 0 ;
bool useMDF = false ;

////////////////////////////////////////////////////////////////////////////////

int getFlagPosition( int argc, char* argv[], const std::string& flag ) ;

void computeCenterAtlasBundleFibers(
                           BundlesData& atlasBundleData,
                           std::vector<float>& medialPointsAtlasBundleFibers ) ;

void computeAverageFiberBundle(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::vector<float>& medialPointAtlasBundle ) ;

void computeGravityCenterAtlasBundle(
                                BundlesData& atlasBundleData,
                                int nbPoints,
                                std::vector<float>& gravityCenterAtlasBundle ) ;

void computeDistancesToCenterBundle(
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        const std::vector<float>& medialPointAtlasBundle,
                        int nbFibersAtlasBundle,
                        std::vector<float>& distancesToCenterAtlasBundle ) ;

void computeNormalVectorFibersAtlasBundle(
                                     BundlesData& atlasBundleData,
                                     std::vector<float>& normalVectorsBundle ) ;

void computeDistancesBetweenMedialPointsBundle(
                      BundlesData& atlasBundleData,
                      const std::vector<float>& medialPointsAtlasBundleFibers,
                      std::vector<float>& distancesBetweenMedialPointsBundle ) ;

void computeAnglesBundle( BundlesData& atlasBundleData,
                          const std::vector<float>& normalVectorsBundle,
                          std::vector<float>& anglesAtlasBundle ) ;


void computeDirectionAnglesBundle(
                               BundlesData& atlasBundleData,
                               const std::vector<float>& normalVectorsBundle,
                               int nbPoints,
                               std::vector<float>& directionAnglesAtlasBundle ) ;

void computeShapeAnglesBundle(
                        BundlesData& atlasBundleData,
                        int nbPoints,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        std::vector<float>& shapeAnglesAtlasBundle ) ;

void computeAverageDisimilarity(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& normalVectorsBundle,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& disimilaritiesAtlasBundle ) ;

void computeAverageDisimilarityMDF(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& normalVectorsBundle,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& disimilaritiesAtlasBundle ) ;
//
