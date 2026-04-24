#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>

#include "bundlesData.h"
#include "ioWrapper.h"

///////////////////////////////////////////////////////////////////////////////
////////// Function to get flag position when parsing arguments ///////////////
///////////////////////////////////////////////////////////////////////////////
int getFlagPosition( int argc, char* argv[], const std::string& flag )
{
    for ( int i = 0 ; i < argc ; i++ )
    {
        std::string arg = argv[ i ] ;
        if ( arg == flag )
        {
            return i ;
        }
    }
    return 0 ;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Main ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{
    const auto start_time = std::chrono::system_clock::now() ;

    int index_input, index_output, index_qbThr, index_verbose, index_help ;

    index_input = getFlagPosition( argc, argv, "-i" ) ;
    index_output = getFlagPosition( argc, argv, "-o" ) ;
    index_qbThr = getFlagPosition( argc, argv, "-qbThr" ) ;
    index_verbose = getFlagPosition( argc, argv, "-v" ) ;
    index_help = getFlagPosition( argc, argv, "-h" ) ;

    if ( index_help > 0 )
    {
        std::cout << "Compute QuickBundles Centroids Command :\n"
                  << "-i : Input tractogram path\n"
                  << "-o : Output centroids tractogram path\n"
                  << "[-qbThr] : Distance threshold (mm) for QuickBundles (default = 10.0)\n"
                  << "[-v] : Set verbosity level at 1\n"
                  << "-h : Show this message" << std::endl ;
        exit( 1 ) ;
    }

    if ( !index_input ) { std::cout << "-i argument required ..." << std::endl ; exit( 1 ) ; }
    if ( !index_output ) { std::cout << "-o argument required ..." << std::endl ; exit( 1 ) ; }

    float qbThr = 10.0f ;
    int verbose = 0 ;

    if ( index_qbThr ) qbThr = std::stof( argv[ index_qbThr + 1 ] ) ;
    if ( index_verbose ) verbose = 1 ;

    std::string inputPath = argv[ index_input + 1 ] ;
    std::string outputPath = argv[ index_output + 1 ] ;

    if (!is_file(inputPath)) { std::cout << "Input file does not exist: " << inputPath << std::endl; exit(1); }

    if (verbose)
    {
        std::cout << "Input Tractogram : " << inputPath << std::endl ;
        std::cout << "Output Path      : " << outputPath << std::endl ;
        std::cout << "qbThreshold      : " << qbThr << " mm" << std::endl ;
    }

    std::cout << "Loading Input tractogram..." << std::endl ;
    BundlesData inputBundles( inputPath.c_str() ) ;

    if (inputBundles.curves_count == 0) {
        std::cout << "Input tractogram has no streamlines. Exiting." << std::endl ;
        exit(1);
    }

    std::cout << "Computing QuickBundles clusters..." << std::endl ;
    std::vector<Cluster> clusters = inputBundles.computeQuickBundles( qbThr ) ;
    std::cout << "Found " << clusters.size() << " clusters." << std::endl ;

    // Create new bundles data for centroids
    int numPoints = inputBundles.pointsPerTrack[0];
    int numClusters = clusters.size();

    std::vector<float> centroidsMatrixTracks;
    centroidsMatrixTracks.reserve(numClusters * 3 * numPoints);
    std::vector<int32_t> pointsPerTrackCentroids(numClusters, numPoints);

    for (const auto& c : clusters) {
        centroidsMatrixTracks.insert(centroidsMatrixTracks.end(), c.centroid.begin(), c.centroid.end());
    }

    BundlesData centroidsBundles;
    centroidsBundles.matrixTracks = std::move(centroidsMatrixTracks);
    centroidsBundles.pointsPerTrack = std::move(pointsPerTrackCentroids);
    centroidsBundles.curves_count = numClusters;
    centroidsBundles.isBundles = inputBundles.isBundles;
    centroidsBundles.isTrk = inputBundles.isTrk;
    centroidsBundles.isTck = inputBundles.isTck;

    // Save output
    std::cout << "Saving centroids tractogram to " << outputPath << std::endl ;
    
    // Copy headers/minf info
    BundlesMinf outputInfo( inputPath.c_str() ) ;
    outputInfo.curves_count = numClusters;
    centroidsBundles.write( outputPath.c_str(), outputInfo ) ;

    const std::chrono::duration< double > duration = std::chrono::system_clock::now() - start_time ;
    std::cout << "Centroids computation complete in " << duration.count() << " seconds." << std::endl ;

    return 0 ;
}
