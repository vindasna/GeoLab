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

    int index_fixed, index_moving, index_output ;
    int index_qbThr, index_nbIter, index_outlierRejectionThr ;
    int index_verbose, index_help ;

    index_fixed = getFlagPosition( argc, argv, "-f" ) ;
    index_moving = getFlagPosition( argc, argv, "-m" ) ;
    index_output = getFlagPosition( argc, argv, "-o" ) ;
    int index_xfm = getFlagPosition( argc, argv, "-xfm" ) ;
    index_qbThr = getFlagPosition( argc, argv, "-qbThr" ) ;
    index_nbIter = getFlagPosition( argc, argv, "-nbIter" ) ;
    index_outlierRejectionThr = getFlagPosition( argc, argv, "-outlierRejectionThr" ) ;
    index_verbose = getFlagPosition( argc, argv, "-v" ) ;
    index_help = getFlagPosition( argc, argv, "-h" ) ;

    if ( index_help > 0 )
    {
        std::cout << "Fast Streamline-based Linear Registration (SLR) Command :\n"
                  << "-f : Fixed tractogram path (Target)\n"
                  << "-m : Moving tractogram path (Source)\n"
                  << "-o : Output directory to save All computed files (matrix, centroids, moved tractogram)\n"
                  << "[-xfm] : Type of transformation (rigid or affine, default = rigid)\n"
                  << "[-qbThr] : Distance threshold (mm) for QuickBundles centroids (default = 10.0)\n"
                  << "[-nbIter] : Number of iterate steps (iterations) for ICP (default = 10)\n"
                  << "[-outlierRejectionThr] : Ratio of worst Centroids matches to reject (0.0 - 0.99, default = 0.20)\n"
                  << "[-v] : Set verbosity level at 1\n"
                  << "[-h] : Show this message" << std::endl ;
        exit( 1 ) ;
    }

    if ( !index_fixed ) { std::cout << "-f argument required ..." << std::endl ; exit( 1 ) ; }
    if ( !index_moving ) { std::cout << "-m argument required ..." << std::endl ; exit( 1 ) ; }
    if ( !index_output ) { std::cout << "-o argument required ..." << std::endl ; exit( 1 ) ; }

    TransformType regType = TransformType::RIGID ;
    float qbThr = 10.0f ;
    int nbIter = 10 ;
    float outlierRejectionThr = 0.20f ;
    int verbose = 0 ;

    if ( index_xfm ) {
        std::string typeStr = argv[ index_xfm + 1 ] ;
        if ( typeStr == "affine" || typeStr == "Affine" ) regType = TransformType::AFFINE ;
        else if ( typeStr == "rigid" || typeStr == "Rigid" ) regType = TransformType::RIGID ;
        else { std::cout << "Invalid -xfm argument: " << typeStr << ", using default Rigid" << std::endl; }
    }
    if ( index_qbThr ) qbThr = std::stof( argv[ index_qbThr + 1 ] ) ;
    if ( index_nbIter ) nbIter = std::stoi( argv[ index_nbIter + 1 ] ) ;
    if ( index_outlierRejectionThr ) outlierRejectionThr = std::stof( argv[ index_outlierRejectionThr + 1 ] ) ;
    if ( index_verbose ) verbose = 1 ;

    std::string fixedPath = argv[ index_fixed + 1 ] ;
    std::string movingPath = argv[ index_moving + 1 ] ;
    std::string outputDir = argv[ index_output + 1 ] ;

    if (!is_file(fixedPath)) { std::cout << "Fixed file does not exist: " << fixedPath << std::endl; exit(1); }
    if (!is_file(movingPath)) { std::cout << "Moving file does not exist: " << movingPath << std::endl; exit(1); }

    // Create output directory if it doesn't exist
    char lastChar = outputDir[ outputDir.size() - 1 ] ;
    if ( lastChar != '/' ) { outputDir = outputDir + '/' ; }
    if ( !is_dir( outputDir ) ) { mkdir( outputDir ) ; }

    if (verbose)
    {
        std::cout << "Fixed Tractogram : " << fixedPath << std::endl ;
        std::cout << "Moving Tractogram: " << movingPath << std::endl ;
        std::cout << "Output Directory : " << outputDir << std::endl ;
        std::cout << "Transformation   : " << (regType == TransformType::RIGID ? "Rigid" : "Affine") << std::endl ;
        std::cout << "qbThreshold      : " << qbThr << " mm" << std::endl ;
        std::cout << "Iterations       : " << nbIter << std::endl ;
        std::cout << "Outlier Rejection: " << outlierRejectionThr * 100.0f << " %" << std::endl ;
    }

    std::cout << "Loading Fixed tractogram..." << std::endl ;
    BundlesData fixedBundles( fixedPath.c_str() ) ;
    std::cout << "Loading Moving tractogram..." << std::endl ;
    BundlesData movingBundles( movingPath.c_str() ) ;

    BundlesData originalMovingCentroids;
    BundlesData movingCentroids;
    BundlesData targetCentroids;

    std::cout << "Starting position registration..." << std::endl ;
    Eigen::Matrix4f transformMatrix = movingBundles.registerFast( fixedBundles, regType, qbThr, nbIter, outlierRejectionThr, &originalMovingCentroids, &movingCentroids, &targetCentroids ) ;
    
    // 1. Save Transformation Matrix
    std::string xfmPath = outputDir + "transform.txt";
    std::ofstream xfmFile(xfmPath);
    if (xfmFile.is_open()) {
        xfmFile << transformMatrix << std::endl;
        xfmFile.close();
        if (verbose) std::cout << "Saved transform to " << xfmPath << std::endl;
    } else {
        std::cerr << "Error writing transform to " << xfmPath << std::endl;
    }

    // Use formats from loaded copies.
    BundlesMinf movingInfo( movingPath.c_str() ) ;
    BundlesMinf fixedInfo( fixedPath.c_str() ) ;

    std::string extMoving = movingInfo.getFormat();
    std::string extFixed = fixedInfo.getFormat();

    // 2. Save Centroids
    std::cout << "Saving Centroids..." << std::endl;
    originalMovingCentroids.isBundles = movingBundles.isBundles;
    originalMovingCentroids.isTrk = movingBundles.isTrk;
    originalMovingCentroids.isTck = movingBundles.isTck;

    movingCentroids.isBundles = movingBundles.isBundles;
    movingCentroids.isTrk = movingBundles.isTrk;
    movingCentroids.isTck = movingBundles.isTck;

    targetCentroids.isBundles = fixedBundles.isBundles;
    targetCentroids.isTrk = fixedBundles.isTrk;
    targetCentroids.isTck = fixedBundles.isTck;

    std::string originalCentroidsPath = outputDir + "original_moving_centroids" + extMoving;
    std::string movedCentroidsPath = outputDir + "moved_centroids" + extMoving;
    std::string targetCentroidsPath = outputDir + "target_centroids" + extFixed;

    BundlesMinf movingCentroidsInfo = movingInfo;
    movingCentroidsInfo.curves_count = movingCentroids.curves_count;

    originalMovingCentroids.write(originalCentroidsPath.c_str(), movingCentroidsInfo);
    movingCentroids.write(movedCentroidsPath.c_str(), movingCentroidsInfo);

    BundlesMinf targetCentroidsInfo = fixedInfo;
    targetCentroidsInfo.curves_count = targetCentroids.curves_count;
    targetCentroids.write(targetCentroidsPath.c_str(), targetCentroidsInfo);

    // 3. Save Registered Moved Tractogram
    std::string movedTractogramPath = outputDir + "moved_tractogram" + extMoving;
    std::cout << "Saving registered tractogram to " << movedTractogramPath << std::endl ;
    movingBundles.write( movedTractogramPath.c_str(), movingInfo ) ;

    const std::chrono::duration< double > duration = std::chrono::system_clock::now() - start_time ;
    std::cout << "Fast registration complete in " << duration.count() << " seconds." << std::endl ;

    return 0 ;
}
