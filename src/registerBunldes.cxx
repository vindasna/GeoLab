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
#include <future>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include <Eigen/Core>
#include <LBFGS.h>

#include "./registerBunldes.h"

using namespace std::chrono_literals ;


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
/////////////////// Compute distance of a fiber to a bundle ////////////////////
////////////////////////////////////////////////////////////////////////////////
float computeDistanceFiberToBundle( BundlesDataFormat& atlasBundleData,
                                    const std::vector<float>& fiber )
{

  int nbPointsAtlasFibers = atlasBundleData.pointsPerTrack[ 0 ] ;
  int nbPointsFiber = ( int )( fiber.size() / 3 ) ;
  if ( nbPointsFiber != nbPointsAtlasFibers )
  {

    std::cout << "ERROR in computeNeighborhood::computeDistanceFiberToBundle : "
              << "the number of points in the fiber and in tha bundle fibers "
              << "must be the same "
              << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPointsFiber ;

  int nbFibersAtlasBundle = atlasBundleData.curves_count ;

  // float distance = 0 ;
  std::vector<float> distances( nbFibersAtlasBundle, 0 ) ;

  #pragma omp parallel for
  for ( int atlasFiberIndex = 0 ; atlasFiberIndex < nbFibersAtlasBundle ;
                                                             atlasFiberIndex++ )
  {

    // distance += atlasBundleData.computeMDADBetweenTwoFibers(
    //                                                atlasBundleData.matrixTracks,
    //                                                fiber,
    //                                                atlasFiberIndex,
    //                                                0,
    //                                                nbPoints ) ;
    distances[ atlasFiberIndex ] =
                             atlasBundleData.computeMDFBetweenTwoFibers(
                                                   atlasBundleData.matrixTracks,
                                                   fiber,
                                                   atlasFiberIndex,
                                                   0,
                                                   nbPoints ) ;

  }

  int minIndex = std::min_element( distances.begin(), distances.end() ) -
                                                             distances.begin() ;

  // distance /= nbFibersAtlasBundle ;
  float distance = distances[ minIndex ] ;
  return distance ;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

  const auto start_time = std::chrono::system_clock::now() ;

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  int index_moving_bundle, index_reference_bundle, index_output_bundle,
                          index_nbFibers, index_transform, index_transform_type,
                                                    index_verbose, index_help ;

  index_moving_bundle = getFlagPosition( argc, argv, "-m") ;
  index_reference_bundle = getFlagPosition( argc, argv, "-r") ;
  index_output_bundle = getFlagPosition( argc, argv, "-o") ;
  index_nbFibers = getFlagPosition( argc, argv, "-nbFibers") ;
  index_transform = getFlagPosition( argc, argv, "-omat") ;
  index_transform_type = getFlagPosition( argc, argv, "-xfmType") ;
  index_verbose = getFlagPosition( argc, argv, "-v") ;
  index_help = getFlagPosition( argc, argv, "-h") ;

  if ( index_help > 0 )
  {

    std::cout << "Function to register bundles using RecoBundles method : \n"
              << "-m : Moving bundle \n"
              << "-r : Reference bundle \n"
              << "-o : Output registered bundle \n"
              << "[-nbFibers] : Select nbFibers in moving and reference "
              << "bundles to do the optimization. If argument not given then "
              << "nbFibers = min ( nbFibresMoving, nbFibersReference ) * 1.5 )."
              << " If -nbFibers = 0 then nbFibers = min ( nbFibresMoving, "
              << "nbFibersReference ) \n"
              << "[-omat] : Output affine transform (.mat) \n"
              << "[-xfmType] : choose between 'affine' and 'rigid' ( default "
              << "= rigid ) \n"
              << "[-v] : Set verbosity level at 1 \n"
              << "[-h] : Show this message " << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_moving_bundle )
  {

    std::cout << "-m argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_reference_bundle )
  {

    std::cout << "-r argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( !index_output_bundle )
  {

    std::cout << "-o argument required ..." << std::endl ;
    exit( 1 ) ;

  }

  if ( index_transform_type )
  {
    if ( argv[ index_transform_type + 1 ] )
    {

      transformType = argv[ index_transform_type + 1 ] ;

    }
    else
    {

      transformType = "rigid" ;

    }

    if ( transformType != "rigid" && transformType!= "affine" )
    {

      std::cout << "ERROR input argument : -xfmType argument must be 'affine' "
                << "or 'rigid' " << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( index_nbFibers )
  {
    if ( argv[ index_nbFibers + 1 ] )
    {

      maxNbFibers = std::stoi( argv[ index_nbFibers + 1 ] ) ;
      useDefaultMaxNbFibers = false ;

    }
    else
    {

      std::cout << "Error : -nbFibers argument must be an integer\n" ;
      exit( 1 ) ;
    }

  }
  else
  {

    useDefaultMaxNbFibers = true ;

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

  /////////////////////////////////////////////////////////////////////////////

  std::string movingFilename( argv[ index_moving_bundle + 1 ] ) ;
  char lastChar = movingFilename[ movingFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    movingFilename = movingFilename.substr( 0, movingFilename.size() - 1 ) ;

  }

  std::string referenceFilename( argv[ index_reference_bundle + 1 ] ) ;
  lastChar = referenceFilename[ referenceFilename.size()  - 1 ] ;
  if ( lastChar == '/' )
  {

    referenceFilename = referenceFilename.substr( 0,
                                                referenceFilename.size() - 1 ) ;

  }

  std::string outputFilename( argv[ index_output_bundle + 1 ] ) ;
  lastChar = outputFilename[ outputFilename.size() - 1 ] ;
  if ( lastChar == '/' )
  {

    outputFilename = outputFilename.substr( 0, outputFilename.size() - 1 ) ;

  }

  std::string transformFilename ;
  if ( index_transform )
  {

    transformFilename = argv[ index_transform + 1 ] ;
    lastChar = transformFilename[ transformFilename.size() - 1 ] ;
    if ( lastChar == '/' )
    {

      transformFilename = transformFilename.substr( 0,
                                                transformFilename.size() - 1 ) ;

    }

  }


  //   xxxxxxxxxxxxxxxxxxxxxxx Reading moving bundle xxxxxxxxxxxxxxxxxxxxxx   //

  std::string movingBundlesFilename ;
  std::string movingBundlesDataFilename ;
  bool isBundlesFormat = false ;

  if ( movingFilename.find( ".bundlesdata" ) != std::string::npos )
  {

    movingBundlesDataFilename = movingFilename ;

    movingBundlesFilename = movingFilename ;
    size_t index = movingFilename.find( ".bundlesdata" ) ;
    movingBundlesFilename.replace( index, 12, ".bundles") ;

    isBundlesFormat = true ;

  }
  else if ( movingFilename.find( ".bundles" ) != std::string::npos )
  {

    movingBundlesFilename = movingFilename ;

    movingBundlesDataFilename = movingFilename ;
    size_t index = movingFilename.find( ".bundles" ) ;
    movingBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

    isBundlesFormat = true ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesDataFormat movingBundlesData ;
  BundlesFormat movingBundleInfo ;

  if ( isBundlesFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << movingBundlesFilename << std::endl ;

    }

    movingBundleInfo.bundlesReading( movingBundlesFilename.c_str(),
                                        verbose ) ;

    if ( verbose )
    {

      std::cout << "Reading : " << movingBundlesDataFilename << std::endl ;

    }
    movingBundlesData.bundlesdataReading( movingBundlesDataFilename.c_str(),
                                        movingBundlesFilename.c_str(),
                                        verbose ) ;

  }

  //   xxxxxxxxxxxxxxxxxxxxx Reading reference bundle xxxxxxxxxxxxxxxxxxxxx   //

  std::string referenceBundlesFilename ;
  std::string referenceBundlesDataFilename ;

  if ( referenceFilename.find( ".bundlesdata" ) != std::string::npos )
  {

    referenceBundlesDataFilename = referenceFilename ;

    referenceBundlesFilename = referenceFilename ;
    size_t index = referenceFilename.find( ".bundlesdata" ) ;
    referenceBundlesFilename.replace( index, 12, ".bundles") ;

    isBundlesFormat = true ;

  }
  else if ( referenceFilename.find( ".bundles" ) != std::string::npos )
  {

    referenceBundlesFilename = referenceFilename ;

    referenceBundlesDataFilename = referenceFilename ;
    size_t index = referenceFilename.find( ".bundles" ) ;
    referenceBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

    isBundlesFormat = true ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }


  BundlesDataFormat referenceBundlesData ;
  BundlesFormat referenceBundleInfo ;

  if ( isBundlesFormat )
  {

    if ( verbose )
    {

      std::cout << "Reading : " << referenceBundlesFilename << std::endl ;

    }

    referenceBundleInfo.bundlesReading( referenceBundlesFilename.c_str(),
                                        verbose ) ;

    if ( verbose )
    {

      std::cout << "Reading : " << referenceBundlesDataFilename << std::endl ;

    }
    referenceBundlesData.bundlesdataReading(
                                        referenceBundlesDataFilename.c_str(),
                                        referenceBundlesFilename.c_str(),
                                        verbose ) ;

  }

  // Sanity checks //
  int nbPointsMoving = movingBundlesData.pointsPerTrack[ 0 ] ;
  int nbCurvesMoving = movingBundlesData.curves_count ;
  for( int fiber = 1 ; fiber < nbCurvesMoving ; fiber++ )
  {

    if ( nbPointsMoving != movingBundlesData.pointsPerTrack[ fiber ] )
    {

      std::cout << "ERROR : All the fibers in the moving bundle must have the "
                << "same number of points, got " << nbPointsMoving
                << " for point 0 and "
                << movingBundlesData.pointsPerTrack[ fiber ]
                << " for point " << fiber << std::endl ;
      exit( 1 ) ;

    }

  }

  int nbPointsReference  = referenceBundlesData.pointsPerTrack[ 0 ] ;
  int nbCurvesReference = referenceBundlesData.curves_count ;
  for( int fiber = 1 ; fiber < nbCurvesReference ; fiber++ )
  {

    if ( nbPointsReference != referenceBundlesData.pointsPerTrack[ fiber ] )
    {

      std::cout << "ERROR : All the fibers in the reference bundle must have "
                << "the same number of points, got " << nbPointsReference
                << " for point 0 and "
                << referenceBundlesData.pointsPerTrack[ fiber ]
                << " for point " << fiber << std::endl ;
      exit( 1 ) ;

    }

  }

  if ( nbPointsMoving != nbPointsReference )
  {

    std::cout << "ERROR : All the fibers in both bundles must have "
              << "the same number of points, got " << nbPointsMoving
              << " for moving fiber and " << nbPointsReference
              << " for reference fiber" << std::endl ;
    exit( 1 ) ;

  }

  int nbPoints = nbPointsReference ;

  if ( verbose )
  {

    std::cout << "Sanity cheks : OK " << std::endl ;

  }


  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  // Getting maxNbFibers


  if ( useDefaultMaxNbFibers )
  {

    if ( nbCurvesMoving < nbCurvesReference )
    {

      maxNbFibers = nbCurvesMoving * 1.5 ;

    }
    else
    {

      maxNbFibers = nbCurvesReference * 1.5 ;

    }

  }
  else if ( !useDefaultMaxNbFibers && maxNbFibers == 0 )
  {

    if ( nbCurvesMoving < nbCurvesReference )
    {

      maxNbFibers = nbCurvesMoving  ;

    }
    else
    {

      maxNbFibers = nbCurvesReference ;

    }

  }




  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  // Selecting maxNbFibers fibers in reference bundle that resemble the most to
  // the reference bundle ( most representative fibers since they are the ones
  // with the smallest average distance to all fibers in the bundle )

  if ( maxNbFibers < nbCurvesReference )
  {

    std::vector<float> selectedReferenceFibers( 3 * nbPoints * maxNbFibers, 0 ) ;
    std::vector<int> pointsPerTrackSelectedFibers( maxNbFibers, 0 ) ;

    if ( verbose )
    {

      std::cout << "Selecting " << maxNbFibers << " fibers of the reference "
                << "bundle which resemble the most to the reference bundle..." ;

    }

    std::vector<float> distancesReferenceFibersToReference( nbCurvesReference,
                                                                           0 ) ;
    #pragma omp parallel for
    for ( int fiberIndex = 0 ; fiberIndex < nbCurvesReference ; fiberIndex++ )
    {

      std::vector<float> fiber( 3 * nbPoints, 0 ) ;
      referenceBundlesData.getFiberFromTractogram(
                                              referenceBundlesData.matrixTracks,
                                              fiberIndex,
                                              nbPoints,
                                              fiber ) ;
      distancesReferenceFibersToReference[ fiberIndex ] =
                   computeDistanceFiberToBundle( referenceBundlesData, fiber ) ;

    }

    std::vector<int> indices( nbCurvesReference, 0 ) ;
    std::iota( indices.begin(), indices.end(), 0 ) ;
    std::partial_sort( indices.begin(), indices.begin() + maxNbFibers,
                       indices.end(),
                       [&distancesReferenceFibersToReference]( int i, int j )
                       { return distancesReferenceFibersToReference[ i ] <
                                distancesReferenceFibersToReference[ j ] ; } ) ;


    for ( int fiberSelectedIndex = 0 ; fiberSelectedIndex < maxNbFibers ;
                                                          fiberSelectedIndex++ )
    {

      int offsetExtractedNeighborhood = 3 * nbPoints * fiberSelectedIndex ;

      // Find minimum element index
      int minIndex = indices[ fiberSelectedIndex ] ;

      std::vector<float> fiber( 3 * nbPoints, 0 ) ;
      // referenceBundlesData.getFiberFromTractogram(
      //                                         referenceBundlesData.matrixTracks,
      //                                         fiberSelectedIndex,
      //                                         nbPoints,
      //                                         fiber ) ;
      referenceBundlesData.getFiberFromTractogram(
                                              referenceBundlesData.matrixTracks,
                                              minIndex,
                                              nbPoints,
                                              fiber ) ;

      std::copy( fiber.begin(),
                 fiber.begin() + 3 * nbPoints,
                 selectedReferenceFibers.begin() + offsetExtractedNeighborhood ) ;

      pointsPerTrackSelectedFibers[ fiberSelectedIndex ] =
                               referenceBundlesData.pointsPerTrack[ minIndex ] ;

    }

    referenceBundlesData.matrixTracks = std::vector<float>() ; // deallocation
    referenceBundlesData.matrixTracks = selectedReferenceFibers ;

    referenceBundlesData.pointsPerTrack = std::vector<int>() ; // deallocation
    referenceBundlesData.pointsPerTrack = pointsPerTrackSelectedFibers ;

    referenceBundlesData.curves_count = maxNbFibers ;


    if ( verbose )
    {

      std::cout << " Done " << std::endl ;

    }

  }


  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  // Selecting maxNbFibers fibers in moving bundle that resemble the most to
  // the reference bundle
  BundlesDataFormat movingBundlesDataCopy( movingBundlesData ) ;
  if ( maxNbFibers < nbCurvesMoving )
  {

    std::vector<float> selectedMovingFibers( 3 * nbPoints * maxNbFibers, 0 ) ;
    std::vector<int> pointsPerTrackSelectedFibers( maxNbFibers, 0 ) ;

    if ( verbose )
    {

      std::cout << "Selecting " << maxNbFibers << " fibers of the moving bundle "
                << "which resemble the most to the reference bundle..." ;

    }

    std::vector<float> distancesMovingFibersToReference( nbCurvesMoving, 0 ) ;
    #pragma omp parallel for
    for ( int fiberIndex = 0 ; fiberIndex < nbCurvesMoving ; fiberIndex++ )
    {

      std::vector<float> fiber( 3 * nbPoints, 0 ) ;
      movingBundlesData.getFiberFromTractogram(
                                              movingBundlesData.matrixTracks,
                                              fiberIndex,
                                              nbPoints,
                                              fiber ) ;
      distancesMovingFibersToReference[ fiberIndex ] =
                   computeDistanceFiberToBundle( referenceBundlesData, fiber ) ;

    }

    std::vector<int> indices( nbCurvesMoving, 0 ) ;
    std::iota( indices.begin(), indices.end(), 0 ) ;
    std::partial_sort( indices.begin(), indices.begin() + maxNbFibers,
                       indices.end(),
                       [&distancesMovingFibersToReference]( int i, int j )
                       { return distancesMovingFibersToReference[ i ] <
                                distancesMovingFibersToReference[ j ] ; } ) ;


    for ( int fiberSelectedIndex = 0 ; fiberSelectedIndex < maxNbFibers ;
                                                          fiberSelectedIndex++ )
    {

      int offsetExtractedNeighborhood = 3 * nbPoints * fiberSelectedIndex ;

      // Find minimum element index
      int minIndex = indices[ fiberSelectedIndex ] ;

      std::vector<float> fiber( 3 * nbPoints, 0 ) ;
      // movingBundlesData.getFiberFromTractogram(
      //                                         movingBundlesData.matrixTracks,
      //                                         fiberSelectedIndex,
      //                                         nbPoints,
      //                                         fiber ) ;
      movingBundlesData.getFiberFromTractogram(
                                              movingBundlesData.matrixTracks,
                                              minIndex,
                                              nbPoints,
                                              fiber ) ;

      std::copy( fiber.begin(),
                 fiber.begin() + 3 * nbPoints,
                 selectedMovingFibers.begin() + offsetExtractedNeighborhood ) ;

      pointsPerTrackSelectedFibers[ fiberSelectedIndex ] =
                                  movingBundlesData.pointsPerTrack[ minIndex ] ;

    }

    movingBundlesData.matrixTracks = std::vector<float>() ; // deallocation
    movingBundlesData.matrixTracks = selectedMovingFibers ;

    movingBundlesData.pointsPerTrack = std::vector<int>() ; // deallocation
    movingBundlesData.pointsPerTrack = pointsPerTrackSelectedFibers ;

    movingBundlesData.curves_count = maxNbFibers ;


    if ( verbose )
    {

      std::cout << " Done " << std::endl ;

    }


  }

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
  if ( verbose )
  {

    std::cout << "Using " << transformType << " transform " << std::endl ;

  }
  int nbTransformCoefficients = 0 ;
  if ( transformType == "rigid" )
  {

    nbTransformCoefficients = 6 ;

  }
  else if ( transformType == "affine" )
  {

    nbTransformCoefficients = 12 ;

  }
  else
  {

    std::cout << "ERROR in registerBundles : invalid transform type \n" ;
    exit( 1 ) ;

  }


  // Creating function to minimize
  BMDDistance fun( referenceBundlesData, movingBundlesData ) ;
  Eigen::VectorXd transformInitial = Eigen::VectorXd::Zero(
                                                     nbTransformCoefficients ) ;
  if ( transformType == "affine" )
  {

    transformInitial[ 0 ] = 1 ;
    transformInitial[ 5 ] = 1 ;
    transformInitial[ 10 ] = 1 ;

  }

  Eigen::VectorXd testGradient = Eigen::VectorXd::Zero(
                                                     nbTransformCoefficients ) ;
  double testDistance = fun( transformInitial, testGradient ) ;
  float normGradient = fun.normGradient( testGradient ) ;

  // Initial distance
  if ( verbose )
  {

    float initialBMD = fun.computeBMD( referenceBundlesData,
                                                           movingBundlesData ) ;
    std::cout << "Initial BMD distance : " << initialBMD << std::endl ;

  }

  bool epsilonFound = false ;
  float epsilon = normGradient * 1e-2 ;
  Eigen::VectorXd transformCoefficients = Eigen::VectorXd::Zero(
                                                     nbTransformCoefficients ) ;
  int nbIterations = 0 ;
  float minimumBmdDistance = 0 ;
  int nbIterationsToFoundEpsilon = 0 ;

  if ( verbose )
  {

    std::cout << "Initial value epsilon : " << epsilon << std::endl ;

  }


  // Optimizing epsilon
  if ( verbose )
  {

    std::cout << "Optimizing value of epsilon... \n" ;

  }

  epsilonFound = false ;
  nbIterationsToFoundEpsilon = 0 ;
  int nbConsecutiveEqualNbOfIterationsOptimizer = 0 ;
  int nbConsecutiveFails = 0 ;
  nbIterations = 0 ;
  while ( !epsilonFound && nbIterationsToFoundEpsilon <
                                                      nbIterationSearchEpsilon )
  {

    if ( verbose )
    {

      printf( " Iteration [ %d / %d ] \r", nbIterationsToFoundEpsilon + 1,
                                                    nbIterationSearchEpsilon ) ;
      fflush( stdout ) ;

    }

    try
    {

      // Set up parameters for Optimization
      LBFGSpp::LBFGSParam<double> param ;
      param.epsilon = epsilon ;
      param.max_iterations = 100 ;
      if ( verbose > 1 )
      {

        std::cout << "Parameters of LBFGS :  { epsilon : " << param.epsilon
                  << ", max_iterations : " << param.max_iterations << " } \n" ;

      }

      // Create solver and function object
      if ( verbose > 1 )
      {

        std::cout << "Creating solver and function object ... " ;

      }
      LBFGSpp::LBFGSSolver<double> solver( param ) ;
      if ( verbose > 1 )
      {

        std::cout << "OK" << std::endl ;

      }


      // Initial guess = identity
      Eigen::VectorXd transformCoefficientsSolution = Eigen::VectorXd::Zero(
                                                     nbTransformCoefficients ) ;

      if ( transformType == "affine" )
      {

        transformCoefficientsSolution[ 0 ] = 1 ;
        transformCoefficientsSolution[ 5 ] = 1 ;
        transformCoefficientsSolution[ 10 ] = 1 ;

      }

      if ( verbose > 1 )
      {

        std::cout << "Initial guess : "
                  << transformCoefficientsSolution.transpose() << std::endl ;

      }

      // transformCoefficientsSolution will be overwritten to be the best one
      if ( verbose > 1 )
      {

        std::cout << "Optimization... "  << std::endl ;

      }

      double bmdDistance = 0 ;

      int niter = solver.minimize( fun, transformCoefficientsSolution,
                                                                 bmdDistance ) ;


      if ( verbose > 1 )
      {

        std::cout << "OK" << std::endl ;

      }


      // Checking stoping criterion

      if ( nbIterations == niter )
      {

        nbConsecutiveEqualNbOfIterationsOptimizer += 1 ;

      }
      else
      {

        nbConsecutiveEqualNbOfIterationsOptimizer = 0 ;

      }

      // if ( nbConsecutiveEqualNbOfIterationsOptimizer == 2 && niter > 10)
      if ( niter > 10 || nbConsecutiveEqualNbOfIterationsOptimizer > 2 )
      {

        epsilonFound = true ;

      }

      // std::cout << "\nNumber iterations optimizer : " << niter ;
      // std::cout << "" << std::flush ;


      // Overwriting results
      epsilon /= 2 ;
      nbIterationsToFoundEpsilon += 1 ;
      transformCoefficients = transformCoefficientsSolution ;
      nbIterations = niter ;
      minimumBmdDistance = bmdDistance ;

      nbConsecutiveFails = 0 ;

    }
    // catch ( const std::logic_error )
    catch ( ... )
    {

      epsilon = ( 2 * epsilon + epsilon ) / 2 ;
      nbIterationsToFoundEpsilon += 1 ;
      nbConsecutiveFails += 1 ;
      if ( nbConsecutiveFails == 3 && nbIterations > 0 )
      {

        epsilonFound = true ;

      }

    }

  }

  if ( verbose > 1 )
  {

    std::cout << "\n epsilon = " << epsilon << " found in "
              << nbIterationsToFoundEpsilon << " iterations \n" ;

  }
  else if ( verbose )
  {

    std::cout << "\n epsilon = " << epsilon << std::endl ;

  }


  std::cout << "Number of iterations LBFGS : " << nbIterations << std::endl ;
  std::cout << "Affinte Transform : \n" << transformCoefficients.transpose()
                                                                  << std::endl ;
  std::cout << "Minimum BMD distance : " << minimumBmdDistance << std::endl ;

  BundlesDataFormat movedBundlesData( movingBundlesDataCopy ) ;
  if ( transformType == "rigid" )
  {

    fun.applyRigidToBundle( movingBundlesDataCopy.matrixTracks,
                             transformCoefficients,
                             movingBundlesDataCopy.curves_count,
                             nbPoints,
                             movedBundlesData.matrixTracks ) ;

  }
  else if ( transformType == "affine" )
  {

    fun.applyAffineToBundle( movingBundlesDataCopy.matrixTracks,
                             transformCoefficients,
                             movingBundlesDataCopy.curves_count,
                             nbPoints,
                             movedBundlesData.matrixTracks ) ;

  }
  else
  {

    std::cout << "ERROR in registerBundles : Invalid transform type \n" ;
    exit( 1 ) ;

  }


  BundlesFormat movedBundleInfo( movingBundleInfo ) ;

  // xxxxxxxxxxxxxxxxxxxxxxxxxxxx Saving results xxxxxxxxxxxxxxxxxxxxxxxxxxxx //

  // Saving .bundles and .bundlesdata
  std::string outputBundlesFilename ;
  std::string outputBundlesDataFilename ;

  if ( outputFilename.find( ".bundlesdata" ) != std::string::npos )
  {

    outputBundlesDataFilename = outputFilename ;

    outputBundlesFilename = outputFilename ;
    size_t index = outputFilename.find( ".bundlesdata" ) ;
    outputBundlesFilename.replace( index, 12, ".bundles") ;

  }
  else if ( outputFilename.find( ".bundles" ) != std::string::npos )
  {

    outputBundlesFilename = outputFilename ;

    outputBundlesDataFilename = outputFilename ;
    size_t index = outputFilename.find( ".bundles" ) ;
    outputBundlesDataFilename.replace( index, 8, ".bundlesdata") ;

  }
  else
  {

    std::cout << "The only tractogram format supported is .bundles"
              << std::endl ;
    exit( 1 ) ;

  }

  movedBundlesData.bundlesdataWriting( outputBundlesDataFilename.c_str(),
                                                                     verbose ) ;
  movedBundleInfo.bundlesWriting( outputBundlesFilename.c_str(), verbose ) ;

  // Saving transform
  if ( index_transform )
  {

    std::fstream file ;
    file.open( transformFilename, std::ios::out ) ;

    if( file.is_open() )
    {
      for ( int i = 0 ; i < nbTransformCoefficients ; i++ )
      {

        file << transformCoefficients[ i ] << "\n" ;
        // file << "\n" ;

      }

      file.close() ;

    }
    else
    {

      std::cout << "Problem saving transform in : " << transformFilename
                                                    << std::endl ;

    }

  }

  const std::chrono::duration< double > duration =
                              std::chrono::system_clock::now() - start_time ;

  if ( verbose )
  {

    std::cout << "Duration : " << duration.count() << std::endl ;

  }

  return( 0 ) ;

}
