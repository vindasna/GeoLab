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

// Optimization
// #include "/home/nv264568/Bureau/Doctorat/Codes/ProjectAtlasFiles/eigen-3.4.0/Eigen/Dense"
// #include "/home/nv264568/Bureau/Doctorat/Codes/ProjectAtlasFiles/yixuan-LBFGSpp-624755d/include/LBFGS.h"
// #include "/home/nv264568/Bureau/Doctorat/Codes/ProjectAtlasFiles/eigen-3.4.0/Eigen/Core"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <LBFGS.h>

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class BundlesDataFormat
{

  public :
  //////////////////////////////// Public Fields ///////////////////////////////
  std::vector<float> matrixTracks ;
  std::vector<int32_t> pointsPerTrack ;
  std::vector<int64_t> fibersWithNans ;
  int curves_count = 0 ;

  ///////////////////////////////// Operators //////////////////////////////////
  float operator[]( int64_t index ) ;

  //////////////////////////////// Constructors ////////////////////////////////
  BundlesDataFormat() ;

  BundlesDataFormat( const char* bundlesdataFilename,
                     const char* bundlesFilename,
                     int verbose ) ;

  BundlesDataFormat( std::vector<float>& matrixTracks,
                     std::vector<int32_t>& pointsPerTrack,
                     int curves_count ) ;

  BundlesDataFormat( const BundlesDataFormat& bundlesData ) ;

  //////////////////////////////// Destructors /////////////////////////////////
  virtual ~BundlesDataFormat() ;

  ////////////////////////////////// Mehtods ///////////////////////////////////
  void bundlesdataReading( const char* bundlesdataFilename,
                           const char* bundlesFilename,
                           int verbose ) ;

  void bundlesdataWriting( const char* bundlesdataFilename,
                           int verbose ) ;

  // void toTRK( const BundlesFormat& bundlesInfo, TrkFormat& trkData ) ;

  float computeLengthFiber( const std::vector<float>& tractogramFibers,
                            int fiberIndex,
                            int nbPoints ) ;
  float computeLengthFiber( const std::vector<float>& fiber ) ;
  float computeLengthFiber( int fiberIndex ) ;


  void resampleFiberEquidistant( const std::vector<float>& inputFiber,
                                 std::vector<float>& outputFiber,
                                 int nbPointsToResample ) ;


  void computeMedialPointFiberTractogram(
                              const std::vector<float>& tractogramFibers,
                              int fiberIndex,
                              int nbPoints,
                              std::vector<float>& medialPointFiberTractogram ) ;
  void computeMedialPointFiberTractogram(
                              int fiberIndex,
                              std::vector<float>& medialPointFiberTractogram ) ;
  void computeMedialPointFiberTractogram(
                              const std::vector<float>& fiber,
                              std::vector<float>& medialPointFiberTractogram ) ;

  int computeMedialPointFiberWithDistance(
                              const std::vector<float>& tractogramFibers,
                              int fiberIndex,
                              int nbPoints,
                              std::vector<float>& medialPointFiberTractogram ) ;
  int computeMedialPointFiberWithDistance(
                              int fiberIndex,
                              std::vector<float>& medialPointFiberTractogram ) ;
  int computeMedialPointFiberWithDistance(
                              const std::vector<float>& fiber,
                              std::vector<float>& medialPointFiberTractogram ) ;

  void computeGravityCenterFiber(
                              const std::vector<float>& tractogramFibers,
                              int fiberIndex,
                              int nbPoints,
                              std::vector<float>& gravityCenterFiber ) ;
  void computeGravityCenterFiber(
                              int fiberIndex,
                              std::vector<float>& gravityCenterFiber ) ;
  void computeGravityCenterFiber(
                              const std::vector<float>& fiber,
                              std::vector<float>& gravityCenterFiber ) ;


  void computeNormalVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& normalVector ) ;
  void computeNormalVectorFiberTractogram( int fiberIndex,
                                           std::vector<float>& normalVector ) ;
  void computeNormalVectorFiberTractogram( const std::vector<float>& fiber,
                                           std::vector<float>& normalVector ) ;

  void computeDirectionVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           const std::vector<float>& normalVector,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& directionVector ) ;
  void computeDirectionVectorFiberTractogram(
                                         int fiberIndex,
                                         const std::vector<float>& normalVector,
                                         std::vector<float>& directionVector ) ;
  void computeDirectionVectorFiberTractogram(
                                         const std::vector<float>& fiber,
                                         const std::vector<float>& normalVector,
                                         std::vector<float>& directionVector ) ;

  float computeMDADBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) ;
  float computeMDADBetweenTwoFibers(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int fiberIndex_1,
                                   int fiberIndex_2,
                                   int nbPoints ) ;

  float computeMDFBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) ;
  float computeMDFBetweenTwoFibers(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int fiberIndex_1,
                                   int fiberIndex_2,
                                   int nbPoints ) ;

  float compareDisimilarityBundles(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int nbFibersTract_1,
                                   int nbFibersTract_2,
                                   int nbPoints1,
                                   int nbPoints2 ) ;
  float compareDisimilarityBundles( const std::vector<float>& tractogramFibers,
                                    int nbFibersTract,
                                    int nbPoints ) ;

  double distanceBetweenBundles( const std::vector<float>& bundle1,
                                 const std::vector<float>& bundle2,
                                 int nbFibersBundle1,
                                 int nbFibersBundle2,
                                 int nbPointsTract_1,
                                 int nbPointsTract_2 ) ;
  double distanceBetweenBundles( const std::vector<float>& bundle,
                                 int nbFibersBundle,
                                 int nbPoints ) ;

  void computeNumberAdjacentFibersBundle1toBundle2(
                                    BundlesDataFormat& bundle1,
                                    BundlesDataFormat& bundle2,
                                    float threshold,
                                    std::vector<int>& nbAdjacentFibersBundle ) ;

  void computeNumberAdjacentFibersBundle1toBundle2(
                                    const std::vector<float>& bundle1,
                                    const std::vector<float>& bundle2,
                                    int nbFibersBundle1,
                                    int nbFibersBundle2,
                                    int nbPoints,
                                    float threshold,
                                    std::vector<int>& nbAdjacentFibersBundle ) ;

  float coverageBundle1toBundle2( BundlesDataFormat& bundle1,
                                  BundlesDataFormat& bundle2,
                                  float threshold,
                                  int verbose ) ;

  float overlapBundle1toBundle2( BundlesDataFormat& bundle1,
                                 BundlesDataFormat& bundle2,
                                 float threshold,
                                 int verbose ) ;


  float bundlesAdjacency( BundlesDataFormat& bundle1,
                          BundlesDataFormat& bundle2,
                          float threshold,
                          int verbose ) ;

  float computeAngleBetweenVectors( const std::vector<float>& vector1,
                                    const std::vector<float>& vector2 ) ;

  float computeAngleBetweenPlanes( const std::vector<float>& vector1,
                                   const std::vector<float>& vector2 ) ;

  float computeAngleBetweenDirections( const std::vector<float>& vector1,
                                       const std::vector<float>& vector2 ) ;


  void projectVectorToPlane( const std::vector<float>& normalToPlane,
                             const std::vector<float>& inputVector,
                             std::vector<float>& projectedVector ) ;

  void computeRotationMatrixFromVectorAndAngle(
                                          const std::vector<float>& vector,
                                          float angle, // in radians
                                          std::vector<float>& rotationMatrix ) ;

  void applyRotationMatrixToVector( const std::vector<float>& vector,
                                    const std::vector<float>& rotationMatrix,
                                    std::vector<float>& rotatedVector ) ;

  void applyRotationMatrixToFiber( const std::vector<float>& fiber,
                                   const std::vector<float>& rotationMatrix,
                                   const std::vector<float>& medialPointFiber,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) ;
  void applyRotationMatrixToFiber( BundlesDataFormat& inputBundlesData,
                                   const std::vector<float>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) ;
  void applyRotationMatrixToFiber( const std::vector<float>& tractogramFibers,
                                   const std::vector<float>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) ;
  void applyRotationMatrixToFiber( const std::vector<float>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) ;
  void applyRotationMatrixToFiber( const std::vector<float>& fiber,
                                   const std::vector<float>& rotationMatrix,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) ;

  void putFibersInSamePlane( const std::vector<float>& normalVector1,
                             const std::vector<float>& normalVector2,
                             const std::vector<float>& tractogramFibers2,
                             int fiberIndex2,
                             int nbPoints,
                             std::vector<float>& fiber2ToPlane1,
                             std::vector<float>& newNormalVectorFiber2 ) ;
  void putFibersInSamePlane( const std::vector<float>& normalVector1,
                             const std::vector<float>& normalVector2,
                             const std::vector<float>& fiber2,
                             int nbPoints,
                             std::vector<float>& fiber2ToPlane1,
                             std::vector<float>& newNormalVectorFiber2 ) ;

  void putFiberInPlaneXY( const std::vector<float>& normalVector,
                          const std::vector<float>& tractogramFibers,
                          int fiberIndex,
                          int nbPoints,
                          std::vector<float>& fiberToPlaneXY,
                          std::vector<float>& newNormalVectorFiber ) ;
  void putFiberInPlaneXY( const std::vector<float>& normalVector,
                          const std::vector<float>& fiber,
                          int nbPoints,
                          std::vector<float>& fiberToPlaneXY,
                          std::vector<float>& newNormalVectorFiber ) ;

  void putFibersInSameDirection( const std::vector<float>& normalVector1,
                                 const std::vector<float>& normalVector2,
                                 const std::vector<float>& directionVector1,
                                 const std::vector<float>& fiber2,
                                 int nbPoints,
                                 std::vector<float>& fiber2ToDirection1 ) ;

  void registerFiber( const std::vector<float>& tractogramFibers2,
                      const std::vector<float>& normalVector1,
                      const std::vector<float>& normalVector2,
                      const std::vector<float>& directionVector1,
                      const std::vector<float>& medialPointFiber1,
                      const std::vector<float>& medialPointFiber2,
                      int fiberIndexTractogram2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2,
                      int verbose ) ;
  void registerFiber( const std::vector<float>& tractogramFibers1,
                      const std::vector<float>& tractogramFibers2,
                      int fiberIndexTractogram1,
                      int fiberIndexTractogram2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2,
                      int verbose ) ;
  void registerFiber( const std::vector<float>& fiber2,
                      const std::vector<float>& normalVector1,
                      const std::vector<float>& normalVector2,
                      const std::vector<float>& directionVector1,
                      const std::vector<float>& medialPointFiber1,
                      const std::vector<float>& medialPointFiber2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2,
                      int verbose ) ;
  void registerFiber( const std::vector<float>& fiber1,
                      const std::vector<float>& fiber2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2,
                      int verbose ) ;


  void registerFiberToPlaneXYAndDirectionX(
                                const std::vector<float>& tractogramFibers,
                                const std::vector<float>& normalVector,
                                const std::vector<float>& medialPointFiber,
                                int fiberIndexTractogram,
                                int nbPoints,
                                std::vector<float>& fiberToPlaneXYAndDirectionX,
                                std::vector<float>& newNormalVectorFiber,
                                int verbose ) ;
  void registerFiberToPlaneXYAndDirectionX(
                                const std::vector<float>& fiber,
                                const std::vector<float>& normalVector,
                                const std::vector<float>& medialPointFiber,
                                int nbPoints,
                                std::vector<float>& fiberToPlaneXYAndDirectionX,
                                std::vector<float>& newNormalVectorFiber,
                                int verbose ) ;

  void getFiberFromTractogram( const std::vector<float>& tractogramFibers,
                               int fiberIndex,
                               int nbPoints,
                               std::vector<float>& fiber ) ;

  void getFiberWithVectors( const std::vector<float>& fiber,
                            const std::vector<float>& referenceFiber,
                            int nbPoints,
                            std::vector<float>& fiberWithVectors ) ;

  void getFiberWithVectors( const std::vector<float>& fiber,
                            int nbPoints,
                            std::vector<float>& fiberWithVectors ) ;



  ///////// Optimization functions /////////
  void getMatrixAndVectorForLeastSquares(
                               const std::vector<float>& tractogramFibers,
                               int fiberIndex,
                               int nbPoints,
                               Eigen::MatrixXf& AMatrix,
                               Eigen::VectorXf& BVector ) ;


  ///////// Algebra functions /////////

  // Point1 is the point where the vector is pointing to
  void vectorFromPoints( const std::vector<float>& point1,
                         const std::vector<float>& point2,
                         std::vector<float>& vector ) ;

  float scalarProduct( const std::vector<float>& vector1,
                       const std::vector<float>& vector2 ) ;

  float normVector( const std::vector<float>& vector ) ;

  void normalizeVector( std::vector<float>& vector ) ;

  void translatePoint( const std::vector<float>& point,
                       const std::vector<float>& unitaryTranslationVector,
                       float translationDistance,
                       std::vector<float>& translatedPoint ) ;
  void translatePoint( const std::vector<float>& point,
                       const std::vector<float>& translationVector,
                       std::vector<float>& translatedPoint ) ;

  void crossProduct( const std::vector<float>& vector1,
                     const std::vector<float>& vector2,
                     std::vector<float>& result ) ;

  int checkIfFiberPointCanBeResampled(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int point,
                                     int nbPoints ) ;
  void resampleFiberWithNan( std::vector<float>& tractogramFibers,
                             int fiberIndex,
                             int nbPoints,
                             int verbose ) ;


} ;
