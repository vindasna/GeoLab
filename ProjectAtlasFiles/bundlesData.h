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

#include <Eigen/Core>
#include <Eigen/Dense>
#include <LBFGS.h>

#include "bundlesMinf.h"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class BundlesData
{

  private :
  float epsilon = 1e-5 ; 

  public :
  /////////////////////////// Common for all formats ///////////////////////////
  std::vector<float> matrixTracks ;
  std::vector<int32_t> pointsPerTrack ;
  std::vector<int64_t> fibersWithNans ;
  int64_t curves_count = 0 ;
  bool isBundles = false ;
  bool isTrk = false ;
  bool isTck = false ;

  ////////////////////////////// For .trk format ///////////////////////////////
  std::vector<std::vector<float>> tracksScalars ;
  std::vector<std::vector<float>> tracksProperties ;


  //////////////////////////////// Constructors ////////////////////////////////
  BundlesData() ;

  BundlesData( const char* bundlesdataFilename ) ;

  BundlesData( std::vector<float>& matrixTracks,
               std::vector<int32_t>& pointsPerTrack,
               std::vector<int64_t>& fibersWithNans,
               std::vector<std::vector<float>>& tracksScalars,
               std::vector<std::vector<float>>& tracksProperties,
               int64_t curves_count,
               bool isBundles,
               bool isTrk,
               bool isTck ) ;

  BundlesData( const BundlesData& bundlesData ) ;

  //////////////////////////////// Destructors /////////////////////////////////
  virtual ~BundlesData() ;

  ///////////////////////////////// Operators //////////////////////////////////
  bool operator == ( const BundlesData& bundleData ) const;
  float operator[]( int64_t index ) const ;

  ////////////////////////////////// Mehtods ///////////////////////////////////
  void read( const char* bundlesFilename ) ;

  void write( const char* bundlesFilename,
              const BundlesMinf& bundleInfo ) const ;

  void write( const char* bundlesFilename,
              const BundlesMinf& bundleInfo,
              const char* referenceAnatomyPath ) const ;

  std::string getFormat() const ;

  bool matrixTracksEquals( std::vector<float>& matrixTracks ) const ;


  float computeLengthFiber( const std::vector<float>& tractogramFibers,
                            int fiberIndex,
                            int nbPoints ) const ;
  float computeLengthFiber( const std::vector<float>& fiber ) const ;
  float computeLengthFiber( int fiberIndex ) const ;


  void resampleFiberEquidistant( const std::vector<float>& inputFiber,
                                 std::vector<float>& outputFiber,
                                 int nbPointsToResample ) const ;


  void computeMedialPointFiberTractogram(
                        const std::vector<float>& tractogramFibers,
                        int fiberIndex,
                        int nbPoints,
                        std::vector<float>& medialPointFiberTractogram ) const ;
  void computeMedialPointFiberTractogram(
                        int fiberIndex,
                        std::vector<float>& medialPointFiberTractogram ) const ;
  void computeMedialPointFiberTractogram(
                        const std::vector<float>& fiber,
                        std::vector<float>& medialPointFiberTractogram ) const ;

  int computeMedialPointFiberWithDistance(
                        const std::vector<float>& tractogramFibers,
                        int fiberIndex,
                        int nbPoints,
                        std::vector<float>& medialPointFiberTractogram ) const ;
  int computeMedialPointFiberWithDistance(
                        int fiberIndex,
                        std::vector<float>& medialPointFiberTractogram ) const ;
  int computeMedialPointFiberWithDistance(
                        const std::vector<float>& fiber,
                        std::vector<float>& medialPointFiberTractogram ) const ;

  void computeGravityCenterFiber(
                              const std::vector<float>& tractogramFibers,
                              int fiberIndex,
                              int nbPoints,
                              std::vector<float>& gravityCenterFiber ) const ;
  void computeGravityCenterFiber(
                              int fiberIndex,
                              std::vector<float>& gravityCenterFiber ) const ;
  void computeGravityCenterFiber(
                              const std::vector<float>& fiber,
                              std::vector<float>& gravityCenterFiber ) const ;


  void computeNormalVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& normalVector ) const ;
  void computeNormalVectorFiberTractogram(
                                      int fiberIndex,
                                      std::vector<float>& normalVector ) const ;
  void computeNormalVectorFiberTractogram(
                                      const std::vector<float>& fiber,
                                      std::vector<float>& normalVector ) const ;

  void computeDirectionVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::vector<float>& medialPointFiberTractogram,
                           const std::vector<float>& normalVector,
                           int fiberIndex,
                           int nbPoints,
                           std::vector<float>& directionVector ) const ;
  void computeDirectionVectorFiberTractogram(
                                   int fiberIndex,
                                   const std::vector<float>& normalVector,
                                   std::vector<float>& directionVector ) const ;
  void computeDirectionVectorFiberTractogram(
                                   const std::vector<float>& fiber,
                                   const std::vector<float>& normalVector,
                                   std::vector<float>& directionVector ) const ;

  float computeMDADBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) const ;
  float computeMDADBetweenTwoFibers(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int fiberIndex_1,
                                   int fiberIndex_2,
                                   int nbPoints ) const ;

  float computeMDFBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::vector<float>& medialPointTractFiber_1,
                              const std::vector<float>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) const ;
  float computeMDFBetweenTwoFibers(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int fiberIndex_1,
                                   int fiberIndex_2,
                                   int nbPoints ) const ;

  float compareDisimilarityBundles(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int nbFibersTract_1,
                                   int nbFibersTract_2,
                                   int nbPoints1,
                                   int nbPoints2 ) const ;
  float compareDisimilarityBundles( const std::vector<float>& tractogramFibers,
                                    int nbFibersTract,
                                    int nbPoints ) const ;

  double distanceBetweenBundles( const std::vector<float>& bundle1,
                                 const std::vector<float>& bundle2,
                                 int nbFibersBundle1,
                                 int nbFibersBundle2,
                                 int nbPointsTract_1,
                                 int nbPointsTract_2 ) const ;
  double distanceBetweenBundles( const std::vector<float>& bundle,
                                 int nbFibersBundle,
                                 int nbPoints ) const ;

  void computeNumberAdjacentFibersBundle1toBundle2(
                              const BundlesData& bundle1,
                              const BundlesData& bundle2,
                              float threshold,
                              std::vector<int>& nbAdjacentFibersBundle ) const ;

  void computeNumberAdjacentFibersBundle1toBundle2(
                              const std::vector<float>& bundle1,
                              const std::vector<float>& bundle2,
                              int nbFibersBundle1,
                              int nbFibersBundle2,
                              int nbPoints,
                              float threshold,
                              std::vector<int>& nbAdjacentFibersBundle ) const ;

  float coverageBundle1toBundle2( const BundlesData& bundle1,
                                  const BundlesData& bundle2,
                                  float threshold ) const ;

  float overlapBundle1toBundle2( const BundlesData& bundle1,
                                 const BundlesData& bundle2,
                                 float threshold ) const ;


  float bundlesAdjacency( const BundlesData& bundle1,
                          const BundlesData& bundle2,
                          float threshold ) const ;

  float computeAngleBetweenVectors( const std::vector<float>& vector1,
                                    const std::vector<float>& vector2 ) const ;

  float computeAngleBetweenPlanes( const std::vector<float>& vector1,
                                   const std::vector<float>& vector2 ) const ;

  float computeAngleBetweenDirections(
                                     const std::vector<float>& vector1,
                                     const std::vector<float>& vector2 ) const ;


  void projectVectorToPlane( const std::vector<float>& normalToPlane,
                             const std::vector<float>& inputVector,
                             std::vector<float>& projectedVector ) const ;

  void computeRotationMatrixFromVectorAndAngle(
                                    const std::vector<float>& vector,
                                    float angle, // in radians
                                    std::vector<float>& rotationMatrix ) const ;

  void applyRotationMatrixToVector( const std::vector<float>& vector,
                                    const std::vector<float>& rotationMatrix,
                                    std::vector<float>& rotatedVector ) const ;

  void applyRotationMatrixToFiber( const std::vector<float>& fiber,
                                   const std::vector<float>& rotationMatrix,
                                   const std::vector<float>& medialPointFiber,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const BundlesData& inputBundlesData,
                                   const std::vector<float>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const std::vector<float>& tractogramFibers,
                                   const std::vector<float>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const std::vector<float>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const std::vector<float>& fiber,
                                   const std::vector<float>& rotationMatrix,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;

  void putFibersInSamePlane( const std::vector<float>& normalVector1,
                             const std::vector<float>& normalVector2,
                             const std::vector<float>& tractogramFibers2,
                             int fiberIndex2,
                             int nbPoints,
                             std::vector<float>& fiber2ToPlane1,
                             std::vector<float>& newNormalVectorFiber2 ) const ;
  void putFibersInSamePlane( const std::vector<float>& normalVector1,
                             const std::vector<float>& normalVector2,
                             const std::vector<float>& fiber2,
                             int nbPoints,
                             std::vector<float>& fiber2ToPlane1,
                             std::vector<float>& newNormalVectorFiber2 ) const ;

  void putFiberInPlaneXY( const std::vector<float>& normalVector,
                          const std::vector<float>& tractogramFibers,
                          int fiberIndex,
                          int nbPoints,
                          std::vector<float>& fiberToPlaneXY,
                          std::vector<float>& newNormalVectorFiber ) const ;
  void putFiberInPlaneXY( const std::vector<float>& normalVector,
                          const std::vector<float>& fiber,
                          int nbPoints,
                          std::vector<float>& fiberToPlaneXY,
                          std::vector<float>& newNormalVectorFiber ) const ;

  void putFibersInSameDirection(
                                const std::vector<float>& normalVector1,
                                const std::vector<float>& normalVector2,
                                const std::vector<float>& directionVector1,
                                const std::vector<float>& fiber2,
                                int nbPoints,
                                std::vector<float>& fiber2ToDirection1 ) const ;

  void registerFiber( const std::vector<float>& tractogramFibers2,
                      const std::vector<float>& normalVector1,
                      const std::vector<float>& normalVector2,
                      const std::vector<float>& directionVector1,
                      const std::vector<float>& medialPointFiber1,
                      const std::vector<float>& medialPointFiber2,
                      int fiberIndexTractogram2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2 ) const ;
  void registerFiber( const std::vector<float>& tractogramFibers1,
                      const std::vector<float>& tractogramFibers2,
                      int fiberIndexTractogram1,
                      int fiberIndexTractogram2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2 ) const ;
  void registerFiber( const std::vector<float>& fiber2,
                      const std::vector<float>& normalVector1,
                      const std::vector<float>& normalVector2,
                      const std::vector<float>& directionVector1,
                      const std::vector<float>& medialPointFiber1,
                      const std::vector<float>& medialPointFiber2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2 ) const ;
  void registerFiber( const std::vector<float>& fiber1,
                      const std::vector<float>& fiber2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::vector<float>& newNormalVectorFiber2 ) const ;


  void registerFiberToPlaneXYAndDirectionX(
                              const std::vector<float>& tractogramFibers,
                              const std::vector<float>& normalVector,
                              const std::vector<float>& medialPointFiber,
                              int fiberIndexTractogram,
                              int nbPoints,
                              std::vector<float>& fiberToPlaneXYAndDirectionX,
                              std::vector<float>& newNormalVectorFiber ) const ;
  void registerFiberToPlaneXYAndDirectionX(
                              const std::vector<float>& fiber,
                              const std::vector<float>& normalVector,
                              const std::vector<float>& medialPointFiber,
                              int nbPoints,
                              std::vector<float>& fiberToPlaneXYAndDirectionX,
                              std::vector<float>& newNormalVectorFiber ) const ;

  void getFiberFromTractogram( const std::vector<float>& tractogramFibers,
                               int fiberIndex,
                               int nbPoints,
                               std::vector<float>& fiber ) const ;

  void getFiberWithVectors( const std::vector<float>& fiber,
                            const std::vector<float>& referenceFiber,
                            int nbPoints,
                            std::vector<float>& fiberWithVectors ) const ;

  void getFiberWithVectors( const std::vector<float>& fiber,
                            int nbPoints,
                            std::vector<float>& fiberWithVectors ) const ;



  ///////// Optimization functions /////////
  void getMatrixAndVectorForLeastSquares(
                               const std::vector<float>& tractogramFibers,
                               int fiberIndex,
                               int nbPoints,
                               Eigen::MatrixXf& AMatrix,
                               Eigen::VectorXf& BVector ) const ;


  ///////// Algebra functions /////////

  // Point1 is the point where the vector is pointing to
  void vectorFromPoints( const std::vector<float>& point1,
                         const std::vector<float>& point2,
                         std::vector<float>& vector ) const ;

  float scalarProduct( const std::vector<float>& vector1,
                       const std::vector<float>& vector2 ) const ;

  float normVector( const std::vector<float>& vector ) const ;

  void normalizeVector( std::vector<float>& vector ) const ;

  void translatePoint( const std::vector<float>& point,
                       const std::vector<float>& unitaryTranslationVector,
                       float translationDistance,
                       std::vector<float>& translatedPoint ) const ;
  void translatePoint( const std::vector<float>& point,
                       const std::vector<float>& translationVector,
                       std::vector<float>& translatedPoint ) const ;

  void crossProduct( const std::vector<float>& vector1,
                     const std::vector<float>& vector2,
                     std::vector<float>& result ) const ;

  int checkIfFiberPointCanBeResampled(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int point,
                                     int nbPoints ) const ;
  void resampleFiberWithNan( std::vector<float>& tractogramFibers,
                             int fiberIndex,
                             int nbPoints ) const ;


} ;
