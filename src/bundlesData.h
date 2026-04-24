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
#include <array>
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

enum class TransformType {
    RIGID,
    AFFINE
};

struct Cluster {
    std::vector<float> centroid; // Flattened array of 3 * numPoints
    int count;                   // Number of streamlines in this cluster
    std::vector<int> indices;    // Original indices of the streamlines
};

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
                        std::array<float, 3>& medialPointFiberTractogram ) const ;
  void computeMedialPointFiberTractogram(
                        int fiberIndex,
                        std::array<float, 3>& medialPointFiberTractogram ) const ;
  void computeMedialPointFiberTractogram(
                        const std::vector<float>& fiber,
                        std::array<float, 3>& medialPointFiberTractogram ) const ;

  int computeMedialPointFiberWithDistance(
                        const std::vector<float>& tractogramFibers,
                        int fiberIndex,
                        int nbPoints,
                        std::array<float, 3>& medialPointFiberTractogram ) const ;
  int computeMedialPointFiberWithDistance(
                        int fiberIndex,
                        std::array<float, 3>& medialPointFiberTractogram ) const ;
  int computeMedialPointFiberWithDistance(
                        const std::vector<float>& fiber,
                        std::array<float, 3>& medialPointFiberTractogram ) const ;

  void computeGravityCenterFiber(
                              const std::vector<float>& tractogramFibers,
                              int fiberIndex,
                              int nbPoints,
                              std::array<float, 3>& gravityCenterFiber ) const ;
  void computeGravityCenterFiber(
                              int fiberIndex,
                              std::array<float, 3>& gravityCenterFiber ) const ;
  void computeGravityCenterFiber(
                              const std::vector<float>& fiber,
                              std::array<float, 3>& gravityCenterFiber ) const ;


  void computeNormalVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::array<float, 3>& medialPointFiberTractogram,
                           int fiberIndex,
                           int nbPoints,
                           std::array<float, 3>& normalVector ) const ;
  void computeNormalVectorFiberTractogram(
                                      int fiberIndex,
                                      std::array<float, 3>& normalVector ) const ;
  void computeNormalVectorFiberTractogram(
                                      const std::vector<float>& fiber,
                                      std::array<float, 3>& normalVector ) const ;

  void computeDirectionVectorFiberTractogram(
                           const std::vector<float>& tractogramFibers,
                           const std::array<float, 3>& medialPointFiberTractogram,
                           const std::array<float, 3>& normalVector,
                           int fiberIndex,
                           int nbPoints,
                           std::array<float, 3>& directionVector ) const ;
  void computeDirectionVectorFiberTractogram(
                                   int fiberIndex,
                                   const std::array<float, 3>& normalVector,
                                   std::array<float, 3>& directionVector ) const ;
  void computeDirectionVectorFiberTractogram(
                                   const std::vector<float>& fiber,
                                   const std::array<float, 3>& normalVector,
                                   std::array<float, 3>& directionVector ) const ;

  float computeMDADBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::array<float, 3>& medialPointTractFiber_1,
                              const std::array<float, 3>& medialPointTractFiber_2,
                              int fiberIndex_1,
                              int fiberIndex_2,
                              int nbPoints ) const ;
  float computeMDADBetweenTwoFibers(
                                   const std::vector<float>& tractogramFibers_1,
                                   const std::vector<float>& tractogramFibers_2,
                                   int fiberIndex_1,
                                   int fiberIndex_2,
                                   int nbPoints ) const ;
  

  
  float computeMDADBetweenTwoFibersAfterAlignement(
              const std::vector<float>& tractogramFibers_1_to_tractogramFibers2,
              const std::vector<float>& tractogramFibers_2,
              int fiberIndex_1,
              int fiberIndex_2,
              bool useMean,
              int nbPoints ) const ;

  float computeMDFBetweenTwoFibers(
                              const std::vector<float>& tractogramFibers_1,
                              const std::vector<float>& tractogramFibers_2,
                              const std::array<float, 3>& medialPointTractFiber_1,
                              const std::array<float, 3>& medialPointTractFiber_2,
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

  float computeAngleBetweenVectors( const std::array<float, 3>& vector1,
                                    const std::array<float, 3>& vector2 ) const ;

  float computeAngleBetweenPlanes( const std::array<float, 3>& vector1,
                                   const std::array<float, 3>& vector2 ) const ;

  float computeAngleBetweenDirections(
                                     const std::array<float, 3>& vector1,
                                     const std::array<float, 3>& vector2 ) const ;


  void projectVectorToPlane( const std::array<float, 3>& normalToPlane,
                             const std::array<float, 3>& inputVector,
                             std::array<float, 3>& projectedVector ) const ;

  void computeRotationMatrixFromVectorAndAngle(
                                    const std::array<float, 3>& vector,
                                    float angle, // in radians
                                    std::array<float, 9>& rotationMatrix ) const ;

  void applyRotationMatrixToVector( const std::array<float, 3>& vector,
                                    const std::array<float, 9>& rotationMatrix,
                                    std::array<float, 3>& rotatedVector ) const ;

  void applyRotationMatrixToFiber( const std::vector<float>& fiber,
                                   const std::array<float, 9>& rotationMatrix,
                                   const std::array<float, 3>& medialPointFiber,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const BundlesData& inputBundlesData,
                                   const std::array<float, 9>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const std::vector<float>& tractogramFibers,
                                   const std::array<float, 9>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const std::array<float, 9>& rotationMatrix,
                                   int fiberIndex,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;
  void applyRotationMatrixToFiber( const std::vector<float>& fiber,
                                   const std::array<float, 9>& rotationMatrix,
                                   int nbPoints,
                                   std::vector<float>& rotatedFiber ) const ;

  void putFibersInSamePlane( const std::array<float, 3>& normalVector1,
                             const std::array<float, 3>& normalVector2,
                             const std::vector<float>& tractogramFibers2,
                             int fiberIndex2,
                             int nbPoints,
                             std::vector<float>& fiber2ToPlane1,
                             std::array<float, 3>& newNormalVectorFiber2 ) const ;
  void putFibersInSamePlane( const std::array<float, 3>& normalVector1,
                             const std::array<float, 3>& normalVector2,
                             const std::vector<float>& fiber2,
                             int nbPoints,
                             std::vector<float>& fiber2ToPlane1,
                             std::array<float, 3>& newNormalVectorFiber2 ) const ;

  void putFiberInPlaneXY( const std::array<float, 3>& normalVector,
                          const std::vector<float>& tractogramFibers,
                          int fiberIndex,
                          int nbPoints,
                          std::vector<float>& fiberToPlaneXY,
                          std::array<float, 3>& newNormalVectorFiber ) const ;
  void putFiberInPlaneXY( const std::array<float, 3>& normalVector,
                          const std::vector<float>& fiber,
                          int nbPoints,
                          std::vector<float>& fiberToPlaneXY,
                          std::array<float, 3>& newNormalVectorFiber ) const ;

  void putFibersInSameDirection(
                                const std::array<float, 3>& normalVector1,
                                const std::array<float, 3>& normalVector2,
                                const std::array<float, 3>& directionVector1,
                                const std::vector<float>& fiber2,
                                int nbPoints,
                                std::vector<float>& fiber2ToDirection1 ) const ;

  void registerFiber( const std::vector<float>& tractogramFibers2,
                      const std::array<float, 3>& normalVector1,
                      const std::array<float, 3>& normalVector2,
                      const std::array<float, 3>& directionVector1,
                      const std::array<float, 3>& medialPointFiber1,
                      const std::array<float, 3>& medialPointFiber2,
                      int fiberIndexTractogram2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::array<float, 3>& newNormalVectorFiber2 ) const ;
  void registerFiber( const std::vector<float>& tractogramFibers1,
                      const std::vector<float>& tractogramFibers2,
                      int fiberIndexTractogram1,
                      int fiberIndexTractogram2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::array<float, 3>& newNormalVectorFiber2 ) const ;
  void registerFiber( const std::vector<float>& fiber2,
                      const std::array<float, 3>& normalVector1,
                      const std::array<float, 3>& normalVector2,
                      const std::array<float, 3>& directionVector1,
                      const std::array<float, 3>& medialPointFiber1,
                      const std::array<float, 3>& medialPointFiber2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::array<float, 3>& newNormalVectorFiber2 ) const ;
  void registerFiber( const std::vector<float>& fiber1,
                      const std::vector<float>& fiber2,
                      int nbPoints,
                      std::vector<float>& fiber2Tofiber1,
                      std::array<float, 3>& newNormalVectorFiber2 ) const ;


  void registerFiberToPlaneXYAndDirectionX(
                              const std::vector<float>& tractogramFibers,
                              const std::array<float, 3>& normalVector,
                              const std::array<float, 3>& medialPointFiber,
                              int fiberIndexTractogram,
                              int nbPoints,
                              std::vector<float>& fiberToPlaneXYAndDirectionX,
                              std::array<float, 3>& newNormalVectorFiber ) const ;
  void registerFiberToPlaneXYAndDirectionX(
                              const std::vector<float>& fiber,
                              const std::array<float, 3>& normalVector,
                              const std::array<float, 3>& medialPointFiber,
                              int nbPoints,
                              std::vector<float>& fiberToPlaneXYAndDirectionX,
                              std::array<float, 3>& newNormalVectorFiber ) const ;

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
  void vectorFromPoints( const std::array<float, 3>& point1,
                         const std::array<float, 3>& point2,
                         std::array<float, 3>& vector ) const ;

  float scalarProduct( const std::array<float, 3>& vector1,
                       const std::array<float, 3>& vector2 ) const ;

  float normVector( const std::array<float, 3>& vector ) const ;

  void normalizeVector( std::array<float, 3>& vector ) const ;

  void translatePoint( const std::array<float, 3>& point,
                       const std::array<float, 3>& unitaryTranslationVector,
                       float translationDistance,
                       std::array<float, 3>& translatedPoint ) const ;
  void translatePoint( const std::array<float, 3>& point,
                       const std::array<float, 3>& translationVector,
                       std::array<float, 3>& translatedPoint ) const ;

  void crossProduct( const std::array<float, 3>& vector1,
                     const std::array<float, 3>& vector2,
                     std::array<float, 3>& result ) const ;

  int checkIfFiberPointCanBeResampled(
                                     const std::vector<float>& tractogramFibers,
                                     int fiberIndex,
                                     int point,
                                     int nbPoints ) const ;
  void resampleFiberWithNan( std::vector<float>& tractogramFibers,
                             int fiberIndex,
                             int nbPoints ) const ;



  void registerTo(const BundlesData& target, int iterations = 10);

  std::vector<Cluster> computeQuickBundles(float distance_threshold);


  void resampleAll(int targetPoints);
  Eigen::Matrix4f registerFast(BundlesData& target, 
                                  TransformType regType = TransformType::RIGID, 
                                  float qbThreshold = 10.0f, 
                                  int iterations = 10, 
                                  float outlierRejectionPct = 0.20f,
                                  BundlesData* outOriginalMovingCentroids = nullptr,
                                  BundlesData* outMovingCentroids = nullptr,
                                  BundlesData* outTargetCentroids = nullptr);

} ;
