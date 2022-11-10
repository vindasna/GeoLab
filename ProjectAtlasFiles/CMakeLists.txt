add_library(BundlesFormat BundlesFormat.cxx BundlesFormat.h)
add_library(BundlesDataFormat BundlesDataFormat.cxx BundlesDataFormat.h)
add_library(TrkFormat TrkFormat.cxx TrkFormat.h)
add_library(AtlasBundles AtlasBundles.cxx AtlasBundles.h)
add_library(RecognizedBundles RecognizedBundles.cxx RecognizedBundles.h)
add_library(BMDDistance BMDDistance.cxx BMDDistance.h)


target_link_libraries(BundlesDataFormat Eigen3::Eigen lbfgspp)
target_link_libraries(TrkFormat BundlesFormat BundlesDataFormat)
target_link_libraries(AtlasBundles TrkFormat stdc++fs)
target_link_libraries(RecognizedBundles AtlasBundles)
target_link_libraries(BMDDistance RecognizedBundles)
