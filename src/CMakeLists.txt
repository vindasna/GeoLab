add_library(ioWrapper ioWrapper.cxx ioWrapper.h)
add_library(NiftiImage niftiImage.cxx niftiImage.h)
add_library(BundlesMinf bundlesMinf.cxx bundlesMinf.h)
add_library(BundlesData bundlesData.cxx bundlesData.h)
add_library(AtlasBundles AtlasBundles.cxx AtlasBundles.h)
add_library(RecognizedBundles RecognizedBundles.cxx RecognizedBundles.h)


target_link_libraries(BundlesData BundlesMinf NiftiImage Eigen3::Eigen lbfgspp)
target_link_libraries(AtlasBundles BundlesData stdc++fs)
target_link_libraries(RecognizedBundles AtlasBundles)
