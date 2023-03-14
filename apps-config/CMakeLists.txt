include_directories ( ${CMAKE_SOURCE_DIR}/src)

file( GLOB LIB_SOURCE *.cxx)

foreach(f ${LIB_SOURCE})
    get_filename_component(sourceName ${f} NAME_WE)
    add_executable(${sourceName} ${f})
    target_link_libraries(${sourceName}
                          stdc++fs
                          ioWrapper
                          RecognizedBundles
                          Eigen3::Eigen
                          lbfgspp
                          ${Boost_LIBRARIES}
                          OpenMP::OpenMP_CXX)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${sourceName}
            DESTINATION bin)
endforeach(f)
configure_file(analyseAtlasBundle.py ${CMAKE_BINARY_DIR}/bin/analyseAtlasBundle.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(dipyServer.py  ${CMAKE_BINARY_DIR}/bin/dipyServer.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(clientRegisterBundles.py ${CMAKE_BINARY_DIR}/bin/clientRegisterBundles.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(clientComputeCentroids.py ${CMAKE_BINARY_DIR}/bin/clientComputeCentroids.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(clientCloseServer.py ${CMAKE_BINARY_DIR}/bin/clientCloseServer.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(scoresPredictedSGT.py ${CMAKE_BINARY_DIR}/bin/scoresPredictedSGT.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(applySupWMA.py ${CMAKE_BINARY_DIR}/bin/applySupWMA.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(piecewiseLinearRegressionDtiMeasures.py ${CMAKE_BINARY_DIR}/bin/piecewiseLinearRegressionDtiMeasures.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(extract_bundles_feat.py ${CMAKE_BINARY_DIR}/bin/extract_bundles_feat.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(recoBundlesLabelsToGeoLab.py ${CMAKE_BINARY_DIR}/bin/recoBundlesLabelsToGeoLab.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(statisticsComparisonWithAtlas.py ${CMAKE_BINARY_DIR}/bin/statisticsComparisonWithAtlas.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(countFoundBundles.py ${CMAKE_BINARY_DIR}/bin/countFoundBundles.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(showColorCodedBundlesByAge.py ${CMAKE_BINARY_DIR}/bin/showColorCodedBundlesByAge.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(plotDistributionBundleByAge.py ${CMAKE_BINARY_DIR}/bin/plotDistributionBundleByAge.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(showColorCodedFibersPerBundle.py ${CMAKE_BINARY_DIR}/bin/showColorCodedFibersPerBundle.py FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
