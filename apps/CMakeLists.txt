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
configure_file(analyseAtlasBundle.py ${CMAKE_BINARY_DIR}/bin/analyseAtlasBundle.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(dipyServer.py  ${CMAKE_BINARY_DIR}/bin/dipyServer.py  COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(clientRegisterBundles.py ${CMAKE_BINARY_DIR}/bin/clientRegisterBundles.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(clientComputeCentroids.py ${CMAKE_BINARY_DIR}/bin/clientComputeCentroids.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(clientCloseServer.py ${CMAKE_BINARY_DIR}/bin/clientCloseServer.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(scoresPredictedSGT.py ${CMAKE_BINARY_DIR}/bin/scoresPredictedSGT.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(applySupWMA.py ${CMAKE_BINARY_DIR}/bin/applySupWMA.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(piecewiseLinearRegressionDtiMeasures.py ${CMAKE_BINARY_DIR}/bin/piecewiseLinearRegressionDtiMeasures.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
configure_file(extract_bundles_feat.py ${CMAKE_BINARY_DIR}/bin/extract_bundles_feat.py COPYONLY FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
