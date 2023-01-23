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
