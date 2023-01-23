include_directories ( ${CMAKE_SOURCE_DIR}/ProjectAtlasFiles)

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
