cmake_minimum_required(VERSION  3.20.2)

set(CMAKE_BUILD_TYPE Release)

# Set the project name and version
project(ProjectAtlas VERSION 1.0 LANGUAGES CXX)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH ${CMAKE_BINARY_DIR} FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Configure a header file to pass the version number to the source code
configure_file(ProjectAtlasConfig.h.in ProjectAtlasConfig.h)

# Set PYTHON_BINARY
set(PYTHON_BINARY "/usr/bin/python3")



######################################## Boost #########################################
find_package(Boost REQUIRED)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    add_definitions( "-DHAS_BOOST" )
endif()
find_package(Boost REQUIRED COMPONENTS system)


######################################### Eigen #########################################
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

######################################### openMP ########################################
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


#################################### Configure files ####################################
# file(GLOB APPS_DIR apps/*)
# foreach(f ${APPS_DIR})
#     get_filename_component(sourceName ${f} NAME)
#     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/apps/${sourceName} ${CMAKE_CURRENT_SOURCE_DIR}/apps-config/${sourceName} FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ )
# endforeach(f)

######################################## LBFGSpp ########################################
add_subdirectory(LBFGSpp-master)

################################# Adding subdirectories #################################
add_subdirectory(src)
# add_subdirectory(apps-config)
add_subdirectory(apps-config)
add_subdirectory(doc)
