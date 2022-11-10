cmake_minimum_required(VERSION  3.20.2)

set(CMAKE_BUILD_TYPE Release)

# Set the project name and version
project(ProjectAtlas VERSION 1.0 LANGUAGES CXX)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH ${CMAKE_BINARY_DIR} FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


# Configure a header file to pass the version number to the source code
configure_file(ProjectAtlasConfig.h.in ProjectAtlasConfig.h)


######### Boost #########
find_package(Boost REQUIRED)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    add_definitions( "-DHAS_BOOST" )
endif()
find_package(Boost REQUIRED COMPONENTS system)
##########################

########## Eigen #########
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

######### openMP #########
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
###########################

######### LBFGSpp #########
add_subdirectory(LBFGSpp-master)
###########################
add_subdirectory(ProjectAtlasFiles)
add_subdirectory(src)
add_subdirectory(doc)
