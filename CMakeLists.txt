# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project(main_test)

option(BUILD_WITH_TBB "build with tbb" OFF)


set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE True)
set(CMAKE_BUILD_TYPE Release)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
endif()

# Boost and its components
find_package( Boost QUIET )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

# include for local directory
# include for local package

add_executable(${PROJECT_NAME} main.cpp)

find_package(Eigen3 REQUIRED NO_MODULE)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

find_package(GMP REQUIRED)

# glog
find_package(glog REQUIRED)
if (glog_FOUND )
  message(STATUS "Found Glog.")
  target_link_libraries(${PROJECT_NAME} glog::glog)
else()
  message(STATUS "No GLog.")
endif()
# include_directories(${GLOG_INCLUDE_DIRS}) 

# tbb
if (BUILD_WITH_TBB)
  find_package(TBB REQUIRED)
  if (TBB_FOUND)
    message(STATUS "Found TBB")
    target_link_libraries(${PROJECT_NAME} tbb::tbb)
  else()
    message(STATUS "No TBB")
  endif()
endif()

target_link_libraries(${PROJECT_NAME}
                      ${gmp_LIBRARIES})


# create_single_source_cgal_program( "main.cpp" )
