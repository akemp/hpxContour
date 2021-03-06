# Require a recent version of cmake
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# This project is C++ based.
project(HPXContour CXX)

# This adds the HPX cmake configuration directory to the search path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    $ENV{HPX_DIR}/share/cmake-2.8/Modules)

# Instruct cmake to find the HPX settings
find_package(HPX)

# Adjust auxiliary build settings for HPX to work properly
if(NOT MSVC)
  set(HPX_RPATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}${HPX_RPATH})
else()
  # enable solution folders for MSVC
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  if(HPX_NO_INSTALL)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/$(Configuration)" CACHE PATH
      "Prefix prepended to install directories." FORCE)
  endif()
endif()

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})

find_package(HPX_NETCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIR})
link_directories(${NETCDF_LIBRARY_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
    ${CMAKE_SOURCE_DIR}/cmake)

find_package(HPX_SHAPELIB REQUIRED)
include_directories(${SHAPELIB_INCLUDE_DIR})
link_directories(${SHAPELIB_LIBRARY_DIR})

	
set (sources 
	contour_lines/contour.cpp
	contour_lines/datastructures.hpp
	contour_lines/gzip.cpp
	contour_lines/io.cpp
	contour_lines/retLine.cpp
	contour_lines/treeHeaders.cpp
	contour_lines/treeHeaders.hpp
	contour_lines/zlib.cpp
	)
# add example executable
add_hpx_executable(contour_lines
             SOURCES ${sources}
             COMPONENT_DEPENDENCIES iostreams netcdf shapelib)