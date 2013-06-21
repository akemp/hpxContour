# Require a recent version of cmake
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# This project is C++ based.
project(HPXContour CXX)

# This adds the HPX cmake configuration directory to the search path.
set(CMAKE_MODULE_PATH
    $ENV{HPX_DIR}/share/cmake-2.8/Modules)

# Instruct cmake to find the HPX settings
find_package(HPX REQUIRED)

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})

find_package(HPX_SHAPELIB)
find_package(HPX_NETCDF)

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
             ${PLAYGROUND_FLAGS} 
             FOLDER 
             COMPONENT_DEPENDENCIES ${hpx_RUNTIME_LIBRARIES}  shapelib netcdf )
