#ifndef HEADER
#define HEADER
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include <fstream>


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

#include <boost/shared_array.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

using std::vector;
using std::cout;
using std::endl;
using std::string;

using namespace boost::geometry;

using namespace boost::program_options;
using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)
typedef model::point<double, 2, cs::cartesian> pointxy;
#endif
